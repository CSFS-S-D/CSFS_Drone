library(lidR)
library(magrittr)
library(data.table)
library(StereoMorph) # Find the distance from a point to a line
library(lidRplugins)
library(ggplot2)
library(rgl)

header = readLASheader("C:/Users/C820392777/Documents/DJI/DJITerra/csfs_drones@colostate.edu/AbbotRidgeHQ_FRWRM_121824/lidars/terra_las/cloud8fa0f164801d49e9.las")
xt=ext(header)
xt[2]-(xt[2]-xt[1])/2
xt[4] -(xt[4]-xt[3])/2
las = readLAS("C:/Users/C820392777/Documents/DJI/DJITerra/csfs_drones@colostate.edu/AbbotRidgeHQ_FRWRM_121824/lidars/terra_las/cloud8fa0f164801d49e9.las",
              filter="-keep_xy 495370 4110460 495410 4110500")


dec_las = decimate_points(las, random_per_voxel(0.25, 2))
plot(dec_las, color="RGB")

norm_las = dec_las %>% classify_ground(csf()) %>%
  normalize_height(knnidw()) %>%
  filter_poi(Classification!=2)
plot(norm_las)

#################################
######      Functions      ######
#################################
# count of returns in a voxel
vox_mt <- function(z)
{
  SViD = length(z)
  
  metrics =list(
    SViD = SViD,
    z_list = list(z)
  )
  return(metrics)}


vox_mt2 <- function(vox, res, max_z)
{
  fullvox <- data.table(vox,
                        id_xy = paste0(vox$X,"-", vox$Y),
                        id_xyz = paste0(vox$X,"-", vox$Y,"-",vox$Z))
  if(is.numeric(max_z) & max_z != max(fullvox$Z)){
    
    z_need <- seq(max(fullvox$Z), max_z, res)[-1]
    new_z <- lapply(z_need, function(x){paste0(unique(fullvox$id_xy),"-", x)})
    new_id <- unlist(new_z)
    empty_z <- data.frame(X = do.call(rbind, lapply(strsplit(new_id, "-"), `[`, 1)),
                          Y = do.call(rbind, lapply(strsplit(new_id, "-"), `[`, 2)),
                          Z = do.call(rbind, lapply(strsplit(new_id, "-"), `[`, 3)),
                          SViD = 0,
                          z_list = NA,
                          id_xy = paste0(do.call(rbind, lapply(strsplit(new_id, "-"), `[`, 1)),
                                         "-",
                                         do.call(rbind, lapply(strsplit(new_id, "-"), `[`, 2))),
                          id_xyz = new_id
    )
    fullvox <- rbind(fullvox, empty_z)
  }
  
  fullvox$SViD[is.na(fullvox$SViD)] <- 0
  fullvox[unlist(lapply(fullvox$z_list, is.null)),]$z_list <- NA
  
  ### P_Di_above is the number of returns above each voxel (Kim et al. 2016 uses returns above as P_Di)
  # calculate points above each voxel
  point_above <- fullvox[,as.list(lapply(unique(Z), function(x){
    return(list(
      m_z = median(unlist(z_list[Z >x]), na.rm = TRUE),
      p_a = sum(SViD[Z > x]),
      p_s = sum(SViD[Z == x]*SViD[Z==x+1])))}
  )
  ),
  by = list(X,Y)]
  
  colnames(point_above)[3:ncol(point_above)] <- paste0(unique(fullvox$Z))
  points_am <- point_above[seq(1, nrow(point_above), 3), ]
  points_am <- suppressWarnings(data.table::as.data.table(lapply(points_am,  function(x) unlist(x))))
  points_am[is.na(points_am)] <- 0
  points_am <- melt(points_am, id.vars = c("X","Y"), variable.name = "Z", value.name = "pa_med")
  points_am$id_xyz = paste0(points_am$X,"-",points_am$Y,"-",points_am$Z)
  
  points_ab <- point_above[seq(2, nrow(point_above), 3), ]
  points_ab <- suppressWarnings(data.table::as.data.table(lapply(points_ab,  function(x) unlist(x))))
  points_ab <- melt(points_ab, id.vars = c("X","Y"), variable.name = "Z", value.name = "npoints_above")
  points_ab$id_xyz = paste0(points_ab$X,"-",points_ab$Y,"-",points_ab$Z)
  
  points_sp <- point_above[seq(3, nrow(point_above), 3), ]
  points_sp <- suppressWarnings(data.table::as.data.table(lapply(points_sp,  function(x) unlist(x))))
  points_sp <- melt(points_sp, id.vars = c("X","Y"), variable.name = "Z", value.name = "spi")
  points_sp$id_xyz = paste0(points_sp$X,"-",points_sp$Y,"-",points_sp$Z)
  
  points_ab <- merge(points_ab,points_am)
  points_ab <- merge(points_ab,points_sp)
  fullvox <- merge(fullvox, points_ab[,c(4,5,6,7)], by = "id_xyz")
  
  metrics = data.table(
    X = as.numeric(fullvox$X),
    Y = as.numeric(fullvox$Y),
    Z = as.numeric(fullvox$Z),
    vox_res = res,
    PDiAbove = fullvox$npoints_above,
    SPI = fullvox$spi
    
  )
  return(metrics)
}


vox_sum <- function(SViD,FRDi, FRSVi, PDiBelow, PDiAbove, SVsum, SVMmdAbove, SVMmd, SVM, SPI, vox_res){
  las_vox <- data.table(SViD,FRDi, FRSVi, PDiBelow, PDiAbove, SVsum, SVMmdAbove, SVMmd, SVM, SPI)
  SVM <- las_vox[las_vox$SVM == 1]
  SVM$SVM <- NULL
  means <- apply(SVM[,6:ncol(SVM)], 2, mean, na.rm = TRUE)
  names(means)<- paste0(names(means), "_", vox_res, "_","mean")
  meds <- apply(SVM[,6:ncol(SVM)], 2, median, na.rm = TRUE)
  names(meds)<- paste0(names(meds), "_", vox_res, "_", "med")
  var <- apply(SVM[,6:ncol(SVM)], 2, var, na.rm = TRUE)
  names(var)<- paste0(names(var), "_",   vox_res, "_","var")
  sd <- apply(SVM[,6:ncol(SVM)], 2, sd, na.rm = TRUE)
  names(sd)<- paste0(names(sd), "_",  vox_res, "_","sd")
  cv <- apply(SVM[,6:ncol(SVM)], 2, function(x) {sd(x, na.rm = T)/mean(x, na.rm = T)})
  names(cv)<- paste0(names(cv), "_",  vox_res, "_","cv")
  IQR <- apply(SVM[,6:ncol(SVM)], 2,IQR, na.rm = T)
  names(IQR)<- paste0(names(IQR), "_",  vox_res, "_","IQR")
  skew <- apply(SVM[,6:ncol(SVM)], 2, function(x) {(sum((x - mean(x, na.rm = T))^3,na.rm = T)
                                                    /length(x))/(sum((x - mean(x, na.rm = T))^2,na.rm = T)
                                                                 /length(x))^(3/2)})
  names(skew)<- paste0(names(skew), "_",  vox_res, "_","skew")
  kurt <- apply(SVM[,6:ncol(SVM)], 2, function(x) {length(x)*sum((x- mean(x, na.rm = T))^4, na.rm = T)/
      (sum((x - mean(x, na.rm = T))^2, na.rm = T)^2)})
  names(kurt)<- paste0(names(kurt), "_", vox_res, "_", "kurt")
  sums <- apply(SVM[,6:ncol(SVM)], 2, sum, na.rm = TRUE)
  names(sums)<- paste0(names(sum), "_", vox_res, "_","sum")
  data_svm <- c(means, meds, var, sd, cv, IQR, skew, kurt)
  
  las_vox <- las_vox[,-c("SVsum", "SVMmd", "SVM", "SVMmdAbove")]
  means <- apply(las_vox, 2, mean, na.rm = TRUE)
  names(means)<- paste0(names(means), "_", vox_res, "_","mean")
  meds <- apply(las_vox, 2, median, na.rm = TRUE)
  names(meds)<- paste0(names(meds), "_", vox_res, "_", "med")
  var <- apply(las_vox, 2, var, na.rm = TRUE)
  names(var)<- paste0(names(var), "_",   vox_res, "_","var")
  sd <- apply(las_vox, 2, sd, na.rm = TRUE)
  names(sd)<- paste0(names(sd), "_",  vox_res, "_","sd")
  cv <- apply(las_vox, 2, function(x) {sd(x, na.rm = T)/mean(x, na.rm = T)})
  names(cv)<- paste0(names(cv), "_",  vox_res, "_","cv")
  IQR <- apply(las_vox, 2,IQR, na.rm = T)
  names(IQR)<- paste0(names(IQR), "_",  vox_res, "_","IQR")
  skew <- apply(las_vox, 2, function(x) {(sum((x - mean(x, na.rm = T))^3,na.rm = T)
                                          /length(x))/(sum((x - mean(x, na.rm = T))^2,na.rm = T)
                                                       /length(x))^(3/2)})
  names(skew)<- paste0(names(skew), "_",  vox_res, "_","skew")
  kurt <- apply(las_vox, 2, function(x) {length(x)*sum((x- mean(x, na.rm = T))^4, na.rm = T)/
      (sum((x - mean(x, na.rm = T))^2, na.rm = T)^2)})
  names(kurt)<- paste0(names(kurt), "_", vox_res, "_", "kurt")
  sums <- apply(las_vox, 2, sum, na.rm = TRUE)
  names(sums)<- paste0(names(sums), "_", vox_res, "_","sum")
  pct_fill_vox <- data.frame(pct_fill_vox = as.numeric(sum(las_vox[,1] > 0 ))/ nrow(las_vox[,1]))
  names(pct_fill_vox) <- paste0(names(pct_fill_vox), "_", vox_res)
  data <- data.frame(c(means, meds, var, sd, cv, IQR, skew, kurt, pct_fill_vox,data_svm))
  
  return(data)
  
}


voxel_summary <- function(las, vox_res = vox_res, max_z ){
  vox <- lidR::voxel_metrics(las, func = vox_mt(Z), res = vox_res, all_voxels = TRUE)
  vox_2 <- vox_mt2(vox, res = vox_res, max_z = max_z)
  out <- suppressWarnings(LAS(vox_2, header = las@header, crs = las@crs, index = las@index ))
  vox_sum(out$SViD, out$FRDi, out$FRSVi, out$PDiBelow, out$PDiAbove, out$SVsum, out$SVMmdAbove,
          out$SVMmd, out$SVM, vox_res[1])
}

vox <- lidR::voxel_metrics(thing, func = vox_mt(Z), res = 0.25, all_voxels = TRUE)
vox_2 = vox_mt2(vox, 0.25, max_z = max(vox$Z))
vox_2$TPDI = ifelse(vox_2$SPI==0,0,1)
out <- suppressWarnings(LAS(vox_2, header = las@header, crs = las@crs, index = las@index ))

vox_rast = pixel_metrics(out, func=~list(spi=sum(SPI), pda=sum(PDiAbove),
                                         tpdi=sum(TPDI)), 0.25)
plot(vox_rast$tpdi)

vals = terra::values(vox_rast$tpdi)
vox_rast$tpdi_norm = (vals-min(vals))/(max(vals)-min(vals))


# Local maxima filter
vox_max = terra::focal(vox_rast, w=3, "sum")
plot(rasterize_canopy(norm_las, res=0.25), type="continuous")

tptrees = locate_trees(vox_max$tpdi, lmf(3))
sptrees = locate_trees(vox_max$spi, lmf(3))
pdtrees = locate_trees(vox_max$pda, lmf(3))
plot(vox_max$tpdi)
plot(tptrees['Z'], col="red", add=T)
ggplot(tptrees)+geom_sf_label(aes(label=treeID))

# 85 is an overtopped tree of interest.
crds = st_coordinates(tptrees[tptrees$treeID==85,])
berf = clip_circle(thing, crds[1,1], crds[1,2], 1)
plot(berf, legend=T)

# merge raster values to point cloud
thing = merge_spatial(norm_las, vox_max$tpdi, attribute="tpdi")
zs_thing = thing$Z
thing$tpdi[is.na(thing$tpdi)] = 0
thing$Z = thing$tpdi
thing = segment_trees(thing, dalponte2016(vox_max$tpdi, tptrees, th_tree=0,
                                          th_seed=0.1, th_cr = 0.01))
thing$Z = zs_thing

# ptrees segmentation
thing = segment_trees(norm_las, ptrees(k=c(85,80,70), hmin=-Inf))
plot(thing, color="treeID")

out = list()
for(i in 1:length(unique(thing$treeID))) {
  zmax = max(na.omit(thing$Z[thing$treeID==i]))
  out[[length(out)+1]]=st_as_sf(thing[thing$treeID==i & thing$Z==zmax,])}
out = do.call(rbind, out)


thing = locate_trees(norm_las, multichm(ws=3, dist_2d=2, dist_3d=1, layer_thickness = 0.02))
plot(vox_max$tpdi)
plot(thing, color="treeID", col="red", add=T)
ggplot(tptrees)+geom_sf_text(aes(label=round(treeID, 0)))


# function to extend line to the top an bottom heights of the point cloud
extend_line = function(p1,p2,height1,height2){
  dvect = c(p2[1]-p1[1], p2[2]-p1[2], p2[3]-p1[3])
  dz1 = abs(p1[3]-height1)
  dz2 = abs(p2[3]-height2)
  sf1 = dz1/dvect[3]
  sf2 = dz2/dvect[3]
  bottom_pt = p1-dvect*sf1
  top_pt = p2+dvect*sf2
  
  return(matrix(c(bottom_pt, top_pt), nrow=2, byrow=T))
}

# 3-d line fitting and ransac


find_stems = function(ttops, las, radius, extend=T){ # radius sort of radius of stem
  
  results = as.data.frame(matrix(nrow=nrow(ttops),ncol=8))
  names(results) = c("treeID","Xb","Yb","Zb","Xt","Yt","Zt","inliers")
  results$treeID = ttops$treeID
  
  
  for(i in 1:nrow(ttops)) {
    center = sf::st_coordinates(ttops[i,])
    center[3] = 0 # Put the point on the ground
    
    tmp_tree = clip_circle(las, center[1],center[2], radius=radius)
    x = tmp_tree$X
    y = tmp_tree$Y
    z = tmp_tree$Z
    
    xyz <- matrix(c(x, y, z), ncol=3)
    xyz_tophalf = matrix(c(x[z>max(z)/2],y[z>max(z)/2],z[z>max(z)/2]),ncol=3)
    xyz_bottomhalf = matrix(c(x[z<=max(z)/2],y[z<=max(z)/2],z[z<=max(z)/2]),ncol=3)
    zmax = max(xyz[,3])
    
    j=0
    k = 500
    in_dist = 0.25
    
    inliers_max = 0
    rando_bmax = c(0,0,0)
    rando_tmax = c(0,0,0)
    while(j<=k) {
      
      rando_top = xyz_tophalf[sample(1:nrow(xyz_tophalf),1),]
      rando_bottom = xyz_bottomhalf[sample(1:nrow(xyz_bottomhalf),1),]
      
      ex_line = extend_line(center, rando_top, 0, zmax)
      
      if(extend==T){
        lineEnds = ex_line
      }else{
        lineEnds = matrix(c(rando_bottom,rando_top), nrow=2, byrow=T)
      }
      
      dists = distancePointToLine(xyz, lineEnds[1,], lineEnds[2,])
      inliers = length(dists[dists<in_dist])
      
      if(inliers > inliers_max) {
        inliers_max = inliers
        rando_bmax = lineEnds[1,]
        rando_tmax = lineEnds[2,]
      }
      
      j=j+1
    }
    results[results$treeID==ttops$treeID[i],2:8] = c(rando_bmax,rando_tmax,inliers_max)
  }
  return(results)
}


stemLines = find_stems(tptrees, thing, radius=0.5, extend=T)
stemLines_p = rbindlist(list(stemLines[,c(1:4)], stemLines[,c(1,5,6,7)]))[order(treeID)]

offs = plot(thing, legend=T, pal=grey.colors(100))
stemLines_p$Xb = stemLines_p$Xb-offs[1]
stemLines_p$Yb = stemLines_p$Yb-offs[2]

segments3d(stemLines_p[,2:4], col="blue", lwd=3)

N <- nrow(xyz) 

mean_xyz <- apply(xyz_tmp, 2, mean)
xyz_pca   <- princomp(xyz_tmp) 
dirVector <- xyz_pca$loadings[, 1]   # PC1

xyz_fit <- matrix(rep(mean_xyz, each = N), ncol=3) + xyz_pca$score[, 1] %*% t(dirVector) 

plot3d(xyz, type="s")
abclines3d(mean_xyz, a = dirVector, col="blue", lwd=2)     # mean + t * direction_vector
for(i in 1:N) segments3d(rbind(xyz[i,], xyz_fit[i,]), col="green3")

thing
