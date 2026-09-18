# Program   snagDetection.R
# Purpose   To detect snags in aerial lidar point clouds
# Person    Andy Whelan
# Date      September 11, 2025
# Modified  August 4, 2026
##################################################################################


################################################################################
#######################    Setup    ############################################
################################################################################
# Libraries
library(lidR)
library(lidRviewer)
library(terra)
library(sf)
library(terrainr)
library(multimode)
library(autothresholdr)


# Environment
setwd("C:/Users/csfsuser/Colostate/Science & Data - GIS Team Server/PROJECTS/Drone/Manitou/")


lasDirs = c("N1_October06_2025/Terra_output/N1_200AGL_20MPH_TFOFF/",
            "N1_October06_2025/Terra_output/N1_300AGL_20MPH_TFOFF/",
            "N1_October06_2025/Terra_output/N1_400AGL_20MPH_TFOFF/",
            "N1_October06_2025/N1-200AGL-20MPH-TF/Terra_output/",
            "N1_October06_2025/N1-200AGL-10MPH-TF/Terra_output/",
            "USGS_lidar/",
            "../Frisco_GRFO/TerraOutput/East0.las")
            
boundsPath = "bounds.shp"
boundsPath = "../Frisco_GRFO/Shapefiles/FriscoBackyard_Units/FriscoBackyard_Units/FriscoBackyard_Units.shp"



################################################################################
##########################    Data    ##########################################
################################################################################
# Wing (2015) suggest to only use first returns, but to normalize intensity, we
# need all returns, so we keep them in for now and remove them after normalization.
# The "-drop_withheld" flag basically drops a bunch
# of points that have been flagged as bad. The drone lidar data probably won't 
# have any withheld, but the USGS lidar data can have a lot.

data = list()
for(i in 1:length(lasDirs)){
  
  lasHeaderEx = readLASheader("../Frisco_GRFO/TerraOutput/East0.las")
  lasHeaderEx = readLASheader(dir(lasDirs[i], full.names=T)[1])
  ctg = catalog(lasDirs[i], filter="-drop_withheld -keep_first")
  bounds = st_read(boundsPath)
  # bounds = subset(bounds, Phase_1==1)
  # bounds = st_union(bounds)
  bounds = st_transform(bounds, st_crs(lasHeaderEx))
  
  # If the simple feature collection has a Z dimension, remove it.
  bounds = st_zm(bounds)
  
  # subset if necessary
  # bounds = bounds[1,]
  
  # clip the catalog to the bounds. This also loads the actual point cloud.
  pc = clip_roi(ctg, st_bbox(bounds))
  pc = clip_roi(pc, bounds)
  # pc$Z = pc$Z*0.3048 # convert feet to meters. XY are already in meters.
  ################################################################################
  ########                      Data Wranglin                             ########
  ################################################################################
  
  # The first thing we need to do is get the Density right. Look at the density 
  # raster and homogenize to lower densities until there aren't obvious (i.e., stripey) 
  # high-density areas. UPDATE: While it turns out you can get rid of the density
  # bands by setting the density value to around 200, it causes the snag segmentation
  # algorithm to take forever. Forever as in maybe days or weeks. I waited half a 
  # day before I threw in the towel. Anyway, I decided to reduce the density to something
  # closer to the example data from the snag detection help. The algorithm was
  # developed using lower density point clouds, I think.
  
  # Thin the point cloud to between 10 and 20 points per square meter. If you see
  # differences it point density in the plot due to flight line overlap (think stripes),
  # thin it further.
  
  pc = decimate_points(pc, homogenize(density = 10, res=5))
  # pc = decimate_points(pc, random(density = 10))
  
  # Normalize height for snag detection
  # pc = classify_ground(pc, csf())
  pc = normalize_height(pc, knnidw())
  
  # Individual tree segmentation
  chm = rasterize_canopy(pc)
  
  # I used a window size function from the cloud2trees package that increases the
  # window size exponentially with height. Adjust as necessary to get the best
  # tree identification.
  ttops = locate_trees(chm, lmf(cloud2trees::itd_ws_functions()$exp_fn))
  
  # segment the point cloud
  pc = segment_trees(pc, dalponte2016(chm, ttops))
  
  # Scale if values are really low
  if(max(pc$Intensity)<10000) pc$Intensity = as.integer(pc$Intensity*((2^16-1)/max(pc$Intensity)))
  
  pc$Intensity = as.integer(pc$Intensity/(2^16-1)*(2^8-1))
  
  data[[length(data)+1]] = pc
}


modes = lapply(data, function(x) {
  
  # We need to isolate woody vegetation, so we'll remove the ground points.
  over_pc = filter_poi(x, Z>=2)

  # finding some modes (high and low points in the histogram) can help us find
  # lower and upper intensity thresholds.
  return(locmodes(over_pc$Intensity, mod0=3, lowsup=0, uppsup=max(x$Intensity))$locations) 
})

Lint_est = function(pc, filt) {
  dif_tb = matrix(nrow=24, ncol=512)
  for(i in 2:25) {
    tmp_dens = density(filter_poi(pc, Z>filt)$Intensity)
    dif_tb[i-1,] = tmp_dens$y
  }
  dim(dif_tb)
  sds = apply(dif_tb, 2, sd)
  xs = density(data[[1]]$Intensity)$x
  plot(sds ~ xs)
  xs[sds == min(sds[xs>25 & xs<55])]
  
  # points that are higher than their neighbors
  maxs = xs[which(c(FALSE, diff(sign(diff(sds))) == -2, FALSE))]
  # points that are lower than their neighbors
  mins = xs[which(diff(sign(diff(sds))) == 2) + 1]
  
  return(list(maxs, mins, max(filter_poi(pc, Z>2)$Intensity)))
  
}

snagMapper = function(pc, Uint, ladj=0) {
  
  # Get the point density requirement (PDR) from the Plot-level point density (PLPD). 
  # We use the entire point cloud for PLPD. If PLPD <= 3, PDR = 3; if 3 < PLPD <= 6, 
  # PDR = 4; if 6 < PLPD <= 12, PDR = 5; if PLPD > 12, PDR = 8.
  if(density(pc)<=3) {
    pdr = 3
  }else if(density(pc)>3 & density(pc)<=6){
    pdr = 4
  }else if(density(pc)>6 & density(pc)<=12){
    pdr = 5
  }else{
    pdr = 8
  }
  
  over_pc = filter_poi(pc, Z>2)
  
  # Set lower and upper intensity thresholds. For Front Range ponderosa pine, a 
  # lower threshold of 35 or a little more seems to work pretty well. I use the modes
  # to help adjust it. I add the value for the first mode to 35. Sometimes that's a
  # little high. We can further adjust it later. For Uint, I just use the last mode -
  # it's close to the max intensity, which works pretty well.
  
  while(!is.na(ladj)) {
      
    Lint = ladj
    
    # from Wing etal. 2015. I don't see much reason to fiddle with these. They're the same for both point clouds.
    bbpr_thresholds <- matrix(c(0.80, 0.80, 0.70,
                                0.85, 0.85, 0.60,
                                0.80, 0.80, 0.60,
                                0.90, 0.90, 0.55),
                              nrow =3, ncol = 4)
    
    
    # Run the stem segmentation algorithm using the appropriate thresholds (Either Lint and Uint
    # or values just higher than the histogram peaks. Might take some trial and error).
    
    over_pc = segment_snags(over_pc, wing2015(BBPRthrsh_mat = bbpr_thresholds, pt_den_req = pdr,
                                    low_int_thrsh = Lint, uppr_int_thrsh = Uint))
    
    # Plot it
    # plot(over_pc, color="snagCls", legend=T)
    # plot(pc, color="snagCls", legend=T)
    # 
    # ladj = readline(prompt=paste0("If this looks crappy, input a new ladj value. Otherwise, just hit return. Current value is: ", ladj, ". "))
    # ladj = as.numeric(ladj)
    ladj = NA
  }
  
  # This gets the point data as a data frame.
  df = over_pc@data
  
  # combine over_pc with the part of the original point cloud <=2
  thing = filter_poi(pc, Z<=2)
  thing = add_lasattribute(thing, 0, "snagCls", "Number identifying a snag class")
  pc = rbind(thing, over_pc)
  
  # This is a custom function we'll use in the aggregate command below. All it does
  # is find the ratio of points classified as snag to all points. The aggregate command
  # will apply it to each tree for us.
  snaggy = function(x){
    length(x[x!=0])/length(x)
  }
  
  # Use "aggregate" to find snag to not snag point ratio for each tree.
  t = aggregate(snagCls~treeID, df, snaggy)
  
  thing = dplyr::left_join(pc@data, t, by="treeID")
  thing$snagCls[is.na(thing$snagCls)] = 0
  pc = add_lasattribute(pc, thing$snagCls, "snagPct", "PercentSnagPoints")

  # tree height
  h = aggregate(Z~treeID, df, max)
  thing = dplyr::left_join(pc@data, h, by="treeID")
  thing$Z.y[is.na(thing$Z.y)] = 0
  pc = add_lasattribute(pc, thing$Z.y, "treeHt", "TreeHeight")
  
  # dbh 
  dbh = pc$treeHt*0.3707240+1.1857200
  pc = add_lasattribute(pc, dbh, "dbh", "TreeDBH")
  
  
  print(paste0("Lint = ", Lint))  
  return(list(pc, Lint))
} 

LintEstOld = function(pc) {  
  over_pc = filter_poi(pc, Z>2)
  bbvf = (length(over_pc$Intensity[over_pc$Intensity<=50 | over_pc$Intensity>=170]))/ # 50 and 170
    (length(over_pc$Intensity[over_pc$Intensity>50 | over_pc$Intensity<170]))
  
  # Calculate these (lower and upper intensity thresholds)
  Lint = 20*bbvf+0.075*(max(over_pc$Intensity)) + 26.5 #26.5
  Uint = 20*bbvf+0.1875*(max(over_pc$Intensity)) + 100.25

    # Some wisdom from Wing etal. 2015 below
  # If Lint<50, Lint=50, If Lint>70, Lint=70
  # If Uint<150, Uint=150 else if Uint>170, Uint=170
  # If Located in SF Area, Lint=Lint+5 else if SF, Uint=Uint-5. SF is "Storrie Fire" 
  # and not relevant to us but they adjusted Lint and Uint for it, so we might read Wing 2015 more closely to see how we might
  # come up with adjustments of our own.
  
  if(Lint<50) Lint = 50 # 50
  if(Lint>70) Lint = 70 # 70
  if(Uint<150) Uint = 150
  if(Uint>170) Uint = 170

  return(c(Lint, Uint))
}
  
snagRast = function(snagLAS, res=1) {
  tmpsc = snagLAS$snagCls
  tmpsc[tmpsc>2] = 0
  snagLAS$snagCls = tmpsc
  rast22off = pixel_metrics(snagLAS, sum(snagCls), res=res)
  rast22off = classify(rast22off, matrix(c(NA, 0), byrow=T, ncol=2))
  return(rast22off)
}
  
modes = matrix(unlist(modes), ncol=5, byrow=T)

zfilt = 2
# Lints200 = Lint_est(data[[1]], zfilt)
# Lints300 = Lint_est(data[[2]], zfilt)
# Lints400 = Lint_est(data[[3]], zfilt)
# Lints220tf = Lint_est(data[[4]], zfilt)
# Lints210tf = Lint_est(data[[5]], zfilt)
# LintsUS = Lint_est(data[[6]], zfilt)

talgs = data.frame(alg = c("IJDefault", "Huang", "Huang2", "Intermodes", "IsoData", "Li", "Mean", "MinErrorI", "Minimum", "Moments", "Otsu", "Percentile", "RenyiEntropy", "Shanbhag", "Triangle"))
talgs$thresh32off = sapply(talgs$alg, function(x) auto_thresh(filter_poi(data[[1]], Z>zfilt)$Intensity, x)[1])
talgs$thresh32off = sapply(talgs$alg, function(x) auto_thresh(filter_poi(data[[2]], Z>zfilt)$Intensity, x)[1])
talgs$thresh22off = sapply(talgs$alg, function(x) auto_thresh(filter_poi(data[[1]], Z>zfilt)$Intensity, x)[1])
talgs$thresh42off = sapply(talgs$alg, function(x) auto_thresh(filter_poi(data[[3]], Z>zfilt)$Intensity, x)[1])
talgs$thresh22on = sapply(talgs$alg, function(x) auto_thresh(filter_poi(data[[4]], Z>zfilt)$Intensity, x)[1])
talgs$thresh21on = sapply(talgs$alg, function(x) auto_thresh(filter_poi(data[[5]], Z>zfilt)$Intensity, x)[1])
talgs$threshUSGS = sapply(talgs$alg, function(x) auto_thresh(filter_poi(data[[6]], Z>zfilt)$Intensity, x)[1])
talgs$threshFrisco = sapply(talgs$alg, function(x) auto_thresh(filter_poi(pc, Z>zfilt)$Intensity, x)[1])


# I changed snagMapper so that it returns the whole over_pc point cloud, not the one thinned with ptr.
snags200 = snagMapper(data[[1]], modes[1,5], ladj=talgs[12,2])
snags300 = snagMapper(data[[2]], modes[2,5], ladj=talgs[12,3])
snags400 = snagMapper(data[[3]], modes[3,5], ladj=talgs[12,4])
snags220tf = snagMapper(data[[4]], modes[4,5], ladj=talgs[12,5])
snags210tf = snagMapper(data[[5]], modes[5,5], ladj=talgs[12,6])
snagsUS = snagMapper(data[[6]], modes[6,5], ladj=talgs[12,7])
snagsFrisco = snagMapper(pc, 86, ladj=36)

plot(snagsFrisco[[1]], color="snagCls", legend=T)
plot(snagsFrisco[[1]], color="RGB")

plttr = function(df, height){
  tmp_snag = df[df$snagCls != 0,]
  tmp_chm = rasterize_canopy(tmp_snag)
  tmp_ttops = locate_trees(tmp_chm, lmf(5))
  plot(tmp_ttops[tmp_ttops$Z > height,]['Z'], pch=20)
}
  
plttr(snags200[[1]], 15)
plttr(snags300[[1]], 15)
plttr(snags400[[1]], 15)
plttr(snagsUS[[1]], 15)

# snag threshold and count table
tholds = matrix(nrow=15,ncol=4)
snagCounts = matrix(nrow=15,ncol=4)
for(i in 1:15){
  for(j in c(1,2,3,6)){
    tmp_snag = snagMapper(data[[j]], modes[j,5], ladj=talgs[i,j+1], ptr=0.33)
    if(j==6) j=4
    tholds[i,j] = tmp_snag[[3]]
    snagCounts[i,j] = tmp_snag[[2]]
  }
}

old200 = LintEstOld(data[[1]])
old300 = LintEstOld(data[[2]])
old400 = LintEstOld(data[[3]])
old220tf = LintEstOld(data[[4]])
old210tf = LintEstOld(data[[5]])
oldUS = LintEstOld(data[[6]])

oldSnags200 = snagMapper(data[[1]], Lint=50.58553, Uint=150)
oldSnags300 = snagMapper(data[[2]], Lint=50.84088, Uint=150)
oldSnags400 = snagMapper(data[[3]], Lint=50.65318, Uint=150)
oldSnags220tf = snagMapper(data[[4]], old220tf[2], ladj=old220tf[1], ptr=0.33)
oldSnags210tf = snagMapper(data[[5]], old210tf[2], ladj=old210tf[1], ptr=0.33)
oldSnagsUS = snagMapper(data[[6]], Lint=50, Uint=150)

ThresholdingComp = data.frame(Altitude = c("200","300","400","200 (20mph TF)",
                                           "200 10mph TF","10,000"),
                              Auto=c(snags200[[2]],
                                     snags300[[2]],
                                     snags400[[2]],
                                     snags220tf[[2]],
                                     snags210tf[[2]],
                                     snagsUS[[2]]),
                              "Wing Et al."=c(oldSnags200[[2]],
                                              oldSnags300[[2]],
                                              oldSnags400[[2]],
                                              oldSnags220tf[[2]],
                                              oldSnags210tf[[2]],
                                              oldSnagsUS[[2]])
                              )

write.csv(ThresholdingComp, "../Github/CSFS_DataHub/snagDetection/Tables/thresholdingComp.csv")

autoSD = sd(ThresholdingComp$Auto)
wingSD = sd(ThresholdingComp$Wing.Et.al.)



t2 = snagRast(snags200[[1]], res=1)
t3 = snagRast(snags300[[1]], res=1)
tj = classify(t2, matrix(c(0,15,0),byrow=T, ncol=3))
tk = classify(t3, matrix(c(0,15,0),byrow=T, ncol=3))
plot(tj)
plot(tk)

plot(thing)
rast200 = snagRast(snags200[[1]], res=1)
rast300 = snagRast(snags300[[1]], res=1)
plot(rast300)
rast400 = snagRast(snags400[[1]], res=1)
rast220tf = snagRast(snags220tf[[1]], res=2)
rast210tf = snagRast(snags210tf[[1]], res=2)
rastUS = snagRast(snagsUS[[1]], res=1)

oldRast200 = snagRast(oldSnags200[[1]], res=2)
oldRast300 = snagRast(oldSnags300[[1]], res=2)
oldRast400 = snagRast(oldSnags400[[1]], res=2)
oldRast220tf = snagRast(oldSnags220tf[[1]], res=2)
oldRast210tf = snagRast(oldSnags210tf[[1]], res=2)
oldRastUS = snagRast(oldSnagsUS[[1]], res=2)

par(mfrow=c(2,2))
  plot(rast400)
  plot(oldRast400)
  plot(rastUS)
  plot(oldRastUS)
dev.off()

# Make some sweet plots and table data ################################################
# Return intensity density histograms
plotDir = "C:/Users/csfsuser/Colostate/Science & Data - GIS Team Server/PROJECTS/Drone/Github/CSFS_DataHub/snagDetection/Graphs/"
png(paste0(plotDir, "ReturnIntensity.png"))
  par(mfrow=c(2,2))
  plot(density(data[[1]]$Intensity), main="UAV 200ft", xlab="Return Intensity")
  plot(density(data[[2]]$Intensity), main="UAV 300ft", xlab="Return Intensity")
  plot(density(data[[3]]$Intensity), main="UAV 400ft", xlab="Return Intensity")
  plot(density(data[[6]]$Intensity), main="USGS 10,000ft", xlab="Return Intensity")
dev.off()

# Lidar scene comparison
bounds3d = extract(dtm, vect(st_cast(bounds, "LINESTRING")), xy=T) |> 
  mutate(Z = Z+30) |> 
  st_as_sf(coords=c("x","y","Z"), crs=crs(dtm))

for(i in c(1,2,3)){
  catalog(lasDirs[i]) |> readLAS(filter="-keep_random_fraction 0.01") |> plot(color="RGB") |> 
    add_flightlines3d(bounds3d, radius=2)
}

ctg = catalog("USGS_lidar/")
las10000 = readLAS(ctg, filter="-keep_random_fraction 0.1 -drop_withheld")
aoi = st_bbox(ctg, crs=crs(ctg)) |> st_as_sfc()
naip = rsi::get_naip_imagery(aoi, start_date="2020-01-01", end_date = "2025-01-01", 
                             output_filename = "naip.tif")
naip = rast(naip)
las10000 = merge_spatial(las10000, naip[[1:3]])

plot(las10000, color="RGB") |> add_flightlines3d(bounds3d, radius=2)
plot(las300, color="RGB") |> add_flightlines3d(bounds3d, radius=2)

rgl::snapshot3d("../Github/CSFS_DataHub/snagDetection/Graphs/las10000.png")
rgl::snapshot3d("../Github/CSFS_DataHub/snagDetection/Graphs/las200.png")
rgl::snapshot3d("../Github/CSFS_DataHub/snagDetection/Graphs/las300.png")
rgl::snapshot3d("../Github/CSFS_DataHub/snagDetection/Graphs/las400.png")
  
png("../Github/CSFS_DataHub/snagDetection/Graphs/footprintComparison.png", width=2875, height=3000)
  par(mfrow=c(2,2))
  for(i in c("las200.png","las300.png","las400.png","las10000.png")) {
    tmpRast = rast(paste0("../Github/CSFS_DataHub/snagDetection/Graphs/",i), noflip=T)
    plotRGB(tmpRast)
  }
dev.off()

# Automatic thresholding thresholds
tableDir = "C:/Users/csfsuser/Colostate/Science & Data - GIS Team Server/PROJECTS/Drone/Github/CSFS_DataHub/snagDetection/Tables/"
table_talgs = talgs[,c(1:4,7)]
names(table_talgs) = c("Algorithm","200ft","300ft","400ft","10,000ft")
write.csv(table_talgs, paste0(tableDir, "Thresholds.csv"))

# snag class vs. RGB point cloud comparison
plot(pc, color="RGB")
plot(snags200[[1]], color="snagCls", legend=T)
rgl::snapshot3d("../Github/CSFS_DataHub/snagDetection/Graphs/snagClassPC.png")
rgl::snapshot3d("../Github/CSFS_DataHub/snagDetection/Graphs/RGB_pc.png")


plot(snags200[[1]], color="snags")
plot(snags300[[1]], color="snags")
plot(snags400[[1]], color="snags")
plot(snags220tf[[1]], color="snags")
plot(snags210tf[[1]], color="snags")
plot(snagsUS[[1]], color="snags")

plot(rast200)
plot(rast300)
plot(rast400)
plot(rast220tf)
plot(rast210tf)
plot(rastUS)

snagPlotter = function(rstr, las, res){
  tmp_chm = rasterize_canopy(las, res=res)
  tmp_rast = project(rstr, tmp_chm)
  extnt = ext(las)
  
  # rastrcl = classify(rstr, matrix(c(0,13.4,0, 13.4,100,1), byrow=T, ncol=3))
  rastrcl = classify(rstr, matrix(c(0,NA), byrow=T, ncol=2))
  poly = vect(st_cast(st_buffer(st_union(st_as_sf(as.polygons(rastrcl))), 0.5), "POLYGON"))
  poly$Z = extract(tmp_chm, poly, mean)[,2]
  poly$snagCls = extract(rstr, poly, mean)[,2]
  poly$area = expanse(poly)
  plot(poly, "Z", xlim=extnt[1:2], ylim=extnt[3:4])
}

snagPlotter(rast200, data[[1]], 2)
snagPlotter(rast300, data[[2]], 2)
snagPlotter(rast400, data[[3]], 2)
snagPlotter(rastUS, data[[6]], 2)

png("../Github/CSFS_DataHub/snagDetection/Graphs/SnagPolys.png", height=480, width=590)
  plot(poly, "snagCls", pax = list(cex.axis=1.5), plg=list(cex=1.5))
dev.off()

plot(poly, "area")
snagpts = centroids(poly)
plot(snagpts, "snagCls")

p = plot(pc, color="RGB")
add_treetops3d(p, st_as_sf(snagpts), z = "Z")







# Estimate intensity from RGB for SFM point clouds
# pc$Intensity = as.integer(pc$R*0.2126+pc$G*0.7152+pc$B*0.0722)

# Intensity values need to be 8-bit. They are often 16-bit. So, we
# need to convert.

# We need to figure out what the high and low threshold intensity values are for
# branch to bowl ratios, which has a strong effect on snag detection.
# over_pc = filter_poi(pc, Classification!=2)
# over_pc = filter_first(over_pc)
# over_pc = decimate_points(over_pc, homogenize(density = 10, res=5))


# Histograms help us find the lower and upper thresholds. Look for two peaks and
# choose thresholds just higher than each. Adjust the histogram as necessary.
hist(over_pc$Intensity, breaks=20)
axis_ticks = seq(floor(min(over_pc$Intensity)), ceiling(max(over_pc$Intensity)), by=10)
axis(side=1, at=axis_ticks, labels=axis_ticks)
plot(density(over_pc$Intensity))


###################################################################################
# The following is an attempt to automatically set the upper and lower intensity  #
# thresholds. Sometimes it doesn't work. In that case, just set the thresholds    #
# to values just higher than the two highest peaks in the histogram.              #
###################################################################################


# Automatically find modes. Sometimes there's crap, uLim tries to get rid of it.
uLim = pc$Intensity[order(pc$Intensity)][round(0.999*length(pc$Intensity))]
modes = locmodes(filter_poi(over_pc, Z>max(over_pc$Z)*0.8 & Z<=max(over_pc$Z))$Intensity, mod0=3, lowsup=0, uppsup=max(pc$Intensity)) 
modes = locmodes(filter_poi(over_pc, Z>24 & Z<=26)$Intensity, mod0=3, lowsup=0, uppsup=max(pc$Intensity)) 
modes = locmodes(over_pc$Intensity, mod0=3, lowsup=0, uppsup=max(pc$Intensity)) 

thing = filter_poi(over_pc, Z>2 & Z<=4)
thing2 = filter_poi(over_pc, Z>max(over_pc$Z)-5 & Z<=max(over_pc$Z))
strf = summary(c(thing$Intensity, thing2$Intensity))

berf = density(thing$Intensity, from=strf[1], to=strf[6])
berf2 = density(thing2$Intensity, from=strf[1], to=strf[6])
plot(berf)
lines(berf2, col="red")

floof = abs((berf$y[1:256]-berf2$y[1:256]))
greeble = order(floof)[1:10]
Lint = berf$x[min(greeble)]
berf2$x[120]

Lint = berf2$x[which.max(berf2$y)]

# modes$locations = modes$locations/(2^16-1)*(2^8-1)

# Calculate bole and branch to foliage intensity ratio. Use the two values just 
# after the peaks from the histogram.
over_pc = filter_poi(data[[3]], Z>2)
bbvf = (length(over_pc$Intensity[over_pc$Intensity<=50 | over_pc$Intensity>=170]))/ # 50 and 170
  (length(over_pc$Intensity[over_pc$Intensity>50 | over_pc$Intensity<170]))

# Calculate these (lower and upper intensity thresholds)
Lint = 20*bbvf+0.075*(max(over_pc$Intensity)) + 26.5 #26.5
Uint = 20*bbvf+0.1875*(max(over_pc$Intensity)) + 100.25
Lint = 35 + modes$locations[1]
# Some wisdom from Wing etal. 2015 below
# If Lint<50, Lint=50, If Lint>70, Lint=70
# If Uint<150, Uint=150 else if Uint>170, Uint=170
# If Located in SF Area, Lint=Lint+5 else if SF, Uint=Uint-5. SF is "Storrie Fire" 
# and not relevant to us but they adjusted Lint and Uint for it, so we might read Wing 2015 more closely to see how we might
# come up with adjustments of our own.

if(Lint<50) Lint = 50 # 50
if(Lint>70) Lint = 70 # 70
if(Uint<150) Uint = 150
if(Uint>170) Uint = 170

# from Wing etal. 2015. I don't see much reason to fiddle with these. They're the same for both point clouds.
bbpr_thresholds <- matrix(c(0.80, 0.80, 0.70,
                            0.85, 0.85, 0.60,
                            0.80, 0.80, 0.60,
                            0.90, 0.90, 0.55),
                          nrow =3, ncol = 4)

# # Lower density to get more snags. Used this for USGS point cloud at Manitou.
# bbpr_thresholds <- matrix(c(0.70, 0.70, 0.60,
#                             0.75, 0.75, 0.50,
#                             0.70, 0.70, 0.50,
#                             0.80, 0.80, 0.45),
#                           nrow =3, ncol = 4)

# Get the point density requirement (PDR) from the Plot-level point density (PLPD). 
# We use the entire point cloud for PLPD. If PLPD <= 3, PDR = 3; if 3 < PLPD <= 6, 
# PDR = 4; if 6 < PLPD <= 12, PDR = 5; if PLPD > 12, PDR = 8.
if(density(pc)<=3) {
  pdr = 3
}else if(density(pc)>3 & density(pc)<=6){
  pdr = 4
}else if(density(pc)>6 & density(pc)<=12){
  pdr = 5
}else{
  pdr = 8
}

# Run the stem segmentation algorithm using the appropriate thresholds (Either Lint and Uint
# or values just higher than the histogram peaks. Might take some trial and error).

pc = segment_snags(over_pc, wing2015(BBPRthrsh_mat = bbpr_thresholds, pt_den_req = pdr,
                                            low_int_thrsh = Lint, uppr_int_thrsh = Uint))

# Plot it
# plot(over_pc, color="snagCls", legend=T)
plot(pc, color="snagCls", legend=T)


# We now have a bunch of points that are classified as snag points, but what we
# really need to know is which trees are snags. So, let's calculate the ratio of 
# snag points in each tree so we can classify trees with a high ratio as snags.

# This gets the point data as a data frame.
df = over_pc@data
df = pc@data

# Take a look at the snag class data.
df$snagCls
table(df$snagCls)

# This is a custom function we'll use in the aggregate command below. All it does
# is find the ratio of points classified as snag to all points. The aggregate command
# will apply it to each tree for us.
snaggy = function(x){
  length(x[x!=0])/length(x)
}

# Use "aggregate" to find snag to not snag point ratio for each tree.
t = aggregate(snagCls~treeID, df, snaggy)


# Choose a threshold percentage of snag points to classify snags (50%, 75%?). 
# This code gives us the treeID numbers of snags.
sngs = t$treeID[t$snagCls>0.33]

# Here we subset the point cloud so we can just look at the snags.
over_snags = over_pc[over_pc$treeID %in% sngs,]
snags = pc[pc$treeID %in% sngs,]

plot(over_snags)
plot(snags)

# This looks alright but there's still some junk/noise we might want to get rid of.
# We're often more interested in taller snags - they represent a lot more biomass,
# they're more dangerous to people, and also better for wildlife. Each identified
# tree has a tree height (Z) associated with it. These data are in the over_ttops data
# frame. We could subset the over_ttops data frame to trees taller than 5-10 meters or
# trees larger than 5 cm dbh and then use the treeIDs from that subset to select 
# only those trees from the over_snags point cloud. 




# dbh 
ttops$dbh = ttops$Z*0.3707240+1.1857200
# over_ttops = st_intersection(over_ttops, bounds);
ttops_dbh5 = ttops[ttops$dbh >= 5,] # leave out really short ones
over_snags_dbh5 = over_snags[over_snags$treeID %in% ttops_dbh5$treeID,]
snags_dbh5 = snags[snags$treeID %in% ttops_dbh5$treeID,]
length(unique(over_snags_dbh5$treeID))
length(unique(snags_dbh5$treeID))
plot(snags_dbh5, legend=T)
view(snags_dbh5)

# There were some weird floating snags, and this gets rid of them. Don't run if
# you don't see a lot of floating snags.
snags_dbh5_3 = aggregate(Z~treeID, over_snags_dbh5@data, min)
snags_dbh5_3 = snags_dbh5_3[snags_dbh5_3$Z < 3,]

over_snags_dbh5_3 = over_snags_dbh5[over_snags_dbh5$treeID %in% snags_dbh5_3$treeID,]
plot(over_snags_dbh5_3, color="treeID")

nrow(snags_dbh5_3)


# make a nice plot
naip = spanner::download_naip_for_las(pc)
naip = rast("naip_imagery.tif")
RGB(naip) = 1:3

sample(naip, matrix(c(NA,0), ncol=2, byrow=T))
pc = merge_spatial(pc, naip[[1:3,]])

p_pc = pc
p_pc$R[p_pc$treeID %in% unique(over_snags_dbh5$treeID)] = as.integer(65535)
p_pc$B[p_pc$treeID %in% unique(over_snags_dbh5$treeID)] = as.integer(0)
p_pc$G[p_pc$treeID %in% unique(over_snags_dbh5$treeID)] = as.integer(0)
plot(p_pc, color="RGB")
plot(pc, color="RGB")
view(pc)

rgl::snapshot3d("exampleGrey.tif")
rgl::snapshot3d("exampleRed.tif")


# Validation
# Ground truthed snags
gt_snags = st_read("Snags_N1_Shapefile")
gt_snags = st_transform(gt_snags, st_crs(lasHeaderEx))
snag_ttops = ttops[ttops$treeID %in% over_snags_dbh5$treeID,]
gt_snags$type = sub("[a-z].*", "", gt_snags$Name)
gt_snags$Z = extract(chm, gt_snags)$Z
gt_snags$Z[is.na(gt_snags$Z)] = 0

p1 = plot(pc, color="snagCls", legend=T)
add_treetops3d(p1, gt_snags)

p2 = plot(over_pc, color="RGB")
add_treetops3d(p2, thing)

plot(ext(bounds))
plot(gt_snags['type'], pch=20, add=T)
plot(snag_ttops['treeID'], pch=20, col="green", add=T)


################################################################################
#######################           Junk              ############################
################################################################################
# filter by tree height
# snags taller than 10m
over_ttops10 = over_ttops[over_ttops$Z > 10,]
over_snags10 = over_snags[over_snags$treeID %in% over_ttops10$treeID,]
length(unique(over_snags10$treeID))
view(over_snags10)


snags_10_3 = aggregate(Z~treeID, over_snags10@data, min)
snags_10_3 = snags_10_3[snags_10_3$Z < 3,]

over_snags10_3 = over_snags10[over_snags10$treeID %in% snags_10_3$treeID,]
plot(over_snags10_3, color="treeID")

nrow(snags_10_3)


# Snags taller than 5m
over_ttops5 = over_ttops[over_ttops$Z > 5,]
over_snags5 = over_snags[over_snags$treeID %in% over_ttops5$treeID,]
length(unique(over_snags5$treeID))
view(over_snags5)

snags_5_3 = aggregate(Z~treeID, over_snags5@data, min)
snags_5_3 = snags_5_3[snags_5_3$Z < 3,]

over_snags5_3 = over_snags5[over_snags5$treeID %in% snags_5_3$treeID,]
plot(over_snags5_3, color="treeID")

nrow(snags_5_3)



################ Normalize lidar return intensity ##############################
# DJI Terra flightlines don't have usable PointsourceID data (they're all 0), so 
# we need to generate it. To do that we run retrieve_flightlines. The dt parameter is 
# important. Try to set it to about how much time it takes the drone to turn around
# at the end of a flight line. Seems like we should skip for UAV lidar.
#################################################################################

# if(all(pc$PointSourceID == 0)) {
#   pc = retrieve_flightlines(pc, dt=10)
#   pc$PointSourceID = pc$flightlineID
# }
# 
# # The number of flightline ID's should match your flight pattern.
# unique(pc$PointSourceID)
# 
# # Now we need some locations for the lidar sensor during the flight. Pmin is the
# # minimum number of points needed to record a sensor location. More points make a
# # more accurate location, but with sparse point clouds, we'll need a lower number.
# sensor = track_sensor(pc, Roussel2020(interval=0.5, pmin=20))
# plot(sensor['PointSourceID'])
# rnge = get_range(pc, sensor)
# summary(rnge)
# 
# # Now we can normalize intensity. Rs is the range of reference. The internet told
# # me to set it to the mean sensor height. For the drone in this case, that's 60m.
# pc1 = normalize_intensity(pc, range_correction(sensor, Rs=mean(rnge), f=2.3))
# summary(pc1$Intensity)
# hist(pc1$Intensity, breaks=50)
# 
# hist(rnge, breaks=50)
# 
# # Now filter out everything but first returns.
# pc = filter_first(pc1)
# 
# 

