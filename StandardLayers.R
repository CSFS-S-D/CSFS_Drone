# Program     StandardLayers.R
# Purpose     Create map layers and graphs from drone data
# Person      Andy Whelan 
# Date        October 31, 2024
# Modified    July 15, 2025
################################################################################


############################ Setup #############################################


# Libraries
library(cloud2trees)
library(tidyverse)
library(terra)
library(sf)
library(lidR)
library(future)
library(viridis)
library(tidyterra)
library(rgl)
library(data.table)


################################ Setup #########################################
# set the working directory to the top level of your project.
setwd("Z:/GIS/Drone/DUFO_GrassyMountain07092024/GrassyMountain_DUFO/OUTPUTS/")

# path to folder with .las point clouds.
las_dir = "Terra_lidar/"

# path to shapefile of project area boundaries.
bounds = st_read("../FlightPlan_GrassyMountain_DUFO/GrassyMountain_Final_Units_dis.shp")

# path to orthomosaic (if needed, otherwise NULL).
ortho = NULL

# DTM and CHM constants for cloud to trees.
dtm_res = 0.25
chm_res = 0.25

# create output file structure
dir.create("Products/")
dir.create("Products/Cloud2Trees_products/")
dir.create("Products/Raster/")
dir.create("Products/Vector/")
dir.create("Products/las/")
dir.create("Products/las/segmented")
dir.create("Products/Raster/temp/")
dir.create("Products/Graphs/")
dir.create("Products/Animations/")

############################### End setup ######################################


################## Tree and crown processing ###################################
# retile if necessary (drone data will likely need to be retiled)
dir.create(paste0(las_dir,"retiled/"))
ctg = catalog(las_dir)
opt_chunk_buffer(ctg) = 0
# Think carefully about the chunk size
opt_chunk_size(ctg) = 300
opt_chunk_alignment(ctg) = c(0,0)
opt_output_files(ctg) = "Z:/GIS/Drone/DUFO_GrassyMountain07092024/GrassyMountain_DUFO/OUTPUTS/Terra_lidar/retiled/retile_{XLEFT}_{YBOTTOM}"
plot(ctg, chunk=T)

catalog_retile(ctg)

# search window size function to detect trees
ws2 = function(x) { y <- dplyr::case_when(is.na(x) ~ 0.001, x < 0 ~ 0.001, x < 
                                      8 ~ 2.5, x > 26.5 ~ 5, TRUE ~ exp(-(3.5 * (1/x)) + (x^0.17)))
return(y)
}

# cloud2Trees
output = cloud2trees(output_dir="Products/Cloud2Trees_products/", 
                     paste0(las_dir,"retiled/"), 
                     dtm_res_m=dtm_res, 
                     chm_res_m=chm_res,
                     estimate_tree_dbh = T,
                     ws = ws2,
                     estimate_dbh_from_cloud=F)
################################################################################


######################## Make map layers #######################################

# load data from disk
chm = rast(paste0("Products/Cloud2Trees_products/point_cloud_processing_delivery/chm_", 
                  as.character(chm_res),"m.tif"))
dtm = rast(paste0("Products/Cloud2Trees_products/point_cloud_processing_delivery/dtm_", 
                  as.character(dtm_res),"m.tif"))
ttops = st_read("Products/Cloud2Trees_products/point_cloud_processing_delivery/final_detected_tree_tops.gpkg")
crowns = st_read("Products/Cloud2Trees_products/point_cloud_processing_delivery/final_detected_crowns.gpkg")


# rename a couple of things
names(chm) = "chm"
names(dtm) = "dtm"

# Fill some holes in the dtm if needed
# dtm = focal(dtm, w=9, fun=mean, na.policy="only", na.rm=T) # fill holes

# ---> crop the rasters to the project boundary with a little buffer
bounds = bounds %>% st_transform(crs = st_crs(chm))
bounds_dis = bounds %>% st_union() %>% st_sf()
bounds_buf = bounds_dis %>% st_buffer(dist= 20)

# orthomosaic
if(!is.null(ortho)){
  ortho = crop(ortho, bounds_buf, mask=T)
  writeRaster(ortho, "Products/Raster/ortho_cropped.tif")
}

# tree_tops
ttops_c = st_intersection(ttops, bounds_buf)
ttops_c$treeID_ref = ttops_c$treeID
ttops_c$treeID = 1:nrow(ttops_c)

# crowns 
crowns = st_intersection(crowns, bounds_buf)

# Make a sort of crown continuity map
clumps <- crowns %>% st_union() %>% st_sf() %>% st_cast("POLYGON")
clumps$area = as.numeric(st_area(clumps))

# slope and aspect 
slope = terrain(dtm, "slope")
aspect = terrain(dtm, "aspect")



# --- Calculate crown Density ----

# ---> functions 
# density
crDensity_f = function(z) {
  vd = length(z[z>2])/length(z)
}

# function to segment point clouds and rasters by tree
crown_cat = function(chunk, chm, ttops, res) {
  las = readLAS(chunk)
  if(lidR::is.empty(las)) return(NULL)
  if(is.null(terra::intersect(ext(las), ext(chm)))) return(NULL)
  xmin=round(ext(las)[1])
  ymin=round(ext(las)[3])
  
  las = classify_noise(las, ivf(res=1, n=10))
  las = filter_poi(las, Classification!=LASNOISE)
  las = classify_ground(las, csf(sloop_smooth=TRUE, rigidness = 2))
  las = normalize_height(las, knnidw(k=10, p=2)) # Need a normalized las for vertical density
  las = filter_poi(las, Z>=0)
  las = segment_trees(las, algorithm = dalponte2016(chm, ttops, max_cr=100))
  fun1 = ~list(treeID = mean(treeID))
  crown_rast = pixel_metrics(las, fun1, res)
  verticalDens1 = pixel_metrics(las, ~crDensity_f(z=Z), res=res)
  
  rStack = c(crown_rast, verticalDens1)
  names(rStack) =c("cr", "vd")
  writeRaster(rStack, paste0("Z:\\GIS\\Drone\\StateForest\\SnowCourses_lidar_102524\\Products\\Raster\\temp\\crown_rasts_",
                             xmin,"_",ymin,".tif"))
  
  las = unnormalize_height(las)
  return(las)
}



# ---- set up the catalog processing ----

# retile to get rid of overlap
ctg = readLAScatalog("Products/Terra_output/")
opt_chunk_size(ctg) = 1000
opt_chunk_buffer(ctg) = 0
opt_filter(ctg) = ""
opt_output_files(ctg) <- "Products/Terra_output/las/retiled/retiled_{XLEFT}_{YBOTTOM}"
ctg@processing_options$progress = T

catalog_retile(ctg)


# Segment point clouds and create density and crown area rasters
ctg = readLAScatalog("Terra_lidar/retiled/")
opt_chunk_size(ctg) = 300
opt_chunk_buffer(ctg) = 10
opt_filter(ctg) = ""
opt_output_files(ctg) <- "Products/las/segmented/segmented_{XLEFT}_{YBOTTOM}"
opt_chunk_alignment(ctg) = c(0,0)
opt = list(raster_alignment=1)
plot(ctg, chunk_pattern=T)

output = catalog_apply(ctg, crown_cat, chm, ttops_c, res=0.25, .options=opt)


# Crop crown area and crown density rasters
rast_list = list.files("Products/Raster/temp/", full.names = T)
rasts = do.call(merge, lapply(rast_list, function(x) rast(x)))
rasts = crop(x=rasts, y=vect(bounds_buf), mask=T)
crown_rast = rasts$cr
density = rasts$vd


# put everything in a raster stack
cropped = lapply(list(chm,dtm,slope,aspect,crown_rast,density), function(x) crop(x,y=bounds_buf, mask=T))
rasters = do.call(c, cropped)

# make a digital surface model
rasters$dsm = classify(rasters$chm, matrix(c(NA,0), ncol=2))+rasters$dtm
varnames(rasters$dsm)

# Make sure everything has names.
names(rasters)[6] = "verticalDensity"

###################################### Save products ###########################################################
# Rasters
for(i in 1:length(names(rasters))) {
  writeRaster(rasters[[i]], paste0("Products/Raster/",names(rasters)[i],".tif"), overwrite=T)
}

# Vector
st_write(ttops_c, "Products/Vector/ttops_c.gpkg", append=F)
st_write(crowns, "Products/Vector/crowns.gpkg", append=F)
st_write(clumps, "Products/Vector/canopy_cont.gpkg", append=F)
st_write(bounds, "Products/Vector/bounds.gpkg", append=F)
st_write(bounds_dis, "Products/Vector/bounds_dis.gpkg", append=F)
st_write(bounds_buf, "Products/Vector/bounds_buf.gpkg", append=F)



###################################### Graphs and other figures ################################################

# dtm
png("Products/Graphs/wholeDTM.png", height=2500, width=4000)
plot(rasters$dtm, pax=list(cex.axis=5), plg=list(cex=5), mar=c(3.1,4.1,2.1,9.1))
dev.off()

png("Products/Graphs/wholeDSM.png", height=2500, width=4000)
plot(rasters$dsm, pax=list(cex.axis=5), plg=list(cex=5), mar=c(3.1,4.1,2.1,9.1))
dev.off()

png("Products/Graphs/wholeCHM.png", height=2500, width=4000)
plot(rasters$chm, pax=list(cex.axis=5), plg=list(cex=5), mar=c(3.1,4.1,2.1,9.1))
dev.off()

png("Products/Graphs/wholeVD.png", height=2500, width=4000)
plot(rasters$vd, pax=list(cex.axis=5), plg=list(cex=5), mar=c(3.1,4.1,2.1,9.1), col=map.pal("byg", n=100))
dev.off()

# Multipanel
png("Products/Graphs/heightProducts.png", height=2500, width=4000)
par(mfrow=c(2,2))
plot(rasters$dsm, pax=list(cex.axis=5),
     plg=list(cex=3), mar=c(5.1,5.1,4.1,6.1))
plot(rasters$dtm, pax=list(cex.axis=5),
     plg=list(cex=3), mar=c(5.1,5.1,4.1,6.1), col=inferno(100))
plot(rasters$slope, pax=list(cex.axis=5), plg=list(cex=3, x="topleft"), mar=c(5.1,5.1,4.1,6.1), col=map.pal("gyr", n=11), 
     breaks=c(0,1.72,3.43,5.71,8.53,11.3,14.04,16.7,21.8,30.96,45,90))
plot(rasters$aspect, pax=list(cex.axis=5),plg=list(cex=3, x="topleft"), mar=c(5.1,5.1,4.1,6.1), 
     col=c("gray","red","orange","yellow","green","skyblue","blue","purple4","magenta","orangered"),
     breaks=c(-1,0,22.5,67.5,112.5,157.5,202.5,247.5,292.5,337.5,360))
dev.off()


#################################### tree map graph ###################################################

png("Products/Graphs/treeMap.png", height=1500, width=4000)
plot(rasters$chm, pax=list(cex.axis=5), plg=list(cex=5), mar=c(3.1,4.1,2.1,9.1))
plot(ttops_c['tree_height_m'], add=T, col="black", pch=20, cex=2)
dev.off()

png("Products/Graphs/treeMap_zoom.png", height=3000, width=4000)
par(mfrow=c(2,1))
plot(rasters$chm, pax=list(cex.axis=5), plg=list(cex=5), mar=c(3.1,4.1,2.1,9.1))
plot(ttops['tree_height_m'], add=T, col="black", pch=20, cex=2)
plot(rasters$chm, pax=list(cex.axis=5), plg=list(cex=5), mar=c(3.1,4.1,2.1,9.1), ext=ext(c(416000,416375, 4494233,4494400)))
plot(ttops['tree_height_m'], add=T, col="black", pch=20, cex=3, ext=ext(c(416000,416375, 4494233,4494400)))
dev.off()


#---> crown density map 
ggplot(crowns) + geom_sf(aes(fill=crown_area_m2, color=crown_area_m2)) + 
  coord_sf(datum=st_crs(32613)) +
  scale_fill_continuous(type="viridis") +
  scale_color_continuous(type="viridis") +
  labs(color=expression(paste("Area (",m^2,")")), fill=expression(paste("Area (",m^2,")"))) +
  theme_bw() +
  theme(axis.text=element_text(size=6),
        axis.text.y=element_text(angle=90),
        legend.title=element_text(size=6),
        legend.text=element_text(size=4),
        legend.position = "inside",
        legend.position.inside=c(0.045,0.24),
        legend.key.size = unit(0.5, "cm"),
        panel.grid.major = element_blank())

ggsave("Products/Graphs/crownArea.png", height=5, width=8, dpi=600)


# clump area map 
ggplot(clumps) + geom_sf(aes(fill=area, color=area)) + 
  coord_sf(datum=st_crs(32613)) +
  scale_fill_continuous(type="viridis") +
  scale_color_continuous(type="viridis") +
  labs(color=expression(paste("Area (",m^2,")")), fill=expression(paste("Area (",m^2,")"))) +
  theme_bw() +
  theme(axis.text=element_text(size=6),
        axis.text.y=element_text(angle=90),
        legend.title=element_text(size=6),
        legend.text=element_text(size=4),
        legend.position = "inside",
        legend.position.inside=c(0.05,0.24),
        legend.key.size = unit(0.5, "cm"),
        panel.grid.major = element_blank())

ggsave("Products/Graphs/clumpArea.png", height=5, width=8, dpi=600)


################################### End ########################################