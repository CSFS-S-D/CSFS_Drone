# Project     ALS_forestry_visualization.R
# Purpose     Visualize forest management with lidar and aerial imagery
# Person      Andy Whelan
# Date        November 5, 2024
# Modified    March 14, 2025
################################################################################


################################################################################
#####                        Setup                                         #####
################################################################################

# R libraries

## install lasR from the r-universe
install.packages("lasR", repos = "https://r-lidar.r-universe.dev")

## install TreeLS from github
remotes::install_github(repo = "tiagodc/TreeLS", upgrade = F)
# install cloud2trees
remotes::install_github(repo = "georgewoolsey/cloud2trees", upgrade = F)
library(reticulate)
library(terra)
library(lidR)
library(tidyverse)
library(rlas)
library(cloud2trees)
library(sf)
library(rgl)

### Get external data for cloud2trees.
get_data(force=T)

### Install Python. Only need to do this once
# install_miniconda(update = T)
# conda_list()
# use_condaenv(condaenv = "r-reticulate")

### Install python modules if needed
# py_install("opencv-python", pip=T)
# py_install("numpy")


# Set working directory
setwd("C:/Users/C820392777/Colostate/Science & Data - GIS Team Server/PROJECTS/Aerial_lidar_visualizations/FCFO_CarterLake_lidar_visulazation/")

# 

# Load python functions from ALS_viz_functions.py
source_python("../../Drone/Github/CSFS_DataHub/ALS_viz_functions.py")
source("../../Drone/Github/CSFS_DataHub/forestVizFunctions.R")



################################################################################
######              Setup                                                 ######
################################################################################
baseDir = "./"
lidarDir = "las/Bundle/"
naipDir = "NAIP/"
projProj = "26913" # epsg code


# Project boundaries. Unit 6 this spring. 1-5 this fall.
bounds = st_read(paste0(baseDir, "Boundaries/"))
# subunits = st_read(paste0(baseDir, "Boundaries/CP_Subunits/"))

# reproject
bounds = st_transform(bounds, crs="EPSG:26913")
# subunits = st_transform(subunits, crs="EPSG:26913")

st_write(st_bbox(bounds) %>% st_as_sfc(), paste0(baseDir, "Boundaries/bbox.shp"))




#reproject the whole dang catalog (I did it on the command line)
# setwd(paste0(baseDir, "Lidar/LiDAR_2025-03-12T14_45_10.353Z"))
# las_list=dir()
# for(i in 1:length(las_list)) {
#   system2(command="las2las64", args=c(" -i ", las_list[i],
#                 " -o ",gsub(".laz","_utm.laz",las_list[i]),
#                 " -proj_epsg ", projProj))
# }
# setwd("../../../../Drone/Github/CSFS_DataHub/")

#### fix the z-values





################################################################################
###########################  process imagery   #################################
################################################################################

# ------------------find orthomosaics------------------------------------------#
imgs = list.files(paste0(baseDir, naipDir),
                  recursive=T,
                  pattern="*.tif",
                  full.names = T)

#------------------------- enhance colors -------------------------------------#


# run the Python color enhancement and reproject to match las catalog if necessary
orthos = list()
for(img in imgs) {
  tmp_rast = rast(img) %>% terra::split(f=c(1,1,1,2)) %>% pluck(1)
  values(tmp_rast) = values(rast(enhance_pic(img, brightness=5, contrast=1.35)))
  if(terra::crs(tmp_rast, proj=T)!=projection(ctg)) project(tmp_rast, projection(ctg))
  orthos[length(orthos)+1] = tmp_rast
}

# # run the python shadow remover
# orthos = list()
# for(img in imgs) {
#   tmp_rast = rast(img) %>% terra::split(f=c(1,1,1,2))
#   values(tmp_rast[[1]])=values(rast(remove_shadow(img)))
#   orthos[length(orthos)+1] = tmp_rast[[1]]
# }

orthos = do.call(merge, orthos)
plotRGB(orthos)


################################################################################
################ Color and segment point clouds ################################
################################################################################

# Load lidar data
las_list=dir(paste0(baseDir, lidarDir), full.names=T)
las_list = las_list[grep("^(?!.*utm)", las_list, perl=T)]
ctg = catalog(las_list)
# writeLAS(readLAS(ctg), "Combined.las")

# function to color and segment point clouds
pc_cleanAndSeg = function(chunk, ortho, chm, ttops){
  pc_clean(chunk=chunk, ortho=ortho) %>% pc_segment(chm=chm, ttops=ttops)
}

# set up the catalog processing
opt_chunk_size(ctg) = 1500
opt_chunk_buffer(ctg) = 0
opt_chunk_alignment(ctg) = c(0,0)
opt = list(raster_alignment=1)
dir.create("las/colorized/")
opt_output_files(ctg) <- "las/colorized/colorized_{XLEFT}_{YBOTTOM}"
ctg@processing_options$progress = T
plot(ctg, chunk_pattern=T)

# Process. pc_clean also colors the point clouds
output = catalog_apply(ctg, pc_clean, ortho=orthos, fix_Z=T, .options=opt)

# index the colorized las files for faster processing
lasFiles = list.files(paste0(baseDir, "/las/colorized/"), full.names = T)
for(f in lasFiles) writelax(f)



################################################################################
#################### classify trees in the point cloud #########################
################################################################################

# load data from cloud2trees from disk
chm = rast(paste0(baseDir, "C2T_products/point_cloud_processing_delivery/chm_1m.tif"))
ttops_list = list.files(paste0(baseDir, "C2T_products/point_cloud_processing_delivery/"), 
                               pattern="*tops*", full.names=T)
ttops = lapply(ttops_list, function(x) st_read(x))
ttops = do.call(rbind, ttops)

ttops$treeIDRef = ttops$treeID
ttops$treeID = 1:nrow(ttops)

st_crs(ttops)==st_crs(chm)
st_crs(ctg) = st_crs(ttops)

# set up the catalog
ctg = catalog(paste0(baseDir, "las/colorized/"))
opt_chunk_size(ctg) = 1500
opt_chunk_buffer(ctg) = 10
opt_chunk_alignment(ctg) = c(0,0)
opt = list(raster_alignment=1)
opt_output_files(ctg) <- paste0(baseDir, "las/segmented/segmented_{XLEFT}_{YBOTTOM}")
ctg@processing_options$progress = T
plot(ctg, chunk_pattern=T)


options(future.globals.maxSize=3000*1024^2)
output = catalog_apply(ctg, pc_segment, chm=chm, ttops=ttops, .options=opt)

# Save ttops with the new treeIDs
dir.create(paste0(baseDir, "Products/Vector"), recursive = T)
st_write(ttops, paste0(baseDir, "Products/Vector/ttops.gpkg"))
st_write(paste0(baseDir, "Products/Vector/crowns.gpkg"))


################################################################################
#####               Virtual thinning                                      ######
################################################################################

take_trees = ttops_thinned[ttops_thinned$is_keep_tree==0,]
las_low = pc_thin(las, take_trees=take_trees, qprob=0.75)
# las_low = pc_thin(las, tree_locs, keep = "(!tree_locs$treeID %in% ttops_thinned$treeID)")

# make an orthomosaic-like image that shows the thinning
low_ortho = ortho_thin(las_low)

# inpaint using the "inpaint" python script sourced above
values(low_ortho) = values(rast(inpaint(as.array(low_ortho))))
writeRaster(low_ortho, "../../../FCFO_CarterLake_lidar_visulazation/Products/Textures/low_ortho_example.png", overwrite=T)

# plot the original and the thinned version
landscape_viz(las, texture_file = "../../../FCFO_CarterLake_lidar_visulazation/Products/Textures/low_ortho_example.png")
landscape_viz(las_low, texture_file = "../../../FCFO_CarterLake_lidar_visulazation/Products/Textures/low_ortho_example.png")

# Pay attention to which window is the Focus window
snapshot3d("../../../FCFO_CarterLake_lidar_visulazation/Products/Images/Current_example.png")
snapshot3d("../../../FCFO_CarterLake_lidar_visulazation/Products/Images/lowThin_example.png")


