# Project     ALS_forestry_visualization.R
# Purpose     Visualize forest management with lidar and aerial imagery
# Person      Andy Whelan
# Date        November 5, 2024
# Modified    November 5, 2024
################################################################################

# Install python modules if needed
# py_install("opencv-python")
# py_install("numpy")


# R libraries
library(reticulate)
library(terra)
library(lidR)
library(tidyverse)
library(rlas)
library(lasR)

# Load python functions from ALS_viz_functions.py
source_python("ALS_viz_functions.py")


########################## Remove shadows ######################################
# find orthomosaics
imgs = list.files("../../../FCFO_CarterLake_lidar_visulazation/naip/", 
                  recursive=T,
                  pattern="*20130716.tif",
                  full.names = T)


# run the Python color enhancement
orthos = list()
for(img in imgs) {
  tmp_rast = rast(img) %>% terra::split(f=c(1,1,1,2))
  values(tmp_rast[[1]]) = values(rast(enhance_pic(img, brightness=5, contrast=1.35)))
  orthos[length(orthos)+1] = tmp_rast[[1]]
}

# # run the python shadow remover
# orthos = list()
# for(img in imgs) {
#   tmp_rast = rast(img) %>% terra::split(f=c(1,1,1,2))
#   values(tmp_rast[[1]])=values(rast(remove_shadow(img)))
#   orthos[length(orthos)+1] = tmp_rast[[1]]
# }

orthos = do.call(merge, orthos)


##################### Add RGB to point clouds ##################################
# set up the catalog processing
ctg = catalog("../../../FCFO_CarterLake_lidar_visulazation/las/Bundle/")
opt_chunk_size(ctg) = 1500
opt_chunk_buffer(ctg) = 0
opt_chunk_alignment(ctg) = c(0,0)
opt = list(raster_alignment=1)
opt_output_files(ctg) <- "../../../FCFO_CarterLake_lidar_visulazation/las/colorized/colorized_{XLEFT}_{YBOTTOM}"
ctg@processing_options$progress = T
plot(ctg, chunk_pattern=T)

cleanLAS = function(chunk, ortho) {
  las = readLAS(chunk)
  if(lidR::is.empty(las)) return(NULL)
  
  las = filter_poi(las, Withheld_flag==F)
  las = classify_noise(las, ivf())
  las = filter_poi(las, Classification!=LASNOISE)
  las = classify_ground(las, csf(sloop_smooth=T))
  las = merge_spatial(las, ortho)
} 

output = catalog_apply(ctg, cleanLAS, ortho=orthos, .options=opt)

# index the colorized las files for faster processing
lasFiles = list.files("../../../FCFO_CarterLake_lidar_visulazation/las/colorized/", full.names = T)
for(f in lasFiles) writelax(f)



#---> run the cloud2trees stuff 
output = cloud2trees(output_dir="C2T_products/", input_las_dir="las/colorized", chm_res_m = 1)


#---> Now classify trees in the point cloud
# load data from disk
chm = rast("C2T_products/point_cloud_processing_delivery/chm_1m.tif")
ttops_list = list.files("C2T_products/point_cloud_processing_delivery/", pattern="*tops*", full.names=T)
ttops = lapply(ttops_list, function(x) st_read(x))
ttops = do.call(rbind, ttops)

ttops$treeIDRef = ttops$treeID
ttops$treeID = 1:nrow(ttops)

# set up the catalog
ctg = catalog("las/colorized/")
opt_chunk_size(ctg) = 1500
opt_chunk_buffer(ctg) = 10
opt_chunk_alignment(ctg) = c(0,0)
opt = list(raster_alignment=1)
opt_output_files(ctg) <- "las/segmented/segmented_{XLEFT}_{YBOTTOM}"
ctg@processing_options$progress = T
plot(ctg, chunk_pattern=T)


pc_segment = function(chunk, chm, ttops) {
  las = readLAS(chunk)
  if(lidR::is.empty(las)) return(NULL)
  
  # tree segmentation takes a normalized point cloud
  las = normalize_height(las, knnidw())
  las = segment_trees(las, dalponte2016(chm, ttops))
  
  # save the tree heights in a separate column
  las = add_lasattribute(las, las$Z, "treeHt", "tree heights")
  
  # restore original Z values
  las = unnormalize_height(las)
}

options(future.globals.maxSize=2400*1024^2)
output = catalog_apply(ctg, pc_segment, chm=chm, ttops=ttops, .options=opt)


#---> thinning
# take tree column
ttops$takeTree = 0 # 0 means not taken, 1 means take

# Low thin
ttops$takeTree[ttops$tree_height_m < 15] = 1 # low thin Rx

# high thin
ttops$takeTree[ttops$tree_height_m > 15] = 1 # low thin Rx

# uniform thin 70% take
ttops$takeTree[ttops$treeID %in% sample(ttops$treeID, length(ttops$treeID)*0.70)] = 1


pc_thin = function(las_files, ttops, qprob=0.75) { # ttops with take/leave column, 
  # qprob = quantile probability to calculate ground color. 0.5 = mean
  las = readLAS(las_files)
  is(lidR::is.empty(las)) return(NULL)
  
  takes = which(las$treeID %in% ttops$treeID[ttops$takeTree==1])
  R_ground = round(quantile(las$R[las$Classification == 2], prob=qprob))
  G_ground = round(quantile(las$G[las$Classification == 2], prob=qprob))
  B_ground = round(quantile(las$B[las$Classification == 2], prob=qprob))
  
  las$Z[takes] = las$Z[takes] - las$treeHt[takes]
  las$R[takes] = as.integer(R_ground)
  las$G[takes] = as.integer(G_ground)
  las$B[takes] = as.integer(B_ground)
  
  return(las)
} 

