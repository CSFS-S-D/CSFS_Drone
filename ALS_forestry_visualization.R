# Project     ALS_forestry_visualization.R
# Purpose     Visualize forest management with lidar and aerial imagery
# Person      Andy Whelan
# Date        November 5, 2024
# Modified    November 5, 2024
################################################################################


################################################################################
#####                        Setup                                         #####
################################################################################

# R libraries
library(reticulate)
library(terra)
library(lidR)
library(tidyverse)
library(rlas)
library(cloud2trees)
library(sf)
library(rgl)

# install_miniconda(update = T)
# conda_list()
use_condaenv(condaenv = "r-reticulate")

# Install python modules if needed
py_install("opencv-python", pip=T)
py_install("numpy")

# Load python functions from ALS_viz_functions.py
source_python("ALS_viz_functions.py")


#---------------- Functions ----------------------------------------------------

# Function to denoise and color point clouds
pc_clean = function(chunk, ortho) {
  las = readLAS(chunk)
  if(lidR::is.empty(las)) return(NULL)
  
  las = filter_poi(las, Withheld_flag==F)
  las = classify_noise(las, ivf())
  las = filter_poi(las, Classification!=LASNOISE)
  las = classify_ground(las, csf(sloop_smooth=T))
  las = merge_spatial(las, ortho)
} 

# Function to segment trees in point clouds
pc_segment = function(chunk, chm, ttops) {
  las = readLAS(chunk)
  if(lidR::is.empty(las)) return(NULL)
  
  # tree segmentation takes a normalized point cloud
  las = normalize_height(las, knnidw())
  las = segment_trees(las, dalponte2016(chm, ttops))
  
  # save the tree heights in a separate field
  las = add_lasattribute(las, las$Z, "treeHt", "tree heights")
  
  # restore original Z values
  las = unnormalize_height(las)
}

# Function to digitally remove trees and replace them with something ground colored.
pc_thin = function(las_files, ttops, take, qprob=0.75) { # ttops with take/leave column, 
  # take = logical statement for trees to take. e.g., height<10
  # qprob = quantile probability to calculate ground color. 0.5 is the mean,
  # 0.75 is the 3rd quartile which results in a brighter color that usually
  # looks a little better.
  
  if(is.list(las_files) | is.character(las_files)) {
    las = readLAS(las_files)
  }else{
    las = las_files
  }
  if(lidR::is.empty(las)) return(NULL)
  
  # take tree column
  ttops$takeTree = 0 # 0 means not taken, 1 means take
  
  # which trees to take
  ttops$takeTree[eval(parse(text=paste0("ttops$",take)))] = 1
  
  takes = which(las$treeID %in% ttops$treeID[ttops$takeTree==1] & las$Classification != 2)
  R_ground = round(quantile(las$R[las$Classification == 2], prob=qprob))
  G_ground = round(quantile(las$G[las$Classification == 2], prob=qprob))
  B_ground = round(quantile(las$B[las$Classification == 2], prob=qprob))
  
  las$Z[takes] = las$Z[takes] - las$treeHt[takes]
  las$R[takes] = as.integer(R_ground)
  las$G[takes] = as.integer(G_ground)
  las$B[takes] = as.integer(B_ground)
  
  las = normalize_height(las, knnidw())
  las = filter_poi(las, !is.na(treeID) | Z<2)
  las = unnormalize_height(las)
  
  return(las)
} 

# Create a 3d boundary and display it on an las plot
pc_add_polygon = function(poly, dtm, offset) {
  # poly is a spatvector
  # dtm is a spatraster
  bounds = poly %>% 
    as.lines() %>% 
    extractAlong(x=dtm, xy=T)
  
  names(bounds) = c("ID","X","Y","Z")
  
  # add some height to points to make them easier to see
  bounds$Z = bounds$Z+offset
  
  # return new, baddass layer
  return(bounds)
}

# make an orthomosaic that reflects virtual thinning. 
ortho_thin = function(las) {
  # creates an image like an orthomosaic from the input las object.
  # Note, there will likely be holes that can be filled with the python script
  # inpaint. The returned file is a 3 channel 8 bit image that works well with
  # inpaint.
  
  fu = ~list(R = mean(R), G=mean(G), B=mean(B))
  low_ortho = pixel_metrics(las, fu, res=1)
  values(low_ortho) = as.integer(values(low_ortho)/257) # convert to 8 bit
  values(low_ortho)[which(is.na(values(low_ortho)))] = 0
  
  return(low_ortho)
}

# Visualize landscapes
landscape_viz = function(las, bg_col="skyblue"){
  dtm = rasterize_terrain(las, res=1, knnidw())
  open3d()
  bg3d(bg_col)
  surface3d(seq(xmin(dtm),xmax(dtm)-1,1), seq(ymax(dtm)-1,ymin(dtm),-1), 
            matrix(values(dtm)[,1],nrow=ncol(dtm)), color="white", specular="black", lit=F,
            texture="../../../FCFO_CarterLake_lidar_visulazation/low_thin_ortho.png")
  points3d(las$X, las$Y, las$Z, 
           color=rgb(las$R, las$G, las$B, maxColorValue = 65535),
           size=4)
}


################################################################################
###########################  process imagery   #################################
################################################################################

#------------------------- enhance colors --------------------------------------
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



################################################################################
##################### Color point clouds #######################################
################################################################################

# set up the catalog processing
ctg = catalog("../../../FCFO_CarterLake_lidar_visulazation/las/Bundle/")
opt_chunk_size(ctg) = 1500
opt_chunk_buffer(ctg) = 0
opt_chunk_alignment(ctg) = c(0,0)
opt = list(raster_alignment=1)
opt_output_files(ctg) <- "../../../FCFO_CarterLake_lidar_visulazation/las/colorized/colorized_{XLEFT}_{YBOTTOM}"
ctg@processing_options$progress = T
plot(ctg, chunk_pattern=T)


output = catalog_apply(ctg, cleanLAS, ortho=orthos, .options=opt)

# index the colorized las files for faster processing
lasFiles = list.files("../../../FCFO_CarterLake_lidar_visulazation/las/colorized/", full.names = T)
for(f in lasFiles) writelax(f)



################################################################################
#################### run the cloud2trees stuff #################################
################################################################################

# first find a good radius function for the locate-trees search window
las = readLAS("../../../FCFO_CarterLake_lidar_visulazation/las/segmented/segmented_480000_4467000.las")
lass = filter_poi(las, X>480500 & X<480600 & Y>4467500 & Y<4467600)
lass = normalize_height(lass, knnidw())
chm = rasterize_canopy(lass, 1, pitfree())

# after messing around with it, this is pretty good for the Carter Lake area.
# There are a bunch of scrubby junipers around in the foothills that result
# in too many trees with a linear function. This probably could be worked on
# more, but a logrhythmic function seemed like a good spot to start.

ws = function(x) {
y <- dplyr::case_when(is.na(x) ~ 0.001, x < 0 ~ 0.001, x < 2 ~
1, x > 30 ~ 5, TRUE ~ 2*log(x/4)+2.5)
return(y)
}

ttops = locate_trees(lass, lmf(ws))
plot(chm)
plot(sf::st_geometry(ttops), add = TRUE, pch = 3)


output = cloud2trees(output_dir="../../../FCFO_CarterLake_lidar_visulazation/C2T_products/", 
                     input_las_dir="../../../FCFO_CarterLake_lidar_visulazation/las/colorized", 
                     chm_res_m = 1,
                     estimate_tree_dbh=T,
                     ws = ws)



################################################################################
#################### classify trees in the point cloud #########################
################################################################################

# load data from disk
chm = rast("../../../FCFO_CarterLake_lidar_visulazation/C2T_products/point_cloud_processing_delivery/chm_1m.tif")
ttops_list = list.files("../../../FCFO_CarterLake_lidar_visulazation/C2T_products/point_cloud_processing_delivery/", pattern="*tops*", full.names=T)
ttops = lapply(ttops_list, function(x) st_read(x))
ttops = do.call(rbind, ttops)

ttops$treeIDRef = ttops$treeID
ttops$treeID = 1:nrow(ttops)

# set up the catalog
ctg = catalog("../../../FCFO_CarterLake_lidar_visulazation/las/colorized/")
opt_chunk_size(ctg) = 1500
opt_chunk_buffer(ctg) = 10
opt_chunk_alignment(ctg) = c(0,0)
opt = list(raster_alignment=1)
opt_output_files(ctg) <- "../../../FCFO_CarterLake_lidar_visulazation/las/segmented/segmented_{XLEFT}_{YBOTTOM}"
ctg@processing_options$progress = T
plot(ctg, chunk_pattern=T)


options(future.globals.maxSize=2400*1024^2)
output = catalog_apply(ctg, pc_segment, chm=chm, ttops=ttops, .options=opt)



################################################################################
#####               Virtual thinning                                      ######
################################################################################

las = readLAS("../../../FCFO_CarterLake_lidar_visulazation/las/segmented/segmented_480000_4459500.las")
las_low = pc_thin(las, ttops, take = "dbh_cm<30")

# make an orthomosaic-like image the shows the thinning
low_ortho = ortho_thin(las_low)

# inpaint using the "inpaint" python script sourced above
values(low_ortho) = values(rast(inpaint(as.array(low_ortho))))

# plot the original and the thinned version
landscape_viz(las)
landscape_viz(las_low)

# Pay attention to which window is the Focus window
snapshot3d("../../../FCFO_CarterLake_lidar_visulazation/Products/Images/Current_example.png")
snapshot3d("../../../FCFO_CarterLake_lidar_visulazation/Products/Images/lowThin_example.png")


