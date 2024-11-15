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


# run the python shadow remover
orthos = list()
for(img in imgs) {
  tmp_rast = rast(img) %>% terra::split(f=c(1,1,1,2))
  values(tmp_rast[[1]])=values(rast(remove_shadow(img)))
  orthos[length(orthos)+1] = tmp_rast[[1]]
}
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


# lasR version to see if it's faster
read = reader_las()
filter1 = delete_points(filter="-drop_withheld")
write1 = write_las(paste0("../../../FCFO_CarterLake_lidar_visulazation/las/cleaned","/*_cleaned.laz"))
pipeline = read + filter1 + classify_with_sor() + delete_noise() + write1
exec(pipeline, on=ctg, with=list(progress=TRUE))

clean_dir = "../../../FCFO_CarterLake_lidar_visulazation/las/cleaned/"
clean_files = dir(clean_dir)
for(i in 1:length(clean_files)) {
  tmp_las = readLAS(paste0(clean_dir,clean_files[i]))
  tmp_las = merge_spatial(tmp_las, orthos)
  writeLAS(tmp_las, paste0("../../../FCFO_CarterLake_lidar_visulazation/las/colorized/",sub("cleaned","colorized",clean_files[i])))
}
