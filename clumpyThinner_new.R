# Program   gappyClumpyThin.R
# Purpose   Virtually remove trees from point clouds to simulate gappy-clumpy
#             thinning ala GTR 373 (Addington etal. 2018)
# Person    Andy Whelan
# Date      February 27, 2025
# Modified  February 27, 2025
################################################################################

#---> Packages. Install 'em if you need 'em.
# install.packages('lidRviewer', repos = c('https://r-lidar.r-universe.dev'))
# install.packages("pak")
# pak::pkg_install("mhahsler/dbscan", upgrade = T)

#---> Libraries
library(lidR)
library(tidyverse) # the tidyverse
library(reticulate)

# spatial analysis
library(terra) # raster
library(sf) # simple features
library(dbscan)

# visualization
library(rgl)

# custom thinning and vizualization functions/
source("gappyClumpyFunctions.R")
source("forestVizFunctions.R")

# custom python scripts to enhance orthomosaics images.
source_python("ALS_viz_functions.py")

################################################################################
######              Setup                                                 ######
################################################################################

#----------------> Load data 
las = readLAS("../../../FCFO_CarterLake_lidar_visulazation/las/segmented/segmented_478500_4459500.las")
tree_locs = st_read("../../../FCFO_CarterLake_lidar_visulazation/Products/Vector/ttops.gpkg")
tree_locs = st_crop(tree_locs, ext(las))


#---------------> Classify forests into low, medium, and high productivity ala
# GTR 373 (sort of).
# 0: inoperable or not forested
# 1: Low density matrix (south aspect, low productivity)
# 2: Medium density matrix (east and west slopes)
# 3: High density matrix (north aspect)

# Classify by above ground point density then reclass to forested and non-forested.  
forRast = forested(las) %>% 
  classify(matrix(c(0,0.20,0,0.20,1,1), ncol=3, byrow=T))


# Get slope and aspect then reclassify to simplify.
dtm = rasterize_terrain(las, res=5)
slope = terrain(dtm, "slope") %>% classify(matrix(c(0,10,1,10,35,2,35,90,NA),
                                                  ncol=3, byrow=T))
aspect = terrain(dtm, "aspect") %>% classify(matrix(c(0,45,0,45,135,90,135,225,180,
                                                      225,315,270,315,360,0),
                                                    ncol=3, byrow=T))


# Classify by forest productivity
forProd = (slope+aspect)*forRast %>% focal(fun="median", w=5)

# One more reclass for the final categories
forProd_rcl = classify(forProd, matrix(c(1,3,
                                         2,3,
                                         91,2,
                                         92,2,
                                         181,1,
                                         182,1,
                                         271,2,
                                         272,2), ncol=2, byrow=T))


# convert to polygon and split into individual polygons then remove some small 
# islands of a different class within larger areas of another class.
forProd_poly = forProd_rcl %>% as.polygons %>% st_as_sf() %>% st_cast("POLYGON") %>% 
  nngeo::st_remove_holes()

# Calculate area
forProd_poly$area = st_area(forProd_poly)



################################################################################
######            Do the virtual thinning                                 ######
################################################################################

# Split by managemement type
prodLow = forProd_poly[forProd_poly$slope==1,] # low
prodMed = forProd_poly[forProd_poly$slope==2,] # medium
prodHigh = forProd_poly[forProd_poly$slope==3,] # high


# Define a dbh limit for overstory trees and a distance in meters to
# separate tree groups.
tree_clump_dist_m <- 6

# Define and overstory cutoff (can also be by height: ostory_ht_m=4)
ostory_dbh_cm <- 1*2.54 # = inches*2.54. 4 because lots of small dbh trees


#-----------> Select target clump proportions for each productivity

# desired BA
target_ba_low = 15 # cannot be > current BA # Low density
target_ba_med = 20 # cannot be > current BA # Medium density
target_ba_high = 25 # cannot be > current BA # High density

# desired QMD
target_qmd_low = 8 # low
target_qmd_med = 9 # medium
target_qmd_high = 10 # high

# check tpa
get_tpa(target_ba_high, target_qmd_high)
# desired proportion (%) of trees in each clump size
# !cannot be create larger proportion of ">25 trees" clump as this would require adding trees...
# c("Individual", "2-4 trees",    "5-9 trees",    "10-15 trees","16-25 trees",">25 trees")
target_pcts_low = c(45, 40, 10, 5, 0, 0) # low density forest
target_pcts_med = c(30, 50, 10, 5, 5, 0) # medium density forest
target_pcts_high = c(10, 15, 50, 10, 10, 5) # high density forest


################################################################################
#####                     processing                                      ######
################################################################################

lowProd = forestThin_gg(prodLow, ttops=tree_locs, tree_clump_dist_m, ostory_dbh_cm,
                      target_ba_low, target_qmd_low, target_pcts_low)
medProd = forestThin_gg(prodMed, ttops=tree_locs, tree_clump_dist_m, ostory_dbh_cm,
                      target_ba_med, target_qmd_med, target_pcts_med)
highProd = forestThin_gg(prodHigh, ttops=tree_locs, tree_clump_dist_m, ostory_dbh_cm,
                      target_ba_high, target_qmd_high, target_pcts_high)

ttops_thinned = do.call(rbind, list(lowProd,medProd,highProd))


################################################################################
#####               Visualization                                         ######
################################################################################

take_trees = ttops_thinned[ttops_thinned$is_keep_tree==0,]
las_thinned = pc_thin(las, take_trees=take_trees, qprob=0.75)

# make an orthomosaic-like image that shows the thinning
thinned_ortho = ortho_thin(las_thinned)

# inpaint using the "inpaint" python script from the python functions script above.
values(thinned_ortho) = values(rast(inpaint(as.array(thinned_ortho))))
writeRaster(thinned_ortho, "../../../FCFO_CarterLake_lidar_visulazation/Products/Textures/thinned_ortho_example.png", overwrite=T)

# plot the original and the thinned version
landscape_viz(las, texture_file = "../../../FCFO_CarterLake_lidar_visulazation/Products/Textures/thinned_ortho_example.png")
landscape_viz(las_thinned, texture_file = "../../../FCFO_CarterLake_lidar_visulazation/Products/Textures/thinned_ortho_example.png")

# Pay attention to which window is the Focus window
snapshot3d("../../../FCFO_CarterLake_lidar_visulazation/Products/Images/Current_example.png")
snapshot3d("../../../FCFO_CarterLake_lidar_visulazation/Products/Images/lowThin_example.png")


