# Program   gappyClumpyThin.R
# Purpose   Virtually remove trees from point clouds to simulate gappy-clumpy
#             thinning ala GTR 373 (Addington etal. 2018)
# Person    Andy Whelan
# Date      February 27, 2025
# Modified  March 14, 2025
################################################################################

#---> Packages. Install 'em if you need 'em.
# install.packages('lidRviewer', repos = c('https://r-lidar.r-universe.dev'))
# install.packages("pak")
# pak::pkg_install("mhahsler/dbscan", upgrade = T)
# remotes::install_github("rspatial/predicts")

#---> Libraries
library(lidR)
library(tidyverse) # the tidyverse
library(reticulate)

# spatial analysis
library(terra) # raster
library(sf) # simple features
library(dbscan)
library(predicts)

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

baseDir = "../../../Aerial_lidar_visualizations/FCFO_Cherokee_Park/"

#----------------> Load data 
bounds = st_read(paste0(baseDir, "Boundaries/CPBoundary"))

# Management actions
bounds$desc_ # 1-6 will be cut

#---> lidar
ctg = catalog(paste0(baseDir, "las/segmented/"))

# reproject bounds to ctg
bounds = st_transform(bounds, projection(ctg))
# bounding box for the area
bb = st_bbox(bounds)

las = clip_roi(ctg, bb)
tree_locs = st_read(paste0(baseDir, "Products/Vector/ttops.gpkg"))
tree_locs = st_transform(tree_locs, projection(ctg))
tree_locs = st_crop(tree_locs, bb)


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

# Split by managemement type. more than about 30000 trees messes stuff up (as.integer has a limit)
prod1 = poly_split(bounds[bounds$desc_==1,], tree_locs) # 
prod2 = poly_split(bounds[bounds$desc_==2,], tree_locs) # 
prod3 = poly_split(bounds[bounds$desc_==3,], tree_locs) # 
prod4 = poly_split(bounds[bounds$desc_==4,], tree_locs) # 
prod5 = poly_split(bounds[bounds$desc_==5,], tree_locs) # 
prod6 = poly_split(bounds[bounds$desc_==6,], tree_locs) # 

# find basal area
trees1 = st_intersection(tree_locs, bounds[bounds$desc_==1,])
trees2 = st_intersection(tree_locs, bounds[bounds$desc_==2,])
trees3 = st_intersection(tree_locs, bounds[bounds$desc_==3,])
trees4 = st_intersection(tree_locs, bounds[bounds$desc_==4,])
trees5 = st_intersection(tree_locs, bounds[bounds$desc_==5,])
trees6 = st_intersection(tree_locs, bounds[bounds$desc_==6,])

# more than about 30000 trees messes stuff up (as.integer has a limit)
nrow(trees1)
nrow(trees2) 
nrow(trees3) 
nrow(trees4) 
nrow(trees5) 
nrow(trees6) 

# ba's
sum(trees1$basal_area_ft2)/trees1$Acres_1[1]
sum(trees2$basal_area_ft2)/trees2$Acres_1[1]
sum(trees3$basal_area_ft2)/trees3$Acres_1[1]
sum(trees4$basal_area_ft2)/trees4$Acres_1[1]
sum(trees5$basal_area_ft2)/trees5$Acres_1[1]
sum(trees6$basal_area_ft2)/trees6$Acres_1[1]


# Per Pete Morin, target BA is 50-70. We'll set it at 60
# ba's after thinning. Use these for target_ba's
# sum(treesLow$basal_area_ft2[treesLow$dbh_cm/2.54<12])/treesLow$acres[1] # 41.8


# Define a dbh limit for overstory trees and a distance in meters to
# separate tree groups.
# ****For this one, we need to use a distance in feet! crs is in feet.****
tree_clump_dist_m <- 6

# Define and overstory cutoff (can also be by height: ostory_ht_m=4)
ostory_dbh_cm <- 2*2.54 # = inches*2.54


#-----------> Select target clump proportions for each productivity

# desired BA
target_ba = 60 # cannot be > current BA # Low density
# target_ba_med = 24.9 # cannot be > current BA # Medium density
# target_ba_high = 28.6 # cannot be > current BA # High density

# desired QMD (from Pete)
target_qmd = 12 #
# target_qmd_med = 10 # medium
# target_qmd_high = 10 # high

# check tpa
get_tpa(target_ba, target_qmd)
# desired proportion (%) of trees in each clump size
# !cannot be create larger proportion of ">25 trees" clump as this would require adding trees...
# c("Individual", "2-4 trees",    "5-9 trees",    "10-15 trees","16-25 trees",">25 trees")
target_pcts = c(10, 20, 40, 20, 10, 0) # low density forest
target_pcts_med = c(30, 50, 10, 5, 5, 0) # medium density forest
target_pcts_high = c(30, 50, 10, 5, 5, 0) # high density forest


################################################################################
#####                     processing                                      ######
################################################################################

# Unit 1 
unit1 = list()
for(i in 1:nrow(prod1)){
  thin_tmp = forestThin_gg(prod1[i,], ttops=trees1, tree_clump_dist_m, ostory_dbh_cm,
                target_ba, target_qmd, target_pcts)
  unit1[[length(unit1)+1]] = thin_tmp
}
unit1 = do.call(rbind, unit1)

# Unit 2
unit2 = list()
for(i in 1:nrow(prod2)){
  thin_tmp = forestThin_gg(prod2[i,], ttops=trees2, tree_clump_dist_m, ostory_dbh_cm,
                target_ba, target_qmd, target_pcts)
  unit2[[length(unit2)+1]] = thin_tmp
}
unit2 = do.call(rbind, unit2)

# Unit 3 has fewer than 30000 trees
unit3 = forestThin_gg(prod3, ttops=tree_locs, tree_clump_dist_m, ostory_dbh_cm,
                      target_ba, target_qmd, target_pcts)

# Unit 4
unit4 = list()
for(i in 1:nrow(prod4)){
  thin_tmp = forestThin_gg(prod4[i,], ttops=trees4, tree_clump_dist_m, ostory_dbh_cm,
                           target_ba, target_qmd, target_pcts)
  unit4[[length(unit4)+1]] = thin_tmp
}
unit4 = do.call(rbind, unit4)

# Unit 5. Didn't run because estimated ba is 57.
# unit5 = list()
# for(i in 1:nrow(prod5)){
#   thin_tmp = forestThin_gg(prod5[i,], ttops=trees5, tree_clump_dist_m, ostory_dbh_cm,
#                            target_ba, target_qmd, target_pcts)
#   unit5[[length(unit5)+1]] = thin_tmp
# }
# unit5 = do.call(rbind, unit5)

# Unit 6
unit6 = list()
for(i in 1:nrow(prod6)){
  thin_tmp = forestThin_gg(prod6[i,], ttops=trees6, tree_clump_dist_m, ostory_dbh_cm,
                           target_ba, target_qmd, target_pcts)
  unit6[[length(unit6)+1]] = thin_tmp
}
unit6 = do.call(rbind, unit6)

ttops_thinned = do.call(rbind, list(unit1,unit2,unit3,unit4,unit6))
st_write(ttops_thinned, paste0(baseDir, "Products/Vector/Units12346_thinned.gpkg"))



################################################################################
#####               Visualization                                         ######
################################################################################

take_trees = 
# take_trees = ttops_thinned[ttops_thinned$is_keep_tree==0,]
las_thinned = pc_thin(las, take_trees=take_trees, qprob=0.75)
dir.create(paste0(baseDir, "Products/las"))
writeLAS(las_thinned, paste0(baseDir, "Products/las/Unit1_thinned.las"))

# make an orthomosaic-like image that shows the thinning
thinned_ortho = ortho_thin(las_thinned)

# inpaint using the "inpaint" python script from the python functions script above.
values(thinned_ortho) = values(rast(inpaint(as.array(thinned_ortho))))
dir.create(paste0(baseDir, "Products/Raster"))
writeRaster(thinned_ortho, paste0(baseDir, "Products/Raster/thinned_ortho.png"), overwrite=T)

# plot the original and the thinned version
landscape_viz(las, texture_file = paste0(baseDir, "Products/Raster/thinned_ortho.png"))
landscape_viz(las_thinned, texture_file = paste0(baseDir, "Products/Raster/thinned_ortho.png"))

# Pay attention to which window is the Focus window
snapshot3d("../../../FCFO_CarterLake_lidar_visulazation/Products/Images/Current_example.png")
snapshot3d("../../../FCFO_CarterLake_lidar_visulazation/Products/Images/lowThin_example.png")


