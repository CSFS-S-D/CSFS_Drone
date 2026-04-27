# Messing around with the Spanner package.

# remotes::install_github("bi0m3trics/spanner")
library(spanner)
library(lidR)
library(sf)

bounds = st_read("../../../Manitou/bounds.shp")

las = readLAS("../../../Manitou/N1_October06_2025/Terra_output/N1_200AGL_20MPH_TFOFF/cloud0.las")#,
              # filter="-keep_xy 467242 4444070 467267 4444095")
bounds = sf::st_transform(bounds, st_crs(las))
las = clip_roi(las, st_bbox(bounds))

bounds = sf::st_zm(bounds)
las = clip_roi(las, st_bbox(bounds))
las = clip_roi(las, bounds)

las = normalize_height(las, knnidw())

trees = get_raster_eigen_treelocs(las,
                                  res=0.2,
                                  pt_spacing = 0.0125,
                                  dens_threshold = 0.05,
                                  neigh_sizes = c(0.333,0.166,0.5),
                                  eigen_threshold = 0.6666,
                                  grid_slice_min = 0.5,
                                  grid_slice_max = 2,
                                  minimum_polygon_area = 0.025,
                                  cylinder_fit_type = "ransac",
                                  max_dia = 0.5,
                                  SDvert = 0.25,
                                  n_best = 25,
                                  n_pts = 20,
                                  inliers = 0.9,
                                  conf = 0.99,
                                  max_angle = 20
                                  )

plot(rasterize_canopy(las, res=0.25))
symbols(sf::st_coordinates(trees)[,1], sf::st_coordinates(trees)[,2],
        circles=trees$Radius^2*3.14, inches=F, add=T, bg="red", fg="red")
lidRviewer::view(las, color="Intensity")

trees
plot(Radius~Error, subset(trees, Error<0.01))


las2 = colorize_las(las, "pcv")
las$R = as.integer(las$R/256)
las$G = as.integer(las$G/256)
las$B = as.integer(las$B/256)

las3 = spanner::merge_las_colors(las2,las, alpha=0.5, method="alpha")
plot(las3, color="RGB")
plot(las, color="RGB")

las = spanner::segment_graph(las, trees)
plot(las, color="treeID")

trees = spanner::process_tree_data(trees, las2)
plot(height~Radius, subset(trees, Error<0.01))

ctg = catalog("../../../Manitou/N1_October06_2025/Terra_output/N1_200AGL_20MPH_TFOFF/temp/")
opt_chunk_buffer(ctg) = 10
opt_chunk_size(ctg) = 50
opt_progress(ctg) = T

func = function(chunk, trees){
  las = readLAS(chunk)
  if(is.empty(las)) return(NULL)
  
  trees = sf::st_crop(trees, st_bbox(las))
  output = spanner::segment_graph(las, trees)
  return(output)
}

thing = catalog_apply(ctg, func, trees=trees)
las = readLAS("../../../Manitou/N1_October06_2025/Terra_output/N1_200AGL_20MPH_TFOFF/temp/temp.las",
              filter="-keep_xy 490000 4330600 490025 4330625")

las = spanner::segment_graph(las, trees)
plot(las, color="Classification")
las@data$treeID
trees$treeID = trees$TreeID

sm_trees = st_crop(trees, st_bbox(las))
las = spanner::segment_graph(las, sm_trees)
plot(las, color="treeID")

for(i in 1:length(thing)) writeLAS(thing[[i]], paste0("tempSpannerLAS/segmented_",i,".las"))
