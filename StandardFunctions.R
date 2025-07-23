standardLayers = function(las_dir, output="./", bounds=NULL, ortho=NULL, dtm_res=1, chm_res=1,
                          ws=cloud2trees::itd_ws_functions()$log_fn) {
  
  library(magrittr)
  
  # For ws function suggest something from itd_ws_functions() in cloud2trees
  
  # create output file structure
  dir.create(paste0(output, "Products/Cloud2Trees_products/"), recursive=T)
  dir.create(paste0(output, "Products/Raster/"))
  dir.create(paste0(output, "Products/Vector/"))
  dir.create(paste0(output, "Products/las/"))
  dir.create(paste0(output, "Products/las/retiled/"))
  dir.create(paste0(output, "Products/las/segmented/"))
  dir.create(paste0(output, "Products/Raster/temp/"))
  dir.create(paste0(output, "Products/Graphs/"))
  dir.create(paste0(output, "Products/Animations/"))
  
  
  # path to shapefile of project area boundaries.
  if(!is.null(bounds)) {
    bounds = sf::st_read(bounds)
  }
  
  # Note: The user should retile, thin/homogenize, and otherwise prep the catalog for consistency 
  # before using this function.
  
  # Run cloud2trees
  ans = cloud2trees::cloud2trees(output_dir=paste0(output, "Products/Cloud2Trees_products/"), 
                       las_dir, 
                       dtm_res_m=dtm_res, 
                       chm_res_m=chm_res,
                       estimate_tree_dbh = T,
                       ws = ws,
                       estimate_dbh_from_cloud=F)
  
  
  
  ######################## Make map layers #######################################
  # load data from disk
  chm = terra::rast(paste0(output, "Products/Cloud2Trees_products/point_cloud_processing_delivery/chm_", 
                    as.character(chm_res),"m.tif"))
  dtm = terra::rast(paste0(output, "Products/Cloud2Trees_products/point_cloud_processing_delivery/dtm_", 
                    as.character(dtm_res),"m.tif"))
  ttops = sf::st_read(paste0(output, "Products/Cloud2Trees_products/point_cloud_processing_delivery/final_detected_tree_tops.gpkg"))
  crowns = sf::st_read(paste0(output, "Products/Cloud2Trees_products/point_cloud_processing_delivery/final_detected_crowns.gpkg"))
  
  
  # rename a couple of things
  names(chm) = "chm"
  names(dtm) = "dtm"
  
  # reproject bounds to match chm
  if(!is.null(bounds)) {
    bounds = sf::st_read(bounds) 
    bounds = bounds %>% sf::st_transform(crs = sf::st_crs(chm))
  }else{
    bounds = sf::st_as_sfc(sf::st_bbox(ttops))
  }
  
  
  # ----------------------> Make layers <------------------------------------
  # orthomosaic
  if(!is.null(ortho)){
    ortho = terra::crop(ortho, bounds, mask=T)
    terra::writeRaster(ortho, paste0(output, "Products/Raster/ortho_cropped.tif"))
  }
  
  # tree_tops
  ttops_c = sf::st_intersection(ttops, bounds)
  ttops_c$treeID_ref = ttops_c$treeID
  ttops_c$treeID = 1:nrow(ttops_c)
  
  # crowns 
  crowns = sf::st_intersection(crowns, bounds)
  
  # Make a sort of crown continuity map
  clumps <- crowns %>% sf::st_union() %>% sf::st_sf() %>% sf::st_cast("POLYGON")
  clumps$area = as.numeric(sf::st_area(clumps))
  
  # slope and aspect 
  slope = terra::terrain(dtm, "slope")
  aspect = terra::terrain(dtm, "aspect")
  
  
  
  # --- Calculate crown Density ----
  
  # function to segment point clouds and rasters by tree
  crown_cat = function(chunk, chm, ttops, res) {
    las = lidR::readLAS(chunk)
    if(lidR::is.empty(las)) return(NULL)
    if(is.null(terra::intersect(lidR::ext(las), terra::ext(chm)))) return(NULL)
    xmin=round(lidR::ext(las)[1])
    ymin=round(lidR::ext(las)[3])
    
    las = lidR::classify_noise(las, lidR::ivf(res=1, n=10))
    las = lidR::filter_poi(las, Classification!=lidR::LASNOISE)
    las = lidR::classify_ground(las, lidR::csf(sloop_smooth=TRUE, rigidness = 2))
    las = lidR::normalize_height(las, lidR::knnidw(k=10, p=2)) # Need a normalized las for vertical density
    las = lidR::filter_poi(las, Z>=0)
    las = lidR::segment_trees(las, algorithm = lidR::dalponte2016(chm, ttops, max_cr=10))
    fun1 = ~list(treeID = mean(treeID))
    crown_rast = lidR::pixel_metrics(las, fun1, res)
    func2 = ~list(vd = length(Z[Z>2])/length(Z))
    verticalDens1 = lidR::pixel_metrics(las, func2, res=res)
    
    rStack = c(crown_rast, verticalDens1)
    names(rStack) =c("cr", "vd")
    terra::writeRaster(rStack, paste0(output, "Products/Raster/temp/crown_rasts_",
                               xmin,"_",ymin,".tif"))
    
    las = lidR::unnormalize_height(las)
    return(las)
  }
  
  
  
  # ---- set up the catalog processing ----
  # Segment point clouds and create density and crown area rasters
  ctg = lidR::readLAScatalog(las_dir)
  lidR::opt_chunk_size(ctg) = 200
  lidR::opt_chunk_buffer(ctg) = 10
  lidR::opt_filter(ctg) = ""
  lidR::opt_output_files(ctg) <- paste0(output, "Products/las/segmented/segmented_{XLEFT}_{YBOTTOM}")
  lidR::opt_chunk_alignment(ctg) = c(0,0)
  opt = list(raster_alignment=1)
  lidR::plot(ctg, chunk_pattern=T)
  
  results = lidR::catalog_apply(ctg, crown_cat, chm, ttops_c, res=chm_res, .options=opt)
  
  
  # Crop crown area and crown density rasters
  rast_list = list.files(paste0(output, "Products/Raster/temp/"), full.names = T)
  rasts = do.call(terra::merge, lapply(rast_list, function(x) terra::rast(x)))
  rasts = terra::crop(x=rasts, y=terra::vect(bounds), mask=T)
  crown_rast = rasts$cr
  density = rasts$vd
  
  
  # put everything in a raster stack
  cropped = lapply(list(chm,dtm,slope,aspect,crown_rast,density), function(x) terra::crop(x,y=bounds, mask=T))
  rasters = do.call(c, cropped)
  
  # make a digital surface model
  rasters$dsm = terra::classify(rasters$chm, matrix(c(NA,0), ncol=2))+rasters$dtm
  terra::varnames(rasters$dsm)
  
  ###################################### Save products ###########################################################
  # Rasters
  for(i in 1:length(names(rasters))) {
    terra::writeRaster(rasters[[i]], paste0(output, "Products/Raster/",names(rasters)[i],".tif"), overwrite=T)
  }
  
  # Vector
  sf::st_write(ttops_c, paste0(output, "Products/Vector/ttops_c.gpkg"), append=F)
  sf::st_write(crowns, paste0(output, "Products/Vector/crowns.gpkg"), append=F)
  sf::st_write(clumps, paste0(output, "Products/Vector/canopy_cont.gpkg"), append=F)
  sf::st_write(bounds, paste0(output, "Products/Vector/bounds.gpkg"), append=F)
  
  print("Congratulations! You did it!")
  
}
dir.create("standardTest/")

standardLayers("../../Pinecliffe_BUFO/Arc/Products/las_segmented/test/", "standardTest/",dtm_res = 1,
               chm_res = 1, ws=3)
