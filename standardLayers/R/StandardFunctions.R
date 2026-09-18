smoothRetile = function(las_dir, output_dir, chunk_size=200, res=1, pct=1, density=NULL) {
  
  # Check if "las_dir" is a catalog or a character path to a directory
  if(is.character(las_dir)){
    ctg = lidR::readLAScatalog(las_dir)
  }else{
    ctg=las_dir
  }
  
  # set density
  if(is.null(density)) {
    dens_list = list()
    for(i in 1:nrow(ctg)){
      dens_list[[length(dens_list)+1]] = lidR::density(ctg[i,])
    }
    density = min(unlist(dens_list))*pct
  }
  
  # Set up catalog
  lidR::opt_chunk_size(ctg) = chunk_size
  lidR::opt_chunk_buffer(ctg) = 0
  lidR::opt_filter(ctg) = ""
  lidR::opt_output_files(ctg) <- paste0(output_dir, "retiled_{XLEFT}_{YBOTTOM}")
  lidR::opt_chunk_alignment(ctg) = c(0,0)
  
  # function to retile and homogenize
  grumpy_booger = function(chunk, boog_den, boog_res){
    las = lidR::readLAS(chunk)
    if(lidR::is.empty(las)) return(NULL)
    
    lidR::decimate_points(las, lidR::homogenize(density=boog_den, res=boog_res))
  }
  
  lidR::catalog_apply(ctg, grumpy_booger, boog_den=density, boog_res=res)
  
  catalog_out = lidR::readLAScatalog(output_dir)
  return(catalog_out)
}


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

  #-------> tree_tops
  ttops_c = sf::st_intersection(ttops, bounds)
  ttops_c$treeID_ref = ttops_c$treeID
  ttops_c$treeID = 1:nrow(ttops_c)
  
  # virtual thin
  ttops_c$lowThin = 1 # 1 means cut
  ttops_c$highThin = 1
  ttops_c$uniformThing = 1
  
  ttops_c$lowThin[ttops_c$dbh_cm > 30] = 0 # 0 means it didn't get cut
  ttops_c$highThin[ttops_c$dbh_cm < 30] = 0
  ttops_c$uniformThin[sample(nrow(ttops_c), nrow(ttops_c)*0.33)] = 0

  #-------> crowns
  crowns = sf::st_intersection(crowns, bounds)

  # virtual thin
  crowns$lowThin = 1
  crowns$highThin = 1
  crowns$uniformThin = 1
  
  crowns$lowThin[crowns$treeID %in% ttops_c$treeID_ref[ttops_c$lowThin == 0]] = 0
  crowns$highThin[crowns$treeID %in% ttops_c$treeID_ref[ttops_c$highThin == 0]] = 0
  crowns$uniformThin[crowns$treeID %in% ttops_c$treeID_ref[ttops_c$uniformThin == 0]] = 0
  
  #-------> Make crown continuity maps
  clumps <- crowns %>% sf::st_union() %>% sf::st_sf() %>% sf::st_cast("POLYGON")
  clumps$area = as.numeric(sf::st_area(clumps))
  
  lowThinClumps <- subset(crowns, lowThin==0) %>% sf::st_union() %>% sf::st_sf() %>% sf::st_cast("POLYGON")
  lowThinClumps$area = as.numeric(sf::st_area(lowThinClumps))
  
  highThinClumps <- subset(crowns, highThin==0) %>% sf::st_union() %>% sf::st_sf() %>% sf::st_cast("POLYGON")
  highThinClumps$area = as.numeric(sf::st_area(highThinClumps))
  
  uniformThinClumps <- subset(crowns, uniformThin==0) %>% sf::st_union() %>% sf::st_sf() %>% sf::st_cast("POLYGON")
  uniformThinClumps$area = as.numeric(sf::st_area(uniformThinClumps))

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
  sf::st_write(lowThinClumps, paste0(output, "Products/Vector/lowThinCanopy_cont.gpkg"), append=F)
  sf::st_write(highThinClumps, paste0(output, "Products/Vector/highThinCanopy_cont.gpkg"), append=F)
  sf::st_write(uniformThinClumps, paste0(output, "Products/Vector/uniformThinCanopy_cont.gpkg"), append=F)
  sf::st_write(bounds, paste0(output, "Products/Vector/bounds.gpkg"), append=F)

  print("Congratulations! You did it!")

}


# functions to help along the way
snagPrep = function(las, boundsPath=NULL, Lmeth="Percentile", Umeth="Minimum") {  
  
  if(is.null(boundsPath)){
    bounds = st_bbox(las)
  }else{
    bounds = st_read(boundsPath)
    bounds = st_transform(bounds, st_crs(las))
    # If the simple feature collection has a Z dimension, remove it.
    bounds = st_zm(bounds)
  }
  
  if(class(las) == "LAS"){
    pc = filter_poi(las, Withheld_flag==F & ReturnNumber==1)
    pc = clip_roi(pc, st_bbox(bounds))
    pc = clip_roi(pc, bounds)
  }else{
    pc = clip_roi(las, st_bbox(bounds))
    pc = clip_roi(pc, bounds)
  }
  
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
  
  
  # Scale if values are really low
  if(max(pc$Intensity)<10000) pc$Intensity = as.integer(pc$Intensity*((2^16-1)/max(pc$Intensity)))
  
  pc$Intensity = as.integer(pc$Intensity/(2^16-1)*(2^8-1))
  
  
  # Automatic thresholding
  pcu = lidR::filter_poi(pc, Z>2)$Intensity
  Lint = autothresholdr::auto_thresh(pcu, method=Lmeth)
  Uint = autothresholdr::auto_thresh(pcu, method=Umeth)
  pc@index$Lint = Lint[1]
  pc@index$Uint = Uint[1]
  
  return(pc)
}

snagMapper = function(pc, Lint, Uint) {
  
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
  
  over_pc = filter_poi(pc, Z>=2)
  
  # Set lower and upper intensity thresholds. For Front Range ponderosa pine, a 
  # lower threshold of 35 or a little more seems to work pretty well. I use the modes
  # to help adjust it. I add the value for the first mode to 35. Sometimes that's a
  # little high. We can further adjust it later. For Uint, I just use the last mode -
  # it's close to the max intensity, which works pretty well.
  
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
  
  
  # This gets the point data as a data frame.
  df = over_pc@data
  
  # combine over_pc with the part of the original point cloud <=2
  thing = filter_poi(pc, Z<2)
  thing = add_lasattribute(thing, 0, "snagCls", "Number identifying a snag class")
  pc = rbind(thing, over_pc)
  
  return(pc)
}

add_LasSnagAttributes = function(pc, ttops) {
  
  # tree height
  thing = dplyr::left_join(pc@data, ttops, by="treeID")
  thing$tree_height_m[is.na(thing$tree_height_m)] = 0
  pc = add_lasattribute(pc, thing$tree_height_m, "treeHt", "TreeHeight")

  # dbh
  thing$dbh_cm[is.na(thing$dbh_cm)] = 0
  pc = add_lasattribute(pc, thing$dbh_cm, "dbh", "TreeDBH")

  
  # This is a custom function we'll use in the aggregate command below. All it does
  # is find the ratio of points classified as snag to all points. The aggregate command
  # will apply it to each tree for us.
  snaggy = function(x){
    length(x[x!=0])/length(x)
  }
  
  # Use "aggregate" to find snag to not snag point ratio for each tree.
  t = aggregate(snagCls~treeID, lidR::filter_poi(pc, Z>2)@data, snaggy)

  thing = dplyr::left_join(pc@data, t, by="treeID")
  thing$snagCls.y[is.na(thing$snagCls.y)] = 0
  pc = add_lasattribute(pc, thing$snagCls.y, "snagPct", "PercentSnagPoints")
  
  return(pc)
}
