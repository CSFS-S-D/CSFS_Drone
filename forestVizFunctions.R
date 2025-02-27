# Project     forestVizFunctions.R
# Purpose     Functions to visualize forestry with point clouds and aerial imagery.
# Person      Andy Whelan
# Date        November 5, 2024
# Modified    February 27, 2025
################################################################################


################################################################################
######            Functions                                                #####
################################################################################

#---> Function to denoise and color point clouds
pc_clean = function(chunk, ortho) {
  las = readLAS(chunk)
  if(lidR::is.empty(las)) return(NULL)

  las = filter_poi(las, Withheld_flag==F)
  las = classify_noise(las, ivf())
  las = filter_poi(las, Classification!=LASNOISE)
  las = classify_ground(las, csf(sloop_smooth=T))
  las = merge_spatial(las, ortho)

  return(las)
}


#---> Function to segment trees in point clouds
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

  return(las)
}


#---> Function to digitally remove trees and replace them with something ground colored.
pc_thin = function(las_files, take_trees=NULL, ttops=NULL, take=NULL, qprob=0.75) { # ttops with take/leave column,
  # take_trees = Optional tree map sfc points object with treeID field of trees to take. 
  # ttops = Optional tree map sfc points object with treeID field from which "keep" will query. 
  # take = If take_trees is null, supply a logical statement for trees to take e.g., height<10
  # qprob = quantile probability to calculate ground color. 0.5 is the mean,
  # 0.75 is the 3rd quartile which results in brighter colors that usually
  # look a little better.

  if(is.list(las_files) | is.character(las_files)) {
    las = readLAS(las_files)
  }else{
    las = las_files
  }
  if(lidR::is.empty(las)) return(NULL)
  
  
  # find a nice ground color
  R_ground = round(quantile(las$R[las$Classification == 2], prob=qprob))
  G_ground = round(quantile(las$G[las$Classification == 2], prob=qprob))
  B_ground = round(quantile(las$B[las$Classification == 2], prob=qprob))

  
  # Is there a treemap?
  if(is.null(take_trees)) {
    # take/keep tree column
    ttops$takeTree = 0 # 0 means take, 1 means keep
  
    # which trees to take
    ttops$takeTree[eval(parse(text=take))] = 1
  
    takes = which(las$treeID %in% ttops$treeID[ttops$takeTree==1] & las$Classification!=2)
    
  }else{
    takes = which(las$treeID %in% take_trees$treeID)# & las$Classification!=2)
  }
  

  las$Z[takes] = las$Z[takes] - las$treeHt[takes]
  las$R[takes] = as.integer(R_ground)
  las$G[takes] = as.integer(G_ground)
  las$B[takes] = as.integer(B_ground)

  las = normalize_height(las, knnidw())
  las = filter_poi(las, !is.na(treeID) | Z<2)
  las = unnormalize_height(las)

  return(las)
}


#---> Function to create a 3d boundary and display it on an las plot
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


#---> make an orthomosaic that reflects virtual thinning.
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


#---> Function to visualize landscapes 
landscape_viz = function(las, bg_col="skyblue", texture_file){
  # texture_file: on disk raster file to use to color the 3D surface underlying 
  # the point cloud.
  dtm = rasterize_terrain(las, res=1, knnidw())
  open3d()
  bg3d(bg_col)
  points3d(las$X, las$Y, las$Z,
           color=rgb(las$R, las$G, las$B, maxColorValue = 65535),
           size=4)
  surface3d(seq(xmin(dtm),xmax(dtm)-1,1), seq(ymax(dtm)-1,ymin(dtm),-1),
            matrix(values(dtm)[,1],nrow=ncol(dtm)), color="white", specular="black", lit=F,
            texture=texture_file)
}
