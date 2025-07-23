#' Segment trees in a point cloud.
#'
#' Segments trees in a point cloud and adds a treeID field. This function
#' normalizes the point cloud by height prior to tree segmentation then
#' unnormalizes heights to return a point cloud with the original Z data. It
#' saves ground-normalized point heights in a field called "normHt." I can't
#' remember why I thought this was useful. Maybe to remove points below a
#' certain height for better visualization or something.
#' @param las A point cloud in .las or .laz format
#' @param chm Spatraster canopy height model.
#' @param ttops SFC MULTIPOINT object with tree locations like what comes out of
#' lidR::locate_trees().
#' @return The segmented point cloud with additional height-above-ground info.
#' @examples
#' # Example code
#' data(CarterLake)
#'
#' las = pc_clean(las = CarterLake, denoise=TRUE, ground=TRUE)
#' chm = lidR::ras
#'
#' segmented_pc = pc_segment(las, chm, ttops)
#' @export
pc_segment = function(las, chm, ttops) {
  las = readLAS(las)
  if(lidR::is.empty(las)) return(NULL)

  # tree segmentation takes a normalized point cloud
  if("normHt" %in% names(las)) {
    zref = las$Z
    las$Z = las$normHt
  }else{
    las = normalize_height(las, knnidw())
  }

  las = segment_trees(las, dalponte2016(chm, ttops))

  # restore original Z values
  las = unnormalize_height(las)

  return(las)
}
