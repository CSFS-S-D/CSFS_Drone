#' Clean up point clouds with lots of noise
#'
#' Classify noise and remove along with points flagged "Withheld" from a point
#' cloud and classify ground if needed. The idea behind this function is to get
#' the point cloud ready for further processing for visualization.
#' @param las A point cloud in .las or .laz format
#' @param denoise Boolean. Set to "TRUE" to remove erroneous points
#' @param ground Boolean. Set to "TRUE" to classify ground points
#' @return The colorized point cloud.
#' @examples
#' # Example code
#' color_las = pc_color(las=grey_las, ortho=ortho)
#' @export
pc_clean = function(las, denoise=TRUE, ground=TRUE) {
  if(lidR::is.empty(las)) return(NULL)

  las = lidR::filter_poi(las, Withheld_flag==F)
  if(denoise==T) {
    las = las %>% lidR::classify_noise(ivf()) %>%
    lidR::filter_poi(Classification!=lidR::LASNOISE)
  }

  # Classify ground
  if(ground==T) {
    las = lidR::classify_ground(las, lidR::csf(sloop_smooth=T))
  }

  # Get and save ground normalized heights and original Z value height (which
  # will be redundant with Z-values after unnormalization but can speed things
  # up later.)
  las = las %>% lidR::normalize_height(lidR::knnidw()) %>%
    lidR::add_lasattribute(las$Z, "normHt", "ground-normilized heights") %>%
    lidR::add_lasattribute(las$Zref, "origZ", "original Z heights") %>%
    lidR::unnormalize_height()

  return(las)
}
