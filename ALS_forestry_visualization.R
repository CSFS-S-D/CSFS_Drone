# Project     ALS_forestry_visualization.R
# Purpose     Visualize forest management with lidar and aerial imagery
# Person      Andy Whelan
# Date        November 5, 2024
# Modified    November 5, 2024
################################################################################

library(reticulate)
library(terra)



# Install python modules if needed
# py_install("opencv-python")
# py_install("numpy")

# Load python functions from ALS_viz_functions.py
source_python("ALS_viz_functions.py")

# Load the ortho image
img_path = "../../../FCFO_CarterLake_lidar_visulazation/naip/m_4010539_sw_13_1_20130716/m_4010539_sw_13_1_20130716.tif"



########################## Remove shadows ######################################
result = remove_shadow(img_path)

# The resulting raster is upside down
thing = rast("../../../FCFO_CarterLake_lidar_visulazation/naip/bananas.tif")
t2 = rast(img_path)
ext(thing) = ext(t2)
crs(thing) = crs(t2)
thing = project(thing, t2)
xmid = xmin(thing)+ncol(thing)/2
ymid = ymin(thing)+nrow(thing)/2
plotRGB(thing, ext=c(xmid,xmid+250,ymax(thing)-250,ymax(thing)))
plotRGB(t2, ext=c(xmid,xmid+250,ymax(thing)-250,ymax(thing)))
plotRGB(thing, axes=T)
plotRGB(t2)

# Display the image
cv2.imshow('Shadow Removed', result)
cv2.imshow('Original', image)
cv2.waitKey(0)
cv2.destroyAllWindows()
#########################################################333

