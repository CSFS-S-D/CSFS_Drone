# Program     GCP_fixer.R
# Purpose     Correct GCPs taken with Emlid.
# Person      Andy Whelan
# Date        October 11, 2024
# Modified    October 11, 2024
################################################################################


######################## Setup ################################################

# GCP data path
gcps_path = "../StateForest/SnowCourses_lidar_102524/Drone_data/stateforest_SOWCOURSES_102524.csv"

# Corrected base station position path
base_path = "../StateForest/SnowCourses_lidar_102524/Drone_data/reach-base-raw_202410251625_RINEX-3_03/reach-base-raw_202410251625 (4).pos"



########################## Correct gcps ########################################
# The base station coordinates have been corrected using OPUS <https://geodesy.noaa.gov/OPUS/> 
# Correct the GCPs based on the difference between the orignial and corrected
# base station coordinates

# Read in data
gcps = read.csv(gcps_path)
base_cor = read.table(base_path, skip=9, header=T)

# Take position data from the base and add it to gcp data.
gcps$base_lat_cor = base_cor$latitude.deg.
gcps$base_lon_cor = base_cor$longitude.deg.
gcps$base_alt_cor = base_cor$height.m.

# calculate differences
gcps$base_lat_diff = gcps$Base.latitude-gcps$base_lat_cor
gcps$base_lon_diff = gcps$Base.longitude-gcps$base_lon_cor
gcps$base_alt_diff = gcps$Base.ellipsoidal.height-gcps$base_alt_cor

# Correct the gcps
gcps$Latitude_cor = gcps$Latitude-gcps$base_lat_diff
gcps$Longitude_cor = gcps$Longitude-gcps$base_lon_diff
gcps$Ellipsoidal.height_cor = gcps$Ellipsoidal.height-gcps$base_alt_diff


##################### Write to file ############################################
write.csv(gcps, sub(".csv", "_cor.csv", gcps_path))
