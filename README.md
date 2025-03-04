<<<<<<< HEAD
# CSFS_DataHub
### GCP_fixer.R
This script came about because getting corrected coordinates for GCPs collected with a paired base station and rover is a multistep process - even with real-time corrections sent from the base. Although the base is stationary, if it's not connected to a correction service (e.g., NTRIP, Trimble RTX), its absolute position cannot be determined without post processing. This also means that although the rover coordinates are accurate within several centimeters of the base station, they're subject to the same errors as the base station. So, after we find the absolute position of the base station using OPUS and additional RINEX data collected by the base station, we need to adjust the rover coordinates the same amount. The GCP_fixer script takes a csv with coordiates of the base station and GCPs, and uses corrected coordinates from OPUS to adjust the GCPs to the best absolute locations possible.  
### StandardLayers.R
StandarLayers takes point cloud products from initial lidar or SFM processing (e.g., DJI Terra, Pix4d) and creates a standard set of raster and vector map layers, and some simple graphs.
### StandardReport.R
Creates a standard report which includes tables and figures and some explanation about how and why they were created.
### ALS_forestry_visualization.R
Currently does a passable job of removing shadows from NAIP imagery. In the future, it will use the shadowless NAIP imagery to map RGB color information to ALS point clouds and display 3D scenes of current and future forest conditions.
### ALS_viz_functions.py
Python OpenCV functions to help remove shadows from NAIP imagery. You need this to run ALS_forestry_visualization.R.
=======
# CSFS_Drone
Scripts to do forestry with drone data
>>>>>>> 64431a52072d4f4cc1466d41e6fc608b1b4a5750
