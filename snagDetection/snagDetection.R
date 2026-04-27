# Program   snagDetection.R
# Purpose   To detect snags in aerial lidar point clouds
# Person    Andy Whelan
# Date      September 11, 2025
# Modified  April 27, 2026
##################################################################################


################################################################################
#######################    Setup    ############################################
################################################################################
# Libraries
library(lidR)
library(lidRviewer)
library(terra)
library(sf)
library(terrainr)

# Environment
lasDir = "../../Frisco_GRFO/LAS/"
boundsPath = "../../Frisco_GRFO/Shapefiles/FriscoBackyard_Units/FriscoBackyard_Units/FriscoBackyard_Units.shp"
lasHeaderEx = readLASheader(dir(lasDir, full.names=T)[2])


################################################################################
##########################    Data    ##########################################
################################################################################
# Wing (2015) suggest to only use first returns, so we can just drop the rest
# when we load the point cloud. The "-drop_withheld" flag basically drops a bunch
# of points that have been flagged as bad. The drone lidar data probably won't 
# have any withheld, but the USGS lidar data can have a lot.
ctg = catalog(lasDir, filter="-drop_withheld -keep_first")
bounds = st_read(boundsPath)
bounds = st_transform(bounds, st_crs(lasHeaderEx))

# subset if necessary
bounds = bounds[3,]

# clip the catalog to the bounds. This also loads the actual point cloud.
pc = clip_roi(ctg, bounds)
pc$Z = pc$Z*0.3048 # convert feet to meters. XY are already in meters.

################################################################################
########                      Data Wranglin                             ########
################################################################################

# The first thing we need to do is get the Density right. Look at the density 
# raster and homogenize to lower densities until there aren't obvious (i.e., stripey) 
# high-density areas. UPDATE: While it turns out you can get rid of the density
# bands by setting the density value to around 200, it causes the snag segmentation
# algrorithm to take forever. Forever as in maybe days or weeks. I waited half a 
# day before I threw in the towel. Anyway, I decided to reduce the density to something
# closer to the example data from the snag detection help. The algorithm was
# developed using lower density point clouds, I think.

plot(rasterize_density(pc)) 
pc = decimate_points(pc, homogenize(density = 1, res=5))

# Normalize height for snag detection
pc = normalize_height(pc, knnidw())

# Intensity values need to be 8-bit. They are often 16-bit. So, we
# need to convert.
summary(pc$Intensity)
pc$Intensity = as.integer(pc$Intensity/(2^16-1)*(2^8-1))

# Histograms help us find the lower and upper thresholds. Look for two peaks and
# choose thresholds just higher than each. Adjust the histogram as necessary.
hist(pc$Intensity)# breaks=200, xlim=c(0,25))

# Calculate bole and branch to foliage intensity ratio. Use the two values just 
# after the peaks from the histogram.
over_pc = filter_poi(pc, Z>=2)
bbvf = (length(over_pc$Intensity[over_pc$Intensity<=3 | over_pc$Intensity>=13]))/
  (length(over_pc$Intensity[over_pc$Intensity>3 | over_pc$Intensity<13]))

# Calculate these (lower and upper intensity thresholds)
Lint = 20*bbvf+0.075*(max(pc$Intensity)) + 26.5
Uint = 20*bbvf+0.1875*(max(pc$Intensity)) + 100.25

# Some wisdom from Wing etal. 2015 below
# If Lint<50, Lint=50, If Lint>70, Lint=70
# If Uint<150, Uint=150 else if Uint>170, Uint=170
# If Located in SF Area, Lint=Lint+5 else if SF, Uint=Uint-5. SF is "Storrie Fire" 
# and not relevant to us but they adjusted Lint and Uint for it, so we might read Wing 2015 more closely to see how we might
# come up with adjustments of our own.

if(Lint<50) Lint = 50
if(Lint>70) Lint = 70
if(Uint<150) Uint = 150
if(Uint>170) Uint = 170

# from Wing etal. 2015. I don't see much reason to fiddle with these. They're the same for both point clouds.
bbpr_thresholds <- matrix(c(0.80, 0.80, 0.70,
                            0.85, 0.85, 0.60,
                            0.80, 0.80, 0.60,
                            0.90, 0.90, 0.55),
                          nrow =3, ncol = 4)

# Get the point density requirement (PDR) from the Plot-level point density (PLPD). 
# We use the entire point cloud for PLPD. If PLPD <= 3, PDR = 3; if 3 < PLPD <= 6, 
# PDR = 4; if 6 < PLPD <= 12, PDR = 5; if PLPD < 12, PDR = 8.
if(density(pc)<=3) {
  pdr = 3
}else if(density(pc)>3 & density(pc)<=6){
  pdr = 4
}else if(density(pc)>6 & density(pc)<=12){
  pdr = 5
}else{
  pdr = 8
}

# Run the stem segmentation algorithm using the appropriate thresholds.
over_pc = segment_snags(over_pc, wing2015(BBPRthrsh_mat = bbpr_thresholds, pt_den_req = 3,
                                            low_int_thrsh = Lint, uppr_int_thrsh = Uint))

# Plot it
plot(over_pc, color="snagCls", legend=T)
plot(filter_poi(over_pc, snagCls>0), color="snagCls", pal=rainbow(5)[-1], legend=T)


# We now have a bunch of points that are classified as snag points, but what we
# really need to know is which trees are snags. So the first thing we'll need to do
# is find the trees. 
over_chm = rasterize_canopy(over_pc)


# I used a window size function from the cloud2trees package that increases the
# window size exponentially with height. Adjust as necessary to get the best
# tree identification.
over_ttops = locate_trees(over_chm, lmf(cloud2trees::itd_ws_functions()$exp_fn))


# segment the point cloud
over_pc = segment_trees(over_pc, dalponte2016(over_chm, over_ttops))
plot(over_pc, color="treeID")

# Calculate the ratio of snag points in each tree so we can classify trees with a
# high ratio as snags.

# This gets the point data as a data frame.
df = over_pc@data

# Take a look at the snag class data.
df$snagCls


# This is a custom function we'll use in the aggregate command below. All it does
# is find the ratio of points classified as snag to all points. The aggregate command
# will apply it to each tree for us.
snaggy = function(x){
  length(x[x!=0])/length(x)
}

# Use "aggregate" to find snag to not snag point ratio for each tree.
t = aggregate(snagCls~treeID, df, snaggy)


# Choose a threshold percentage of snag points to classify snags (50%, 75%?). 
# This code gives us the treeID numbers of snags.
sngs = t$treeID[t$snagCls>0.50]

# Here we subset the point cloud so we can just look at the snags.
over_snags = over_pc[over_pc$treeID %in% sngs,]
plot(over_snags)
view(over_snags)


# This looks alright but there's still some junk/noise we might want to get rid of.
# We're often more interested in taller snags - they represent a lot more biomass,
# they're more dangerous to people, and also better for wildlife. Each identified
# tree has a tree height (Z) associated with it. These data are in the over_ttops data
# frame. We could subset the over_ttops data frame to trees taller than 5-10 meters or
# trees larger than 5 cm dbh and then use the treeIDs from that subset to select 
# only those trees from the over_snags point cloud. 



# filter by tree height
# snags taller than 10m
over_ttops10 = over_ttops[over_ttops$Z > 10,]
over_snags10 = over_snags[over_snags$treeID %in% over_ttops10$treeID,]
nrow(over_snags10)
view(over_snags10)


snags_10_3 = aggregate(Z~treeID, over_snags10@data, min)
snags_10_3 = snags_10_3[snags_10_3$Z < 3,]

over_snags10_3 = over_snags10[over_snags10$treeID %in% snags_10_3$treeID,]
plot(over_snags10_3, color="treeID")

nrow(snags_10_3)


# Snags taller than 5m
over_ttops5 = over_ttops[over_ttops$Z > 5,]
over_snags5 = over_snags[over_snags$treeID %in% over_ttops5$treeID,]
nrow(over_snags5)
view(over_snags5)

snags_5_3 = aggregate(Z~treeID, over_snags5@data, min)
snags_5_3 = snags_5_3[snags_5_3$Z < 3,]

over_snags5_3 = over_snags5[over_snags5$treeID %in% snags_5_3$treeID,]
plot(over_snags5_3, color="treeID")

nrow(snags_5_3)

# dbh 
over_ttops$dbh = over_ttops$Z*0.3707240+1.1857200
# over_ttops = st_intersection(over_ttops, bounds);
over_ttops_dbh5 = over_ttops[over_ttops$dbh >= 5,]
over_snags_dbh5 = over_snags[over_snags$treeID %in% over_ttops_dbh5$treeID,]
nrow(over_snags_dbh5)
view(over_snags_dbh5)

snags_dbh5_3 = aggregate(Z~treeID, over_snags_dbh5@data, min)
snags_dbh5_3 = snags_dbh5_3[snags_dbh5_3$Z < 3,]

over_snags_dbh5_3 = over_snags_dbh5[over_snags_dbh5$treeID %in% snags_dbh5_3$treeID,]
plot(over_snags_dbh5_3, color="treeID")

nrow(snags_dbh5_3)

