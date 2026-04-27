# Code to do snag detection
# February 9, 2026
################################################################################


# Setup
library(lidR)
library(lidRviewer)
library(terra)
library(sf)
library(terrainr)


# Data (hopefully your computer can handle this)
# Wing (2015) suggest to only use first returns, so we can just drop the rest
# when we load the point cloud. The "-drop_withheld" flag basically drops a bunch
# of points that have been flagged as bad. The drone lidar data probably won't 
# have any withheld, but the USGS lidar data can have a lot.
las = readLAS("../../../Manitou/N1_October06_2025/Terra_output/N1_200AGL_20MPH_TFOFF/cloud0.las", filter="-drop_withheld -keep_first")
ctg = catalog("../../../Frisco_GRFO/LAS", filter="-drop_withheld -keep_first")
readLASheader("../../../Frisco_GRFO/LAS/USGS_LPC_CO_Central_and_WesternCO_2016_A16_LD28261617_colorized_1_1_1.las")
bounds = st_read("../../../Manitou/bounds.shp")
bounds = st_transform(bounds, st_crs(las))

# Use "view" from the lidRviewer package for large point clouds. It's way better
# about memory usage.
view(las, color="RGB")

# clip the catalog to the bounds. This also loads the actual point cloud.
pc2 = clip_roi(ctg, st_bbox(c(xmin=406000, xmax=406500, ymin=4379800, ymax=4380300)))


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
plot(rasterize_density(las))
las = decimate_points(las, homogenize(density=20, res=5))

plot(rasterize_density(pc2)) # This one looks good.
pc2 = decimate_points(pc2, homogenize(density = 1, res=5))

# Normalize height for snag detection
las = normalize_height(las, knnidw())
pc2 = normalize_height(pc2, knnidw())

# Intensity values need to be 8-bit. From the drone, they are 16-bit. So, we
# need to convert.
las$Intensity = as.integer(las$Intensity/(2^16-1)*(2^8-1))
# pc2$Intensity = as.integer(pc2$Intensity/(2^13-1)*(2^8-1))
pc2$Intensity = as.integer(pc2$Intensity/(max(pc2$Intensity)/255))

# Histograms help us find the lower and upper thresholds. Look for two peaks and
# choose thresholds just higher than each.
hist(las$Intensity)
hist(pc2$Intensity)

# Calculate bole and branch to foliage intensity ratio for drone data.
over_las = filter_poi(las, Z>=2)
bbvf = (length(over_las$Intensity[over_las$Intensity<=50 | over_las$Intensity>=170]))/
  (length(over_las$Intensity[over_las$Intensity>50 | over_las$Intensity<170]))

# Calculate these (lower and upper intensity thresholds)
Lint = 20*bbvf+0.075*(max(las$Intensity)) + 26.5
uint = 20*bbvf+0.1875*(max(las$Intensity)) + 100.25

# Calculate bole and branch to foliage intensity ratio for aerial data.
over_pc2 = filter_poi(pc2, Z>=2)
bbvf_pc2 = (length(over_pc2$Intensity[over_pc2$Intensity<=30 | over_pc2$Intensity>=170]))/
  (length(over_pc2$Intensity[over_pc2$Intensity>30 | over_pc2$Intensity<170]))

# Calculate these (lower and upper intensity thresholds)
Lint_pc2 = 20*bbvf_pc2+0.075*(max(pc2$Intensity)) + 26.5
uint_pc2 = 20*bbvf_pc2+0.1875*(max(pc2$Intensity)) + 100.25

# Some wisdom from Wing etal. 2015 below
# If Lint<50, Lint=50, If Lint>70, Lint=70
# If Uint<150, Uint=150 else if Uint>170, Uint=170
# If Located in SF Area, Lint=Lint+5 else if SF, Uint=Uint-5. SF is "Storrie Fire" 
# and not relevant to us but they adjusted Lint and Uint for it, so we might read Wing 2015 more closely to see how we might
# come up with adjustments of our own.


# from Wing etal. 2015. I don't see much reason to fiddle with these. They're the same for both point clouds.
bbpr_thresholds <- matrix(c(0.80, 0.80, 0.70,
                            0.85, 0.85, 0.60,
                            0.80, 0.80, 0.60,
                            0.90, 0.90, 0.55),
                          nrow =3, ncol = 4)


# Run the stem segmentation algorithm using the appropriate thresholds.
over_las = segment_snags(over_las, wing2015(BBPRthrsh_mat = bbpr_thresholds, pt_den_req = 8,
                                            low_int_thrsh = Lint, uppr_int_thrsh = uint))

over_pc2 = segment_snags(over_pc2, wing2015(BBPRthrsh_mat = bbpr_thresholds, pt_den_req = 1,
                                            low_int_thrsh = 40, uppr_int_thrsh = 140))


plot(over_las, color="snagCls", legend=T)
plot(over_pc2, color="snagCls", legend=T)
plot(filter_poi(over_las, snagCls>0), color="snagCls", pal=rainbow(5)[-1], legend=T)
plot(filter_poi(over_pc2, snagCls>0), color="snagCls", pal=rainbow(5)[-1], legend=T)


# We now have a bunch of points that are classified as snag points, but what we
# really need to know is which trees are snags. So the first thing we'll need to do
# is find the trees. I realize you're probably well aquainted with most of this.
# find individual trees
over_chm = rasterize_canopy(over_las)
over_pc2_chm = rasterize_canopy(over_pc2)

# I used a window size function from the cloud2trees package that increases the
# window size exponentially with height. It seemed appropriate because The cloud2trees
# window size functions were actually developed largely on data from the Manitou
# experimental forest.
over_ttops = locate_trees(over_chm, lmf(cloud2trees::itd_ws_functions()$exp_fn))
over_pc2_ttops = locate_trees(over_pc2_chm, lmf(cloud2trees::itd_ws_functions()$exp_fn))

# segment the point cloud
over_las = segment_trees(over_las, dalponte2016(over_chm, over_ttops))
over_pc2 = segment_trees(over_pc2, dalponte2016(over_pc2_chm, over_pc2_ttops))
plot(over_las, color="treeID")
plot(over_pc2, color="treeID")

# Calculate the ratio of snag points in each tree so we can classify trees with a
# high ratio as snags.

# This gets the point data as a data frame.
df = over_las@data
df2 = over_pc2@data

# Take a look at the snag class data.
df$snagCls
df2$snagCls

# This is a custom function we'll use in the aggregate command below. All it does
# is find the ratio of points classified as snag to all points. The aggregate command
# will apply it to each tree for us.
snaggy = function(x){
  length(x[x!=0])/length(x)
}

# Use "aggregate" to find snag to not snag point ratio for each tree.
t = aggregate(snagCls~treeID, df, snaggy)
t2 = aggregate(snagCls~treeID, df2, snaggy)

# Let's just say trees with 75% snag points are snags. This code gives us the
# treeID numbers of snags.
sngs = t$treeID[t$snagCls>0.50]
sngs2 = t2$treeID[t2$snagCls>0.50]

# Here we subset the point cloud so we can just look at the snags.
over_snags = over_las[over_las$treeID %in% sngs,]
over_pc2_snags = over_pc2[over_pc2$treeID %in% sngs2,]
plot(over_snags)
view(over_snags)
view(over_pc2_snags)

# This looks alright but there's still some junk/noise we might want to get rid of.
# We're often more interested in taller snags - they represent a lot more biomass,
# they're more dangerous to people, and also better for wildlife. Each identified
# tree has a tree height (Z) associated with it. These data are in the over_ttops data
# frame. We could subset the over_ttops data frame to trees taller than 5-10 meters 
# and then use the treeIDs from that subset to select only those trees from the 
# over_snags point cloud. Do you want to figure out how to do it?


# The other thing you could do, if you're itching for more R fun, is to apply this 
# code to USGS lidar data for the same area. You'll have to think carefully about
# the settings for the segment_snag and wing2015 functions.


# filter by tree height
# snags taller than 10m

over_ttops10 = over_ttops[over_ttops$Z > 10,]
over_snags10 = over_snags[over_snags$treeID %in% over_ttops10$treeID,]

snags_10_3 = aggregate(Z~treeID, over_snags10@data, min)
snags_10_3 = snags_10_3[snags_10_3$Z < 3,]

over_snags10_3 = over_snags10[over_snags10$treeID %in% snags_10_3$treeID,]
plot(over_snags10_3, color="treeID")

nrow(snags_10_3)


# Snags taller than 5m
over_ttops5 = over_ttops[over_ttops$Z > 5,]
over_snags5 = over_snags[over_snags$treeID %in% over_ttops5$treeID,]

snags_5_3 = aggregate(Z~treeID, over_snags5@data, min)
snags_5_3 = snags_5_3[snags_5_3$Z < 3,]

over_snags5_3 = over_snags5[over_snags5$treeID %in% snags_5_3$treeID,]
plot(over_snags5_3, color="treeID")

nrow(snags_5_3)

# dbh 
over_ttops$dbh = over_ttops$Z*0.3707240+1.1857200
over_ttops = st_intersection(over_ttops, st_transform(bounds, st_crs(las)));
over_ttops_dbh5 = over_ttops[over_ttops$dbh >= 5,]
over_snags_dbh5 = over_snags[over_snags$treeID %in% over_ttops_dbh5$treeID,]

snags_dbh5_3 = aggregate(Z~treeID, over_snags_dbh5@data, min)
snags_dbh5_3 = snags_dbh5_3[snags_dbh5_3$Z < 3,]

over_snags_dbh5_3 = over_snags_dbh5[over_snags_dbh5$treeID %in% snags_dbh5_3$treeID,]
plot(over_snags_dbh5_3, color="treeID")

nrow(snags_dbh5_3)
nrow(over_ttops_dbh5)

# dbh for aerial point cloud 
over_pc2_ttops$dbh = over_pc2_ttops$Z*0.3707240+1.1857200
over_pc2_ttops = st_intersection(over_pc2_ttops, st_transform(bounds, st_crs(pc2)))
over_pc2_ttops_dbh5 = over_pc2_ttops[over_pc2_ttops$dbh >= 5,]
over_pc2_snags_dbh5 = over_pc2_snags[over_pc2_snags$treeID %in% over_pc2_ttops_dbh5$treeID,]

snags_pc2_dbh5_3 = aggregate(Z~treeID, over_pc2_snags_dbh5@data, min)
snags_pc2_dbh5_3 = snags_pc2_dbh5_3[snags_pc2_dbh5_3$Z < 3,]

over_pc2_snags_dbh5_3 = over_pc2_snags_dbh5[over_pc2_snags_dbh5$treeID %in% snags_pc2_dbh5_3$treeID,]
plot(over_pc2_snags_dbh5_3, color="treeID")

nrow(snags_pc2_dbh5_3)
nrow(over_pc2_ttops_dbh5)
