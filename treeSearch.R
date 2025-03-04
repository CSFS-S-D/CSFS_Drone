# Program   treeSearch.R
# Purpose   find the optimal search-window size to find trees using the lmf
#           function in lidR
# Person    Andy Whelan
# Date      January 29, 2025
# modified  January 29, 2025
################################################################################


################################################################################
#####                       Setup                                       ########
################################################################################

# Libraries. This downloads the latesst versions and loads them
devtools::install_github("NEONScience/NEON-geolocation/geoNEON", force=T)

packages = c("neonUtilities","neonOS","terra","tidyverse","sf","lidR",
             "geoNEON","likelihood")
for(package in packages){
  if(suppressMessages(!require(package,character.only=T))){
    install.packages(package)
    suppressMessages(library(package,character.only=T))
  }
}

options(stringsAsFactors=F)

# Load data from NEON for Niwot Ridge
veglist <- loadByProduct(dpID="DP1.10098.001", 
                         site="NIWO", 
                         package="basic", 
                         release="RELEASE-2023",
                         check.size = FALSE)


# Load NEON lidar data for Niwot Ridge
las = catalog("C:/Users/csfsuser/Downloads/NEON_lidar-point-cloud-line (1)/NEON_lidar-point-cloud-line/NEON.D13.NIWO.DP1.30003.001.2020-08.basic.20250130T162219Z.RELEASE-2025/") 


################################################################################
########                Processing                                    ##########
################################################################################

# Get plant locations
vegmap <- getLocTOS(veglist$vst_mappingandtagging, 
                    "vst_mappingandtagging")

# Add measurement data
veg <- joinTableNEON(veglist$vst_apparentindividual, 
                     vegmap, 
                     name1="vst_apparentindividual",
                     name2="vst_mappingandtagging")

# Narrow it down to trees taller than 2m that have coordinates 
veg = veg[which(veg$height>=2),]
veg = veg[which(!is.na(veg$adjEasting)),]
veg = veg[which(veg$eventID.y=="vst_NIWO_2019"),]
nrow(veg)

# Take only the tallest of multi-bole trees. Not sure if I want to do this.
# veg <- veg %>%
#   group_by(adjEasting, adjNorthing) %>%
#   filter(row_number() == which.max(height)) %>%
#   ungroup()



# Map with uncertainty
symbols(veg$adjEasting[which(veg$plotID=="NIWO_057")], 
        veg$adjNorthing[which(veg$plotID=="NIWO_057")], 
        circles=veg$stemDiameter[which(veg$plotID=="NIWO_057")]/100/2, 
        inches=F, xlab="Easting", ylab="Northing")

symbols(veg$adjEasting[which(veg$plotID=="NIWO_057")], 
        veg$adjNorthing[which(veg$plotID=="NIWO_057")], 
        circles=veg$adjCoordinateUncertainty[which(veg$plotID=="NIWO_057")], 
        inches=F, add=T, fg="lightblue")

# Get locations of plots with trees
vegplot <- veglist$vst_perplotperyear[which(
  veglist$vst_perplotperyear$treesPresent=="Y"),]

# remove duplicates
vegplot = vegplot[which(!duplicated(vegplot$plotID)),]

# Keep only plots that are in "veg"
vegplot = vegplot[vegplot$plotID %in% veg$plotID,]

# trees per plot
tpp = veg %>% group_by(plotID) %>% count() %>% as.data.frame()

# sort
vegplot = vegplot[order(vegplot$plotID),]

# check a tautology for some reason. The universe is sane, I guess.
vegplot$plotID==vegplot$plotID


# Clip out point clouds for each plot
pcs = list()
tpi = list() 
for(i in 1:nrow(vegplot)) {
  # clip with buffer so TPI covers the whole plot
  # different size tower and distributed plots. Add 1 m for buffer.
  size=11
  
  tmplas = clip_rectangle(las, vegplot$easting[i]-size, vegplot$northing[i]-size, 
                          vegplot$easting[i]+size, vegplot$northing[i]+size) %>% 
    classify_noise(sor()) %>% 
    filter_poi(Classification!=LASNOISE) %>% 
    normalize_height(knnidw()) %>% 
    filter_poi(Z>0)
  
  # TPI 
  # tmpTRI = tmplas %>% rasterize_canopy() %>% terrain(v="TRI")
  # tmplas = merge_spatial(tmplas, tmpTRI, attribute="TRI")

  # clip the buffer off
  clippedlas = clip_rectangle(tmplas, vegplot$easting[i]-(size-1), vegplot$northing[i]-(size-1), 
                              vegplot$easting[i]+(size-1), vegplot$northing[i]+(size-1))
  
  # save
  pcs[[length(pcs)+1]] = clippedlas
}
plot(pcs[[1]])
plot(pcs[[1]], color="TRI")

plot(rasterize_canopy(pcs[[1]], res=1))
ttops = locate_trees(pcs[[1]], lmf(2))
# SF object to make it easier to plot with the chm
vegsf = data.frame(plotID = veg$plotID, 
                  plantID = veg$individualID,
                  X=veg$adjEasting,
                  Y=veg$adjNorthing) %>% 
  st_as_sf(coords=c("X","Y"))
plot(vegsf[vegsf$plotID == vegplot$plotID[1],]['plotID'], pch=20, col="red", add=T)
plot(ttops['treeID'], pch=20, col="blue", add=T)
nrow(vegsf[vegsf$plotID==vegplot$plotID[2],])
nrow(ttops)
plot(vegsf[vegsf$plotID==vegplot$plotID[2],]['plotID'])


# Model dataframe
model_df = data.frame(ntrees=tpp$n)
for(i in 1:length(pcs)) {pcs[[i]]$TRI[pcs[[i]]$TRI==0 | is.na(pcs[[i]]$TRI)] = 0.1}
model_df$pcs = unlist(pcs)


################################################################################
###########                  Modeling                                ###########
################################################################################


#------------------ estimate parameter. ########################################
# set output directory
# setwd("R:/foreco4/projects/documents/ltm_lidar/Output/Chen")


  
log_model <- function(a, c, k, pcs) {
  logws = function(x) {(c)/(1+a*exp(-k*x))}
  result = data.frame(treeEst = matrix(0, nrow=length(pcs)))
  for(i in 1:length(pcs)){
    treestmp = locate_trees(pcs[[i]], lmf(logws))
    result$treeEst[i] = nrow(treestmp)
  }
  return(result$treeEst)
}

par <- list(a=20, c=1, k=1)
par_lo <- list(a=1, c=0.001, k=0)
par_hi <- list(a=50, c=10, k=20)


var <- list(mean = "predicted", x ="ntrees", log = TRUE,
            pcs="pcs")


# model -------------------------------------------------------------------#

results_log_model <-anneal(model = log_model, par = par, var = var, 
                           source_data = model_df, 
                           par_lo = par_lo,
                           par_hi = par_hi, 
                           pdf = dnorm, 
                           dep_var = "ntrees", 
                           initial_temp = 6, 
                           temp_red = 0.6, 
                           max_iter = 500)


(best <- results_log_model$best_pars)
(lower <- results_log_model$lower_limits)
(upper <- results_log_model$upper_limits)



thing = readLAS("C:/Users/csfsuser/Downloads/NEON_lidar-point-cloud-line SOAP/NEON_lidar-point-cloud-line/NEON.D17.SOAP.DP1.30003.001.2019-06.basic.20250205T223241Z.RELEASE-2025/NEON_D17_SOAP_DP1_297000_4100000_classified_point_cloud_colorized.laz") 
thing = thing %>%  normalize_height(knnidw()) %>% filter_poi(Z<50 & Z>=0)
plot(thing)
ws = function(x) {(best$c)/(1+best$a*exp(-best$k*x))}
treesws = locate_trees(thing, lmf(ws))
treesfx = locate_trees(thing, lmf(best$c))


# write_results(results_g_model, paste("Model_",e,".txt", sep=""), data=T, print_whole_hist=F)
# save(results_g_model, file=paste("Model_",e,".Rdata", sep=""))


parms_fit.lm = ("ntrees")
parms_fit.lm
model_df$ntrees
