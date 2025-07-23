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
# ## install lasR from the r-univers
# install.packages("lasR", repos = "https://r-lidar.r-universe.dev")
# 
# ## install TreeLS from github
# remotes::install_github(repo = "tiagodc/TreeLS", upgrade = F)
# 
# ## install LadderFuelsR
# remotes::install_github(repo = "olgaviedma/LadderFuelsR", upgrade = F)
# ## install leafR
# remotes::install_github(repo = "DRAAlmeida/leafR", upgrade = F)
# 
# ## get cloud2trees
# remotes::install_github(repo = "georgewoolsey/cloud2trees", upgrade = F)
# 
# # Libraries. This downloads the latest versions and loads them
# devtools::install_github("NEONScience/NEON-geolocation/geoNEON", force=T)

packages = c("neonUtilities","neonOS","terra","tidyverse","sf","lidR",
             "geoNEON","likelihood", "rcompanion", "cloud2trees")
for(package in packages){
  if(suppressMessages(!require(package,character.only=T))){
    install.packages(package)
    suppressMessages(library(package,character.only=T))
  }
}

options(stringsAsFactors=F)

# Load NEON data
veglist <- loadByProduct(dpID="DP1.10098.001", 
                         site="NIWO", 
                         package="basic", 
                         release="RELEASE-2023",
                         check.size = FALSE)

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


# Narrow it down to trees taller than 2m  
veg = veg[which(veg$height>=2),]

# remove all but the tallest duplicates
veg = veg %>% group_by(plotID, individualID) %>% filter(height==max(height))
veg = veg[!duplicated(veg$individualID),]

# remove dead trees
veg = veg[grep("Live", veg$plantStatus),]

# only trees with species information
veg = subset(veg, !is.na(veg$taxonID))


# Take only the tallest of multi-bole trees. Not sure if I want to do this.
# veg <- veg %>%
#   group_by(adjEasting, adjNorthing) %>%
#   filter(row_number() == which.max(height)) %>%
#   ungroup()


# Get locations of plots with trees
vegplot <- veglist$vst_perplotperyear[which(
  veglist$vst_perplotperyear$treesPresent=="Y"),]

# remove duplicates
vegplot = vegplot[which(!duplicated(vegplot$plotID)),]

# Keep only plots that are in "veg"
vegplot = vegplot[vegplot$plotID %in% veg$plotID,]

# figure out which subplots were sampled
# first narrow veg down to plots with coordinates
vegjunior = veg[which(!is.na(veg$adjEasting)),]
vegjunior = vegjunior %>% filter(plotID %in% vegplot$plotID[vegplot$plotType=="tower"])
subs = aggregate(subplotID~plotID, vegjunior, unique) %>% 
  {lapply(.$subplotID, `length<-`, 4)} %>% 
  unlist() %>% 
  matrix(ncol=4, byrow=T) %>% 
  as.data.frame() %>% 
  mutate(plotID=unique(vegjunior$plotID)[order(unique(vegjunior$plotID))]) %>% 
  data.table::melt(variable.name="subplot", id.vars=ncol(.), measure.vars=1:4) %>% 
  na.omit()

vegplot = merge(vegplot, subs, all=T)
vegplot$subplotID = vegplot$value

# sort
vegplot = vegplot[order(vegplot$plotID),]

# check a tautology for some reason. The universe is sane, I guess.
vegplot$plotID==vegplot$plotID



# Find the most common tree species by plot...
plotTypes = veg %>% group_by(plotID, taxonID) %>% 
  summarize(count=length(taxonID)) %>% 
  left_join(veg %>% group_by(plotID) %>% 
              summarize(totplot=length(plotID))) %>% 
  mutate(plotPct=round((count/totplot)*100))

plotTypes$forType_plot = "NA"
for(i in unique(plotTypes$plotID)) {
  temp_df = subset(plotTypes, plotID== i)
  temp_df = temp_df[order(temp_df$plotPct, decreasing=T),]
  if(any(temp_df$plotPct>=80)){
    plotTypes$forType_plot[plotTypes$plotID==i] = temp_df$taxonID[temp_df$plotPct>=80]
  }else if(sum(temp_df$plotPct[1:2])>=80) {
    plotTypes$forType_plot[plotTypes$plotID==i] = paste0(temp_df$taxonID[1],temp_df$taxonID[2])
  }else{
    plotTypes$forType_plot[plotTypes$plotID==i] = "MIXCON"
  }
}

plotTypes = plotTypes %>% group_by(plotID) %>% summarise(forType_plot=max(forType_plot))
table(plotTypes$forType_plot)


# ...and subplot
subplotTypes = veg %>% group_by(plotID, subplotID, taxonID) %>% 
  summarize(count=length(taxonID)) %>% 
  left_join(veg %>% group_by(plotID, subplotID) %>% 
              summarize(totplot=length(subplotID))) %>% 
  mutate(subplotPct=round((count/totplot)*100), uniqueID=paste0(plotID,"_",subplotID))


subplotTypes$forType_subplot = "NA"
for(i in unique(subplotTypes$uniqueID)) {
  temp_df = subset(subplotTypes, uniqueID== i)
  temp_df = temp_df[order(temp_df$subplotPct, decreasing=T),]
  if(any(temp_df$subplotPct>=80)){
    subplotTypes$forType_subplot[subplotTypes$uniqueID==i] = temp_df$taxonID[temp_df$subplotPct>=80]
  }else if(sum(temp_df$subplotPct[1:2])>=80) {
    subplotTypes$forType_subplot[subplotTypes$uniqueID==i] = paste0(temp_df$taxonID[1],temp_df$taxonID[2])
  }else{
    subplotTypes$forType_subplot[subplotTypes$uniqueID==i] = "MIXCON"
  }
}

subplotTypes = subplotTypes %>% group_by(plotID, subplotID) %>% summarise(forType_subplot=max(forType_subplot))
table(subplotTypes$forType_subplot)



veg = left_join(veg, plotTypes, by=join_by(plotID))
veg = left_join(veg, subplotTypes, by=join_by(plotID, subplotID))
vegplot = left_join(vegplot, plotTypes, by=join_by(plotID))
vegplot = left_join(vegplot, subplotTypes, by=join_by(plotID, subplotID))


# Mean crown diameter = maxCrownDiameter+ninetyCrownDiameter/2
veg$meanCrownDiameter = (veg$maxCrownDiameter+veg$ninetyCrownDiameter)/2

# find a good predictive model for ws for each forest type
mdl_df = subset(veg, !is.na(ninetyCrownDiameter))
plot(ninetyCrownDiameter~height, mdl_df)

# remove outliers
# mdl_df = subset(mdl_df, maxCrownDiameter<15)
mdl_df = subset(mdl_df, height<35)

# Separate tower plots from distributed. They look like different forests, different elevation, slope, and aspect.
mdl_df = left_join(mdl_df, vegplot[!duplicated(vegplot$plotID),c("plotID","plotType")], by="plotID")
tower_mdl_df = subset(mdl_df, plotType=="tower")
dist_mdl_df = subset(mdl_df, plotType=="distributed")


# Functions
run_models = function(data, response) {
  output = list()
  null_mod = lm(data[[response]]~1, data)
  lin_mod = lm(data[[response]]~height, data)
  exp_mod = tryCatch({nls(data[[response]] ~ a*b^height, data=data, start=list(a=1,b=1))},
                     error = function(e) {NA})
  log_mod = tryCatch({nls(data[[response]] ~ exp(-(a*(1/height))+(height^b)), data=data, start=list(a=1,b=1))},
                     error = function(e) {NA})
  mono_mod = tryCatch({nls(data[[response]] ~ a*(1-exp(-b*height)), data=data, start=list(a=5,b=1))},
                      error = function(e) {NA})
  logit_mod = tryCatch({nls(data[[response]]~c/(1+a*exp(-b*height)), data=data, start=list(a=5, b=0.4, c=4))},
                       error = function(e) {NA})
  mm_mod = tryCatch({nls(data[[response]]~a*height/(b+height), data=data, start=list(a=1, b=2))},
                    error = function(e) {NA})
  h3_mod = tryCatch({nls(data[[response]]~a*height^2/(b^2+height^2), data=data, start=list(a=1,b=1))},
                    error = function(e) {NA})
  
  mods = list(null_mod, lin_mod,exp_mod,log_mod,mono_mod,logit_mod,mm_mod,h3_mod)
  names(mods) = c("null","linear","exponential","logarithmic","monomolecular","logistic",
                  "MichaelisMenten","Holling3")
  return(mods)
}

summarize_models = function(model_list) {
  row_tmp = list()
  for(i in 1:length(model_list)) {
    if(all(is.na(model_list[[i]]))) next
    summary_tmp = summary(model_list[[i]]) %>% 
      coef() %>% as.data.frame() %>% mutate(md_name=names(model_list)[i],
                                            AIC_tmp = AIC(model_list[[i]]))
    r2_tmp = nagelkerke(model_list[[i]], model_list[[1]])
    summary_tmp$R2 = r2_tmp$Pseudo.R.squared[3,] # Nagelkerke
    
    row_tmp[[length(row_tmp)+1]] = summary_tmp
  }
  row_tmp=do.call(rbind, row_tmp)
  row_tmp$deltaAIC = row_tmp$AIC - min(row_tmp$AIC)
  return(row_tmp)
}   

model_tables = function(model_data, best=T, response) {
  type_tmp = na.omit(c(unique(model_data$forType_plot), NULL))

  if(best) {
    func = function(x) {subset(model_data, forType_plot==x) %>% 
                            run_models(response) %>% 
                            summarize_models() %>% 
                            subset(deltaAIC==0) %>% 
                            mutate(forType=x, Response=response)
    }
  }else{
    func = function(x) {subset(model_data, forType_plot==x) %>% 
                            run_models(response) %>% 
                            summarize_models() %>%  
                            mutate(forType=x, Response=response)
    }
  }
    
  best_mod_tmp = lapply(type_tmp, func)
  
  allFor_mod_tmp = run_models(model_data, response) %>% 
    summarize_models() %>% 
    mutate(forType="ALLTYPES", Response=response)
  
  if(best) allFor_mod_tmp = subset(allFor_mod_tmp, deltaAIC==0) 
  
  best_mod_tmp[[length(best_mod_tmp)+1]] = allFor_mod_tmp  
  
  best_mod_tmp = do.call(rbind, best_mod_tmp)
}

plotter = function(mdl_data, ws, height_list, response, ...) {
  # args = list(...)
  # argNames = names(list(...))
  
  if(length(unique(mdl_data$forType_plot))>1){
    plot_name = "All Forest Types"
  }else{
    plot_name = mdl_data$forType_plot[1]
  }
  
  plot(eval(parse(text=response))~height, data=mdl_data, cex=2, cex.main=3, 
       ylab=response, cex.lab=3, cex.axis=2,
       main=plot_name, ...)
  diams = ws(height_list)
  lines(diams~height_list, col="green", lwd=2)
} 


# Model summary data
tower_mdl_summary_long = model_tables(tower_mdl_df, response="ninetyCrownDiameter",best=F)
tower_mdl_summary = model_tables(tower_mdl_df, "ninetyCrownDiameter", best=T)
dist_mdl_summary_long = model_tables(dist_mdl_df, response="ninetyCrownDiameter",best=F)
dist_mdl_summary = model_tables(dist_mdl_df, "ninetyCrownDiameter", best=T)

# Let's just do the forest types that had a bunch of trees (>30?). Also, combine
# pienablal and ablalpien into "sprfir."

table(tower_mdl_df$forType_plot)
tower_ablapien_df = subset(tower_mdl_df, forType_plot=="ABLALPIEN")
tower_pienabla_df = subset(tower_mdl_df, forType_plot=="PIENABLAL")
tower_sprfir_df = subset(tower_mdl_df, forType_plot=="PIENABLAL" | forType_plot=="ABLALPIEN")
  tower_sprfir_df$forType_plot = "sprfir"


table(dist_mdl_df$forType_plot)
dist_pico_df = subset(dist_mdl_df, forType_plot=="PICOL")
dist_abla_df = subset(dist_mdl_df, forType_plot=="ABLAL")
dist_ablapien_df = subset(dist_mdl_df, forType_plot=="ABLALPIEN")
dist_pienabla_df = subset(dist_mdl_df, forType_plot=="PIENABLAL")
dist_mixcon_df = subset(dist_mdl_df, forType_plot=="MIXCON")
dist_sprfir_df = subset(dist_mdl_df, forType_plot=="PIENABLAL" | forType_plot=="ABLALPIEN")
  dist_sprfir_df$forType_plot = "sprfir"

  
# Run models for mixed forest types
tower_sprfir_summary_long = model_tables(tower_sprfir_df, "ninetyCrownDiameter", best=F)
dist_sprfir_summary_long = model_tables(dist_sprfir_df, "ninetyCrownDiameter", best=F)


# ws functions -----------------------------------------------------------------------
# tower
tower_ablapien_ws = function(x) {dplyr::case_when(is.na(x) ~ 0.001,
                                       TRUE ~ 0.3857049*1.3255355^x)}
tower_pienabla_ws = function(x) {dplyr::case_when(is.na(x) ~ 0.001, x <= 0 ~ 0.001,
                                       TRUE ~ 4.1031350*x/(7.3406803+x))}
tower_sprfir_ws = function(x) {dplyr::case_when(is.na(x) ~ 0.001,
                                       TRUE ~ 3.68152014/(1+8.98927765*exp(-0.38743418*x)))}
tower_alltyp_ws = function(x) {dplyr::case_when(is.na(x) ~ 0.001,
                                       TRUE ~ 3.7018386/(1+8.5178998*exp(-0.3751006*x)))}

# distributed
dist_pico_ws = function(x) {dplyr::case_when(is.na(x) ~ 0.001, x <= 0 ~ 0.001,
                                       TRUE ~ 2.70109769*x^2/(3.92763459^2+x^2))}
dist_abla_ws = function(x) {dplyr::case_when(is.na(x) ~ 0.001, x <= 0 ~ 0.001,
                                       TRUE ~ 4.1364199*(1-exp(-0.1256920*x)))}
dist_ablapien_ws = function(x) {dplyr::case_when(is.na(x) ~ 0.001, x <= 0 ~ 0.001,
                                       TRUE ~ 4.22466232*(1-exp(-0.09069298*x)))}
dist_pienabla_ws = function(x) {dplyr::case_when(is.na(x) ~ 0.001,
                                       TRUE ~ 1.07996268*1.08388971^x)}
dist_sprfir_ws = function(x) {dplyr::case_when(is.na(x) ~ 0.001,
                                       TRUE ~ 3.72351435/(1+4.22109271*exp(-0.21718661*x)))}
dist_mixcon_ws = function(x) {dplyr::case_when(is.na(x) ~ 0.001, x <= 0 ~ 0.001,
                                       TRUE ~ 2.69508856/(1+18.21764821*exp(-0.81247595*x)))}
dist_alltyp_ws = function(x) {dplyr::case_when(is.na(x) ~ 0.001, x <= 0 ~ 0.001,
                                       TRUE ~ 3.33336298*(1-exp(-0.15498102*x)))}

# tower curves
plot(ninetyCrownDiameter~height, tower_mdl_df, pch=20)
lines(dist_ablapien_ws(heights)~heights, col="yellow", lwd=2)
lines(dist_pienabla_ws(heights)~heights, col="green", lwd=2)
lines(dist_sprfir_ws(heights)~heights, col="blue", lwd=2)
lines(dist_alltyp_ws(heights)~heights, col="black", lwd=2)

# distributed curves
plot(ninetyCrownDiameter~height, dist_mdl_df, pch=20)
lines(dist_pico_ws(heights)~heights, col="red", lwd=2)
lines(dist_abla_ws(heights)~heights, col="orange", lwd=2)
lines(dist_ablapien_ws(heights)~heights, col="yellow", lwd=2)
lines(dist_pienabla_ws(heights)~heights, col="green", lwd=2)
lines(dist_sprfir_ws(heights)~heights, col="blue", lwd=2)
lines(dist_mixcon_ws(heights)~heights, col="purple", lwd=2)
lines(dist_alltyp_ws(heights)~heights, col="black", lwd=2)

# piflpsme_ws = function(x) {dplyr::case_when(is.na(x) ~ 0.001, x <= 0 ~ 0.001,
#                                        TRUE ~ 3.740385)}



#--------------------------------------------> Save outputs
# dir.create("treeSearch/outputs", recursive=T)
dir.create("treeSearch/outputs/graphs")
dir.create("treeSearch/outputs/tables")
dir.create("treeSearch/outputs/ws_functions")

# tables
write.csv(tower_mdl_summary_long, "treeSearch/outputs/tables/TowerModelSummary_NIWO.csv")
write.csv(tower_mdl_summary, "treeSearch/outputs/tables/TowerBestModelsSummary_NIWO.csv")
write.csv(dist_mdl_summary_long, "treeSearch/outputs/tables/DistModelSummary_NIWO.csv")
write.csv(dist_mdl_summary, "treeSearch/outputs/tables/distBestModelsSummary_NIWO.csv")

# ws functions
save(tower_ablapien_ws, file="treeSearch/outputs/ws_functions/tower_ablapien_NIWO_ws.RData")
save(tower_pienabla_ws, file="treeSearch/outputs/ws_functions/tower_pienabla_NIWO_ws.RData")
save(tower_sprfir_ws, file="treeSearch/outputs/ws_functions/tower_sprfir_NIWO_ws.RData")
save(tower_alltyp_ws, file="treeSearch/outputs/ws_functions/tower_alltyp_NIWO_ws.RData")
save(dist_mixcon_ws, file="treeSearch/outputs/ws_functions/dist_mixcon_NIWO_ws.RData")
save(dist_pico_ws, file="treeSearch/outputs/ws_functions/dist_pico_NIWO_ws.RData")
save(dist_abla_ws, file="treeSearch/outputs/ws_functions/dist_abla_NIWO_ws.RData")
save(dist_ablapien_ws, file="treeSearch/outputs/ws_functions/dist_ablapien_NIWO_ws.RData")
save(dist_pienabla_ws, file="treeSearch/outputs/ws_functions/dist_pienabla_NIWO_ws.RData")
save(dist_sprfir_ws, file="treeSearch/outputs/ws_functions/dist_sprfir_NIWO_ws.RData")
save(dist_alltyp_ws, file="treeSearch/outputs/ws_functions/dist_alltyp_NIWO_ws.RData")


# tower graphs
png("treeSearch/outputs/graphs/tower_bestFits_NIWO.png", height=2500, width=2500)
par(mfrow=c(2,2))
  # ablapien
  par(mar=c(5.1,5.1,4.1,6.1))
  plotter(tower_ablapien_df, ws=tower_ablapien_ws, 0:35, "ninetyCrownDiameter", xlim=c(0,20), ylim=c(0,6))
  text(x=6.8,y=5.8, label="y=0.386*1.326^x",
       font=2, col="red", cex=3)
  # pienabla
  par(mar=c(5.1,5.1,4.1,6.1))
  plotter(tower_pienabla_df, ws=tower_pienabla_ws, 0:35, "ninetyCrownDiameter", xlim=c(0,20), ylim=c(0,6))
  text(x=6.8,y=5.8, label="y=4.103*x/(7.341+x)",
       font=2, col="red", cex=3)
  # spruce-fir
  par(mar=c(5.1,5.1,4.1,6.1))
  plotter(tower_sprfir_df, ws=tower_sprfir_ws, 0:35, "ninetyCrownDiameter", xlim=c(0,20), ylim=c(0,6))
  text(x=6.8,y=5.8, label="y=3.682/(1+8.989*exp(-0.387*x))",
       font=2, col="red", cex=3)
  # All Forest types
  par(mar=c(5.1,5.1,4.1,6.1))
  plot(ninetyCrownDiameter~height, data=tower_mdl_df, cex=2, cex.main=3, 
       cex.lab=3, cex.axis=2,
       main="All Forest Types", xlim=c(0,20), ylim=c(0,6))
  heights =0:35
  diams = tower_alltyp_ws(heights)
  lines(diams~heights, col="green", lwd=2)
  text(x=6.8,y=5.8, label="y=3.702/(1+8.518*exp(-0.375*x))",
       font=2, col="red", cex=3)
dev.off()


# distributed plots graphs
png("treeSearch/outputs/graphs/dist_bestFits_PicoPienAbla_NIWO.png", height=2500, width=2500)
par(mfrow=c(2,2))
  # Pinus contorta
  par(mar=c(5.1,5.1,4.1,6.1))
  plotter(dist_pico_df, ws=dist_pico_ws, 0:35, "ninetyCrownDiameter", xlim=c(0,20), ylim=c(0,7))
  text(x=6.8,y=6.8, label="y=2.701*x^2/(3.928^2+x^2)",
       font=2, col="red", cex=3)
  # abla
  par(mar=c(5.1,5.1,4.1,6.1))
  plotter(dist_abla_df, ws=dist_abla_ws, 0:35, "ninetyCrownDiameter", xlim=c(0,20), ylim=c(0,7))
  text(x=6.8,y=6.8, label="y=4.136*(1-exp(-0.126*x))",
       font=2, col="red", cex=3)
  # ablapien
  par(mar=c(5.1,5.1,4.1,6.1))
  plotter(dist_ablapien_df, ws=dist_ablapien_ws, 0:35, "ninetyCrownDiameter", xlim=c(0,20), ylim=c(0,7))
  text(6.8,y=6.8, label="y=4.225*(1-exp(-0.091*x))",
       font=2, col="red", cex=3)
  # pienabla
  par(mar=c(5.1,5.1,4.1,6.1))
  plotter(dist_pienabla_df, ws=dist_pienabla_ws, 0:35, "ninetyCrownDiameter", xlim=c(0,20), ylim=c(0,7))
  text(6.8,y=6.8, label="y=1.080*1.084^x",
       font=2, col="red", cex=3)
dev.off()


# distributed spruce-fir and all forest types
png("treeSearch/outputs/graphs/dist_bestFits_sprfirAllTyp_NIWO.png", height=1250, width=2500)
par(mfrow=c(1,2))
# sprfir
par(mar=c(5.1,5.1,4.1,6.1))
plotter(dist_sprfir_df, ws=dist_sprfir_ws, 0:35, "ninetyCrownDiameter", xlim=c(0,20), ylim=c(0,7))
text(x=6.8,y=6.8, label="y=3.724/(1+4.221*exp(-0.217*x))",
     font=2, col="red", cex=3)
# alltyp
par(mar=c(5.1,5.1,4.1,6.1))
plotter(dist_mdl_df, ws=dist_alltyp_ws, 0:35, "ninetyCrownDiameter", xlim=c(0,20), ylim=c(0,7))
text(x=6.8,y=6.8, label="y=4.136*(1-exp(-0.126*x))",
     font=2, col="red", cex=3)
dev.off()



# Load NEON lidar data for cloud2trees dbh and basal area comparison with NEON
ctg = catalog("../../Data/NIWO_lidar/")

# create index files. Speeds things up.
lidR:::catalog_laxindex(ctg)


# Canopy height model/treetop location plots for visual assessment.
# SF object to make it easier to plot with the chm
vegplot_sf = data.frame(plotID = vegplot$plotID,
                   subplotID = vegplot$subplotID,
                   plotType = vegplot$plotType,
                   X=vegplot$easting,
                   Y=vegplot$northing) %>% 
  st_as_sf(coords=c("X","Y"))
st_crs(vegplot_sf) = "EPSG:32616"


veg_sf = subset(veg, !is.na(adjEasting)) %>% 
  data.frame(plotID = .$plotID,
             subplotID = .$subplotID,
             X=.$adjEasting,
             Y=.$adjNorthing) %>% 
  st_as_sf(coords=c("X","Y"))
st_crs(veg_sf) = "EPSG:32616"


# Pinus contorta areas
vegplot_sf$plotID[vegplot$forType=="PICOL"] # NIWO_001 and NIWO_002
crds001 = st_coordinates(vegplot_sf[vegplot_sf$plotID=="NIWO_001",][1,])
crds002 = st_coordinates(vegplot_sf[vegplot_sf$plotID=="NIWO_002",][1,])

# Get the point cloud for NIWO_001 and normalize.
las001 = readLAS(ctg, filter = paste0("-keep_xy ", crds001[1]-20," ",crds001[2]-20," ",
                                      crds001[1]+20," ",crds001[2]+20))
las001 = las001 %>% classify_noise(sor()) %>% filter_poi(Classification!=LASNOISE) %>% 
  normalize_height(knnidw())
chm001 = rasterize_canopy(las001, res=0.5)

ttops001 = locate_trees(las001, lmf(pico_ws))

# Plot
png("treeSearch/outputs/graphs/chmTtops_NIWO_001.png", height=1000, width=1000)
  zoom(chm001, e=ext(c(crds001[1]-20,crds001[1]+20,crds001[2]-20,crds001[2]+20)),
       pax=list(cex.axis=2), plg=list(cex=2), main="Pinus contorta NIWO_001", cex.main=2)
  plot(ttops001['Z'], pch=20, add=T, col="red", cex=2)
dev.off()
  
# NIWO_002
las002 = readLAS(ctg, filter = paste0("-keep_xy ", crds002[1]-20," ",crds002[2]-20," ",
                                      crds002[1]+20," ",crds002[2]+20))
las002 = las002 %>% classify_noise(sor()) %>% filter_poi(Classification!=LASNOISE) %>% 
  normalize_height(knnidw())
chm002 = rasterize_canopy(las002, res=0.5)

ttops002 = locate_trees(las002, lmf(pico_ws))

# plot
png("treeSearch/outputs/graphs/chmTtops_NIWO_002.png", height=1000, width=1000)
  zoom(chm002, e=ext(c(crds002[1]-20,crds002[1]+20,crds002[2]-20,crds002[2]+20)),
       pax=list(cex.axis=2), plg=list(cex=2), main="Pinus contorta NIWO_002", cex.main=2)
  plot(ttops002['Z'], pch=20, add=T, col="red", cex=2)
dev.off()


# spruce_fir areas
vegplot_sf$plotID[vegplot$forType=="ABLALPIEN"] 
crds005 = st_coordinates(vegplot_sf[vegplot_sf$plotID=="NIWO_005",][1,])
crds057 = st_coordinates(vegplot_sf[vegplot_sf$plotID=="NIWO_057",][1,])

# NIWO_005
las005 = readLAS(ctg, filter = paste0("-keep_xy ", crds005[1]-20," ",crds005[2]-20," ",
                                      crds005[1]+20," ",crds005[2]+20))
las005 = las005 %>% classify_noise(sor()) %>% filter_poi(Classification!=LASNOISE) %>% 
  normalize_height(knnidw())
chm005 = rasterize_canopy(las005, res=0.5)
ttops005 = locate_trees(las005, lmf(sprfir_ws))

# plot
png("treeSearch/outputs/graphs/chmTtops_NIWO_005.png", height=1000, width=1000)
  zoom(chm005, e=ext(c(crds005[1]-20,crds005[1]+20,crds005[2]-20,crds005[2]+20)),
       pax=list(cex.axis=2), plg=list(cex=2), main="Spruce/Fir NIWO_005", cex.main=2)
  plot(ttops005['Z'], pch=20, add=T, col="red", cex=2)
dev.off()


# NIWO_057
las057 = readLAS(ctg, filter = paste0("-keep_xy ", crds057[1]-20," ",crds057[2]-20," ",
                                      crds057[1]+20," ",crds057[2]+20))
las057 = las057 %>% classify_noise(sor()) %>% filter_poi(Classification!=LASNOISE) %>% 
  normalize_height(knnidw())
chm057 = rasterize_canopy(las057, res=0.5)
ttops057 = locate_trees(las057, lmf(sprfir_ws))

# plot
png("treeSearch/outputs/graphs/chmTtops_NIWO_057.png", height=1000, width=1000)
  zoom(chm057, e=ext(c(crds057[1]-20,crds057[1]+20,crds057[2]-20,crds057[2]+20)),
       pax=list(cex.axis=2), plg=list(cex=2), main="Spruce/Fir NIWO_057", cex.main=2)
  plot(ttops057['Z'], pch=20, add=T, col="red", cex=2)
dev.off()


# Picea Engelmannii areas
vegplot_sf$plotID[vegplot$forType=="PIEN"] 
crds004 = st_coordinates(vegplot_sf[vegplot_sf$plotID=="NIWO_004",][1,])
crds067 = st_coordinates(vegplot_sf[vegplot_sf$plotID=="NIWO_067",][1,])

# Plot NIWO_004
las004 = readLAS(ctg, filter = paste0("-keep_xy ", crds004[1]-20," ",crds004[2]-20," ",
                                      crds004[1]+20," ",crds004[2]+20))
las004 = las004 %>% classify_noise(sor()) %>% filter_poi(Classification!=LASNOISE) %>% 
  normalize_height(knnidw())
chm004 = rasterize_canopy(las004, res=0.5)
ttops004 = locate_trees(las004, lmf(pien_ws))

# graph
png("treeSearch/outputs/graphs/chmTtops_NIWO_004.png", height=1000, width=1000)
zoom(chm004, e=ext(c(crds004[1]-20,crds004[1]+20,crds004[2]-20,crds004[2]+20)),
     pax=list(cex.axis=2), plg=list(cex=2), main="Picea Engelmannii NIWO_004", cex.main=2)
plot(ttops004['Z'], pch=20, add=T, col="red", cex=2)
dev.off()


# Plot NIWO_067
las067 = readLAS(ctg, filter = paste0("-keep_xy ", crds067[1]-20," ",crds067[2]-20," ",
                                      crds067[1]+20," ",crds067[2]+20))
las067 = las067 %>% classify_noise(sor()) %>% filter_poi(Classification!=LASNOISE) %>% 
  normalize_height(knnidw())
chm067 = rasterize_canopy(las067, res=0.5)
ttops067 = locate_trees(las067, lmf(pien_ws))

# plot
png("treeSearch/outputs/graphs/chmTtops_NIWO_067.png", height=1000, width=1000)
zoom(chm067, e=ext(c(crds067[1]-20,crds067[1]+20,crds067[2]-20,crds067[2]+20)),
     pax=list(cex.axis=2), plg=list(cex=2), main="Picea Engelmannii NIWO_067", cex.main=2)
plot(ttops067['Z'], pch=20, add=T, col="red", cex=2)
dev.off()




# cloud2trees comparisons
dir.create("treeSearch/cloud2trees/NIWO_001/", recursive = T)
dir.create("treeSearch/cloud2trees/NIWO_002/")
dir.create("treeSearch/cloud2trees/NIWO_004/")
dir.create("treeSearch/cloud2trees/NIWO_005/")
dir.create("treeSearch/cloud2trees/NIWO_057/")
dir.create("treeSearch/cloud2trees/NIWO_067/")
writeLAS(las001, "treeSearch/cloud2trees/NIWO_001/las001.las")
writeLAS(las002, "treeSearch/cloud2trees/NIWO_002/las002.las")
writeLAS(las004, "treeSearch/cloud2trees/NIWO_004/las004.las")
writeLAS(las005, "treeSearch/cloud2trees/NIWO_005/las005.las")
writeLAS(las057, "treeSearch/cloud2trees/NIWO_057/las057.las")
writeLAS(las067, "treeSearch/cloud2trees/NIWO_067/las067.las")

# run cloud2trees
output = cloud2trees::cloud2trees("treeSearch/cloud2trees/NIWO_001/",
                                  input_las_dir = "treeSearch/cloud2trees/NIWO_001/",
                                  chm_res_m = 0.5,
                                  ws = pico_ws,
                                  estimate_tree_dbh = T)

output = cloud2trees::cloud2trees("treeSearch/cloud2trees/NIWO_002/",
                                  input_las_dir = "treeSearch/cloud2trees/NIWO_002/",
                                  chm_res_m = 0.5,
                                  ws = pico_ws,
                                  estimate_tree_dbh = T)

output = cloud2trees::cloud2trees("treeSearch/cloud2trees/NIWO_004/",
                                  input_las_dir = "treeSearch/cloud2trees/NIWO_004/",
                                  chm_res_m = 0.5,
                                  ws = pien_ws,
                                  estimate_tree_dbh = T)

output = cloud2trees::cloud2trees("treeSearch/cloud2trees/NIWO_067/",
                                  input_las_dir = "treeSearch/cloud2trees/NIWO_067/",
                                  chm_res_m = 0.5,
                                  ws = pien_ws,
                                  estimate_tree_dbh = T)

output = cloud2trees::cloud2trees("treeSearch/cloud2trees/NIWO_005/",
                                  input_las_dir = "treeSearch/cloud2trees/NIWO_005/",
                                  chm_res_m = 0.5,
                                  ws = sprfir_ws,
                                  estimate_tree_dbh = T)

output = cloud2trees::cloud2trees("treeSearch/cloud2trees/NIWO_057/",
                                  input_las_dir = "treeSearch/cloud2trees/NIWO_057/",
                                  chm_res_m = 0.5,
                                  ws = sprfir_ws,
                                  estimate_tree_dbh = T)

# NEON
sum((thing$stemDiameter/2)^2*pi)/1600
trees001 = st_read("treeSearch/cloud2trees/NIWO_001/point_cloud_processing_delivery/final_detected_tree_tops.gpkg")
sum((trees001$dbh_m))/6400*10000


# Clip out point clouds for each plot
pcs = list()
plot_centers = list() 
distSum = list()
trees_mapped =list()
for(i in 1:nrow(vegplot)) {
  # clip with buffer so TPI covers the whole plot
  # different size tower and distributed plots. Add 1 m for buffer.
  
  vegx = vegplot$easting[i]
  vegy = vegplot$northing[i]
  
  xleft = vegx-10
  xright = vegx+10
  ybottom = vegy-10
  ytop = vegy+10
  
  # Plot sizes:
  # Tower: 40x40 (1600 m^2)
  # Distributed: 20x20 (400 m^2)
  # for trees, smallest plot size 20x20? (I think so).
  # For distributed take the whole plot, for tower, find the subplots
  
  # get plot/subplot size(s) and locations
  subplotID = vegplot$value[i]
  if(vegplot$value[i] %in% c("21","23","39","41")) {
    size=20
    if(subplotID==21 | subplotID==31){
      xleft=vegx-size;ybottom=vegy-size;xright=vegx;ytop=vegy
    }else if(subplotID==39 | subplotID==40) {
      xleft=vegx-size;ybottom=vegy;xright=vegx;ytop=vegy+size
    }else if(subplotID==41) {
      xleft=vegx;ybottom=vegy;xright=vegx+size;ytop=vegy+size
    }else if(subplotID==23 | subplotID==32) {
      xleft=vegx;ybottom=vegy-size;xright=vegx+size;ytop=vegy
    }
  }
  
  
  # Find plot center and add 0 for height above ground
  plotCenter_tmp = c(vegx,vegy, 0)
  
  # Find tree locations
  if(is.na(vegplot$value[i])) {
    treesx = na.omit(veg$adjEasting[veg$plotID == vegplot$plotID[i]])
    treesy = na.omit(veg$adjNorthing[veg$plotID == vegplot$plotID[i]])
    treesz = veg$height[veg$plotID == vegplot$plotID[i]]
    if(!is.null(attr(treesx, "na.action"))) {treesz = treesz[-c(attr(treesx, "na.action"))]} 
  }else{
    treesx = na.omit(veg$adjEasting[veg$plotID == vegplot$plotID[i] & 
                                      veg$subplotID==vegplot$value[i]]) 
    treesy = na.omit(veg$adjNorthing[veg$plotID == vegplot$plotID[i] & 
                                       veg$subplotID==vegplot$value[i]])
    treesz = veg$height[veg$plotID == vegplot$plotID[i] & 
                          veg$subplotID==vegplot$value[i]]
    
    if(!is.null(attr(treesx, "na.action"))) {treesz = treesz[-c(attr(treesx, "na.action"))]} 
  }
  
  trees_tmp = cbind(treesx,treesy,treesz)
  if(nrow(trees_tmp)==0) {
    pcs[[length(pcs)+1]] = NA
    plot_centers[[length(plot_centers)+1]] = NA
    distSum[[length(distSum)+1]] = NA
    trees_mapped[[length(trees_mapped)+1]] = NA
    next
  }
  
  distsSum_tmp = sum(proxy::dist(matrix(plotCenter_tmp, ncol=3,byrow=T), trees_tmp))
  
  trees_mapped[[length(trees_mapped)+1]] = trees_tmp
  
  # cut out the point cloud
  tmplas = clip_rectangle(ctg, xleft, ybottom, xright, ytop) %>% 
    classify_noise(sor()) %>% 
    filter_poi(Classification!=LASNOISE) %>% 
    normalize_height(knnidw())
    
    
  plot_centers[[length(plot_centers)+1]] = plotCenter_tmp
  distSum[[length(distSum)+1]] = distsSum_tmp
  pcs[[length(pcs)+1]] = tmplas
}
plot_centers = do.call(rbind, plot_centers)
distSum = do.call(rbind, distSum)



# Get the boundaries so we can download the right lidar tiles.
summary(veg$adjNorthing)
summary(veg$adjEasting)


plot_example = unique(veg$plotID)[2]
# Map with uncertainty
symbols(veg$adjEasting[which(veg$plotID==plot_example)], 
        veg$adjNorthing[which(veg$plotID==plot_example)], 
        circles=veg$stemDiameter[which(veg$plotID==plot_example)]/100/2, 
        inches=F, xlab="Easting", ylab="Northing")

symbols(veg$adjEasting[which(veg$plotID==plot_example)], 
        veg$adjNorthing[which(veg$plotID==plot_example)], 
        circles=veg$adjCoordinateUncertainty[which(veg$plotID==plot_example)], 
        inches=F, add=T, fg="lightblue")


# trees per plot
tpp = veg %>% group_by(plotID, subplotID) %>% count() %>% as.data.frame() %>% 
  na.omit()








# plots to see if it looks like I did it right
# add a column for the ws model for each forest type
unique(vegplot$forType_plot)
vegplot$forType_model = vegplot$forType_subplot
vegplot$forType_model[vegplot$plotType=="distributed"] = vegplot$forType_plot[vegplot$plotType=="distributed"]
vegplot$model = "pico_ws"
vegplot$model[vegplot$forType_model=="PIFL2"] = "pifl_ws"
vegplot$model[vegplot$forType_model=="ABLAL"] = "abla_ws"
vegplot$model[vegplot$forType_model=="ABLALPIEN" | vegplot$forType_model=="PIENABLAL"] = "sprfir_ws"
vegplot$model[vegplot$forType_model=="MIXCON"] = "alltyp_ws"


comparison_table = function(plt) {
  
  dir.create(paste0("treeSearch/cloud2trees/",vegplot$plotID[plt],"/",vegplot$value[plt],"/"), recursive=T)
  writeLAS(pcs[[plt]], paste0("treeSearch/cloud2trees/",vegplot$plotID[plt],"/",vegplot$value[plt],"/",vegplot$value[plt],".las"))
  
  if(vegplot$plotType[plt]=="tower") {
    ws_tmp = paste0("tower_",vegplot$model[plt])
  }else{
    ws_tmp = paste0("dist_",vegplot$model[plt])
  }
  
  tmp = cloud2trees(paste0("treeSearch/cloud2trees/",vegplot$plotID[plt],"/",vegplot$value[plt],"/"),
                    input_las_dir = paste0("treeSearch/cloud2trees/",vegplot$plotID[plt],"/",vegplot$value[plt],"/"),
                    chm_res_m = 0.5,
                    ws = eval(parse(text=ws_tmp)))
  
  out = trees_dbh(tmp$treetops)
  
  veg1 = subset(veg_sf, plotID==vegplot$plotID[plt] &
                  (subplotID == vegplot$value[plt] | is.na(vegplot$value[plt])))
  neonBA = sum((veg1$stemDiameter/2)^2*pi)/400
  modelBA = sum((out$dbh_cm/2)^2*pi)/400
  
  return = data.frame(plotID=vegplot$plotID[plt],
                      plotType =vegplot$plotType[plt],
                     forType = vegplot$forType_plot[plt],
                     model = ws_tmp,
                     neonBA = neonBA,
                     modelBA = modelBA)
}

comps = list()
for(i in 1:length(pcs)) {
  plt=i
  comps[[length(comps)+1]] = try(comparison_table(plt))
}

splitem = function(x, splitter) {if(class(x)==splitter){return(x)}}
comps_tmp = lapply(comps, function(x) splitem(x, "data.frame"))
comps_tmp = do.call(rbind, comps_tmp)


plot(modelBA~neonBA, subset(comps_tmp, model=="dist_pico_ws"))
t = lm(modelBA~neonBA, subset(comps_tmp, model=="dist_pico_ws"))
summary(t)

plot(modelBA~neonBA, subset(comps_tmp, forType=="PIENABLAL" & plotType=="tower"))
t = lm(modelBA~neonBA, subset(comps_tmp, forType=="PIENABLAL" & plotType=="tower"))
summary(t)



  dir.create(paste0("treeSearch/cloud2trees/",vegplot$plotID[plt],"/",vegplot$value[plt],"/"), recursive=T)
  writeLAS(pcs[[plt]], paste0("treeSearch/cloud2trees/",vegplot$plotID[plt],"/",vegplot$value[plt],"/",vegplot$value[plt],".las"))
  
  tmp = cloud2trees(paste0("treeSearch/cloud2trees/",vegplot$plotID[plt],"/",vegplot$value[plt],"/"),
                    input_las_dir = paste0("treeSearch/cloud2trees/",vegplot$plotID[plt],"/",vegplot$value[plt],"/"),
                    chm_res_m = 0.5,
                    ws = eval(parse(text=vegplot$model[plt])))  
}

plt=6
vegplot$forType[plt]
chm = rasterize_canopy(pcs[[plt]], res=0.5)
plot(chm)
plot(veg_sf[veg_sf$plotID == vegplot$plotID[plt] & (veg_sf$subplotID == vegplot$value[plt] | is.na(vegplot$value[plt])),]['height'], 
     pch=20, col="red", axes=T, add=T)
ttops = locate_trees(pcs[[plt]], lmf(sprfir_ws))
plot(ttops['treeID'], pch=20, col="blue", add=T)







# Model dataframe
model_df = data.frame(plot = vegplot$plotID, subplot=vegplot$value)
model_df$distSum = distSum
model_df = cbind(model_df, plot_centers)
names(model_df)[4:6] = c("plCentX", "plCentY","plCentZ")
model_df$pcs = unlist(pcs)
model_df$trees_mapped = trees_mapped
model_df$errors = 0
model_df = model_df[-26,] # was NA



pipa2_df = model_df[model_df$plot %in% plotTypes$plotID[plotTypes$taxonID=="PIPA2"],]
ablal_df = model_df[model_df$plot %in% plotTypes$plotID[plotTypes$forType=="ABLAL"],]
picol_df = model_df[model_df$plot %in% plotTypes$plotID[plotTypes$taxonID=="PICOL"],]
pien_df = model_df[model_df$plot %in% plotTypes$plotID[plotTypes$taxonID=="PIEN"],]
pifl2_df = model_df[model_df$plot %in% plotTypes$plotID[plotTypes$taxonID=="PIFL2"],]
picol_df$plot
plotTypes
# for(i in 1:length(pcs)) {pcs[[i]]$TRI[pcs[[i]]$TRI==0 | is.na(pcs[[i]]$TRI)] = 0.1}

################################################################################
###########                  Modeling                                ###########
################################################################################


#------------------ estimate parameter. ########################################
# set output directory
# setwd("R:/foreco4/projects/documents/ltm_lidar/Output/Chen")

# ws_steps = function(x) {
#   y <- dplyr::case_when(is.na(x) ~ 0.001, x < 0 ~ 0.001, x < 2 ~
#                           1, x > c ~ 5, TRUE ~ b + (x * m))
#   return(y)
# }

#-----------------------------------> null model
null_model = function(int, pcs, trees_mapped, distSum) {
  nullws = function(x) {dplyr::case_when(is.na(x) ~ 0.001, x < 0 ~ 0.001, TRUE ~ int)}
  result = data.frame(errors = matrix(0, nrow=length(pcs)))
  for(i in 1:length(pcs)){
    treestmp = locate_trees(pcs[[i]], lmf(nullws))
    
    # Do the two datasets have the same # of trees? 
    crds_lmf = st_coordinates(treestmp)
    
    if(nrow(crds_lmf) > nrow(trees_mapped[[i]])) {
      crds1=trees_mapped[[i]]; crds2=crds_lmf; multiplier=1
    }else if(nrow(crds_lmf) < nrow(trees_mapped[[i]])) {
      crds1=crds_lmf; crds2=trees_mapped[[i]]; multiplier=-1
    }else if(nrow(crds_lmf) == nrow(trees_mapped[[i]])) {
      crds1=crds_lmf; crds2=trees_mapped[[i]]; multiplier=1
    }else if(nrow(crds_lmf)==0) {
      crds1=matrix(c(0,0,0),ncol=3,byrow=T); crds2=trees_mapped[[i]]; multiplier=-1
    }
    
  
    errors_tmp = rep(0,nrow(crds2))
    for(tree in 1:nrow(crds1)) {
      errors_tmp[tree] = min(proxy::dist(matrix(crds1[tree,],ncol=3,byrow=T), crds2))
    }
    
    # Penalty for finding the wrong number of trees
    if(nrow(crds1)!=nrow(crds2)) {
      for(tree in (nrow(crds1)+1):nrow(crds2)) {
        errors_tmp[tree] = mean(proxy::dist(rbind(crds1,crds2)))
      }
    }
    
    
    errorSum = distSum[i] + sum(errors_tmp) * multiplier
    
    result$errors[i] = errorSum
  }
  return(result$errors)
}


par <- list(int=1)
par_lo <- list(int=0.1)
par_hi <- list(int=10)


var <- list(mean = "predicted", x ="distSum", log = TRUE,
            pcs="pcs", trees_mapped="trees_mapped", distSum="distSum")


# model
results_null_pipa2 <-anneal(model = null_model, par = par, var = var, 
                           source_data = pipa2_df, 
                           par_lo = par_lo,
                           par_hi = par_hi, 
                           pdf = dnorm, 
                           dep_var = "distSum", 
                           initial_temp = 4, 
                           temp_red = 0.6, 
                           max_iter = 500)
results_null_picol <-anneal(model = null_model, par = par, var = var, 
                           source_data = picol_df, 
                           par_lo = par_lo,
                           par_hi = par_hi, 
                           pdf = dnorm, 
                           dep_var = "distSum", 
                           initial_temp = 4, 
                           temp_red = 0.6, 
                           max_iter = 500)
results_null_pien <-anneal(model = null_model, par = par, var = var, 
                           source_data = pien_df, 
                           par_lo = par_lo,
                           par_hi = par_hi, 
                           pdf = dnorm, 
                           dep_var = "distSum", 
                           initial_temp = 4, 
                           temp_red = 0.6, 
                           max_iter = 500)
results_null_pifl2 <-anneal(model = null_model, par = par, var = var, 
                           source_data = pifl2_df, 
                           par_lo = par_lo,
                           par_hi = par_hi, 
                           pdf = dnorm, 
                           dep_var = "distSum", 
                           initial_temp = 4, 
                           temp_red = 0.6, 
                           max_iter = 500)


(best <- results_null_model$best_pars)
(lower <- results_null_model$lower_limits)
(upper <- results_null_model$upper_limits)
results_null_pifl2$aic_corr



#-----------------------------------> Michaelis-Menten function
MM_model <- function(a, b, c, pcs, trees_mapped, distSum) {
  # logws = function(x) {(c)/(1+a*exp(-k*x))}
  MMws = function(x) {dplyr::case_when(is.na(x) ~ 0.001, x < 0 ~ 0.001,
                                       c >= x ~ 0.001,
                                        TRUE ~ a/(b/(x-c)+1))} 
  #x > 26.5 ~ 5, 
  result = data.frame(errors = matrix(0, nrow=length(pcs)))
  for(i in 1:length(pcs)){
    treestmp = locate_trees(pcs[[i]], lmf(MMws))
    
    
    # Do the two datasets have the same # of trees? 
    crds_lmf = st_coordinates(treestmp)
    
    if(nrow(crds_lmf) > nrow(trees_mapped[[i]])) {
      crds1=trees_mapped[[i]]; crds2=crds_lmf; multiplier=1
    }else if(nrow(crds_lmf) < nrow(trees_mapped[[i]])) {
      crds1=crds_lmf; crds2=trees_mapped[[i]]; multiplier=-1
    }else if(nrow(crds_lmf) == nrow(trees_mapped[[i]])) {
      crds1=crds_lmf; crds2=trees_mapped[[i]]; multiplier=1
    }else if(nrow(crds_lmf)==0) {
      crds1=matrix(c(0,0,0),ncol=3,byrow=T); crds2=trees_mapped[[i]]; multiplier=-1
    }
    
    
    errors_tmp = rep(0,nrow(crds2))
    for(tree in 1:nrow(crds1)) {
      errors_tmp[tree] = min(proxy::dist(matrix(crds1[tree,],ncol=3,byrow=T), crds2))
    }
    
    # Penalty for finding the wrong number of trees
    if(nrow(crds1)!=nrow(crds2)) {
      for(tree in (nrow(crds1)+1):nrow(crds2)) {
        errors_tmp[tree] = mean(proxy::dist(rbind(crds1,crds2)))
      }
    }
    
    
    errorSum = distSum[i] + sum(errors_tmp) * multiplier
    
    result$errors[i] = errorSum
  }
  return(result$errors)
}

par <- list(a=1,b=2,c=1)
par_lo <- list(a=0.1,b=0, c=0)
par_hi <- list(a=50, b=50,c=10)


var <- list(mean = "predicted", x ="distSum", log = TRUE,
            pcs="pcs", trees_mapped="trees_mapped", distSum="distSum")


# model
results_MM_pipa2 <-anneal(model = MM_model, par = par, var = var, 
                           source_data = pipa2_df, 
                           par_lo = par_lo,
                           par_hi = par_hi, 
                           pdf = dnorm, 
                           dep_var = "distSum", 
                           initial_temp = 4, 
                           temp_red = 0.6, 
                           max_iter = 500)
results_MM_picol <-anneal(model = MM_model, par = par, var = var, 
                           source_data = picol_df, 
                           par_lo = par_lo,
                           par_hi = par_hi, 
                           pdf = dnorm, 
                           dep_var = "distSum", 
                           initial_temp = 4, 
                           temp_red = 0.6, 
                           max_iter = 500)
results_MM_model <-anneal(model = MM_model, par = par, var = var, 
                           source_data = picol_df, 
                           par_lo = par_lo,
                           par_hi = par_hi, 
                           pdf = dnorm, 
                           dep_var = "distSum", 
                           initial_temp = 4, 
                           temp_red = 0.6, 
                           max_iter = 500)
results_MM_pifl2 <-anneal(model = MM_model, par = par, var = var, 
                           source_data = pifl2_df, 
                           par_lo = par_lo,
                           par_hi = par_hi, 
                           pdf = dnorm, 
                           dep_var = "distSum", 
                           initial_temp = 4, 
                           temp_red = 0.6, 
                           max_iter = 500)


(best <- results_MM_pipa2$best_pars)
(lower <- results_MM_model$lower_limits)
(upper <- results_MM_model$upper_limits)
results_MM_model$aic_corr



#-----------------------------------> Logistic function
logistic_model <- function(a, c, k, pcs, trees_mapped, distSum) {
  logisticws = function(x) {dplyr::case_when(is.na(x) ~ 0.001, x < 0 ~ 0.001, 
                                        TRUE ~ c/(1+a*exp(-k*x)))} 
  #x > 26.5 ~ 5, 
  result = data.frame(errors = matrix(0, nrow=length(pcs)))
  for(i in 1:length(pcs)){
    treestmp = locate_trees(pcs[[i]], lmf(logisticws))
    
    
    # Do the two datasets have the same # of trees? 
    crds_lmf = st_coordinates(treestmp)
    
    if(nrow(crds_lmf) > nrow(trees_mapped[[i]])) {
      crds1=trees_mapped[[i]]; crds2=crds_lmf; multiplier=1
    }else if(nrow(crds_lmf) < nrow(trees_mapped[[i]])) {
      crds1=crds_lmf; crds2=trees_mapped[[i]]; multiplier=-1
    }else if(nrow(crds_lmf) == nrow(trees_mapped[[i]])) {
      crds1=crds_lmf; crds2=trees_mapped[[i]]; multiplier=1
    }else if(nrow(crds_lmf)==0) {
      crds1=matrix(c(0,0,0),ncol=3,byrow=T); crds2=trees_mapped[[i]]; multiplier=-1
    }
    
    
    errors_tmp = rep(0,nrow(crds2))
    for(tree in 1:nrow(crds1)) {
      errors_tmp[tree] = min(proxy::dist(matrix(crds1[tree,],ncol=3,byrow=T), crds2))
    }
    
    # Penalty for finding the wrong number of trees
    if(nrow(crds1)!=nrow(crds2)) {
      for(tree in (nrow(crds1)+1):nrow(crds2)) {
        errors_tmp[tree] = mean(proxy::dist(rbind(crds1,crds2)))
      }
    }
    
    
    errorSum = distSum[i] + sum(errors_tmp) * multiplier
    
    result$errors[i] = errorSum
  }
  return(result$errors)
}

par <- list(a=10, c=5,k=1)
par_lo <- list(a=0, c=0.1, k=0)
par_hi <- list(a=50, c=15, k=10)


var <- list(mean = "predicted", x ="distSum", log = TRUE,
            pcs="pcs", trees_mapped="trees_mapped", distSum="distSum")


# model
results_logistic_model <-anneal(model = logistic_model, par = par, var = var, 
                           source_data = pipa2_df, 
                           par_lo = par_lo,
                           par_hi = par_hi, 
                           pdf = dnorm, 
                           dep_var = "distSum", 
                           initial_temp = 6, 
                           temp_red = 0.6, 
                           max_iter = 500)


(best <- results_logistic_model$best_pars)
(lower <- results_logistic_model$lower_limits)
(upper <- results_logistic_model$upper_limits)
results_logistic_model$aic_corr



#-----------------------------------> logarithmic function
log_model <- function(a, c, k, pcs, trees_mapped, distSum) {
  # logws = function(x) {(c)/(1+a*exp(-k*x))}
  logws = function(x) {dplyr::case_when(is.na(x) ~ 0.001, x < 0 ~ 0.001, 
                                        TRUE ~ log((x+c)^k)+a)} 
  #x > 26.5 ~ 5, 
  result = data.frame(errors = matrix(0, nrow=length(pcs)))
  for(i in 1:length(pcs)){
    treestmp = locate_trees(pcs[[i]], lmf(logws))
    
    # Do the two datasets have the same # of trees? 
    crds_lmf = st_coordinates(treestmp)
    
    if(nrow(crds_lmf) > nrow(trees_mapped[[i]])) {
      crds1=trees_mapped[[i]]; crds2=crds_lmf; multiplier=1
    }else if(nrow(crds_lmf) < nrow(trees_mapped[[i]])) {
      crds1=crds_lmf; crds2=trees_mapped[[i]]; multiplier=-1
    }else if(nrow(crds_lmf) == nrow(trees_mapped[[i]])) {
      crds1=crds_lmf; crds2=trees_mapped[[i]]; multiplier=1
    }else if(nrow(crds_lmf)==0) {
      crds1=matrix(c(0,0,0),ncol=3,byrow=T); crds2=trees_mapped[[i]]; multiplier=-1
    }
    
    
    errors_tmp = rep(0,nrow(crds2))
    for(tree in 1:nrow(crds1)) {
      errors_tmp[tree] = min(dist(matrix(crds1[tree,],ncol=3,byrow=T), crds2))
    }
    
    # Penalty for finding the wrong number of trees
    if(nrow(crds1)!=nrow(crds2)) {
      for(tree in (nrow(crds1)+1):nrow(crds2)) {
        errors_tmp[tree] = mean(dist(rbind(crds1,crds2)))
      }
    }
    
    
    errorSum = distSum[i] + sum(errors_tmp) * multiplier
    
    result$errors[i] = errorSum
  }
  return(result$errors)
}

par <- list(a=1, c=2, k=1)
par_lo <- list(a=0.1, c=1, k=0)
par_hi <- list(a=10, c=10, k=10)


var <- list(mean = "predicted", x ="distSum", log = TRUE,
            pcs="pcs", trees_mapped="trees_mapped", distSum="distSum")


# model
results_log_model <-anneal(model = log_model, par = par, var = var, 
                           source_data = picol_df, 
                           par_lo = par_lo,
                           par_hi = par_hi, 
                           pdf = dnorm, 
                           dep_var = "distSum", 
                           initial_temp = 6, 
                           temp_red = 0.6, 
                           max_iter = 500)


(best <- results_log_model$best_pars)
(lower <- results_log_model$lower_limits)
(upper <- results_log_model$upper_limits)

 
results_null_model$best_pars

# See how good we did. 
best = summary(lin_mod)
mm_best = results_MM_ablal$best_pars
MMws = function(x) {dplyr::case_when(is.na(x) ~ 0.001, x < 0 ~ 0.001, 
                                     TRUE ~ 4.9774009*x/(5.3992996+x))}

mdlThing_ests = subset(mdl_dfThing_summary, deltaAIC==0, select="Estimate")
mono_ws = function(x) {dplyr::case_when(is.na(x) ~ 0.001, x <= 0 ~ 0.001, 
                                     TRUE ~ mdlThing_ests[1,]*(1-exp(-mdlThing_ests[2,]*x)))}

lin_ws = function(x) {dplyr::case_when(is.na(x) ~ 0.001, x <= 0 ~ 0.001, 
                                     TRUE ~ best$coefficients[1]+best$coefficients[2]*x)}

logit_ws = function(x) {dplyr::case_when(is.na(x) ~ 0.001, x <= 0 ~ 0.001,
                                     TRUE ~ 3.68398/(1+5.10323*exp(-0.37438*x)))}

exp_ws = function(x) {dplyr::case_when(is.na(x) ~ 0.001, x <= 0 ~ 0.001,
                                     TRUE ~ 3.788005*1.042941^x)}

tpp = list()
for(i in 1:nrow(model_df)) {
  tppe_tmp = nrow(locate_trees(model_df$pcs[[i]], lmf(MMws)))
  
  tpp[[length(tpp)+1]] = tppe_tmp
}
tpp
tppo = lapply(model_df$trees_mapped, function(x) nrow(x))
# tppo = vegsf %>% group_by(plotID) %>% count() %>% select(n)
# tppo$n

plot(unlist(tpp)~unlist(tppo))
t = lm(unlist(tpp) ~ unlist(tppo))
summary(t)
abline(a=5.8393, b=0.8194)
abline(a=0,b=1,lty=2)

# Mean absolute percentage error
MM_mape = sum((abs(unlist(tppo)-unlist(tpp))/unlist(tppo)))*100/length(tpp)



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


max(unlist(lapply(pcs, function(x) max(x$Z))))



las = readLAS("../../Bofo/CalWood/Polygon2/Terra_Output/cloud0.las")
las = las %>% classify_noise(sor()) %>% 
  filter_poi(Classification!=LASNOISE) %>%
  normalize_height(knnidw()) %>% 
  filter_poi(Z>0)
plot(las, color="RGB")

# # lasR style. This shit seems to hang up and never finish. Faster my ass.
# library(lasR)
# del = triangulate(filter = keep_ground())
# norm = transform_with(del, "-")
# write = write_las(ofile="../../Bofo/CalWood/Polygon2/Terra_Output/cloud0Norm.laz")
# pipeline = del + norm + write
# ans1 = exec(pipeline, on = "../../Bofo/CalWood/Polygon2/Terra_Output/cloud0.las")
# ans1


del = triangulate(filter = keep_first())
chm = rasterize(0.25, del)
chm2 = pit_fill(chm)
seed2 = local_maximum_raster(chm2, 2)
seed3 = local_maximum_raster(chm2, 3)
seed4 = local_maximum_raster(chm2, 4)
tree2 = region_growing(chm2, seed2)
tree3 = region_growing(chm2, seed3)
tree4 = region_growing(chm2, seed4)
pipeline = del + chm + chm2 +  seed2 + seed3 + seed4 + tree2 + tree3 + tree4
ans = exec(pipeline, on = "../../Bofo/CalWood/Polygon2/Terra_Output/cloud0.las")


# canopy height model
chm = rasterize_canopy(las,res=0.25)

# tree top height to crown diameter data function
ttopZD = function(las, chm, ws) {
  ttops_tmp = locate_trees(chm, lmf(ws))
  las_tmp = segment_trees(las, dalponte2016(chm, ttops_tmp))
  crowns_tmp = as.data.frame(crown_metrics(las_tmp, func=.stdtreemetrics, geom="convex")) %>% 
    filter(convhull_area>0) %>% 
    mutate(radius = sqrt(convhull_area/pi)) %>% 
    mutate(diam = radius*2)
  
  return(crowns_tmp)
}


df1 = ttopZD(las, chm, 1)
df1.5 = ttopZD(las, chm, 1.5)
df2 = ttopZD(las, chm, 2)
df3 = ttopZD(las, chm, 3)
df3.5 = ttopZD(las, chm, 3.5)
df4 = ttopZD(las, chm, 4)
df5 = ttopZD(las, chm, 5)

plot(radius~Z, df1)
plot(radius~Z, df1.5)
plot(radius~Z, df2)
plot(radius~Z, df3)
plot(radius~Z, df4)
plot(radius~Z, df5)
plot(diam~Z, df1)
plot(diam~Z, df1.5)
plot(diam~Z, df2)
plot(diam~Z, df3)
plot(diam~Z, df4)
plot(diam~Z, df5)
dev.off()

run_models = function(data, ...)
#-----------------------------------> Michaelis-Menten function
null_model = function(int){
  int
}
par_null <- list(int=1)
par_nulllo <- list(int=0)
par_nullhi <- list(int=20)

var_null <- list(mean = "predicted", x ="diam",log = TRUE)

results_null_model <-anneal(model = null_model, par = par_null, var = var_null, 
                           source_data = df1.5, 
                           par_lo = par_nulllo,
                           par_hi = par_nullhi, 
                           pdf = dnorm, 
                           dep_var = "diam", 
                           initial_temp = 6, 
                           temp_red = 0.6, 
                           max_iter = 10000)

MM_model <- function(a, b, c, Z) {
  # a*Z/(b+Z)
  a/(b/(Z-c)+1)
}

h3_model <- function(a, b, c, Z) {
  a*Z^2/(b^2+Z^2)
}

lin_model <- function(a,b,Z) {
  a*Z+b
}

logistic_model = function(a,b,c,Z) {             
  #c/(1+a*exp(-b*Z))
  (exp(a+b*Z)/(1+exp(a+b*Z)))+c
} 

# Monomolecular function
mono_model = function(a,b,c,Z) {             
  a*(1-exp(-b*Z+b*c))
} 

par <- list(a=1,b=1, c=1)
par_lo <- list(a=0.1,b=0, c=-10)
par_hi <- list(a=10, b=10, c=10)

par <- list(a=1,b=1)
par_lo <- list(a=0.1,b=0)
par_hi <- list(a=10, b=10)

par <- list(a=1,b=1, c=1)
par_lo <- list(a=0.1,b=0, c=0.001)
par_hi <- list(a=10, b=10, c=10)


var <- list(mean = "predicted", x ="radius", Z="Z",log = TRUE)

results_MMG_model <-anneal(model = MM_model, par = par, var = var, 
                                source_data = df4, 
                                par_lo = par_lo,
                                par_hi = par_hi, 
                                pdf = dnorm, 
                                dep_var = "radius", 
                                initial_temp = 6, 
                                temp_red = 0.6, 
                                max_iter = 10000)

results2 = list(results_MMG_model$best_pars, results_MMG_model$R2, results_MMG_model$aic_corr)
# results3 = list(results_MMG_model$best_pars, results_MMG_model$R2, results_MMG_model$aic_corr)
# results4 = list(results_MMG_model$best_pars, results_MMG_model$R2, results_MMG_model$aic_corr)
# results5 = list(results_MMG_model$best_pars, results_MMG_model$R2, results_MMG_model$aic_corr)
# results8 = list(results_MMG_model$best_pars, results_MMG_model$R2, results_MMG_model$aic_corr)
# results10 = list(results_MMG_model$best_pars, results_MMG_model$R2, results_MMG_model$aic_corr)
# results12 = list(results_MMG_model$best_pars, results_MMG_model$R2, results_MMG_model$aic_corr)
# results16 = list(results_MMG_model$best_pars, results_MMG_model$R2, results_MMG_model$aic_corr)

# See how good we did. 
mm_best = results2
MMws = function(x) {dplyr::case_when(is.na(x) ~ 0.001, x <= 0 ~ 0.001, x < 4 ~ 1,
                                     TRUE ~ mm_best[[1]]$a/(mm_best[[1]]$b/(x-mm_best[[1]]$c)+1))}
monows = function(x) {dplyr::case_when(is.na(x) ~ 0.001, x < 2 ~ 0.001, 
                                     TRUE ~ mm_best[[1]]$a*(1-exp(-mm_best[[1]]$b*x))+mm_best[[1]]$c)}
linws = function(x) {dplyr::case_when(is.na(x) ~ 0.001, x <= 0 ~ 0.001, 
                                     TRUE ~ mm_best[[1]]$a*(1-exp(-mm_best[[1]]$b*x))+mm_best[[1]]$c)}



ttops = ttopZD(las, chm, monows)
plot(radius~Z, df2)
t = lm(diam~Z, ttops)
summary(t)
ttops = locate_trees(chm, lmf(MMws))

plot(chm)
zoom(chm)
plot(ttops['treeID'], col="red", pch=20, add=T)
dev.off()
