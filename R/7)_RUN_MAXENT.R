#########################################################################################################################
#############################################  FIT MAXENT TO GBIF DATA ################################################## 
#########################################################################################################################


#########################################################################################################################
## This code takes a table of all species occurrences (rows) and environmental values (columns), and runs SDMs for a list
## of taxa.


## The HIA brief again:

## The first module will focus on x plant species identified in the project’s Target Species List, and will develop maps 
## that demonstrate each species’ suitability to both current and future climates across Australia.

## These maps will be used to demonstrate how well or poorly a particular species will be able to tolerate future conditions 
## in urban centres across Australia as the climate changes, based on our current understanding of species’ climatic 
## requirements.



#########################################################################################################################
## Load packages, functions and data
#source('./R/HIA_LIST_MATCHING.R')
rasterTmpFile()


#########################################################################################################################
## 1). PREPARE DATA TABLE FOR SDM ANALYSIS
#########################################################################################################################


#########################################################################################################################
## Create a table with all the variables
COMBO.RASTER.ALL  <- dplyr::select(CLEAN.TRUE, searchTaxon, lon, lat,
                                   
                                   Annual_mean_temp,     Mean_diurnal_range,  Isothermality,     Temp_seasonality, 
                                   Max_temp_warm_month,  Min_temp_cold_month, Temp_annual_range, Mean_temp_wet_qu,
                                   Mean_temp_dry_qu,     Mean_temp_warm_qu,   Mean_temp_cold_qu, 
                                   
                                   Annual_precip,        Precip_wet_month,    Precip_dry_month,  Precip_seasonality,   
                                   Precip_wet_qu,        Precip_dry_qu,       Precip_warm_qu,    Precip_col_qu)


#########################################################################################################################
## Create a spatial points object, and change to a projected system to calculate distance more accurately 
coordinates(COMBO.RASTER.ALL)    <- ~lon+lat
proj4string(COMBO.RASTER.ALL)    <- '+init=epsg:4326'
COMBO.RASTER.ALL                 <- spTransform(COMBO.RASTER.ALL, CRS('+init=ESRI:54009'))


## Now split using the data using the species column, and get the unique occurrence cells
COMBO.RASTER.SPLIT.ALL <- split(COMBO.RASTER.ALL, COMBO.RASTER.ALL$searchTaxon)
occurrence_cells_all   <- lapply(COMBO.RASTER.SPLIT.ALL, function(x) cellFromXY(template.raster, x))
length(occurrence_cells_all)   ## this is a list of dataframes, where the number of rows for each being the species table


#########################################################################################################################
## Now get just one record within each 10*10km cell.
SDM.DATA.ALL <- mapply(function(x, cells) {
  x[!duplicated(cells), ]
}, COMBO.RASTER.SPLIT.ALL, occurrence_cells_all, SIMPLIFY = FALSE) %>% do.call(rbind, .)


## Check to see we have 19 variables + the species for the standard predictors, and 19 for all predictors
## Check data :: template, data table and species 
dim(template.raster)
dim(SDM.DATA.ALL)
length(unique(SDM.DATA.ALL$searchTaxon))
class(Koppen_1975)


## What resolution is the template raster at?
xres(template.raster);yres(template.raster)


#########################################################################################################################
## Now create the variables needed to access current environmental conditions + their names in the functions
sdm.predictors <- c("Annual_mean_temp",    "Mean_diurnal_range",  "Isothermality",      "Temp_seasonality",  
                    "Max_temp_warm_month", "Min_temp_cold_month", "Temp_annual_range",  
                    "Mean_temp_wet_qu",    "Mean_temp_dry_qu",    
                    "Mean_temp_warm_qu",   "Mean_temp_cold_qu",   
                    
                    "Annual_precip",       "Precip_wet_month",   "Precip_dry_month",    "Precip_seasonality",  
                    "Precip_wet_qu",       "Precip_dry_qu",      
                    "Precip_warm_qu",      "Precip_col_qu")


## A-priori worldclim predictors
sdm.select     <- c("Annual_mean_temp", "Temp_seasonality",    "Max_temp_warm_month", "Min_temp_cold_month",
                    "Annual_precip",    "Precip_seasonality",  
                    "Precip_wet_month", "Precip_dry_month")      


## Create a raster stack of current environmental conditions if needed
i  <- match(sdm.predictors, sdm.predictors)
ff <- file.path('./data/base/worldclim/world/0.5/bio/current',
                sprintf('bio_%02d.tif', i))


## Name the grids :: these should be indentical
env.grids.current = stack(sub('0.5', '1km', ff))
names(env.grids.current) <- sdm.predictors[i]
identical(names(env.grids.current),sdm.predictors)


#########################################################################################################################
## Loop over the species list and plot the occurrence data for each to check the data bias
# TAXA = as.list(sort(unique(intersect(SDM.DATA.ALL$searchTaxon, spp.mile))))
# 
# for (i in 1:length(TAXA)) {
# 
#   ## Create points for each species
#   spp.points <- SDM.DATA.ALL[SDM.DATA.ALL$searchTaxon == TAXA[i], ] %>%
#     spTransform(CRS('+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'))
# 
#   ## Print to file
#   save_name = gsub(' ', '_', TAXA[i])
#   save_dir  = "output/maxent/summary"
#   png(sprintf('%s/%s_%s.png', save_dir,
#               save_name, "Australian_points"),
#       3236, 2000, units = 'px', res = 300)
# 
#   ## set margins
#   par(mar   = c(3, 3, 5, 3),  ## b, l, t, r
#       #mgp   = c(9.8, 2.5, 0),
#       oma   = c(1.5, 1.5, 1.5, 1.5))
# 
#   ## Plot just the Australian points
#   plot(aus, main = TAXA[i])
#   points(spp.points, col = "red", cex = .3, pch = 19)
# 
#   ## Finish the device
#   dev.off()
# 
# }




#########################################################################################################################
## 2). RUN SDMs USING A-PRIORI VARIABLES
#########################################################################################################################


#########################################################################################################################
## Run Maxent using a random selection of background points. 
## synonyms = c("Waterhousea floribunda", "Sannantha virgata", "Callistemon citrinus")
## 'Angophora costata' %in% SDM.DATA.ALL$searchTaxon  
projection(template.raster);projection(SDM.DATA.ALL);projection(Koppen_1975)


## Create new variable for analysis species - hence SDM.DATA.ALL needs to be update each time...........................
analysis.spp = unique(SDM.DATA.ALL$searchTaxon)   ## could use this in future


## Loop over all the species spp = spp.combo[31]
lapply(analysis.spp, function(spp){ 
  
  ## Skip the species if the directory already exists, before the loop
  outdir <- 'output/maxent/SET_VAR_KOPPEN'
  if(dir.exists(file.path(outdir, gsub(' ', '_', spp)))) {
    message('Skipping ', spp, ' - already run.')
    invisible(return(NULL))
    
  }
  
  ## Print the taxa being processed to screen
  if(spp %in% SDM.DATA.ALL$searchTaxon) {
    message('Doing ', spp) 
    
    ## Subset the records to only the taxa being processed
    occurrence <- subset(SDM.DATA.ALL, searchTaxon == spp)
    
    ## Now get the background points. These can come from any spp, other than the modelled species.
    background <- subset(SDM.DATA.ALL, searchTaxon != spp)
    
    ## Finally fit the models using FIT_MAXENT_TARG_BG. Also use tryCatch to skip any exceptions
    tryCatch(
      FIT_MAXENT_TARG_BG(occ                     = occurrence, 
                         bg                      = background, 
                         sdm.predictors          = sdm.select, 
                         name                    = spp, 
                         outdir                  = outdir, 
                         template.raster,
                         min_n                   = 20,     ## This should be higher...
                         max_bg_size             = 100000, 
                         Koppen                  = Koppen_1975,
                         background_buffer_width = 200000,
                         shapefiles              = TRUE,
                         features                = 'lpq',
                         replicates              = 5,
                         responsecurves          = TRUE),
      
      ## https://stackoverflow.com/questions/19394886/trycatch-in-r-not-working-properly
      #function(e) message('Species skipped ', spp)) ## skip any species for which the function fails
      error = function(cond) {
        
        message(paste('Species skipped ', spp))
        
      })
    
  } else {
    
    message(spp, ' skipped - no data.')         ## This condition ignores species which have no data...
    
  }  
  
})



# #########################################################################################################################
# ## Run analyses using parallel processing
# 
# ## To run the function in a cluster, we need to export all the objects which are not defined in the function call
# cl <- makeCluster(4)
# clusterExport(cl, c('template.raster', 'SDM.DATA.ALL', 'FIT_MAXENT', 'FIT_MAXENT_TARG_BG'))
# clusterEvalQ(cl, {
#   
#   ## These packages are required :: check they are updated
#   require(ff)
#   require(rgeos)
#   require(sp)
#   require(raster)
#   require(rJava)
#   require(dismo)
#   require(things)
#   
# })
# 
# 
# ## Loop over all the species spp = spp.combo[31]
# parLapply(cl, spp.combo, function(spp) { # for serial, lapply(spp.combo, function(spp)
#   
#   ## Skip the species if the directory already exists, before the loop
#   outdir <- 'output/maxent/SET_VAR_KOPPEN'
#   if(dir.exists(file.path(outdir, gsub(' ', '_', spp)))) {
#     message('Skipping ', spp, ' - already run.')
#     invisible(return(NULL))
#     
#   }
#   
#   ## Print the taxa being processed to screen
#   if(spp %in% SDM.DATA.ALL$searchTaxon) {
#     message('Doing ', spp) 
#     
#     ## Subset the records to only the taxa being processed
#     occurrence <- subset(SDM.DATA.ALL, searchTaxon == spp)
#     
#     ## Now get the background points. These can come from any spp, other than the modelled species.
#     background <- subset(SDM.DATA.ALL, searchTaxon != spp)
#     
#     ## Finally fit the models using FIT_MAXENT_TARG_BG. Also use tryCatch to skip any exceptions
#     tryCatch(
#       FIT_MAXENT_TARG_BG(occ                     = occurrence, 
#                          bg                      = background, 
#                          sdm.predictors          = sdm.select, 
#                          name                    = spp, 
#                          outdir                  = outdir, 
#                          template.raster,
#                          min_n                   = 20,     ## This should be higher...
#                          max_bg_size             = 100000, 
#                          Koppen                  = Koppen_1975,
#                          background_buffer_width = 200000,
#                          shapefiles              = TRUE,
#                          features                = 'lpq',
#                          replicates              = 5,
#                          responsecurves          = TRUE),
#       
#       ## https://stackoverflow.com/questions/19394886/trycatch-in-r-not-working-properly
#       #function(e) message('Species skipped ', spp)) ## skip any species for which the function fails
#       error = function(cond) {
#         
#         message(paste('Species skipped ', spp))
#         
#       })
#     
#   } else {
#     
#     message(spp, ' skipped - no data.')         ## This condition ignores species which have no data...
#     
#   }  
#   
# })
# 
# 
# stopCluster(cl)
# 
# 
# #########################################################################################################################
# ## Look through the output directory for species where the code didn't finish. This is usually species with < 27 files
# ## in the top "species_name" directory
# dd <- list.dirs('H:/green_cities_sdm/output/maxent/SET_VAR_KOPPEN', recursive = F)
# 
# 
# ## Loop over all the directories in the maxent output folder
# ff <- sapply(dd, function(d) {
#   
#   list.files(d, full.names = TRUE)
#   
# })
# 
# 
# ## Now remove all the files in each directory, if that directory has < 27 files in the top level
# n <- lengths(ff)
# 
# lapply(ff[n < 27], function(x) {
#   
#   file.remove(x)
#   unlink(dirname(x[1]), recursive = TRUE)
#   
# })





#########################################################################################################################
## OUTSTANDING SDM TASKS:
#########################################################################################################################

## 1). Make list of exotic spp on the top 200 sellers, which have no data and/or produce bad maps

## 2). Re-download ALA data. Add exotic urban inventory data that Ale is compiling (create diagram of data integration, to highlight knowlegde gaps)
##     Four columns ::
##     Species
##     Common name
##     LAT/LONG
##     SOURCE

##     Extra sources:
##     iNaturalist  - Got the data for plants, 20k records from Oceania
##     Flickr       - looks hard....

## 3). Try to thin records for ~100 spp with boundary bias: random sampling of those species records, by state and environment

## 4). Use more forgiving thresholds (10%) for all species, OR just those with bad maps:
##    "Maximum.training.sensitivity.plus.specificity.Logistic.threshold"
##    "X10.percentile.training.presence.Logistic.threshold"
##    "X10.percentile.training.presence.training.omission"
  
## 5). Calculate the TSS for all species (and MESS maps for a few species)

## 7). Decide on gamma diversity for this article - 200 spp - could it go to GEB? Plan figures and tables for MS





#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################