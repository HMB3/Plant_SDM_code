#########################################################################################################################
#############################################  FIT MAXENT TO GBIF DATA ################################################## 
#########################################################################################################################


## This code takes a table of all species occurrences (rows) and environmental values (columns), and runs SDMs for a list
## of taxa


## The brief again:

## The first module will focus on x plant species identified in the project’s Target Species List, and will develop maps 
## that demonstrate each species’ suitability to both current and future climates across Australia.

## These maps will be used to demonstrate how well or poorly a particular species will be able to tolerate future conditions 
## in urban centres across Australia as the climate changes, based on our current understanding of species’ climatic 
## requirements.


#########################################################################################################################
## Load packages, functions and data
load("./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT.RData")
load("./data/base/HIA_LIST/COMBO/HIA_SDM_DATA_ALL_VAR.RData")
load("./data/base/HIA_LIST/COMBO/SDM_TEMPLATE_RASTER.RData")

source('./R/HIA_LIST_MATCHING.R')
source('./R/HIA_CLEAN_MATCHING.R')
source('./R/GREEN_CITIES_FUNCTIONS.R')
source('./R/SDM_FUNCTIONS.R')


## Check data 
str(template.raster)
str(SDM.DATA.ALL)





#########################################################################################################################
## 3). RUN SDM FOR SELECTED VARIABLES
#########################################################################################################################


#########################################################################################################################
## Chose a-priori worldclim predictors
sdm.predictors <- c("Annual_mean_temp",    "Mean_diurnal_range",  "Isothermality",      "Temp_seasonality",  
                    "Max_temp_warm_month", "Min_temp_cold_month", "Temp_annual_range",  "Mean_temp_wet_qu",    
                    "Mean_temp_dry_qu",    "Mean_temp_warm_qu",   "Mean_temp_cold_qu",  "Annual_precip", 
                    "Precip_wet_month",    "Precip_dry_month",    "Precip_seasonality", "Precip_wet_qu",     
                    "Precip_dry_qu",       "Precip_warm_qu",      "Precip_col_qu")


## Are the latest experimental species in there?
head(test.spp, 10)


########################################################################################################################
## We can run Maxent from a cluster of cores on the local computer. Here we send (i.e. export) all the necessary 
## ingredients to the cluster. So that's the:


## template.raster raster, 
## data frame and 
## maxent function

## Use these functions to debug
## debugonce(FIT_MAXENT)
## undebug(FIT_MAXENT)

## Use this function to find all functions relating to a search term/topic (e.g. 'debug')
## apropos('debug')
## apropos('read')


## 100 species takes about 4 hours...
cl <- makeCluster(4)
clusterExport(cl, c('template.raster', 'SDM.DATA.ALL', 'FIT_MAXENT', 'COR_VARIABLES'))
clusterEvalQ(cl, {
  
  require(ff)
  require(rgeos)
  require(sp)
  require(raster)
  require(rJava)
  require(dismo)
  require(things)
  
})


########################################################################################################################
## Now use 'lapply' to run maxent for multiple species
lapply(spp.all[1:length(spp.all)], function(x)  { # for serial, parLapply(cl, species[1:8], function(x) { # for parallel 
  
  ## Print the taxa being processed to screen
  if(x %in% SDM.DATA.ALL$searchTaxon) {
    message('Doing ', x) 
    
    ## Subset the records to only the taxa being processed
    occurrence         <- subset(SDM.DATA.ALL, searchTaxon == x)
    occurrence$species <- x
    
    ## Now get the background points. These can come from anywhere in the whole dataset,other than the species used.
    background         <- subset(SDM.DATA.ALL, searchTaxon != x)
    background$species <- x
    
    ## The create a vector of the sdm.predictors used. 
    ## This should be based on an ecological framework! 
    sdm.predictors <- sdm.predictors # vector of used sdm.predictors
    min_n          = 20
    
    ## Fit the models using FIT_MAXENT. Would be good to make skipping exisitng outputs an argument
    FIT_MAXENT_SIMP(occ                     = occurrence, 
                    bg                      = background, 
                    sdm.predictors          = sdm.predictors, 
                    name                    = x, 
                    outdir                  = 'output/maxent/STD_VAR_ALL', 
                    template.raster,
                    min_n                   = 20,   ## This should be higher...
                    max_bg_size             = 100000,
                    background_buffer_width = 200000,
                    shapefiles              = TRUE,
                    features                = 'lpq',
                    replicates              = 5,
                    responsecurves          = TRUE)
    
  } else {
    
    message(x, ' skipped - no data.')         ## This condition ignores species which have no data
    
  }  
  
})
  


stopCluster(cl)





#########################################################################################################################
## OUTSTANDING SDM TASKS:
#########################################################################################################################


## Further mapping and cleaning of GBIF data needed for the important species: Use worldview?

## Need a condition in the loop to skip folders in the list, can this be an argument in the loop set to "TRUE" OR "FALSE"? Can be time-stamped

## Store the list of skipped species...



#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################