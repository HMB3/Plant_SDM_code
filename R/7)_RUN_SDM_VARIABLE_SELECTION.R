#########################################################################################################################
#############################################  FIT MAXENT TO GBIF DATA ################################################## 
#########################################################################################################################


#########################################################################################################################
## This code takes a table of all species occurrences (rows) and environmental values (columns), and runs SDMs for a list
## of taxa


## The HIA brief again:

## The first module will focus on x plant species identified in the project’s Target Species List, and will develop maps 
## that demonstrate each species’ suitability to both current and future climates across Australia.

## These maps will be used to demonstrate how well or poorly a particular species will be able to tolerate future conditions 
## in urban centres across Australia as the climate changes, based on our current understanding of species’ climatic 
## requirements.



#########################################################################################################################
## Load packages, functions and data
#load("./data/base/HIA_LIST/COMBO/COMBO_RASTER_CONTEXT.RData")
#load("./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT.RData")
load("./data/base/HIA_LIST/COMBO/SDM_TEMPLATE_RASTER.RData")

load("./data/base/HIA_LIST/COMBO/HIA_SDM_DATA_ALL_VAR.RData")
#load("./data/base/HIA_LIST/COMBO/HIA_SDM_DATA_TEST_SPP.RData")
#load("./data/base/HIA_LIST/COMBO/SDM_DATA_TEST_CLEAN.RData")

source('./R/HIA_LIST_MATCHING.R')
source('./R/MAXENT_FUNCTIONS.R')


## Check data 
dim(template.raster)
dim(SDM.DATA.ALL)
length(unique(SDM.DATA.ALL$searchTaxon))


## What resolution is the template raster at? This will affect the density sampling...
xres(template.raster);yres(template.raster)


#########################################################################################################################
## Chose a-priori worldclim predictors: 8, 9, 18 and 19 are suspect
sdm.predictors <- c("Annual_mean_temp",    "Mean_diurnal_range",  "Isothermality",      "Temp_seasonality",  
                    "Max_temp_warm_month", "Min_temp_cold_month", "Temp_annual_range",  
                    "Mean_temp_wet_qu",    "Mean_temp_dry_qu",    
                    "Mean_temp_warm_qu",   "Mean_temp_cold_qu",   
                    
                    "Annual_precip",       "Precip_wet_month",   "Precip_dry_month",    "Precip_seasonality",  
                    "Precip_wet_qu",       "Precip_dry_qu",      
                    "Precip_warm_qu",      "Precip_col_qu")


sdm.select     <- c("Annual_mean_temp", "Temp_seasonality",    "Max_temp_warm_month", "Min_temp_cold_month",
                    "Annual_precip",    "Precip_seasonality",  
                    "Precip_wet_month", "Precip_dry_month")      




#########################################################################################################################
## 3). RUN SDMs USING A-PRIORI VARIABLES
#########################################################################################################################


## Add kernel density to reduce the bias, and add minimum bg size...


## spp = kop.spp[8]
## Run without cluster
lapply(kop.spp, function(spp) { # for serial, parLapply(cl, species[1:8], function(x) { # for parallel 
  
  ## Print the taxa being processed to screen
  if(spp %in% SDM.DATA.ALL$searchTaxon) {
    message('Doing ', spp) 
    
    ## Subset the records to only the taxa being processed
    occurrence <- subset(SDM.DATA.ALL, searchTaxon == spp)
    
    ## Now get the background points. These can come from anywhere in the whole dataset,
    ## other than the species used.
    background <- subset(SDM.DATA.ALL, searchTaxon != spp)
    
    ## The create a vector of the sdm.predictors used. 
    sdm.predictors <- sdm.select # vector of used sdm.predictors 
    
    ## Finally fit the models using FIT_MAXENT. Use tryCatch to skip any exceptions
    tryCatch(
      FIT_MAXENT(occ                     = occurrence, 
                 bg                      = background, 
                 sdm.predictors          = sdm.select, 
                 name                    = spp, 
                 outdir                  = 'output/maxent/SET_VAR_DENSITY', 
                 template.raster,
                 min_n                   = 20,   ## This should be higher...
                 max_bg_size             = 100000, ## need a min bg size?
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





#########################################################################################################################
## 4). RUN SDMs USING BACKWARDS SELECTION
#########################################################################################################################


#########################################################################################################################
## Are the latest experimental species in there?
head(test.spp, 10)
head(spp.all,  10)
head(kop.spp,  10)


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
## traceback()


## Check the old variables are gone :
names(SDM.DATA.ALL)


## 100 species takes about 4 hours...
# cl <- makeCluster(4)
# clusterExport(cl, c('template.raster', 'SDM.DATA.ALL', 'FIT_MAXENT_SELECTION'))
# clusterEvalQ(cl, {
# 
#   require(ff)
#   require(rgeos)
#   require(sp)
#   require(raster)
#   require(rJava)
#   require(dismo)
#   require(things)
# 
# })


## Regarding maxent settings ...........................................................................................


## Also note that for some species there are only a few variables with < 0.7 correlation. So to reduce this problem :


## Increase cor_thr to 0.85, 
## Skip the species with < 5 uncorrelated variables
## Increase k_thr to avoid commission error


## Just project onto current and we have a look at how many species are falling out. If too many are failing, 
## then we just have to change our approach which means you then have to reproject onto future.


## spp = "Acacia implexa"
lapply(kop.spp, function(spp)  { # for serial, parLapply(cl, species[1:8], function(spp) { # for parallel 
  
  ## Print the taxa being processed to screen
  if(spp %in% SDM.DATA.ALL$searchTaxon) {
    message('Doing ', spp) 
    
    ## Subset the records to only the taxa being processed
    occurrence         <- subset(SDM.DATA.ALL, searchTaxon == spp)
    occurrence$species <- spp
    
    ## Now get the background points. These can come from anywhere in the whole dataset,other than the species used.
    background         <- subset(SDM.DATA.ALL, searchTaxon != spp)
    background$species <- spp
    
    ## Then create a vector of the sdm.predictors used: all bioclim variables
    sdm.predictors     <- sdm.predictors # vector of used sdm.predictors
    
    ## Fit the models using FIT_MAXENT. Would be good to make skipping exisitng outputs an argument
    tryCatch(  ## catch any error in the following code, print a message, skipping to next species 
    FIT_MAXENT_SELECTION(occ                     = occurrence, 
                         bg                      = background, 
                         sdm.predictors          = sdm.predictors, 
                         name                    = spp, 
                         outdir                  = 'output/maxent/SEL_VAR',    ## Change the outdir on the new run
                         template.raster,
                         min_n                   = 20,   ## This should be higher...
                         max_bg_size             = 100000,
                         background_buffer_width = 200000,
                         shapefiles              = FALSE,
                         features                = 'lpq',
                         replicates              = 5,
                         cor_thr                 = 0.85, 
                         pct_thr                 = 5, 
                         k_thr                   = 5, 
                         responsecurves          = FALSE), 
    
    ## https://stackoverflow.com/questions/19394886/trycatch-in-r-not-working-properly
    #function(e) message('Species skipped ', spp)) ## skip any species for which the function fails
    error = function(cond) {
      
      message(paste('Species skipped ', spp))
      
    })
    
  } else {
    
    message(spp, ' skipped - no data.')         ## This condition ignores species which have no data...
    
  }  
  
})


stopCluster(cl)





#########################################################################################################################
## OUTSTANDING SDM TASKS:
#########################################################################################################################


## What is being used for background records? How is background selection being done? 
  
## What is the cut-off for correlation and what is the minimum number of variables? 
  
## Which variables are showing up as important

## Further mapping and cleaning of GBIF data needed for the important species: Use worldview?

## Store the list of skipped species...



#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################