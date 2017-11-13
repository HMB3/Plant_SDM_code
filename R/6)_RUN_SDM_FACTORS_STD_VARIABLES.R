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
source('./R/HIA_LIST_MATCHING.R')
source('./R/HIA_CLEAN_MATCHING.R')
source('./R/GREEN_CITIES_FUNCTIONS.R')
source('./R/SDM_FUNCTIONS.R')

load("./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT.RData")
load("./data/base/HIA_LIST/COMBO/HIA_SDM_DATA_STD_VAR.RData")
load("./data/base/HIA_LIST/COMBO/SDM_TEMPLATE_RASTER.RData")


## Check data 
str(template.raster)
str(SDM.DATA)


## Require packages
p <- c('ff',    'things',         'raster',    'dismo',        'sp',           'latticeExtra', 'data.table', 
       'rgdal', 'rgeos',          'gdalUtils', 'rmaxent',      'readr',        'dplyr',        'tidyr',
       'readr', 'rnaturalearth',  'rasterVis', 'RColorBrewer', 'latticeExtra', 'parallel')
sapply(p, require, character.only = TRUE)



#########################################################################################################################
## 3). RUN SDM FOR SELECTED VARIABLES
#########################################################################################################################


#########################################################################################################################
## Create species subsets for analysis
## All species
spp.all  <- unique(HIA.SPP$Binomial)
str(spp.all)                 ## 605


## The trial species
test.spp = sort(unique(c(renee.full$Species, "Betula pendula", "Fraxinus excelsior", "Quercus robur", "Fagus sylvatica")))
test.spp 


## Combine with 45 from the main list
HIA.SAMPLE = head(COMBO.NICHE.CONTEXT, 53)[, c("searchTaxon")]
test.spp   = sort(unique(c(test.spp, HIA.SAMPLE)))


## Have a look at the breakdown:
TEST.CONTEXT  = COMBO.NICHE.CONTEXT[COMBO.NICHE.CONTEXT$searchTaxon %in% test.spp, ][, c(2:18)]
kable(TEST.CONTEXT)


## Now reverse the order, so we can start another R session from the other end
test.reverse = sort(test.spp, decreasing = TRUE)


## Chose a-priori worldclim predictors
sdm.predictors    <- c("Annual_mean_temp", "Temp_seasonality",    "Max_temp_warm_month", "Min_temp_cold_month",
                       "Annual_precip",    "Precip_seasonality",  "Precip_wet_month",    "Precip_dry_month")


########################################################################################################################
## We can run Maxent from a cluster of cores on the local computer. Here we send (i.e. export) all the necessary ingredients 
## to the cluster. So that's the:


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
clusterExport(cl, c('template.raster', 'SDM.DATA', 'FIT_MAXENT'))
clusterEvalQ(cl, {
  
  require(ff)
  require(rgeos)
  require(sp)
  require(raster)
  require(rJava)
  require(dismo)
  require(things)
  
})


## first check if the outdir exists. Could wrap in a function to create true/false skipping condition
#check.outdir   = 'output/maxent/STD_VAR_ALL'


########################################################################################################################
## Now use 'lapply' to run maxent for multiple species
lapply(spp.all[1:length(spp.all)], function(x) { # for serial, parLapply(cl, species[1:8], function(x) { # for parallel 
  
  ## Print the taxa being processed to screen
  if(x %in% SDM.DATA$searchTaxon) {
    message('Doing ', x)
    
    ## Subset the records to only the taxa being processed
    occurrence <- subset(SDM.DATA, searchTaxon == x)
    
    ## Now get the background points. These can come from anywhere in the whole dataset,
    ## other than the species used.
    background <- subset(SDM.DATA, searchTaxon != x)
    
    ## The create a vector of the sdm.predictors used. 
    ## This should be based on an ecological framework! 
    sdm.predictors <- sdm.predictors # vector of used sdm.predictors 
    
    ## Finally fit the models using FIT_MAXENT
    ## There is no switch in the function to skip outputs that exist.
    ## Given all the changes likely to be made to the models, this could be wise...
    FIT_MAXENT(occ                     = occurrence, 
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


## Potential errors...

## This is because you don't have maxent correctly installed yet 
# Error in .local(x, p, ...) : args not understood:
#   replicates = 5, responsecurves = TRUE, threshold = FALSE, hinge = FALSE

## Ignore this error, doesn't mean anything
# In addition: Warning messages:
#   1: In writeOGR(swd_occ, outdir_sp, "occ_swd", "ESRI Shapefile", overwrite_layer = TRUE) :
#   Field names abbreviated for ESRI Shapefile driver
# 2: In writeOGR(swd_bg, outdir_sp, "bg_swd", "ESRI Shapefile", overwrite_layer = TRUE) :
#   Field names abbreviated for ESRI Shapefile driver


## Look at the output...





#########################################################################################################################
## 4). SDM FOR SELECTED VARIABLES USING CULTIVATED RECORDS
#########################################################################################################################


## Fit maxent function needs to be update to do the variable selection


########################################################################################################################
## We can run Maxent from a cluster of cores on the local computer. Here we send (i.e. export) all the necessary ingredients 
## to the cluster. So that's the:


## template.raster raster, 
## data frame and 
## maxent function


## 100 species takes about 4 hours...
cl <- makeCluster(6)
clusterExport(cl, c('template.raster', 'SDM.DATA', 'FIT_MAXENT'))
clusterEvalQ(cl, {
  
  require(ff)
  require(rgeos)
  require(sp)
  require(raster)
  require(rJava)
  require(dismo)
  require(things)
  
})


## Now use 'lapply' to run maxent for multiple species


## This code would be exactly the same, except for the cultivated/non cultivated output
lapply(test.spp[1:length(test.spp)], function(x) { # for serial, parLapply(cl, species[1:8], function(x) { # for parallel 
  
  ## Print the taxa being processed to screen
  if(x %in% SDM.DATA$searchTaxon) {
    message('Doing ', x)
  
  ## Subset the records to only the taxa being processed
  ## Here, add condition for cultivated
  ## occurrence <- subset(SDM.DATA, searchTaxon == x & CULTIVATED == "CULTIVATED")
  occurrence <- subset(SDM.DATA, searchTaxon == x & CULTIVATED == "CULTIVATED")
  
  ## Now get the background points. These can come from anywhere in the whole dataset,
  ## other than the species used.
  background <- subset(SDM.DATA, searchTaxon != x & CULTIVATED == "CULTIVATED")
  
  ## The create a vector of the sdm.predictors used. 
  ## This should be based on an ecological framework! 
  sdm.predictors <- sdm.predictors # vector of used sdm.predictors 
  
  ## Finally fit the models using FIT_MAXENT
  ## There is no switch in the function to skip outputs that exist.
  ## Given all the changes likely to be made to the models, this could be wise...
  FIT_MAXENT(occ                     = occurrence, 
             bg                      = background, 
             sdm.predictors          = sdm.predictors, 
             name                    = x, 
             outdir                  = 'output/maxent/STD_VAR_CULT', 
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
## 5). RUN SDM FOR SELECTED VARIABLES USING UN-CULTIVATED RECORDS
#########################################################################################################################


## Fit maxent function needs to be updated to do the variable selection


########################################################################################################################
## We can run Maxent from a cluster of cores on the local computer. Here we send (i.e. export) all the necessary ingredients 
## to the cluster. So that's the:


## template.raster raster, 
## data frame and 
## maxent function


## 100 species takes about 4 hours...
cl <- makeCluster(6)
clusterExport(cl, c('template.raster', 'SDM.DATA', 'FIT_MAXENT'))
clusterEvalQ(cl, {
  
  require(ff)
  require(rgeos)
  require(sp)
  require(raster)
  require(rJava)
  require(dismo)
  require(things)
  
})


## Now use 'lapply' to run maxent for multiple species


## This code would be exactly the same, except for the cultivated/non cultivated output
lapply(test.spp[1:length(test.spp)], function(x) { # for serial, parLapply(cl, species[1:8], function(x) { # for parallel 
  
  ## Print the taxa being processed to screen
  if(x %in% SDM.DATA$searchTaxon) {
    message('Doing ', x)
    
    ## Subset the records to only the taxa being processed
    ## Here, add condition for cultivated
    ## occurrence <- subset(SDM.DATA, searchTaxon == x & CULTIVATED == "CULTIVATED")
    occurrence <- subset(SDM.DATA, searchTaxon == x & CULTIVATED == "UNKNOWN")
    
    ## Now get the background points. These can come from anywhere in the whole dataset,
    ## other than the species used.
    background <- subset(SDM.DATA, searchTaxon != x & CULTIVATED == "UNKNOWN")
    
    ## The create a vector of the sdm.predictors used. 
    ## This should be based on an ecological framework! 
    sdm.predictors <- sdm.predictors # vector of used sdm.predictors 
    
    ## Finally fit the models using FIT_MAXENT
    ## There is no switch in the function to skip outputs that exist.
    ## Given all the changes likely to be made to the models, this could be wise...
    FIT_MAXENT(occ                     = occurrence, 
               bg                      = background, 
               sdm.predictors          = sdm.predictors, 
               name                    = x, 
               outdir                  = 'output/maxent/STD_VAR_UNCULT', 
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