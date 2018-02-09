#########################################################################################################################
#############################################  FIT MAXENT TO GBIF DATA ################################################## 
#########################################################################################################################


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
#load("./data/base/HIA_LIST/COMBO/COMBO_RASTER_CONTEXT.RData")
load("./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT.RData")
load("./data/base/HIA_LIST/COMBO/HIA_SDM_DATA_ALL_VAR.RData")
load("./data/base/HIA_LIST/COMBO/SDM_TEMPLATE_RASTER.RData")

source('./R/HIA_LIST_MATCHING.R')
source('./R/HIA_CLEAN_MATCHING.R')
source('./R/SDM_FUNCTIONS.R')


## Check data 
str(template.raster)
str(SDM.DATA.ALL)


## Check experimental taxa again
'Swainsona formosa'  %in% SDM.DATA.ALL$searchTaxon
'Templetonia retusa' %in% SDM.DATA.ALL$searchTaxon 
'Dodonaea baueri'    %in% SDM.DATA.ALL$searchTaxon 
'Platanus hispanica' %in% SDM.DATA.ALL$searchTaxon 
'Kennedia beckxiana' %in% SDM.DATA.ALL$searchTaxon 
exp.spp  = c('Swainsona formosa', 'Templetonia retusa', 'Dodonaea baueri', 'Platanus hispanica', 'Kennedia beckxiana')  
exp.rev  = sort(exp.spp, decreasing = TRUE)
miss.spp = c('Corymbia tessellaris', 'Metrosideros excelsa')
exp.eg   = c('Ficus brachypoda', 'Flindersia australis', 'Xanthastemon paradoxus')




#########################################################################################################################
## 3). RUN SDM FOR SELECTED VARIABLES
#########################################################################################################################


#########################################################################################################################
## Chose a-priori worldclim predictors
sdm.predictors <- c("Annual_mean_temp",    "Mean_diurnal_range",  "Isothermality",      "Temp_seasonality",  
                    "Max_temp_warm_month", "Min_temp_cold_month", "Temp_annual_range",  
                    #"Mean_temp_wet_qu",    "Mean_temp_dry_qu",    
                    "Mean_temp_warm_qu",   "Mean_temp_cold_qu",  "Annual_precip", 
                    "Precip_wet_month",    "Precip_dry_month",   "Precip_seasonality", "Precip_wet_qu",     
                    "Precip_dry_qu")       
#"Precip_warm_qu",     "Precip_col_qu")



## Are the latest experimental species in there?
head(test.spp, 10)
head(spp.all, 10)


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
cl <- makeCluster(4)
clusterExport(cl, c('template.raster', 'SDM.DATA.ALL', 'FIT_MAXENT_SELECTION'))
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
## Run for all species
lapply(test.spp, function(spp)  { # for serial, parLapply(cl, species[1:8], function(spp) { # for parallel 
  
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
    min_n               = 20
    
    ## Fit the models using FIT_MAXENT. Would be good to make skipping exisitng outputs an argument
    #browser()
    FIT_MAXENT_SELECTION(occ                     = occurrence, 
                         bg                      = background, 
                         sdm.predictors          = sdm.predictors, 
                         name                    = spp, 
                         outdir                  = 'output/maxent/STD_VAR_ALL',    ## Change the outdir on the new run
                         template.raster,
                         min_n                   = 20,   ## This should be higher...
                         max_bg_size             = 100000,
                         background_buffer_width = 200000,
                         shapefiles              = TRUE,
                         features                = 'lpq',
                         replicates              = 5,
                         cor_thr = 0.7, 
                         pct_thr = 5, 
                         k_thr = 2, 
                         responsecurves          = TRUE)
    
  } else {
    
    message(spp, ' skipped - no data.')         ## This condition ignores species which have no data
    
  }  
  
})


stopCluster(cl)





########################################################################################################################
## 100 species takes about 4 hours...
cl <- makeCluster(2)
clusterExport(cl, c('template.raster', 'SDM.DATA.ALL', 'FIT_MAXENT_SELECTION'))
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
## Run for all species
lapply(test.reverse, function(spp)  { # for serial, parLapply(cl, species[1:8], function(spp) { # for parallel 
  
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
    min_n               = 20
    
    ## Fit the models using FIT_MAXENT. Would be good to make skipping exisitng outputs an argument
    #browser()
    FIT_MAXENT_SELECTION(occ                     = occurrence, 
                         bg                      = background, 
                         sdm.predictors          = sdm.predictors, 
                         name                    = spp, 
                         outdir                  = 'output/maxent/STD_VAR_ALL',    ## Change the outdir on the new run
                         template.raster,
                         min_n                   = 20,   ## This should be higher...
                         max_bg_size             = 100000,
                         background_buffer_width = 200000,
                         shapefiles              = TRUE,
                         features                = 'lpq',
                         replicates              = 5,
                         cor_thr = 0.7, 
                         pct_thr = 5, 
                         k_thr = 2, 
                         responsecurves          = TRUE)
    
  } else {
    
    message(spp, ' skipped - no data.')         ## This condition ignores species which have no data
    
  }  
  
})

  
##
stopCluster(cl)


## Now save .RData file for the next session...
save.image("STEP_7_SDM_BACKWARDS_SELECTION.RData")
save.session(file = 'STEP_7.Rda')





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