#########################################################################################################################
#############################################  FIT MAXENT TO GBIF DATA ################################################## 
#########################################################################################################################


#########################################################################################################################
## This code takes a table of all species occurrences (rows) and environmental values (columns), and runs SDMs for a list
## of taxa. These taxa will be modelled in two sets:

## Species with Australian boundary bias (e.g. Lomandra Longifolia) will be modelled using random selection of background
## points.

## Species without Australian boundary bias (e.g. Ficus Brachypoda) will be modelled using targetted selection of 
## background points.


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
#load("./data/base/HIA_LIST/COMBO/SDM_TEMPLATE_RASTER.RData")
#load('data/hasData_cells.rds')
#load('data/template_hasData.tif')
source('./R/HIA_LIST_MATCHING.R')

## Load SDM data
template.raster = raster("./data/template_hasData.tif")
template.cells  = readRDS("./data/hasData_cells.rds")
SDM.DATA.ALL = readRDS("./data/base/HIA_LIST/COMBO/SDM_DATA_CLEAN_042018.rds")


## Check data :: template, data table and species 
dim(template.raster)
dim(SDM.DATA.ALL)
length(unique(SDM.DATA.ALL$searchTaxon))


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


## Create raster stack of current environmental conditions if needed
i  <- match(sdm.predictors, sdm.predictors)
ff <- file.path('./data/base/worldclim/world/0.5/bio/current',
                sprintf('bio_%02d.tif', i))


## Name the grids :: these should be indentical
env.grids.current = stack(sub('0.5', '1km', ff))
names(env.grids.current) <- sdm.predictors[i]
identical(names(env.grids.current),sdm.predictors)


## Create polygon of land surface
LAND       = readRDS("F:/green_cities_sdm/data/base/CONTEXTUAL/LAND_world.rds")
aus = readRDS("F:/green_cities_sdm/data/base/CONTEXTUAL/aus_states.rds") %>%
  spTransform(CRS('+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'))


#########################################################################################################################
## Loop over the species list and plot the occurrence data for each to check the data bias
TAXA = as.list(sort(unique(intersect(SDM.DATA.ALL$searchTaxon, extra.spp))))
  
for (i in 1:length(TAXA)) {
  
  ## Create points for each species
  spp.points <- SDM.DATA.ALL[SDM.DATA.ALL$searchTaxon == TAXA[i], ] %>%
    spTransform(CRS('+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'))
  
  ## Print to file
  save_name = gsub(' ', '_', TAXA[i])
  save_dir  = "output/maxent/summary"
  png(sprintf('%s/%s_%s.png', save_dir,
              save_name, "Australian_points"),
      3236, 2000, units = 'px', res = 300)
  
  ## set margins
  par(mar   = c(3, 3, 5, 3),  ## b, l, t, r
      #mgp   = c(9.8, 2.5, 0),
      oma   = c(1.5, 1.5, 1.5, 1.5))
  
  ## Plot just the Australian points
  plot(aus, main = TAXA[i])
  points(spp.points, col = "red", cex = .3, pch = 19)
  
  ## Finish the device
  dev.off()
  
}





#########################################################################################################################
## 3). RUN SDMs USING A-PRIORI VARIABLES
#########################################################################################################################


#########################################################################################################################
## Run Maxent using a random selection of background points. 
projection(template.raster);projection(SDM.DATA.ALL)


## Use both targetted selection and random selection


## spp = kop.spp[1]
## Run without cluster
lapply(kop.spp, function(spp) { # for serial, parLapply(cl, species[1:8], function(x) { # for parallel 
  
  ## Print the taxa being processed to screen
  if(spp %in% SDM.DATA.ALL$searchTaxon) {
    message('Doing ', spp) 
    
    ## Subset the records to only the taxa being processed
    occurrence <- subset(SDM.DATA.ALL, searchTaxon == spp)
    
    ## Now get the background points. These can come from any spp, other than the modelled species.
    background <- subset(SDM.DATA.ALL, searchTaxon != spp)
    
    ## The create a vector of the sdm.predictors used. 
    #sdm.predictors <- sdm.select # vector of used sdm.predictors 
    
    ## Finally fit the models using FIT_MAXENT. Use tryCatch to skip any exceptions
    tryCatch(
      FIT_MAXENT_TARG_BG(occ                     = occurrence, 
                         bg                      = background, 
                         sdm.predictors          = sdm.select, 
                         name                    = spp, 
                         outdir                  = 'output/maxent/SET_VAR_COORDCLEAN', 
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
## Run Maxent using a targetted selection of background points


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
    #sdm.predictors <- sdm.select # vector of used sdm.predictors 
    
    ## Finally fit the models using FIT_MAXENT. Use tryCatch to skip any exceptions
    tryCatch(
      FIT_MAXENT_RAND_BG(occ                     = occurrence, 
                         #bg                      = background, 
                         name                    = spp, 
                         outdir                  = 'output/maxent/SET_VAR_COORDCLEAN',
                         sdm.predictors          = sdm.select,
                         env.grids               = env.grids.current,
                         template.raster         = template.raster,
                         template.cells          = template.cells,
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


## Save image
save.image("STEP_7_RUN_SDM.RData")





#########################################################################################################################
## OUTSTANDING SDM TASKS:
#########################################################################################################################


## 1). Create a list of species with boundary bias, and without boundary bias (done). Only model those species with > 20 AUS records.
##     Create list of species with > 20 records.
##     Track the functional coverage of the modelled spp (done).

## 3). Re-process the niches for extra species :: we have niches for 6700 taxa, including all but 200 of the Kachenko taxa (done).

## 4). Get the random background points maxent function working (need John's help to code).
  
## 5). Add the koppen zone constraint to both maxent functions  (need John's help to code).

## 6). Set a minium no. of background points, as well as a maximum (John to advise).

## 7). Summarise all maxent output, check species thresholds :: maxent tables (AIC), predicted maps, occ/bg points, response curves, etc.
##     Choose spp (Hugh, Linda, Rach to review each species, and an indepdendent expert?).

## 8). Get the Koppen summary idea working for species that are not modelled :: where will the koppens with records be in 2030/50/70? 
##     Use Darren Kriticos's grid of changing Koppen with decades (but this only uses two GCMs)

## 9). Model extra species :: take another ~300 spp from the intersection of any grown spp and the risky/innovative species.

## 10). Decide format to present ~600 species to HIA





#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################