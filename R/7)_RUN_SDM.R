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


## Alex bush :: can we do that for our species - native/exotic, native/urban range


#########################################################################################################################
## Load packages, functions and data
#load("./data/base/HIA_LIST/COMBO/COMBO_RASTER_CONTEXT.RData")
#load("./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT.RData")
#load("./data/base/HIA_LIST/COMBO/SDM_TEMPLATE_RASTER.RData")
#load('data/hasData_cells.rds')
#load('data/template_hasData.tif')
source('./R/HIA_LIST_MATCHING.R')
 

## Load SDM data :: template rasters, point data and koppen zones
template.raster = raster("./data/template_hasData.tif")
template.cells  = readRDS("./data/hasData_cells.rds")
SDM.DATA.ALL    = readRDS("./data/base/HIA_LIST/COMBO/SDM_DATA_CLEAN_052018.rds")
Koppen_1975     = raster('data/Koppen_1000m_Mollweide54009.tif')
Koppen_zones    = unique(readOGR('data/base/CONTEXTUAL/WC05_1975H_Koppen_Shapefile/WC05_1975H_Koppen_Kriticos_2012.shp')@data[, 1:2])


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


## Create polygon of land surface
LAND       = readRDS("F:/green_cities_sdm/data/base/CONTEXTUAL/LAND_world.rds")
aus = readRDS("F:/green_cities_sdm/data/base/CONTEXTUAL/aus_states.rds") %>%
  spTransform(ALB.CONICAL)


#########################################################################################################################
## Loop over the species list and plot the occurrence data for each to check the data bias
TAXA = as.list(sort(unique(intersect(SDM.DATA.ALL$searchTaxon, spp.mile))))

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
projection(template.raster);projection(SDM.DATA.ALL);projection(Koppen_1975)


## Loop over all the species spp = spp.combo[31]
lapply(spp.combo, function(spp){ 
  
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



#########################################################################################################################
## Run analyses using parallel processing

## To run the function in a cluster, we need to export all the objects which are not defined in the function call
cl <- makeCluster(4)
clusterExport(cl, c('template.raster', 'SDM.DATA.ALL', 'FIT_MAXENT', 'FIT_MAXENT_TARG_BG'))
clusterEvalQ(cl, {
  
  ## These packages are required :: check they are updated
  require(ff)
  require(rgeos)
  require(sp)
  require(raster)
  require(rJava)
  require(dismo)
  require(things)
  
})


## Loop over all the species spp = spp.combo[31]
parLapply(cl, spp.combo, function(spp) { # for serial, lapply(spp.combo, function(spp)
  
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


stopCluster(cl)


#########################################################################################################################
## Look through the output directory for species where the code didn't finish. This is usually species with < 27 files
## in the top "species_name" directory
dd <- list.dirs('H:/green_cities_sdm/output/maxent/SET_VAR_KOPPEN', recursive = F)


## Loop over all the directories in the maxent output folder
ff <- sapply(dd, function(d) {
  
  list.files(d, full.names = TRUE)
  
})


## Now remove all the files in each directory, if that directory has < 27 files in the top level
n <- lengths(ff)

lapply(ff[n < 27], function(x) {
  
  file.remove(x)
  unlink(dirname(x[1]), recursive = TRUE)
  
})





#########################################################################################################################
## OUTSTANDING SDM TASKS:
#########################################################################################################################


## 1). Created a list of species with boundary bias, and without boundary bias (done). Only model those species with > 20 AUS records.
##     Create list of species with > 20 records.

## 3). Re-processed the niches for extra species :: we have niches for 6700 taxa, including all but 200 of the Kachenko taxa (done).

## 4). Don't model species with strong boundary bias :: these models are unreliable, best to get the data from QLD and Victoria.
  
## 5). Add the koppen zone constraint to the maxent functions  (done).


## 6). Summarise all maxent output, check species thresholds :: maxent tables (AIC), predicted maps, occ/bg points, response curves, etc.
##     Choose final species using all output (Hugh, Linda, Rach to review each species, and an indepdendent expert?).
##     Track the functional coverage of the modelled spp

## 7). Get the Koppen summary idea working for species that are not modelled :: where will the koppens with records be in 2030/50/70? 
##     Use Darren Kriticos's grid of changing Koppen with decades (but this only uses two GCMs).

## 8). Model extra species :: take another ~n spp from the intersection of any grown spp and the risky/innovative species.

## 9). Decide format to present ~n species to HIA





#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################