#########################################################################################################################
###################### PREDICT MAXENT TO FUTURE CLIMATES AND SUMARISE THE RESULTS ####################################### 
#########################################################################################################################


#########################################################################################################################
## This code takes the output of the Maxent species distribution models using current conditions, and generates 
## a prediction of habitat suitability for current and future environmental conditions. 


## The predictions from all 6 GCMs are then combined into a single habitat suitability layer, and the total area of habitat 
## gained, lost or remaining stable is calculated (i.e. for AUS). 


## Using this combined layer, the area occupied by species within areal units (e.g. significant urban areas or postal areas), 
## under each projection (2030, 2050 and 2070) is also calculated. 





#########################################################################################################################
## 1). PROJECT MAXENT MODELS UNDER MULTIPLE CLIMATE SCEANARIOS FOR 2030, 2050 AND 2070
#########################################################################################################################



#########################################################################################################################
## Use a list of GCM scenarios to create maps of habitat suitability. See CSIRO for details ::
## https://www.climatechangeinaustralia.gov.au/en/support-and-guidance/faqs/eight-climate-models-data/
head(gcms.50) ; head(gcms.70) ; head(gcms.30)


#########################################################################################################################
## MESS maps measure the similarity between the new environments, and those in the training sample.
## When model predictions are projected into regions, times or spatial resolutions not analysed in the training data, 
## it may be important to measure the similarity between the new environments and those in the training sample 
## (Elith et al. 2010), as models are not so reliable when predicting outside their domain (Barbosa et al. 2009). 
## The Multivariate Environmental Similarity Surfaces (MESS) analysis measures the similarity in the analysed variables 
## between any given locality in the projection dataset and the localities in the reference (training) dataset 
## (Elith et al. 2010).


## This error needs de-bugging in project_maxent_grids_mess......................................

## Converting raster MESS maps to polygons under mc85bi30 scenario for Moringa_oleifera
# H:\green_cities_sdm
# Warning messages:
#   1: In ogrFIDs(dsn = dsn, layer = layer) : no features found
# 2: In value[[3L]](cond) : no features found

## Also, probably best to remove the plotting from the mapping function, and create a separate function that does this
## It was tricky to get the plotting working on Linux (HPC) and Windows, so maybe it's best to create the files on the
## HPC, and then create the plots locally.


#########################################################################################################################
## Create 2030 maps :: can the try catch be looped over the top
tryCatch(
  project_maxent_grids_mess(shp_path      = "./data/base/CONTEXTUAL/", ## Path for shapefile
                            aus_shp       = "aus_states.rds",          ## Shapefile, e.g. Australian states
                            world_shp     = "LAND_world.rds",          ## World shapefile          
                            
                            scen_list     = scen_2030,                 ## List of climate scenarios
                            species_list  = map_spp,                   ## List of species folders with maxent models
                            maxent_path   = bs_path,                   ## Output folder
                            climate_path  = "./data/base/worldclim/aus/1km/bio", ## climate data
                            
                            grid_names    = grid.names,                ## This must include the soil variables, 
                            time_slice    = 30,                        ## Time period
                            current_grids = aus.grids.current,         ## predictor grids - this must include soil variables too
                            create_mess   = "TRUE",
                            nclust        = 1),
  
  ## If the species fails, write a fail message to file. 
  error = function(cond) {
    
    ## This will write the error message inside the text file, 
    ## but it won't include the species
    file.create(file.path(bs_path, "mapping_failed_2030.txt"))
    cat(cond$message, file=file.path(bs_path, "mapping_failed_2030.txt"))
    warning(cond$message)
    
  })




#########################################################################################################################
## Create 2050 maps
tryCatch(
  project_maxent_grids_mess(shp_path      = "./data/base/CONTEXTUAL/", ## Path for shapefile
                            aus_shp       = "aus_states.rds",          ## Shapefile, e.g. Australian states
                            world_shp     = "LAND_world.rds",          ## World shapefile
                            
                            scen_list     = scen_2050,                 ## List of climate scenarios
                            species_list  = map_spp,                   ## List of species folders with maxent models
                            maxent_path   = bs_path,                   ## Output folder
                            climate_path  = "./data/base/worldclim/aus/1km/bio", ## climate data
                            
                            grid_names    = grid.names,
                            time_slice    = 50,                        ## Time period
                            current_grids = aus.grids.current,         ## predictor grids
                            create_mess   = "TRUE",
                            nclust        = 1),
  
  ## If the species fails, write a fail message to file.
  error = function(cond) {
    
    ## This will write the error message inside the text file, 
    ## but it won't include the species
    file.create(file.path(bs_path, "mapping_failed_2050.txt"))
    cat(cond$message, file=file.path(bs_path, "mapping_failed_2050.txt"))
    warning(cond$message)
    
  })


#########################################################################################################################
## Create 2070 maps
## Then, loop over the species folders and climate scenarios
tryCatch(
  project_maxent_grids_mess(shp_path      = "./data/base/CONTEXTUAL/", ## Path for shapefile
                            aus_shp       = "aus_states.rds",          ## Shapefile, e.g. Australian states
                            world_shp     = "LAND_world.rds",          ## World shapefile
                            
                            scen_list     = scen_2070,                 ## List of climate scenarios
                            species_list  = map_spp,                   ## List of species folders with maxent models
                            maxent_path   = bs_path,                   ## Output folder
                            climate_path  = "./data/base/worldclim/aus/1km/bio", ## climate data
                            
                            grid_names    = grid.names,
                            time_slice    = 70,                        ## Time period
                            current_grids = aus.grids.current,         ## predictor grids
                            create_mess   = "TRUE",
                            nclust        = 1),
  
  ## If the species fails, write a fail message to file.
  error = function(cond) {
    
    ## This will write the error message inside the text file, 
    ## but it won't include the species
    file.create(file.path(bs_path, "mapping_failed_2070.txt"))
    cat(cond$message, file=file.path(bs_path, "mapping_failed_2070.txt"))
    warning(cond$message)
    
  })





#########################################################################################################################
## 2). SUMARIZE MAXENT RESULTS FOR EACH SPECIES ACROSS MULTIPLE GCMs, INSIDE SPATIAL UNITS (E.G SUAs)
#########################################################################################################################


#########################################################################################################################
## So now use the 10th% Logistic threshold for each species from the maxent models to threhsold the rasters of habitat 
## suitability (0-1) For each GCM. For each species, summ the 6 GCMS to create a binary raster with cell values between 0-6. 
## These cell values represent the number of GCMs where that cell had a suitability value above the threshold determined 
## by maxent. We classify a cell has suitable if it met the threshold in > 4 GCMs, and use this combined raster to compare 
## current and future suitability, measuring if the suitability of each cell is changing over time, remaining stable or was 
## never suitable.


#########################################################################################################################
## Combine GCM predictions and calculate gain and loss for 2030
## Here we can add the mask of novel environments to SUA aggregation


## Also, add ncores to this function....................................................................................


## Then loop over the species folders and climate scenarios
tryCatch(mapply(SUA_cell_count,                                  ## Function aggreagating GCM predictions by spatial unit
                unit_path     = "./data/base/CONTEXTUAL/SUA/",   ## Data path for the spatial unit of analysis
                unit_shp      = "SUA_2016_AUST.rds",             ## Spatial unit of analysis - E.G. SUAs
                unit_vec      = "SUA_2016_VEC.rds",              ## Vector of rasterized unit cells
                world_shp     = "LAND_world.rds",                ## Polygon for AUS maps
                aus_shp       = "aus_states.rds",                ## Polygon for World maps

                DIR_list      = SDM.RESULTS.DIR,                 ## List of directories with rasters
                species_list  = map_spp,                         ## List of species' directories
                number_gcms   = 6,                               ## The number of GCMs used (could be determined from object)
                maxent_path   = bs_path,                         ## Directory of maxent results
                thresholds    = percent.10.log,                  ## List of maxent thresholds
                time_slice    = 30,                              ## Time period, eg 2030
                write_rasters = TRUE),

         ## If the species fails, write a fail message to file.
         error = function(cond) {

           ## This will write the error message inside the text file,
           ## but it won't include the species
           file.create(file.path(bs_path, "sua_count_failed_2030.txt"))
           cat(cond$message, file=file.path(bs_path, "sua_count_failed_2030.txt"))
           warning(cond$message)

         })




#########################################################################################################################
## Combine GCM output for 2050
tryCatch(mapply(SUA_cell_count,                                  ## Function aggreagating GCM predictions by spatial unit
                unit_path     = "./data/base/CONTEXTUAL/SUA/",   ## Data path for the spatial unit of analysis
                unit_shp      = "SUA_2016_AUST.rds",             ## Spatial unit of analysis - E.G. SUAs
                unit_vec      = "SUA_2016_VEC.rds",              ## Vector of rasterized unit cells
                world_shp     = "LAND_world.rds",                ## Polygon for AUS maps
                aus_shp       = "aus_states.rds",                ## Polygon for World maps

                DIR_list      = SDM.RESULTS.DIR,                 ## List of directories with rasters
                species_list  = map_spp,                         ## List of species' directories
                number_gcms   = 6,                               ## The number of GCMs used (could be determined from object)
                maxent_path   = bs_path,                         ## Directory of maxent results
                thresholds    = percent.10.log,                  ## List of maxent thresholds
                time_slice    = 50,                              ## Time period, eg 2030
                write_rasters = TRUE),


         ## If the species fails, write a fail message to file.
         error = function(cond) {

           ## This will write the error message inside the text file,
           ## but it won't include the species
           file.create(file.path(bs_path, "sua_count_failed_2050.txt"))
           cat(cond$message, file=file.path(bs_path, "sua_count_failed_2050.txt"))
           warning(cond$message)

         })






#########################################################################################################################
## Combine GCM output for 2070
## Then loop over the species folders and climate scenarios
tryCatch(mapply(SUA_cell_count,                                  ## Function aggreagating GCM predictions by spatial unit
                unit_path     = "./data/base/CONTEXTUAL/SUA/",   ## Data path for the spatial unit of analysis
                unit_shp      = "SUA_2016_AUST.rds",             ## Spatial unit of analysis - E.G. SUAs
                unit_vec      = "SUA_2016_VEC.rds",              ## Vector of rasterized unit cells
                world_shp     = "LAND_world.rds",                ## Polygon for AUS maps
                aus_shp       = "aus_states.rds",                ## Polygon for World maps

                DIR_list      = SDM.RESULTS.DIR,                 ## List of directories with rasters
                species_list  = map_spp,                         ## List of species' directories
                number_gcms   = 6,                               ## The number of GCMs used (could be determined from object)
                maxent_path   = bs_path,                         ## Directory of maxent results
                thresholds    = percent.10.log,                  ## List of maxent thresholds
                time_slice    = 70,                              ## Time period, eg 2030
                write_rasters = TRUE),

         ## If the species fails, write a fail message to file.
         error = function(cond) {

           ## This will write the error message inside the text file,
           ## but it won't include the species
           file.create(file.path(bs_path, "sua_count_failed_2070.txt"))
           cat(cond$message, file=file.path(bs_path, "sua_count_failed_2070.txt"))
           warning(cond$message)

         })





#########################################################################################################################
## OUTSTANDING MAPPING TASKS:


## 1). Will we aggregate by SUA, or postcode? The aggregation step could be run in a separate script, rather than
##     than inside the function





#########################################################################################################################
#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################