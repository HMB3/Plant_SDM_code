#########################################################################################################################
###################### PREDICT MAXENT TO FUTURE CLIMATES AND SUMARISE THE RESULTS ####################################### 
#########################################################################################################################


#########################################################################################################################
## This code takes the output of the Maxent species distribution models using current conditions, and generates 
## a prediction of habitat suitability for current and future environmental conditions. 


## The predictions from all 6 GCMs are then combined into a single habitat suitability layer, and the total area of habitat 
## gained, lost or remaining stable is calculated (i.e. for AUS). 


## Using this combined layer, the area occupied by species within areal units (e.g. significant urban areas or SUAs), 
## under each projection (2030, 2050 and 2070) is also calculated. 





#########################################################################################################################
## 1). PROJECT MAXENT MODELS UNDER MULTIPLE CLIMATE SCEANARIOS FOR 2030, 2050 AND 2070
#########################################################################################################################



#########################################################################################################################
## Use a list of GCM scenarios to create maps of habitat suitability. See CSIRO for details ::
## https://www.climatechangeinaustralia.gov.au/en/support-and-guidance/faqs/eight-climate-models-data/
head(gcms.50) ; head(gcms.70) ; head(gcms.30)


## This function uses MESS maps to constrain the climate suitability surface


#########################################################################################################################
## MESS maps measure the similarity between the new environments, and those in the training sample.
## When model predictions are projected into regions, times or spatial resolutions not analysed in the training data, 
## it may be important to measure the similarity between the new environments and those in the training sample 
## (Elith et al. 2010), as models are not so reliable when predicting outside their domain (Barbosa et al. 2009). 
## The Multivariate Environmental Similarity Surfaces (MESS) analysis measures the similarity in the analysed variables 
## between any given locality in the projection dataset and the localities in the reference (training) dataset 
## (Elith et al. 2010).


## Jonh on MESS :

# MESS describes whether/how novel the predictors' values at a grid cell are, relative to the environmental 
# space captured by the model-fitting (background and occurrence locations) data. You can't use the current 
# climate MESS mask to describe future novelty - you really need to recalculate MESS for each climate you're 
# projecting to (or ensemble thereof). 

# Instead, I think the model predictions should be shown, but with an indication of where these predictions
# are uncertain due to extrapolation to novel climates. Essentially, we're just less confident about model
# behaviour when one or more predictors' values are novel, and that means we're uncertain about both suitability
# and unsuitability. If we just label it as unsuitable, it kind of implies that we're as confident that it's
# unsuitable as we are for non-novel areas that are predicted to be unsuitable. One way to present it is to
# make a map of model predictions and overlay a pattern (e.g. hatching) that indicates where climate is novel.


## The MESS map function needs to be modified so that it can handle species with no novel areas..........................
## Why are these species producing strange results? Could be the multiple lists


#########################################################################################################################
## Create 2030 maps
project_maxent_grids_mess(shp_path      = "./data/base/CONTEXTUAL/", ## Path for shapefile
                          aus_shp       = "aus_states.rds",          ## Shapefile, e.g. Australian states
                          world_shp     = "LAND_world.rds",          ## World shapefile          
                          
                          scen_list     = scen_2030,              ## List of climate scenarios
                          species_list  = map_spp,                   ## List of species folders with maxent models
                          maxent_path   = bs_path,                   ## Output folder
                          climate_path  = "./data/base/worldclim/aus/1km/bio", ## climate data
                          grid_names    = bs.predictors,             ## names of the predictor grids
                          time_slice    = 30,                        ## Time period
                          current_grids = aus.grids.current,         ## predictor grids
                          create_mess   = "TRUE")


#########################################################################################################################
## Create 2050 maps
# env.grids.2050 = project_maxent_grids_mess(shp_path      = "./data/base/CONTEXTUAL/", ## shapefile, e.g. Australia
#                                            shp           = "aus_states.rds",          ## Australian states
#                                            scen_list     = scen_2050,                 ## List of climate scenarios
#                                            species_list  = map_spp,                   ## List of species folders with maxent models
#                                            maxent_path   = maxent_path,               ## Output folder
#                                            climate_path  = "./data/base/worldclim/aus/1km/bio", ## climate data
#                                            grid_names    = grid.names,                ## names of the predictor grids
#                                            time_slice    = 50,                        ## Time period
#                                            current_grids = aus.grids.current,         ## predictor grids
#                                            create_mess   = "TRUE")


#########################################################################################################################
## Create 2070 maps
# env.grids.2070 = project_maxent_grids_mess(shp_path      = "./data/base/CONTEXTUAL/", ## shapefile, e.g. Australia
#                                            shp           = "aus_states.rds",          ## Australian states
#                                            scen_list     = scen_2070,                 ## List of climate scenarios
#                                            species_list  = map_spp,                   ## List of species folders with maxent models
#                                            maxent_path   = maxent_path,               ## Output folder
#                                            climate_path  = "./data/base/worldclim/aus/1km/bio", ## climate data
#                                            grid_names    = grid.names,                ## names of the predictor grids
#                                            time_slice    = 70,                        ## Time period
#                                            current_grids = aus.grids.current,         ## predictor grids
#                                            create_mess   = "TRUE")





#########################################################################################################################
## 2). SUMARIZE MAXENT RESULTS FOR EACH SPECIES ACROSS MULTIPLE GCMs, INSIDE SIGNIFCANT URBAN AREAS
#########################################################################################################################


#########################################################################################################################
## So we can use the values in these columns to threhsold the rasters of habitat suitability (0-1) when combining them.
## For each species, we will create a binary raster with cell values between 0-6. These cell values represent the number of GCMs 
## where that cell had a suitability value above the threshold determined by maxent (e.g. the max training sensitivity + specif).
## We then count the number of cells lost, gained, stable and unchanged, both across Australia, and within SUAs


## This function now needs to use the suitability rasters which have novel environments masked out. So that's :: 
## hs_current_not_novel
## hs_future_not_novel - for each scenario


#########################################################################################################################
## Combine GCM predictions and calculate gain and loss for 2030
## Here we can add the mask of novel environments to SUA aggregation
suitability.2030 = mapply(SUA_cell_count,                                  ## Function aggreagating GCM predictions by spatial unit
                          unit_path     = "./data/base/CONTEXTUAL/SUA/",   ## Data path for the spatial unit of analysis 
                          unit_shp      = "SUA_2016_AUST.rds",             ## Spatial unit of analysis - E.G. SUAs
                          unit_vec      = "SUA_2016_VEC.rds",              ## Vector of rasterized unit cells
                          world_shp     = "LAND_world.rds",                ## Polygon for AUS maps           
                          aus_shp       = "aus_states.rds",                ## Polygon for World maps  
                          
                          DIR_list      = SDM.RESULTS.DIR,                 ## List of directories with rasters
                          species_list  = map_spp,                         ## List of species' directories
                          maxent_path   = bs_path,                         ## Directory of maxent results
                          thresholds    = percent.10.log,                  ## List of maxent thresholds
                          time_slice    = 30,                              ## Time period, eg 2030
                          write_rasters = TRUE)


#########################################################################################################################
## Combine GCM output for 2050 
# suitability.2050 = mapply(SUA_cell_count, 
#                           unit_path     = "./data/base/CONTEXTUAL/SUA/",   ## Data path for the spatial unit of analysis 
#                           unit_file     = "SUA_2016_AUST.rds",             
#                           unit_vec      = "SUA_2016_VEC.rds", 
#                           DIR_list      = SDM.RESULTS.DIR,
#                           species_list  = map_spp,
#                           maxent_path   = maxent_path,
#                           thresholds    = percent.10.log,
#                           percentiles   = percent.10.om,
#                           time_slice    = 50,
#                           write_rasters = TRUE)


#########################################################################################################################
## Combine GCM output for 2070 
# suitability.2070 = mapply(SUA_cell_count,
#                           unit_path     = "./data/base/CONTEXTUAL/SUA/",   ## Data path for the spatial unit of analysis 
#                           unit_file     = "SUA_2016_AUST.rds",             ## Spatial unitt of analysis - E.G. SUAs
#                           unit_vec      = "SUA_2016_VEC.rds", 
#                           DIR_list      = SDM.RESULTS.DIR,
#                           species_list  = map_spp,
#                           maxent_path   = maxent_path,
#                           thresholds    = percent.10.log,
#                           percentiles   = percent.10.om,
#                           time_slice    = 70,
#                           write_rasters = TRUE)





#########################################################################################################################
## OUTSTANDING MAPPING TASKS:


## 1). Figure out how to make the try catch work





#########################################################################################################################
#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################