#########################################################################################################################
###################### PREDICT MAXENT TO FUTURE CLIMATES AND SUMARISE THE RESULTS ####################################### 
#########################################################################################################################


#########################################################################################################################
## This code takes the output of the Maxent species distribution models (i.e. using current conditions), and generates 
## a prediction of habitat suitability for current and future environmental conditions. 


## The predictions from all 6 GCMs are then combined into a single habitat suitability layer, and the total area of habitat 
## gained, lost or remaining stable is calculated (i.e. for AUS). 


## Using this combined layer, the area occupied by species within areal units (e.g. significant urban areas or SUAs), 
## under each projection (2030, 2050 and 2070) is also calculated. 





#########################################################################################################################
## 1). PROJECT MAXENT MODELS UNDER MULTIPLE CLIMATE SCEANARIOS FOR 2030, 2050 AND 2070
#########################################################################################################################



#########################################################################################################################
## Use a list of GCM scenarios to create maps of habitat suitability 
## Eight of the 40 CMIP5 models assessed in this project have been selected for use in provision of application-ready data. 
## This facilitates efficient exploration of climate projections for Australia.
## https://www.climatechangeinaustralia.gov.au/en/support-and-guidance/faqs/eight-climate-models-data/
head(gcms.50) ; head(gcms.70) ; head(gcms.30)


#########################################################################################################################
## For each species, use a function to create raster files and maps under all six GCMs at each time step


#########################################################################################################################
## MESS maps measure the similarity between the new environments and those in the training sample.
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


#########################################################################################################################
## E.G. arguments to run the algorithm inside the function 
shp_path      = "./data/base/CONTEXTUAL/"   
aus_shp       = "aus_states.rds"
world_shp     = "LAND_world.rds"
x             = scen_2030[3]    
species       = map_spp[1]   ## "Eucalyptus_camaldulensis" it breaks on Euc. camuldulensis
maxent_path   = maxent_path  
climate_path  = "./data/base/worldclim/aus/1km/bio" 
grid_names    = grid.names   
current_grids = aus.grids.current
time_slice    = 30
create_mess   = "TRUE"
MESS_folder   = "MESS_output"


## Reverse sort the list because it is very slow.........
map_spp_rev  = sort(map_spp, decreasing = TRUE) 


## The MESS map function needs to be modified so that it can handle species with no novel areas..........................
## 
## Also add a try catch which writes a failed file as per step 7.........................................................


#########################################################################################################################
## Create 2030 maps
lapply(map_spp, function(spp){ 
  
  ## Create a file for the species being mapped
  outdir <- maxent_path
  dir_name = file.path(maxent_path, gsub(' ', '_', spp))
  write.csv(data.frame(), file.path(dir_name, "map_in_progress.txt"))
  
  ## Print the taxa being processed to screen
  if(spp %in% SDM.SPAT.OCC.BG$searchTaxon) {
    message('Doing ', spp) 
    
    tryCatch(project_maxent_grids_mess(shp_path      = "./data/base/CONTEXTUAL/", ## Path for shapefile
                                                        aus_shp       = "aus_states.rds",          ## Shapefile, e.g. Australian states
                                                        world_shp     = "LAND_world.rds",          ## Polygon for AUS maps           
                                                        
                                                        scen_list     = scen_2030[1],                 ## List of climate scenarios
                                                        species_list  = map_spp,                   ## List of species folders with maxent models
                                                        maxent_path   = maxent_path,               ## Output folder
                                                        climate_path  = "./data/base/worldclim/aus/1km/bio", ## climate data
                                                        grid_names    = grid.names,                ## names of the predictor grids
                                                        time_slice    = 30,                        ## Time period
                                                        current_grids = aus.grids.current,         ## predictor grids
                                                        create_mess   = "TRUE"),  ## write mess maps?
                              
                              error = function(cond) {
                                
                                message(spp, ' failed')  
                                write.csv(data.frame(), file.path(dir_name, "failed.txt"))
                                
                              })
    
  } else {
    
    message(spp, ' skipped - no data.')         ## This condition ignores species which have no data...
    write.csv(data.frame(), file.path(dir_name, "completed.txt"))
    
  }  
  
  ## now add a file to the dir to denote that it has completed
  write.csv(data.frame(), file.path(dir_name, "completed.txt"))
  
})


#########################################################################################################################
## Create 2050 maps
env.grids.2050 = tryCatch(project_maxent_grids_mess(shp_path      = "./data/base/CONTEXTUAL/", ## shapefile, e.g. Australia
                                                    shp           = "aus_states.rds",          ## Australian states
                                                    scen_list     = scen_2050,                 ## List of climate scenarios
                                                    species_list  = map_spp,                   ## List of species folders with maxent models
                                                    maxent_path   = maxent_path,               ## Output folder
                                                    climate_path  = "./data/base/worldclim/aus/1km/bio", ## climate data
                                                    grid_names    = grid.names,                ## names of the predictor grids
                                                    time_slice    = 50,                        ## Time period
                                                    current_grids = aus.grids.current,         ## predictor grids
                                                    create_mess   = "TRUE"),  ## write mess maps?
                          
                          ## Skip species
                          error = function(cond) {
                            
                            message(paste('Species skipped - check', spp))
                            
                          })



#########################################################################################################################
## Create 2070 maps
env.grids.2070 = tryCatch(project_maxent_grids_mess(shp_path      = "./data/base/CONTEXTUAL/", ## shapefile, e.g. Australia
                                                    shp           = "aus_states.rds",          ## Australian states
                                                    scen_list     = scen_2070,                 ## List of climate scenarios
                                                    species_list  = map_spp,                   ## List of species folders with maxent models
                                                    maxent_path   = maxent_path,               ## Output folder
                                                    climate_path  = "./data/base/worldclim/aus/1km/bio", ## climate data
                                                    grid_names    = grid.names,                ## names of the predictor grids
                                                    time_slice    = 70,                        ## Time period
                                                    current_grids = aus.grids.current,         ## predictor grids
                                                    create_mess   = "TRUE"),  ## write mess maps?
                          
                          ## Skip species
                          error = function(cond) {
                            
                            message(paste('Species skipped - check', spp))
                            
                          })






#########################################################################################################################
## 2). SUMARIZE MAXENT RESULTS FOR EACH SPECIES ACROSS MULTIPLE GCMs, INSIDE SIGNIFCANT URBAN AREAS
#########################################################################################################################


#########################################################################################################################
## So we can use the values in these columns to threhsold the rasters of habitat suitability (0-1) when combining them.
## For each species, we will create a binary raster with cell values between 0-6. These cell values represent the number of GCMs 
## where that cell had a suitability value above the threshold determined by maxent (e.g. the max training sensitivity + specif).
## We then count the number of cells lost, gained, stable and unchanged, both across Australia, and within SUAs


## Test problematic species by entering the values and running the function manually
unit_path     = "./data/base/CONTEXTUAL/"   ## Data path for the spatial unit of analysis
unit_file     = "./SUA/SUA_2016_AUST.rds"   ## Spatial unit of analysis - E.G. SUAs
unit_vec      = "./SUA/SUA_2016_VEC.rds"    ## Vector of rasterized unit cells
world_shp     = "LAND_world.rds"                              
aus_shp       = "aus_states.rds"

DIR           = SDM.RESULTS.DIR[1]              ## List of directories with rasters
species       = map_spp[1]                      ## List of species' directories
maxent_path   = maxent_path                     ## Directory of maxent results
thresh        = percent.10.log[1]               ## List of maxent thresholds
time_slice    = 30                              ## Time period, eg 2030
write_rasters = "TRUE"                          ## Save the rasters?



#########################################################################################################################
## Reverse sort the lists
SDM.DIR.REV     = sort(SDM.RESULTS.DIR, decreasing = TRUE)
map_spp_rev     = sort(map_spp,         decreasing = TRUE) 
percent.log.rev = percent.10.log[sort(order(percent.10.log), decreasing = TRUE)]
percent.om.rev  = percent.10.om[sort(order(percent.10.om),   decreasing = TRUE)]


## Check the length matches - order should be correct, as well as the length
length(SDM.RESULTS.DIR);length(map_spp);length(percent.10.log);length(percent.10.om)
length(SDM.DIR.REV);length(map_spp_rev);length(percent.log.rev);length(percent.om.rev)


#########################################################################################################################
## Combine GCM predictions and calculate gain and loss for 2030
## Here we can add the mask of novel environments to SUA aggregation
suitability.2030 = tryCatch(mapply(SUA_cell_count,                                  ## Function aggreagating GCM predictions by spatial unit
                                   unit_path     = "./data/base/CONTEXTUAL/SUA/",   ## Data path for the spatial unit of analysis 
                                   unit_shp      = "SUA_2016_AUST.rds",             ## Spatial unit of analysis - E.G. SUAs
                                   unit_vec      = "SUA_2016_VEC.rds",              ## Vector of rasterized unit cells
                                   world_shp     = "LAND_world.rds",                ## Polygon for AUS maps           
                                   aus_shp       = "aus_states.rds",                ## Polygon for World maps  
                                   
                                   DIR_list      = SDM.RESULTS.DIR,                 ## List of directories with rasters
                                   species_list  = map_spp,                         ## List of species' directories
                                   maxent_path   = maxent_path,                     ## Directory of maxent results
                                   thresholds    = percent.10.log,                  ## List of maxent thresholds
                                   time_slice    = 30,                              ## Time period, eg 2030
                                   write_rasters = TRUE),                           ## Save the combined rasters?
                            
                            error = function(cond) {
                              
                              message(paste('Species skipped - check inputs', spp))
                              
                            })


#########################################################################################################################
## Combine GCM output for 2050 
suitability.2050 = tryCatch(mapply(SUA_cell_count, 
                                   unit_path     = "./data/base/CONTEXTUAL/SUA/",   ## Data path for the spatial unit of analysis 
                                   unit_file     = "SUA_2016_AUST.rds",             
                                   unit_vec      = "SUA_2016_VEC.rds", 
                                   DIR_list      = SDM.RESULTS.DIR,
                                   species_list  = map_spp,
                                   maxent_path   = maxent_path,
                                   thresholds    = percent.10.log,
                                   percentiles   = percent.10.om,
                                   time_slice    = 50,
                                   write_rasters = TRUE),
                            
                            error = function(cond) {
                              
                              message(paste('Species skipped - check inputs', spp))
                              
                            })


#########################################################################################################################
## Combine GCM output for 2070 
suitability.2070 = tryCatch(mapply(SUA_cell_count,
                                   unit_path     = "./data/base/CONTEXTUAL/SUA/",   ## Data path for the spatial unit of analysis 
                                   unit_file     = "SUA_2016_AUST.rds",             ## Spatial unitt of analysis - E.G. SUAs
                                   unit_vec      = "SUA_2016_VEC.rds", 
                                   DIR_list      = SDM.RESULTS.DIR,
                                   species_list  = map_spp,
                                   maxent_path   = maxent_path,
                                   thresholds    = percent.10.log,
                                   percentiles   = percent.10.om,
                                   time_slice    = 70,
                                   write_rasters = TRUE),
                            
                            error = function(cond) {
                              
                              message(paste('Species skipped - check inputs', spp))
                              
                            })


#########################################################################################################################
## Reverse the order of the lists to speed up the output
# mapply(SUA_cell_count,
#        unit_path     = "./data/base/CONTEXTUAL/SUA/",   ## Data path for the spatial unit of analysis 
#        unit_file     = "SUA_2016_AUST.rds",             ## Spatial unitt of analysis - E.G. SUAs
#        unit_vec      = "SUA_2016_VEC.rds", 
#        DIR_list      = SDM.DIR.REV,
#        species_list  = map_spp_rev,
#        maxent_path   = maxent_path,
#        thresholds    = percent.log.rev,
#        percentiles   = percent.om.rev,
#        time_slice    = 70,  ## 50, 70
#        write_rasters = FALSE)







#########################################################################################################################
## 3). MOVE MAXENT FOLDERS TO NEW LOCATION FOR EACH RUN
#########################################################################################################################


# ## Create a list of folders for this run of species:EG hollow bearing species 
# run_path         <- "./output/maxent/HOLLOW_SPP"
# run_pat          <- map_spp
# run_pat          <- paste(run_pat, sep = "", collapse = "|")
# maxent_run_list  <- list.files(maxent_path, pattern = run_pat, full.names = TRUE)


# ## Save the list of maxent directories to file
# lapply(maxent_run_list, write, sprintf("%s%s_maxent_folders.txt", maxent_path, save_run), append = TRUE)
# maxent_file_list = sprintf("%s_maxent_folders.txt", save_run)
# 
# 
# ## cd H:/green_cities_sdm/output/maxent/SUA_TREES_ANALYSIS
# ## cat maxent_file_list | xargs -J % cp % run_path
# ## cat HOLLOW_SPP_maxent_folders.txt | xargs -J % cp % HOLLOW_SPP
# 
# 
# ## Copy or move these files to a specific folder
# ## Then you search for a file pattern in that directory
# ## This is very slow, would be better done in unix, etc. But 
# file.copy(from      = maxent_run_list, 
#           to        = run_path, 
#           overwrite = FALSE, 
#           recursive = TRUE, 
#           copy.mode = TRUE)


#########################################################################################################################
## OUTSTANDING MAPPING TASKS:
#########################################################################################################################


## What can be done with the MESS maps? Can they be used to create novelty masks for every GCM, then apply these masks 
## to the combine step? I.e. remove the novel predictions.............................................................  


## Fix the species mapping and combine loops over two lists at once : mapply




#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################