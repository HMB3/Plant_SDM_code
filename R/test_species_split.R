PET      = raster("./data/base/worldclim/world/1km/pet_he_yr1.tif")
AI       = raster("./data/base/worldclim/world/1km/ai_yr1.tif")
twi      = raster("./data/base/ACLEP/TWI.tif")
tpi      = raster("./data/base/ACLEP/TPI.tif")


## Test species
GBIF.spp = unique(outstanding.spp[1501:length(outstanding.spp)]) 


## 
GBIF.spp      = as_utf8(GBIF.spp, normalize = TRUE)  ## Check the species names have the right characters
save_run      = "EXTRA_2000"                     ## a variable to append the run name to the output files
                                                     ## ALL_EVREGREEN_MAY_2018 is the latest version of the niche data

read_data     = "FALSE"                              ## Read intermediary data between the steps?
save_data     = "TRUE"                               ## Save data?
check_maps    = "FALSE"                              ## Create maps, shapefiles and histograms of each speices?

GBIF_path     = "./data/base/HIA_LIST/GBIF_UPDATE/"  ## The path where GBIF data is stored
ALA_path      = "./data/base/HIA_LIST/ALA_UPDATE/"   ## The path where ALA data is stored place
DATA_path     = "./data/ANALYSIS/"                   ## The path where the final data from/for analyses are stored. 
OCC_SOURCE    = c("ALA", "GBIF", "INVENTORY")
#OCC_SOURCE    = c("ALA", "GBIF")                     ## Data source :: for trees, ALL. For non trees, ALA/GBIF, for occ and bg
calc_niche    = "TRUE"

maxent_path   = './output/maxent/SPLIT_TEST_INV/'      ## The directory where maxent files are saved               
maxent_dir    = 'output/maxent/SPLIT_TEST_INV'         ## Another version of the path needed to run maxent loop

bs_path       = './output/maxent/SPLIT_TEST_INV_BS/'   ## Backwards selection directory
bs_dir        = 'output/maxent/SPLIT_TEST_INV_BS'      ## the directory to harvest results : BS dir?
results_dir   = bs_dir 


#########################################################################################################################
## Now source each step in the workflow. 
message('Running SDMs and maps for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")


## Step 1 :: download GBIF and ALA data
#source('./R/1_GBIF_ALA_DOWNLOAD.R', echo = TRUE)


# Step 2 :: combine GBIF occurrence data with ALA data and filter to records > 1950
source('./R/2_ALA_DATA_FILTER_TAXO_SCIENTIFIC_NAME.R', echo = TRUE)
source('./R/3_GBIF_DATA_TAXO_SCIENTIFIC_NAME.R',     echo = TRUE)


# Step 4 :: combine GBIF, ALA and tree inventory data into a single table,
# clean taxonomy,
# extract environmental condtions
source('./R/4_ALA_GBIF_TAXO_COMBINE.R',   echo = TRUE)
source('./R/INVENTORY_RASTER.R',          echo = TRUE)


# Step 5 :: clean the occurrence data using the 'CleanCoordinates' function in the CoordinateCleaner package, to remove
# records near herbaria, duplicates, etc. & add contextual info for each record (taxonomic and horticultural).
# Then estimate the environmental and geographic range of each species
source('./R/5_GEO_CLEAN_DATA.R',         echo = TRUE)
source('./R/CALCULATE_1KM_NICHES.R',     echo = TRUE)


# Step 6 :: Prepare the SDM table and clean the spatial outliers at 1km resolution.
source('./R/6_PREPARE_SDM_TABLE_1KM.R',  echo = TRUE)


# ## Step 7 :: Run maxent on a table of all species, using targetted background selection, then backwards selection
# ## Step 8 :: Create habitat suitability maps for each species using six GCMs and three time slices (2030/50/70).
#source('./R/7_RUN_MAXENT.R',      echo = TRUE)
#source('./R/8_MAP_SDM_COMBINE.R', echo = TRUE) 



tryCatch(
  project_maxent_grids_mess(shp_path      = "./data/base/CONTEXTUAL/", ## Path for shapefile
                            aus_shp       = "aus_states.rds",          ## Shapefile, e.g. Australian states
                            world_shp     = "LAND_world.rds",          ## World shapefile          
                            
                            scen_list     = scen_2030,                 ## List of climate scenarios
                            species_list  = map_spp[39:42],                   ## List of species folders with maxent models
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
  
  
## Then loop over the species folders and climate scenarios
tryCatch(mapply(SUA_cell_count,                                  ## Function aggreagating GCM predictions by spatial unit
                unit_path     = "./data/base/CONTEXTUAL/SUA/",   ## Data path for the spatial unit of analysis
                unit_shp      = "SUA_2016_AUST.rds",             ## Spatial unit of analysis - E.G. SUAs
                unit_vec      = "SUA_2016_VEC.rds",              ## Vector of rasterized unit cells
                world_shp     = "LAND_world.rds",                ## Polygon for AUS maps
                aus_shp       = "aus_states.rds",                ## Polygon for World maps

                DIR_list      = SDM.RESULTS.DIR[39:42],            ## List of directories with rasters
                species_list  = map_spp[39:42],                    ## List of species' directories
                number_gcms   = 6,                               ## The number of GCMs used (could be determined from object)
                maxent_path   = bs_path,                         ## Directory of maxent results
                thresholds    = percent.10.log[39:42],                  ## List of maxent thresholds
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



