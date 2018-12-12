#########################################################################################################################
#######################################  TRY RUNNING ALL SDM CODE AT ONCE ############################################### 
#########################################################################################################################


## This code runs the whole SDM workflow for the HIA project, for a subset of species (e.g. whichever you supply)
## In order to run this from katana, we need to use only the necessary packages. 
on_windows = switch(Sys.info()[['sysname']], Windows = TRUE, FALSE)


#########################################################################################################################
## Read in all data to run the SDM code :: species lists, shapefile, rasters & tables
#source('./R/HIA_TREE_LIST.R')
if (on_windows) {
  
    load("H:/green_cities_sdm/TEST_RUN.RData")
  
} else {
  
    cmd_args <- commandArgs(TRUE)
    rdata_file = cmd_args[1]
    if (is.na(rdata_file)) {rdata_file = "../TEST_RUN.RData"}
    message (rdata_file)
    load(rdata_file)
    
}


#########################################################################################################################
## 1). SETUP
#########################################################################################################################


#########################################################################################################################
## Load only the packages needed for the analysis
p <- c('ff',    'things',    'raster',        'dismo',        'sp',           'latticeExtra', 'data.table',
       'rgdal', 'rgeos',     'gdalUtils',     'rmaxent',      'readr',        'plyr',         'dplyr',
       'tidyr', 'readr',     'rnaturalearth', 'rasterVis',    'RColorBrewer', 'latticeExtra', 'parallel',
       'ALA4R', 'stringr',   'Taxonstand',    'CoordinateCleaner', 'gsubfn',  'PerformanceAnalytics',
       'rvest', 'magrittr',  'devtools',      'ggplot2',      'reshape2',     'rmarkdown', 'flexdashboard', 'shiny', 'rgbif',
       'ENMeval', 'tibble',  'ncdf4',         'Cairo', 'taxonlookup', 'kgc', 'maptools', 'DataCombine', 'mgcv', 'rsq')


## Require packages
sapply(p, require, character.only = TRUE)
source_gist('26e8091f082f2b3dd279',             filename = 'polygonizer.R')
source_gist('c6a1cb61b8b6616143538950e6ec34aa', filename = 'hatch.R')


## Source functions, and set temporary directory (for both raster files and just generally)
source('./R/GREEN_CITIES_FUNCTIONS.R')
source('./R/MAXENT_FUNCTIONS.R')
source('./R/MAPPING_FUNCTIONS.R')
source('./R/HIA_CLEAN_MATCHING.R')
rasterOptions(tmpdir = file.path('./RTEMP'))


#########################################################################################################################
##  If on Katana, override the raster objects so we use the local versions
if (on_windows == "FALSE") {

message ("Loading world raster stack")
world.grids.current = stack(
  file.path('./data/base/worldclim/world/0.5/bio/current',
            sprintf('bio_%02d', 1:19)))

## Create a raster stack of current Australian environmental conditions
## This is used for exrtacting worldclim data for the tree inventories
## And running SDMs
message ("Loading inventory raster stack")
inventory.grids.current = stack(
  file.path('./data/base/worldclim/aus/1km/bio/current/WGS/',
            sprintf('bio_%02d.tif', 1:19)))


## Create a raster stack of current Australian environmental conditions
## This is used for projecting SDMs onto Australia
message ("Loading aus raster stack ./data/base/worldclim/aus/1km/bio/current")
aus.grids.current <- stack(
  file.path('./data/base/worldclim/aus/1km/bio/current',   ## ./data/base/worldclim/aus/1km/bio
            sprintf('bio_%02d.tif', 1:19)))


## Divid the current Australian bioclim values by 10
## The future values are divided nby 10 within the mapping function
message ("Iterating over aus current grids")
for(i in 1:11) {
  
  ## simple loop
  message(i)
  aus.grids.current[[i]] <- aus.grids.current[[i]]/10
  
}


#########################################################################################################################
##  Also include the template data
template.raster = raster("./data/template_hasData.tif")
template.cells  = readRDS("./data/hasData_cells.rds")
Koppen_1975     = raster('data/Koppen_1000m_Mollweide54009.tif')
PET             = raster("./data/base/worldclim/world/1km/pet_he_yr1.tif")


## If running locally on the windows machine, just use the data from the Rdata file
} else {
  
  message ("Raster data already loaded in Rdata object")
  
}


#########################################################################################################################
## Now set global analysis variables
GBIF.spp      = WPW.spp
GBIF.spp.rev  = sort(GBIF.spp, decreasing = TRUE)                ## the list reversed - only needed for a big list

save_run      = "WPW_TEST"                                       ## a variable to append the run name to the output files
map_spp_list  = gsub(" ", "_", GBIF.spp)                         ## species list with "_" for mapping
map_spp_rev   = sort(map_spp_list, decreasing = TRUE)            ## reversed, so we can run two at once

GBIF_path     = "./data/base/HIA_LIST/GBIF/"                     ## The path where GBIF data is stored
ALA_path      = "./data/base/HIA_LIST/ALA/"                      ## The path where ALA data is stored  place

maxent_path   = './output/maxent/SUA_TREES_ANALYSIS/'            ## The directory where files are saved               
maxent_dir    = 'output/maxent/SUA_TREES_ANALYSIS'               ## Another version of the path that John's coding needs to run a loop
save_data     = 'TRUE'                                           ## Arguments for saving the intermediary output - i.e. niches
read_data     = 'FALSE'                                          ## Leave these the same - saves data, but doesn't read back in
save_path     = 'data/base/HIA_LIST/COMBO'





#########################################################################################################################
## 2). RUN SDM WORKFLOW
#########################################################################################################################


## Run the SDM workflow for the species set


#########################################################################################################################
## Now source each step in the workflow. 
message('Running SDM workflow for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")


## Step 3 :: combine GBIF occurrence data with ALA data and filter to records > 1950
source('./R/1_GBIF_ALA_DOWNLOAD.R', echo = TRUE)


## Step 3 :: combine GBIF occurrence data with ALA data and filter to records > 1950
source('./R/ALA_DATA_FILTER_TAXO_SCIENTIFIC_NAME.R', echo = TRUE)
source('./R/3_GBIF_DATA_TAXO_SCIENTIFIC_NAME.R',     echo = TRUE)


## Step 4 :: combine GBIF, ALA and urban occurrence data into a single table, extract environmental condtions
## INVENTORY_RASTER will be a problem for species that are not in Alessandro's dataset.
## Also, consider constructing the niche dataset separately to SDMs, so we can model one species at a time
## This is dealt with by save_data     = 'FALSE'
source('./R/4_ALA_GBIF_URBAN_COMBINE.R',  echo = TRUE)
source('./R/INVENTORY_RASTER.R',          echo = TRUE)


## Step 5 :: clean the occurrence data using the 'CleanCoordinates' function in the CoordinateCleaner package to remove
## records near herbaria, duplicates, etc. & add contextual info for each record (taxonomic and horticultural) 
## Then prepare the SDM table
## Then clean the spatial outliers
source('./R/5_GBIF_ALA_CLEAN_NICHES.R',  echo = TRUE)
source('./R/6_PREPARE_SDM_TABLE_1KM.R',  echo = TRUE)


## Step 7 :: Run maxent on a table of all species
source('./R/7_RUN_MAXENT.R', echo = TRUE)


## Step 8 :: Create habitat suitability maps for each species using six GCMs and three time slices (2030/50/70). 
## Then summarise maxent results and estimate species presences in significant urban areas under climate change
## Takes awhile, so probably run different time slices (2030, 2050, 2070) in separate R sessions
source('./R/8_MAP_SDM_COMBINE.R', echo = TRUE)





#########################################################################################################################
## 3). COMBINE THE NICHE RUNS TOGETHER
#########################################################################################################################


## Create a list of all dataframes with the extension from this run
# COMBO.NICHE.list  = list.files(save_path, pattern = 'COMBO_NICHE_CONTEXT_EVERGREEN',  full.names = TRUE, recursive = TRUE)
# COMBO.RASTER.list = list.files(save_path, pattern = 'COMBO_RASTER_CONTEXT_EVERGREEN', full.names = TRUE, recursive = TRUE)


# #########################################################################################################################
# ## Now combine the niche tables for each species into one table
# COMBO.NICHE.ALL <- COMBO.NICHE.list %>%
# 
#   ## pipe the list into lapply
#   lapply(function(x) {
# 
#     ## create the character string
#     f <- paste0(x)
# 
#     ## load each .csv file
#     d <- readRDS(f)
#     d
# 
#   }) %>%
# 
#   ## finally, bind all the rows together
#   bind_rows
# 
# 
# ## Update this
# str(COMBO.NICHE.ALL)
# dim(COMBO.NICHE.ALL)
# 
# 
# ## Make sure the Species are unique
# COMBO.NICHE.ALL = COMBO.NICHE.ALL[!duplicated(COMBO.NICHE.ALL[,c('searchTaxon')]),]
# dim(COMBO.NICHE.ALL)
# length(unique(COMBO.NICHE.ALL$searchTaxon))
# length(setdiff(Manuel$Species, COMBO.NICHE.ALL$searchTaxon))
# 
# 
# #########################################################################################################################
# ## Now combine the raster tables for each species into one table
# COMBO.RASTER.ALL <- COMBO.RASTER.list %>%
# 
#   ## pipe the list into lapply
#   lapply(function(x) {
# 
#     ## create the character string
#     f <- paste0(x)
# 
#     ## load each .csv file
#     d <- readRDS(f)
#     d
# 
#   }) %>%
# 
#   ## finally, bind all the rows together
#   bind_rows
# 
# 
# ## This is a summary of maxent output for current conditions
# dim(COMBO.RASTER.ALL)
# names(COMBO.RASTER.ALL)[1:10]
# 
# 
# #########################################################################################################################
# ## Save the niche and raster data
# saveRDS(COMBO.NICHE.ALL,  paste0('data/base/HIA_LIST/COMBO/COMBO_NICHE_ALL_',  save_run, '.rds'))
# saveRDS(COMBO.RASTER.ALL, paste0('data/base/HIA_LIST/COMBO/COMBO_RASTER_ALL_', save_run, '.rds'))





#########################################################################################################################
## OUTSTANDING WORKFLOW TASKS:
#########################################################################################################################


## Find points that make the code not reproducible - look for all hard coded paths

## Flatten the structure for an example species

## Improve the raster extract step

## Figure out how to make step 8 parallel - help from Shawn


#########################################################################################################################
###################################################### TBC ############################################################## 
#########################################################################################################################
