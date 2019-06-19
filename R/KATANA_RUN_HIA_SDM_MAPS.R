#########################################################################################################################
#######################################  TRY RUNNING ALL SDM CODE AT ONCE ############################################### 
#########################################################################################################################



## This code runs the whole SDM workflow for the HIA project, for a subset of species (e.g. whichever you supply)





#########################################################################################################################
## 1). SETUP DATA & PACKAGES
#########################################################################################################################


## First, check which operating system we're on
on_windows = switch(Sys.info()[['sysname']], Windows = TRUE, FALSE)


#########################################################################################################################
## Then read in all data to run the SDM code :: species lists, shapefile, rasters & tables
if (on_windows) {
  
  ## Is up to date
  load("KATANA_RUN_DATA.RData")
  
} else {
  
  cmd_args <- commandArgs(TRUE)
  rdata_file = cmd_args[1]
  if (is.na(rdata_file)) {rdata_file = "KATANA_RUN_DATA.RData"}
  message (rdata_file)
  load(rdata_file)
  
}

#  reassert after loading object
on_windows = switch(Sys.info()[['sysname']], Windows = TRUE, FALSE)


#########################################################################################################################
## Load only the packages needed for the analysis
p <- c('ff',      'things',    'raster',        'dismo',             'sp',           'latticeExtra',          'data.table',
       'rgdal',   'rgeos',     'gdalUtils',     'rmaxent',           'readr',        'plyr',                  'dplyr',
       'tidyr',   'readr',     'rnaturalearth', 'rasterVis',         'RColorBrewer', 'latticeExtra',          'parallel',
       'ALA4R',   'stringr',   'Taxonstand',    'CoordinateCleaner', 'gsubfn',       'PerformanceAnalytics',  'utf8',
       'rvest',   'magrittr',  'devtools',      'ggplot2',           'reshape2',     'rmarkdown',             'flexdashboard', 
       'shiny',   'rgbif',     'ENMeval',       'tibble',            'ncdf4',        'Cairo',                 'taxonlookup',  
       'kgc',     'betareg',   'gridExtra',     'grid',              'lattice',      'ConR',                  'writexl',
       'sf',      'ggmap',     'DataCombine',   'exactextractr',     'mgcv', 'doSNOW', 'tidyverse')


## Require packages
sapply(p, require, character.only = TRUE)
#devtools::source_gist('26e8091f082f2b3dd279')
# source_gist('26e8091f082f2b3dd279',             filename = 'polygonizer.R')
# source_gist('c6a1cb61b8b6616143538950e6ec34aa', filename = 'hatch.R')
#  try to avoid github rate limiting
# devtools::source_gist('306e4b7e69c87b1826db',   filename = 'diverge0.R')


## Source functions, and set temporary directory (for both raster files and generally)
source('./R/WPW_GENERAL_FUNCTIONS.R')
source('./R/WPW_MAXENT_FUNCTIONS.R')
source('./R/WPW_MAPPING_FUNCTIONS.R')
rasterOptions(tmpdir = './RTEMP')


#########################################################################################################################
##  If on Katana, override the raster objects so we use the local versions
if (!on_windows) {
  message ("Loading world raster stack")
  world.grids.current = stack(
    file.path('./data/base/worldclim/world/0.5/bio/current',
              sprintf('bio_%02d', 1:19)))
  
  ## Create a raster stack of current Australian environmental conditions
  ## This is used for exrtacting worldclim data for the tree inventories and running SDMs
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
  
  ## Divide the current Australian bioclim values by 10
  ## The future values are divided by 10 within the mapping function
  message ("Iterating over aus current grids")
  for(i in 1:11) {
    
    ## simple loop
    message(i)
    aus.grids.current[[i]] <- aus.grids.current[[i]]/10
    
  }
  
  
  #########################################################################################################################
  ##  Also include the template data
  template.raster.1km     = raster("./data/world_koppen/template_hasData.tif")
  template.raster.1km.84  = raster("./data/world_koppen/template_1km_WGS84.tif")
  Koppen_1975_1km         = raster("./data/world_koppen/Koppen_1000m_Mollweide54009.tif")
  PET                     = raster("./data/base/worldclim/world/1km/pet_he_yr1.tif")
  
  
  ## If running locally on the windows machine, 
  ## just use the data from the .Rdata file
} else {
  
  message ("Raster data already loaded in Rdata object")
  
}





#########################################################################################################################
## 2). SET GLOBAL VARIABLES
#########################################################################################################################


#########################################################################################################################
## Now set global analysis variables :: these assume you are using an R project folder structure
# least_mapped <- readRDS("./data/ANALYSIS/mapped_least.rds")
# less_mapped  <- readRDS("./data/ANALYSIS/mapped_less.rds")
# out_spp      <- readRDS("./data/ANALYSIS/outstanding_spp.rds")


## Run the species 500 or 1000 at a time
#GBIF.spp = unique(WPW.spp)               ## your list of species.
#GBIF.spp = unique(WPW.non.tree) 
GBIF.spp = unique(WPW.non.tree)
#GBIF.spp = unique(WPW.non.tree)
#GBIF.spp   = TI.HIA    ##
#GBIF.spp = c("Acacia falcata")
#GBIF.spp = unique(gsub("_", " ", unique(out_spp$searchTaxon)))


## Subset for PBS array jobs ::
if (Sys.getenv("PBS_ARRAYID") != "") {
  
  max_i = as.integer(Sys.getenv("WPW_MAX_I"))
  i = as.integer (Sys.getenv("PBS_ARRAYID"))
  GBIF.spp = GBIF.spp[1:length(GBIF.spp) %% 100 == (i-1)]
  
}


#########################################################################################################################
## The required folders must be created on katana
GBIF.spp      = as_utf8(GBIF.spp, normalize = TRUE)   ## Check the species names have the right characters
save_run      = "ALL_SHRUB_JULY_2018"                 ## a variable to append the run name to the output files
                                                      ## Need to include tree or not for HIA list


## If running the trees, use all three data sources
if (grepl("TREE", save_run)) {
  OCC_SOURCE    = c("ALA", "GBIF", "INVENTORY")
  
  ## If running the non-trees, use ALA and GBIF only
} else {
  ## Get the background records from any source
  OCC_SOURCE    = c("ALA", "GBIF") 
  
}


## Create directories to store the maxent output
dir.create('./output/maxent/SPLIT_TEST_INV/')
dir.create('./output/maxent/SPLIT_TEST_INV_BS/')


#########################################################################################################################
## Reading and writing?
read_data     = "FALSE"                              ## Read intermediary data between the steps?
save_data     = "TRUE"                               ## Save data?
check_maps    = "FALSE"                              ## Create maps, shapefiles and histograms of each speices?


#########################################################################################################################
## Data paths
GBIF_path     = "./data/base/HIA_LIST/GBIF_UPDATE/"  ## The path where GBIF data is stored
ALA_path      = "./data/base/HIA_LIST/ALA_UPDATE/"   ## The path where ALA data is stored place
DATA_path     = "./data/ANALYSIS/"                   ## The path where the final data from/for analyses are stored. 
calc_niche    = "TRUE"                               ## Calculate global niches for each species?


#########################################################################################################################
## Analysis/results paths
maxent_path   = './output/maxent/SPLIT_TEST_INV/'      ## The directory where maxent files are saved               
maxent_dir    = 'output/maxent/SPLIT_TEST_INV'         ## Another version of the path needed to run maxent loop

bs_path       = './output/maxent/SPLIT_TEST_INV_BS/'   ## Backwards selection directory
bs_dir        = 'output/maxent/SPLIT_TEST_INV_BS'      ## As above
results_dir   = bs_dir                                 ## The directory where the maxent results are harvested :: backwards selection





#########################################################################################################################
## 3). RUN ANALYSIS
#########################################################################################################################


## Run the SDM work flow for the species set -
## Steps 1 - 5 are hard to complete on the cluster across a huge set of binomials. 
## So I processed all species in groups of 500 through steps 1-6 (run_500, _100..._3500, etc.)
## Only steps 7 and 8 were run on Katana.

## This makes the work flow less repeatable, see Outstanding tasks ........................................................


#########################################################################################################################
## Now source each step in the workflow. 
message('Running SDMs and maps for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")


## Step 1 :: download GBIF and ALA data
# source('./R/1_GBIF_ALA_DOWNLOAD.R', echo = TRUE)


## Step 2 :: combine GBIF occurrence data with ALA data and filter to records > 1950
# source('./R/2_ALA_DATA_FILTER_TAXO_SCIENTIFIC_NAME.R', echo = TRUE)
# source('./R/3_GBIF_DATA_TAXO_SCIENTIFIC_NAME.R',       echo = TRUE)


## Step 4 :: combine GBIF, ALA and tree inventory data into a single table,
## clean taxonomy, then extract raster values
## extract environmental conditions
# source('./R/4_ALA_GBIF_TAXO_COMBINE.R',   echo = TRUE)
# source('./R/INVENTORY_RASTER.R',          echo = TRUE)


## Step 5 :: clean the occurrence data using the 'CleanCoordinates' function in the CoordinateCleaner package, to remove
## records near herbaria, duplicates, etc. & add contextual info for each record (taxonomic and horticultural).
## Then estimate the environmental and geographic range of each species
# source('./R/5_GEO_CLEAN_DATA.R',         echo = TRUE)
# source('./R/CALCULATE_1KM_NICHES.R',     echo = TRUE)


## Step 6 :: Prepare the SDM table and clean the spatial outliers at 1km resolution.
# source('./R/6_PREPARE_SDM_TABLE_1KM.R',  echo = TRUE)
# 
# 
## Step 7 :: Run maxent on a table of all species, using targeted background selection, then backwards selection
## Step 8 :: Create habitat suitability maps for each species using six GCMs and three time slices (2030/50/70).
source('./R/7_RUN_MAXENT.R',      echo = TRUE)
source('./R/8_MAP_SDM_COMBINE.R', echo = TRUE)


## COLLATE SDM RESULTS & PLOT MESS MAPS
#source('./R/COLLATE_MAXENT_KATANA_RESULTS.R', echo = TRUE)





#########################################################################################################################
##     To-do list


############################################################################
##     1). Katana checklist :

##         Trees have inventory data
##         Non-trees have no inventory data

##         Direct transfer of species folders :: how big will they be?         


############################################################################
##     2). Further problems to solve in steps 7 and 8 ::

##     - Summarise the Results table, and also create search for PNG files :: MESS_panel & occ_source

##     - Calculate the environmental ranges and add them to MAXENT results  

##     - Omission rate caculation

##     - Store error messages in step 7 and 8 if possible (didn't work)

##     - Change the mapping function to handle species with no  novel areas :: Use Eucalyptus Camuldelunsis

##     - Barcharts for Australia and the whole world, histograms and convex hull (R users group)

##     - Don't aggregate to postal areas first, write a separate function to post-process (modify SUA_cell_count)





#########################################################################################################################
## OUTSTANDING WORKFLOW TASKS:
#########################################################################################################################


## 1).  Introduce more switches for operating system, etc.

## 2).  Tidy up all the code (using piping, etc.)

## 3).  Make the code more modular, less monolithic

## 3).  Iron out some of the points which are not as reproducible (e.g the background points in step 6)

## 4).  Improve the reading of objects and files - e.g. the RData object is cumbersome.... 





#########################################################################################################################
## OUTSTANDING WORKFLOW TASKS:
#########################################################################################################################


## 1). Table of all the species modelled/screened for the evergreen list... 

##     - Check the latest output. Species in CLEAN.GROW should match species in COMBO.NICHE
##     - They do, yay!
##     - CSV output


## 2).  All directories for evergreen species modelled, zipped up (using WG drives) ::
## 
##      - SDMs (two folders, all variables, and Backwards selected)
##      - Check the output for a few species


##      - Maps for 2030/50/70 for 6 GCMs with MESS 
##      - Gain/loss/stable rasters + tables for 2030/50/70
##      - Local mess function is working. 
##        Need a folder for checking all species 
##      - MESS.png + global occ + 2070 gain/loss.png + Maxent table in one folder
##      - This will need to be stored on the G:drive : check_evergreen_species


## 3).  H:drive (2TB SSD that I've been working off) and G:drive (8TB HDD for backup).

##      - Clean up files on H: ..............................................................................
##      - Put species to check in the results
##      - Store the full evergreen results in G:
##      - Run the code from inside G to test. Keep the full folders on the temporary drive


## 4). Latest code for running SDMs and maps on windows and linux. This include the PBS script


## 5). Mark down file of example species - a 'read me' file


## 6). If possible, environmental ranges, histograms and convex hulls for all species...





#########################################################################################################################
###################################################### TBC ############################################################## 
#########################################################################################################################
