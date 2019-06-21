#########################################################################################################################
#######################################  TRY RUNNING ALL SDM CODE AT ONCE ############################################### 
#########################################################################################################################


#########################################################################################################################
##     To-do list

##     1). Trees have inventory data
##         Non-trees have no inventory data

##     2). Re-run steps 1-5, 1000 at a time
##     -   Make the raster extraction on just unique cells, not all cells. Currently on 

##     3). Problems to solve in steps 7 and 8 ::

##     - Omission rate caculation?    

##     - Store error messages in step 7 and 8 if possible.

##     - Change the mapping function to handle species with no  novel areas :: Use Eucalyptus Camuldelunsis

##     - Add 'ncores' argument to both functions (use a few locally, but 1 on Katana)

##     - Barcharts for Australia and the whole world (I can probably figure this out)


## This code runs the whole SDM workflow for the HIA project, for a subset of species (e.g. whichever you supply)
## In order to run this from katana, we need to use only the necessary packages.
## save.image("TEST_RUN.RData") 
on_windows = switch(Sys.info()[['sysname']], Windows = TRUE, FALSE)


#########################################################################################################################
## Read in all data to run the SDM code :: species lists, shapefile, rasters & tables
#source('./R/HIA_TREE_LIST.R')
if (on_windows) {
  
  ## Is up to date
  load("TEST_RUN.RData")
  
} else {
  
  cmd_args <- commandArgs(TRUE)
  rdata_file = cmd_args[1]
  if (is.na(rdata_file)) {rdata_file = "../UPDATE_DATA.RData"}
  message (rdata_file)
  load(rdata_file)
  
}





#########################################################################################################################
## 1). SETUP
#########################################################################################################################


#########################################################################################################################
## Load only the packages needed for the analysis
p <- c('ff',      'things',    'raster',        'dismo',             'sp',           'latticeExtra',          'data.table',
       'rgdal',   'rgeos',     'gdalUtils',     'rmaxent',           'readr',        'plyr',                  'dplyr',
       'tidyr',   'readr',     'rnaturalearth', 'rasterVis',         'RColorBrewer', 'latticeExtra',          'parallel',
       'ALA4R',   'stringr',   'Taxonstand',    'CoordinateCleaner', 'gsubfn',       'PerformanceAnalytics',  'utf8',
       'rvest',   'magrittr',  'devtools',      'ggplot2',           'reshape2',     'rmarkdown',             'flexdashboard', 
       'shiny',   'rgbif',     'ENMeval',       'tibble',            'ncdf4',        'Cairo',                 'taxonlookup',  
       'kgc',     'betareg',   'gridExtra',     'grid',              'lattice',      'ConR',                  'writexl',
       'ggmap',   'DataCombine')


## Require packages
sapply(p, require, character.only = TRUE)
source_gist('26e8091f082f2b3dd279',             filename = 'polygonizer.R')
source_gist('c6a1cb61b8b6616143538950e6ec34aa', filename = 'hatch.R')
devtools::source_gist('306e4b7e69c87b1826db',   filename = 'diverge0.R')


## Source functions, and set temporary directory (for both raster files and generally)
source('./R/WPW_GENERAL_FUNCTIONS.R')
source('./R/WPW_MAXENT_FUNCTIONS.R')
source('./R/WPW_MAPPING_FUNCTIONS.R')
rasterOptions(tmpdir = './RTEMP')


#########################################################################################################################
##  If on Katana, override the raster objects so we use the local versions
if (on_windows == "FALSE") {
  
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
  ## The future values are divided nby 10 within the mapping function
  message ("Iterating over aus current grids")
  for(i in 1:11) {
    
    ## simple loop
    message(i)
    aus.grids.current[[i]] <- aus.grids.current[[i]]/10
    
  }
  
  
  #########################################################################################################################
  ##  Also include the template data
  template.raster.1km     = raster("./data/template_hasData.tif")
  template.raster.1km.84  = raster("./data/template_hasData.tif")
  Koppen_1975_1km         = raster('data/Koppen_1000m_Mollweide54009.tif')
  PET                     = raster("./data/base/worldclim/world/1km/pet_he_yr1.tif")
  
  
  ## If running locally on the windows machine, just use the data from the Rdata file
} else {
  
  message ("Raster data already loaded in Rdata object")
  
}


#########################################################################################################################
## Now set global analysis variables :: these assume you are using an R project folder structure
## Replace any non - UTF8 character in the string with UTF
## WPW.tree[1:500] ## WPW.tree[501:1110] 
# ## WPW.non.tree[1:1000] ## WPW.tree[1001:2000] ## WPW.tree[2001:2907]
# Sys.setenv("PBS_ARRAYID" = 1)
# 
# 
# ## Run the species 500 or 1000 at a time
# GBIF.spp = WPW.tree[1:500]
# 
# 
# ##  Subset for PBS array jobs
# if (length(Sys.getenv("PBS_ARRAYID"))) {
#   
#     i = as.integer (Sys.getenv("PBS_ARRAYID"))
#     GBIF.spp = GBIF.spp[1:length(GBIF.spp) %% 200 == i]
#     
#}


#########################################################################################################################
## These folders must be created on katana, or uploaded
GBIF.spp      = WPW.spp[2000:2500]                        ## Your list of species  
GBIF.spp      = as_utf8(GBIF.spp, normalize = TRUE)  ## Check the species names have the right characters
save_run      = "EVERGREEN_2500"                     ## a variable to append the run name to the output files
read_data     = "FALSE"                              ## Read intermediary data between the steps?
save_data     = "TRUE"                               ## Save data?
check_maps    = "FALSE"                              ## Create maps, shapefiles and histograms of each speices?

GBIF_path     = "./data/base/HIA_LIST/GBIF_UPDATE/"  ## The path where GBIF data is stored
ALA_path      = "./data/base/HIA_LIST/ALA_UPDATE/"   ## The path where ALA data is stored place
DATA_path     = "./data/ANALYSIS/"                   ## The path where the final data from/for analyses are stored 
OCC_SOURCE    = c("ALA", "GBIF", "INVENTORY")        ## Data source :: for trees, ALL. For non trees, ALA/GBIF, for occ and bg
#OCC_SOURCE    = c("ALA", "GBIF")
calc_niche    = "TRUE"

maxent_path   = './output/maxent/HIA_TEST_INV/'      ## The directory where maxent files are saved               
maxent_dir    = 'output/maxent/HIA_TEST_INV'         ## Another version of the path needed to run maxent loop

bs_path       = './output/maxent/HIA_TEST_INV_BS/'   ## Backwards selection directory
bs_dir        = 'output/maxent/HIA_TEST_INV_BS'      ## the directory to harvest results : BS dir?
results_dir   = bs_dir   





#########################################################################################################################
## 2). RUN SDMS, MAPS and COMBINE
#########################################################################################################################


## Run the SDM workflow for the species set.
## Steps 1 - 5 are hard to complete across a huges set og binomials. So try just using the data I created locally to run
## only steps 7 and 8.


#########################################################################################################################
## Now source each step in the workflow. 
message('Running SDMs and maps for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")


## Step 1 :: download GBIF and ALA data
#source('./R/1_GBIF_ALA_DOWNLOAD.R', echo = TRUE)


## Step 2 :: combine GBIF occurrence data with ALA data and filter to records > 1950
source('./R/ALA_DATA_FILTER_TAXO_SCIENTIFIC_NAME.R', echo = TRUE)
source('./R/3_GBIF_DATA_TAXO_SCIENTIFIC_NAME.R',     echo = TRUE)


## Step 4 :: combine GBIF, ALA and tree inventory data into a single table, extract environmental condtions
source('./R/4_ALA_GBIF_TAXO_COMBINE.R',   echo = TRUE)
source('./R/INVENTORY_RASTER.R',          echo = TRUE)


## Step 5 :: clean the occurrence data using the 'CleanCoordinates' function in the CoordinateCleaner package, to remove
## records near herbaria, duplicates, etc. & add contextual info for each record (taxonomic and horticultural). 
## Then prepare the SDM table
## Then clean the spatial outliers
source('./R/5_GEO_CLEAN_DATA.R',         echo = TRUE)
source('./R/6_PREPARE_SDM_TABLE_1KM.R',  echo = TRUE)


## Step 7 :: Run maxent on a table of all species, using targetted background selection, then backwards selection
## Step 8 :: Create habitat suitability maps for each species using six GCMs and three time slices (2030/50/70). 
source('./R/7_RUN_MAXENT.R',      echo = TRUE)
source('./R/8_MAP_SDM_COMBINE.R', echo = TRUE)





#########################################################################################################################
## OUTSTANDING WORKFLOW TASKS:
#########################################################################################################################




## 1) Figure out how to make the whole process parallel - help from Shawn.




#########################################################################################################################
###################################################### TBC ############################################################## 
#########################################################################################################################