#########################################################################################################################
#######################################  TRY RUNNING ALL SDM CODE AT ONCE ############################################### 
#########################################################################################################################


#########################################################################################################################
## Shawn's advice is to make the code more modular, rather than monolithic.
## Small examples of how to do this could be applied across the different steps 



## 1). Check that the paralled approach will work for a few species, using the MESS approach - EG for 10 species.


## 2). Check the spatial outlier step can be improved, before re-running all the species through these steps.
##     Use Alex Ziska's example to calcualte
##     
##     Check the Area of occupancy calculations work - could be causing problems.





## This code runs the whole SDM workflow for the HIA project, for a subset of species (e.g. whichever you supply)
## In order to run this from katana, we need to use only the necessary packages.
## save.image("TEST_RUN.RData") 
on_windows = switch(Sys.info()[['sysname']], Windows = TRUE, FALSE)


#########################################################################################################################
## Read in all data to run the SDM code :: species lists, shapefile, rasters & tables
#source('./R/HIA_TREE_LIST.R')
if (on_windows) {
  
    load("UPDATE_DATA.RData")
  
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
       'ALA4R',   'stringr',   'Taxonstand',    'CoordinateCleaner', 'gsubfn',       'PerformanceAnalytics',
       'rvest',   'magrittr',  'devtools',      'ggplot2',           'reshape2',     'rmarkdown',             'flexdashboard', 'shiny', 'rgbif',
       'ENMeval', 'tibble',    'ncdf4',         'Cairo',             'taxonlookup',  'kgc',        'maptools', 'DataCombine', 'mgcv', 'rsq', 'utf8',
       'betareg', 'hydroTSM',  'bomrang',       'gridExtra',         'grid',         'lattice', 'ConR', 'writexl')


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
template.raster.1km  = raster("./data/template_hasData.tif")
Koppen_1975_1km      = raster('data/Koppen_1000m_Mollweide54009.tif')
PET                  = raster("./data/base/worldclim/world/1km/pet_he_yr1.tif")


## If running locally on the windows machine, just use the data from the Rdata file
} else {
  
  message ("Raster data already loaded in Rdata object")
  
}


#########################################################################################################################
## Now set global analysis variables :: these assume you are using an R project folder structure
## Replace any non - UTF8 character in the string with UTF
Sys.setenv("PBS_ARRAYID" = 1)
GBIF.spp = HOLLOW.SPP[4] ## your list of species


##  Subset for PBS array jobs
if (length(Sys.getenv("PBS_ARRAYID"))) {
  
    i = as.integer (Sys.getenv("PBS_ARRAYID"))
    GBIF.spp = GBIF.spp[1:length(GBIF.spp) %% 200 == i]
    
}


#########################################################################################################################
## These folders must be created on katana, or uploaded
GBIF.spp      = as_utf8(GBIF.spp, normalize = TRUE)  ## Check the species names have the right characters
save_run      = "WPW_TEST"                           ## a variable to append the run name to the output files
read_data     = "FALSE"                              ## only needed to load the SDM table in
save_data     = "FALSE"                              ## only needed to load the SDM table in

GBIF_path     = "./data/base/HIA_LIST/GBIF/"  ## The path where GBIF data is stored
ALA_path      = "./data/base/HIA_LIST/ALA/"   ## The path where ALA data is stored place
DATA_path     = "./data/ANALYSIS/"            ## The path where the final data for analyses are stored 

maxent_path   = './output/maxent/HOLLOW_SPP_PROP_SAMPLE_ALL_VARS/'    ## The directory where files are saved               
maxent_dir    = 'output/maxent/HOLLOW_SPP_PROP_SAMPLE_ALL_VARS'       ## Another version of the path needed to run maxent loop

bs_path       = './output/maxent/HOLLOW_SPP_PROP_SAMPLE_ALL_VARS_BS/'
bs_dir        = 'output/maxent/HOLLOW_SPP_PROP_SAMPLE_ALL_VARS_BS'    ## Backwards selection directory
results_dir   = bs_dir                                                ## the directory to harvest results : BS dir?
   




#########################################################################################################################
## 2). RUN SDMS, MAPS and COMBINE
#########################################################################################################################


## Run the SDM workflow for the species set.
## Steps 1 - 5 are hard to complete across a huges set og binomials. So try just using the data I created locally to run
## only steps 7 and 8.


#########################################################################################################################
## Now source each step in the workflow. 
message('Running SDMs and maps for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")


## Step 3 :: combine GBIF occurrence data with ALA data and filter to records > 1950
source('./R/1_GBIF_ALA_DOWNLOAD.R', echo = TRUE)


## Step 3 :: combine GBIF occurrence data with ALA data and filter to records > 1950
source('./R/ALA_DATA_FILTER_TAXO_SCIENTIFIC_NAME.R', echo = TRUE)
source('./R/3_GBIF_DATA_TAXO_SCIENTIFIC_NAME.R',     echo = TRUE)


## Step 4 :: combine GBIF, ALA and tree inventory data into a single table, extract environmental condtions
source('./R/4_ALA_GBIF_URBAN_COMBINE.R',  echo = TRUE)
source('./R/INVENTORY_RASTER.R',          echo = TRUE)


## Step 5 :: clean the occurrence data using the 'CleanCoordinates' function in the CoordinateCleaner package, to remove
## records near herbaria, duplicates, etc. & add contextual info for each record (taxonomic and horticultural). 
## Then prepare the SDM table
## Then clean the spatial outliers
source('./R/5_GBIF_ALA_CLEAN_NICHES.R',  echo = TRUE)
source('./R/6_PREPARE_SDM_TABLE_1KM.R',  echo = TRUE)



## Step 7 :: Run maxent on a table of all species, using targetted background selection, then backwards selection
source('./R/7_RUN_MAXENT.R', echo = TRUE)


## Step 8 :: Create habitat suitability maps for each species using six GCMs and three time slices (2030/50/70). 
## Then summarise maxent results and estimate species presences in significant urban areas under climate change
## Then run MESS maps. Or, the MESS maps could be run before the combination step
source('./R/8_MAP_SDM_COMBINE.R', echo = TRUE)





#########################################################################################################################
## OUTSTANDING WORKFLOW TASKS:
#########################################################################################################################




## 1) Figure out how to make the whole process parallel - help from Shawn.




#########################################################################################################################
###################################################### TBC ############################################################## 
#########################################################################################################################