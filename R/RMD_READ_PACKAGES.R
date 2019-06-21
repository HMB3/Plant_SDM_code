#########################################################################################################################
## Load only the packages needed for the analysis
p <- c('ff',      'things',    'raster',        'dismo',             'sp',           'latticeExtra',          'data.table',
       'rgdal',   'rgeos',     'gdalUtils',     'rmaxent',           'readr',        'plyr',                  'dplyr',
       'tidyr',   'readr',     'rnaturalearth', 'rasterVis',         'RColorBrewer', 'latticeExtra',          'parallel',
       'ALA4R',   'stringr',   'Taxonstand',    'CoordinateCleaner', 'gsubfn',       'PerformanceAnalytics',  'utf8',
       'rvest',   'magrittr',  'devtools',      'ggplot2',           'reshape2',     'rmarkdown',             'flexdashboard', 
       'shiny',   'rgbif',     'ENMeval',       'tibble',            'ncdf4',        'Cairo',                 'taxonlookup',  
       'kgc',     'betareg',   'gridExtra',     'grid',              'lattice',      'ConR',                  'writexl',
       'sf',      'ggmap',     'DataCombine',   'exactextractr',     'mgcv', 'doSNOW', 'tidyverse', 'knitr')


## Require packages
sapply(p, require, character.only = TRUE)
#devtools::source_gist('26e8091f082f2b3dd279')
# source_gist('26e8091f082f2b3dd279',             filename = 'polygonizer.R')
# source_gist('c6a1cb61b8b6616143538950e6ec34aa', filename = 'hatch.R')
#  try to avoid github rate limiting
# devtools::source_gist('306e4b7e69c87b1826db',   filename = 'diverge0.R')


## Source functions, and set temporary directory (for both raster files and generally)
load("KATANA_RUN_DATA.RData")
source('./R/WPW_GENERAL_FUNCTIONS.R')
source('./R/WPW_MAXENT_FUNCTIONS.R')
source('./R/WPW_MAPPING_FUNCTIONS.R')
rasterOptions(tmpdir = './RTEMP')


## Set variables
#########################################################################################################################
## Reading and writing?
GBIF.spp = unique(WPW.spp)
read_data     = "TRUE"                              ## Read intermediary data between the steps?
save_data     = "FALSE"                              ## Save data?
check_maps    = "FALSE"                              ## Create maps, shapefiles and histograms of each speices?
save_run      = "ALL_EVREGREEN_JUNE_2018"  

#########################################################################################################################
## Data paths
GBIF_path     = "./data/base/HIA_LIST/GBIF_UPDATE/"  ## The path where GBIF data is stored
ALA_path      = "./data/base/HIA_LIST/ALA_UPDATE/"   ## The path where ALA data is stored place
DATA_path     = "./data/ANALYSIS/"                   ## The path where the final data from/for analyses are stored. 
calc_niche    = "TRUE"


#########################################################################################################################
## Analysis/results paths
maxent_path   = './output/maxent/HIA_TEST_INV/'      ## The directory where maxent files are saved               
maxent_dir    = 'output/maxent/HIA_TEST_INV'         ## Another version of the path needed to run maxent loop

bs_path       = './output/maxent/HIA_TEST_INV_BS/'   ## Backwards selection directory
bs_dir        = 'output/maxent/HIA_TEST_INV_BS'      ## the directory to harvest results : BS dir?
results_dir   = bs_dir