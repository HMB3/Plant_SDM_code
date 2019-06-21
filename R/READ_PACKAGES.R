## Load only the packages needed for the analysis
p <- c('ff',    'things',    'raster',        'dismo',        'sp',           'latticeExtra', 'data.table', 
       'rgdal', 'rgeos',     'gdalUtils',     'rmaxent',      'readr',        'plyr',         'dplyr',        
       'tidyr', 'readr',     'rnaturalearth', 'rasterVis',    'RColorBrewer', 'latticeExtra', 'parallel',     
       'ALA4R', 'stringr',   'Taxonstand',    'CoordinateCleaner', 'gsubfn',  'PerformanceAnalytics',
       'rvest', 'magrittr',  'devtools',      'ggplot2',      'reshape2',     'rmarkdown', 'flexdashboard', 'shiny', 'rgbif',
       'ENMeval', 'tibble',  'ncdf4',         'Cairo', 'velox', 'taxonlookup', 'kgc', 'maptools', 'DataCombine', 'mgcv', 'knitr', 'utf8')


## Require packages
sapply(p, require, character.only = TRUE)
source_gist('26e8091f082f2b3dd279',             filename = 'polygonizer.R')
source_gist('c6a1cb61b8b6616143538950e6ec34aa', filename = 'hatch.R')


## Source functions, and set temporary directory (for both raster files and just generally)
source('./R/GREEN_CITIES_FUNCTIONS.R')
source('./R/MAXENT_FUNCTIONS.R')
source('./R/MAPPING_FUNCTIONS.R')
source('./R/HIA_CLEAN_MATCHING.R')
rasterOptions(tmpdir = file.path('/green_cities_sdm/RTEMP'))


## Set variables
GBIF.spp      = native.good.models[1:5]                  ## A list of species
GBIF.spp      = as_utf8(GBIF.spp, normalize = TRUE)
GBIF.spp.rev  = sort(GBIF.spp,   decreasing = TRUE)              ## the list reversed - only needed for a big list

save_run      = "SUA_RMD"                                        ## a variable to append the run name to the output file
map_spp_list  = gsub(" ", "_", GBIF.spp)                         ## species list with "_" for mapping
map_spp_rev   = sort(map_spp_list, decreasing = TRUE)            ## reversed, so we can run two at once

GBIF_path     = "./data/base/HIA_LIST/GBIF/"                     ## The path where GBIF data is stored
ALA_path      = "./data/base/HIA_LIST/ALA/"                      ## The path where ALA data is stored place
DATA_path     = "./data/ANALYSIS/"                               ## The path where the final data for analyses are stored 
SHP_path      = "./data/ANALYSIS"                                ## The data path for readOGR dsn  
unit_path     = "./data/base/CONTEXTUAL/SUA/"                    ## The data path for the spatial unti of analysis 
unit_file     = "SUA_2016_AUST.rds"                              ## The spatial untit of analysis - E.G. SUAs
unit_vec      = "SUA_2016_VEC.rds"                               ## A vector of the spatial units


maxent_path   = './output/maxent/SUA_TREES_ANALYSIS/'            ## The directory where files are saved               
maxent_dir    = 'output/maxent/SUA_TREES_ANALYSIS'               ## Another version of the path needed to run maxent loop
save_data     = 'FALSE'                                          ## Arguments for saving the intermediary output - i.e. niches
read_data     = 'FALSE'                                          ## Leave these the same - saves data, but doesn't read back in
calc_niche    = 'TRUE'                                           ## Calculate niches?
OCC_SOURCE    = 'ALL'                                            ## Create niches using ALA, GBIF and inventory data 