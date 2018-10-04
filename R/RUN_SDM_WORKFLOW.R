#########################################################################################################################
#######################################  TRY RUNNING ALL SDM CODE AT ONCE ############################################### 
#########################################################################################################################


## This code runs the whole SDM workflow for the HIA project, for a subset of species (e.g. whichever you supply)


#########################################################################################################################
## 1). RUN THE DRAFT CODE
#########################################################################################################################


## In order to run this from katana, we need to use only the necessary packages. Then figure out the order we need to 
## load the data in


#########################################################################################################################
## Read in all data to run the SDM code :: species lists, shapefile, rasters & tables
source('./R/HIA_TREE_LIST.R')


## Set global species variables here : species lists, and saving directories
GBIF.spp      = TPL.SPP # setdiff(TPL.SPP, intersect(TPL.SPP, ala.download))  ## your list of species
GBIF.spp.rev  = sort(GBIF.spp, decreasing = TRUE)                   ## the list reversed - only needed for a big list

save_run      = "COMBINED_SUA_400_SPP_300km_SPAT"     ## a variable to append the run name to the output files
map_spp_list  = gsub(" ", "_", GBIF.spp)                            ## species list with "_" for mapping
map_spp_rev   = sort(map_spp_list, decreasing = TRUE)               ## reversed

GBIF_path     = "./data/base/HIA_LIST/GBIF/OCC_SEARCH/"             ## The path where GBIF data is stored
ALA_path      = "./data/base/HIA_LIST/ALA/TREES_TEST/"              ## The path where ALA data is stored - duplicated if in the same place

maxent_path   = 'output/maxent/SUA_TREES_OLD_ALA/'                  ## The directory where files are saved               


#########################################################################################################################
## Now source each step in the workflow. 


## Step 3 :: combine GBIF occurrence data with ALA data and filter to records > 1950
source('./R/ALA_DATA_FILTER_TAXO_SCIENTIFIC_NAME.R', echo = TRUE)
source('./R/3)_GBIF_DATA_TAXO_SCIENTIFIC_NAME.R',    echo = TRUE)


## Step 4 :: combine GBIF, ALA and urban occurrence data into a single table, extract environmental condtions 
source('./R/4)_ALA_GBIF_URBAN_COMBINE.R', echo = TRUE)
source('./R/INVENTORY_RASTER.R',          echo = TRUE)


## Step 5 :: clean the occurrence data using the 'CleanCoordinates' function in the CoordinateCleaner package to remove
## records near herbaria, duplicates, etc. & add contextual info for each record (taxonomic and horticultural) 
## Then prepare the SDM table
## Then clean the spatial outliers
source('./R/5)_GBIF_ALA_CLEAN_NICHES.R',  echo = TRUE)
source('./R/6)_PREPARE_SDM_TABLE.R',      echo = TRUE)
source('./R/CLEAN_SPATIAL_OUTLIERS.R',    echo = TRUE)


## Step 7 :: Run maxent on a table of all species
source('./R/7)_RUN_MAXENT.R', echo = TRUE)


## Step 8 :: Create habitat suitability maps for each species using six GCMs and three time slices (2030/50/70). 
## Then summarise maxent results and estimate species presences in significant urban areas under climate change
## Takes awhile, so probably run different time slices (2030, 2050, 2070) in separate R sessions
source('./R/8)_PREDICT_SDM_SCENARIOS.R', echo = TRUE)





#########################################################################################################################
## OUTSTANDING WORKFLOW TASKS:
#########################################################################################################################


## Find points that make the code not reproducible
## Figure out how to make step 8 parallel


#########################################################################################################################
###################################################### TBC ############################################################## 
#########################################################################################################################