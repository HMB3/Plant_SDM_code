#########################################################################################################################
#######################################  TRY RUNNING ALL SDM CODE AT ONCE ############################################### 
#########################################################################################################################


## This code runs the whole SDM workflow for the HIA project, for a subset of species (e.g. those checked by Linda)


#########################################################################################################################
## 1). RUN THE DRAFT CODE
#########################################################################################################################


#########################################################################################################################
## Read in all data to run the SDM code :: species lists, shapefile, rasters & tables
source('./R/HIA_LIST_MATCHING.R')


## Set global species list variables here...............................................................................
## Now running more tree species
GBIF.spp     = exotic.trees
map_spp_list = new_trees


#########################################################################################################################
## Now source each step in the workflow 
## Step 3 :: combine GBIF occurrence data with ALA data (and hopefully urban data) and filter to records > 1950
source('./R/3)_GBIF_DATA_FILTER.R', echo = TRUE)


## Step 4 :: combine GBIF, ALA and urban occurrence data into a single table, extracts environmental condtions 
## & add contextual info for each record (taxonomic and horticultural) 
source('./R/4)_ALA_GBIF_URBAN_COMBINE.R', echo = TRUE)


## Step 5 :: clean the occurrence data using the 'CleanCoordinates' function in the CoordinateCleaner package to remove
## records near herbaria, duplicates, etc.
source('./R/5)_GBIF_ALA_CLEAN.R', echo = TRUE)


## Step 7 :: Run maxent on a table of all species
source('./R/7)_RUN_MAXENT.R', echo = TRUE)


## Step 8 :: Create habitat suitability maps for each species using six GCMs and three time slices (2030/50/70). 
## Then summary maxent results and estimate species presences in significant urban areas under climate change
## Takes awhile, so probably run different time slices in separate R sessions
source('./R/8)_PREDICT_SDM_SCENARIOS.R', echo = TRUE)





#########################################################################################################################
###################################################### TBC ############################################################## 
#########################################################################################################################