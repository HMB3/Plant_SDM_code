#########################################################################################################################
#######################################  TRY RUNNING ALL SDM CODE AT ONCE ############################################### 
#########################################################################################################################


## This code runs the whole SDM workflow for the HIA project, for a subset of species (e.g. those checked by Linda)


#########################################################################################################################
## 1). RUN THE DRAFT CODE
#########################################################################################################################


## In order to run this from katana, we need to use only the necessary packages. Then figure out the order we need to 
## load the data in


#########################################################################################################################
## Read in all data to run the SDM code :: species lists, shapefile, rasters & tables
source('./R/HIA_TREE_LIST.R')


## Set global species list variables here
GBIF.spp      = camp.spp[1:3] #intersect(TPL.SPP, TREE.200.SPP) # workaround for ALA problem use TPL checked species
map_spp_list  = camp_spp[1:3]

GBIF_path     = "./data/base/HIA_LIST/GBIF/OCC_SEARCH/"
ALA_path      = "./data/base/HIA_LIST/ALA/TREE_SPECIES/"

save_dir      = 'output/maxent/CAMPBELLTOWN/'
out_dir       = 'output/maxent/CAMPBELLTOWN'
maxent_path   = 'output/maxent/CAMPBELLTOWN/'


#########################################################################################################################
## Now source each step in the workflow. 
## Step 3 :: combine GBIF occurrence data with ALA data (and hopefully urban data) and filter to records > 1950
source('./R/3)_GBIF_DATA_TAXO_SCIENTIFIC_NAME.R',    echo = TRUE)
source('./R/ALA_DATA_FILTER_TAXO_SCIENTIFIC_NAME.R', echo = TRUE)


## Step 4 :: combine GBIF, ALA and urban occurrence data into a single table, extract environmental condtions 
## & add contextual info for each record (taxonomic and horticultural) 
source('./R/4)_ALA_GBIF_URBAN_COMBINE.R', echo = TRUE)


## Step 5 :: clean the occurrence data using the 'CleanCoordinates' function in the CoordinateCleaner package to remove
## records near herbaria, duplicates, etc.
source('./R/5)_GBIF_ALA_CLEAN_NICHES.R', echo = TRUE)


## Save .RData file at this point, then load in another file to run steps 7 and 8
## This would need all n species processed, then a script that splits


## Step 7 :: Run maxent on a table of all species
source('./R/7)_RUN_MAXENT.R', echo = TRUE)


## Step 8 :: Create habitat suitability maps for each species using six GCMs and three time slices (2030/50/70). 
## Then summary maxent results and estimate species presences in significant urban areas under climate change
## Takes awhile, so probably run different time slices in separate R sessions
source('./R/8)_PREDICT_SDM_SCENARIOS.R', echo = TRUE)





#########################################################################################################################
###################################################### TBC ############################################################## 
#########################################################################################################################