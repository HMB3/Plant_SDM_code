#########################################################################################################################
#######################################  TRY RUNNING ALL SDM CODE AT ONCE ############################################### 
#########################################################################################################################


## This code runs the whole SDM workflow for the HIA project, for a subset of species (e.g. whichever you supply)


#########################################################################################################################
## 1). RUN THE DRAFT CODE
#########################################################################################################################


## In order to run this from katana, we need to use only the necessary packages. 


#########################################################################################################################
## Read in all data to run the SDM code :: species lists, shapefile, rasters & tables
source('./R/HIA_TREE_LIST.R')
## load("H:/green_cities_sdm/TEST_RUN.RData")

## Next step is to model the differences between the SUA species list, and the HIA species list
## That's 377 species currently. Then there are another 1000-odd species on the larger clean list
## It's time to start thinking strategically as to how to structure the analyses. So for one species, how much data is
## taken up? Need to back this up, then conisder what the wesbite will eventually need - make sure all the necessary 
## files are create by the functions, or can be created by running another script/loop. Two directory structures:


## 1). All files, for Macquarie science IT backup

## 2). Reduced file set, for the web developers


## Do we need the rasters themselves, to be read somehow by another function, or do we need already formatted images?
## Depends on how the website works - what language, how they handle objects, etc. Just saving the rasters is better,
## because then style decisons could be made after - that's a separate job.


## One of the outstanding tasks is how to tune the spatial cleaning, and also run the mess map creation. Discuss this with
## Rony, see if he has time to play around with different settings.


#########################################################################################################################
## Step one is to create the same taxonomy for the next lot of species - GBIF, TPL, etc
## Then crunch the data
## Then check the file system flattening systems Shawn sent. Will it work the same way?
## Then Try running the next lot of species remotely. Probably no time to figure out how to make the mapping and summary parallel


#########################################################################################################################
## %%&*). For rapid assessment of species, how could we mine the data for potentially useful species?
## Could search both the overall loss/gain table, and also the SUA table, for the summary of the cells gained and lost.
## So top and tail the list - check maps for the biggest losers and gainers overall. Prioritise checking these species


## Set global species variables here : species lists, and saving directories
## GBIF.spp = sort(trimws(unique(c(MOD.3.SPP$searchTaxon, trait.spp))))
GBIF.spp      = TPL.HIA                                          ## your list of species
GBIF.spp.rev  = sort(GBIF.spp, decreasing = TRUE)                ## the list reversed - only needed for a big list

save_run      = "HIA_SPP_25_GROW"                                ## a variable to append the run name to the output files
map_spp_list  = gsub(" ", "_", GBIF.spp)                         ## species list with "_" for mapping
map_spp_rev   = sort(map_spp_list, decreasing = TRUE)            ## reversed, so we can run two at once

GBIF_path     = "./data/base/HIA_LIST/GBIF/OCC_SEARCH/"          ## The path where GBIF data is stored
ALA_path      = "./data/base/HIA_LIST/ALA/TREES_TEST/"           ## The path where ALA data is stored - duplicated if in the same place

maxent_path   = './output/maxent/SUA_TREES_ANALYSIS/'            ## The directory where files are saved               
maxent_dir    = 'output/maxent/SUA_TREES_ANALYSIS'               ## Another version of the path that John's coding needs to run a loop
save_data     = 'FALSE'
save_data     = "TRUE"


#########################################################################################################################
## Now source each step in the workflow. 
message('Running SDM workflow for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")


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
source('./R/6)_PREPARE_SDM_TABLE_1KM.R',  echo = TRUE)


## Step 7 :: Run maxent on a table of all species
source('./R/7)_RUN_MAXENT.R', echo = TRUE)


## Step 8 :: Create habitat suitability maps for each species using six GCMs and three time slices (2030/50/70). 
## Then summarise maxent results and estimate species presences in significant urban areas under climate change
## Takes awhile, so probably run different time slices (2030, 2050, 2070) in separate R sessions
source('./R/8)_MAP_SDM_COMBINE.R', echo = TRUE)





#########################################################################################################################
## OUTSTANDING WORKFLOW TASKS:
#########################################################################################################################


## Find points that make the code not reproducible
## Figure out how to make step 8 parallel


#########################################################################################################################
###################################################### TBC ############################################################## 
#########################################################################################################################