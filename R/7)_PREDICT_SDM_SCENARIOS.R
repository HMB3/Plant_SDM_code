#########################################################################################################################
######################################### PREDICT MAXENT TO FUTURE CLIMATES ############################################# 
#########################################################################################################################


## This code takes the output of the Maxent species distribution models and generates a prediction of habitat suitability 
## for current and future environmental conditions. The data table is the format of all species occurrences (rows) and 
## environmental variables (columns)


#########################################################################################################################
## 1). READ IN THE CURRENT AND FUTURE ENVIRONMENTAL DATA
#########################################################################################################################


#########################################################################################################################
## Load packages
load("./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT.RData")
source('./R/HIA_LIST_MATCHING.R')


#########################################################################################################################
## create a list of GCM scenarios (below is from CSIRO): 


## Eight of the 40 CMIP5 models assessed in this project have been selected for use in provision of application-ready data. 
## This facilitates efficient exploration of climate projections for Australia.
## https://www.climatechangeinaustralia.gov.au/en/support-and-guidance/faqs/eight-climate-models-data/
# scen = c("ip85bi50", "mc85bi50", "mg85bi50", "mi85bi50", "mp85bi50", 
#          "mr85bi50", "no85bi50", "ac85bi50", "bc85bi50", "cc85bi50", 
#          "cn85bi50", "gf85bi50", "gs85bi50", "hd85bi50", "he85bi50",
#          "hg85bi50", "in85bi50")


## Create a lookup table of GCMs
h <- read_html('http://www.worldclim.org/cmip5_30s') 
gcms <- h %>% 
  html_node('table') %>% 
  html_table(header=TRUE) %>% 
  filter(rcp85 != '')

id.50 <- h %>% 
  html_nodes(xpath = "//a[text()='bi']/@href") %>% 
  as.character %>% 
  grep('85bi50', ., value=TRUE) %>% 
  basename %>% 
  sub('\\.zip.', '', .)

id.70 <- h %>% 
  html_nodes(xpath = "//a[text()='bi']/@href") %>% 
  as.character %>% 
  grep('85bi70', ., value=TRUE) %>% 
  basename %>% 
  sub('\\.zip.', '', .)

gcms.50 <- cbind(gcms, id.50)
gcms.50$GCM = sub(" \\(#\\)", "", gcms$GCM)  ## sub replaces first instance in a string, gsub = global

gcms.70 <- cbind(gcms, id.70)
gcms.70$GCM = sub(" \\(#\\)", "", gcms$GCM) 


## Now create the scenario lists across which to loop
gcms.50 ; gcms.70


## Just get the 8 models picked by CSIRO for Australia, for 2050 and 2070
## Also will be a better way to get at this...
scen_2050 = c("mc85bi50", "no85bi50", "ac85bi50", "cn85bi50", "gf85bi50", "hg85bi50")
scen_2070 = c("mc85bi70", "no85bi70", "ac85bi70", "cn85bi70", "gf85bi70", "hg85bi70")


#########################################################################################################################
## Then create a stack of current environmental conditions, and an Australia shapefile for mapping
aus <- ne_states(country = 'Australia') %>% 
  subset(!grepl('Island', name))

env.grids.current <- stack(
  file.path('./data/base/worldclim/aus/0.5/bio/current',
            sprintf('bio_%02d.tif', 1:19)))

## Now divide the current temperature grids by 10
for(i in 1:11) {
  
  message(i)
  env.grids.current[[i]] <- env.grids.current[[i]]/10
  
}


## What do the grids look like? 11,177,684 NAs???!!!!!
summary(env.grids.current[[1]])
summary(env.grids.current[[5]])
summary(env.grids.current[[11]])





#########################################################################################################################
## 2). CREATE SPECIES LISTS FOR MODEL RUNS 
#########################################################################################################################


## All species on the growers list, and also a test list which includes Renee's species
species_list  <- basename(list.dirs('./output/maxent/STD_VAR_ALL',   recursive = FALSE))
species_rev   = sort(species_list, decreasing = TRUE)


## Now make the test species directory names
test_spp = gsub(" ", "_", test.spp)
# test_spp = intersect(test_spp, species_list)
# test_rev = sort(test_spp, decreasing = TRUE)




#########################################################################################################################
## 3). PROJECT MAXNET MODELS FOR MULTIPLE CLIMATE SCEANARIOS, 2050 AND 2070
#########################################################################################################################


#########################################################################################################################
## For each species, use a function to create raster files and maps of all six climate scenarios (GCMs)
## Note that some of the experimental species - e.g. Kennedia_beckxiana - still have to be processed
env.grids.2050 = project.grids.2050(scen_2050, test_spp)
env.grids.2070 = project.grids.2070(scen_2070, test_spp)


#########################################################################################################################
## First list all the directories containing the rasters
SDM.RESULTS.DIR <- test_spp[c(1:length(test_spp))] %>%
  
  ## Pipe the list into lapply
  lapply(function(species) {
    
    ## Create the character string
    m <-   sprintf('./output/maxent/STD_VAR_ALL/%s/full/', species)
    m 
    
  }) %>%
  
  ## Bind the list together
  c()


#########################################################################################################################
## Then for each species, use a function to combine the raster files for each climate scenarios into one layer
suitability.2050 = mapply(ensemble.2050, SDM.RESULTS.DIR, test_spp)
suitability.2070 = mapply(ensemble.2070, SDM.RESULTS.DIR, test_spp)



#########################################################################################################################
## OUTSTANDING PREDICTION TASKS:
#########################################################################################################################


## We are missing two scenarios recommended for Australia: CanESM2 & CESM1-CAM5. These two are not on the worldclim list:

## https://www.climatechangeinaustralia.gov.au/en/support-and-guidance/faqs/eight-climate-models-data/

## http://www.worldclim.org/cmip5_30s
  

## Consider the final format needed: which files: table, plots/maps, files. Disk space important for both local and web... 


## What are the options for presenting consensus layers of all scenarios? EG show the occurrences, the average and then 
## the combined layers? Need the format first, as the details are tricky (e.g. scale bars, etc.)





#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################