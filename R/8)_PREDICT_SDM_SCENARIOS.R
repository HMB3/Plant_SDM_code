#########################################################################################################################
###################### PREDICT MAXENT TO FUTURE CLIMATES AND SUMMARISE THE RESULTS ###################################### 
#########################################################################################################################


#########################################################################################################################
## This code takes the output of the Maxent species distribution models (i.e. using current conditions), and generates 
## a prediction of habitat suitability for current and future environmental conditions. The input data table is in the 
## format of all species occurrences (rows) and environmental variables (columns).

## Then the predictions from all GCMs (currently using six) are then combined into a single habitat suitability layer.
## Using this combined layer, the loss or gain of species within areal units (e.g. significant urban areas or LGAs), 
## between time periods (current, 2030, 2070) can be calculated. 


#########################################################################################################################
## 1). READ IN THE CURRENT AND FUTURE ENVIRONMENTAL DATA FOR MULTIPLE GCMs
#########################################################################################################################


#########################################################################################################################
## Load packages
load("./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT.RData")
source('./R/HIA_LIST_MATCHING.R')


#########################################################################################################################
## Create a list of GCM scenarios, which are used to create maps of habitat suitability 


## Eight of the 40 CMIP5 models assessed in this project have been selected for use in provision of application-ready data. 
## This facilitates efficient exploration of climate projections for Australia.
## https://www.climatechangeinaustralia.gov.au/en/support-and-guidance/faqs/eight-climate-models-data/
# scen = c("ip85bi50", "mc85bi50", "mg85bi50", "mi85bi50", "mp85bi50", 
#          "mr85bi50", "no85bi50", "ac85bi50", "bc85bi50", "cc85bi50", 
#          "cn85bi50", "gf85bi50", "gs85bi50", "hd85bi50", "he85bi50",
#          "hg85bi50", "in85bi50")


## Create a lookup table of GCMs using the WORLDCLIM website
h <- read_html('http://www.worldclim.org/cmip5_30s') 
gcms <- h %>% 
  html_node('table') %>% 
  html_table(header = TRUE) %>% 
  filter(rcp85 != '')

id.50 <- h %>% 
  html_nodes(xpath = "//a[text()='bi']/@href") %>% 
  as.character %>% 
  grep('85bi50', ., value = TRUE) %>% 
  basename %>% 
  sub('\\.zip.', '', .)

id.70 <- h %>% 
  html_nodes(xpath = "//a[text()='bi']/@href") %>% 
  as.character %>% 
  grep('85bi70', ., value = TRUE) %>% 
  basename %>% 
  sub('\\.zip.', '', .)

gcms.50 <- cbind(gcms, id.50)
gcms.50$GCM = sub(" \\(#\\)", "", gcms$GCM)  ## sub replaces first instance in a string, gsub = global

gcms.70 <- cbind(gcms, id.70)
gcms.70$GCM = sub(" \\(#\\)", "", gcms$GCM) 


## Now create the scenario lists across which to loop
gcms.50 ; gcms.70


## Just get the 6 models picked by CSIRO for Australia, for 2050 and 2070
## Also will be a better way to get at this...
scen_2050 = c("mc85bi50", "no85bi50", "ac85bi50", "cn85bi50", "gf85bi50", "hg85bi50")
scen_2070 = c("mc85bi70", "no85bi70", "ac85bi70", "cn85bi70", "gf85bi70", "hg85bi70")


#########################################################################################################################
## Then create a stack of current environmental conditions, and an Australia shapefile for the mapping later...
aus <- ne_states(country = 'Australia') %>% 
  subset(!grepl('Island', name))

env.grids.current <- stack(
  file.path('./data/base/worldclim/aus/0.5/bio/current',
            sprintf('bio_%02d.tif', 1:19)))

## Now divide the current temperature grids by 10
for(i in 1:11) {
  
  message(i)
  env.grids.current[[i]] <- env.grids.current[[ i]]/10
  
}


## What do the grids look like? The 11,177,684 NAs values could be the ocean
summary(env.grids.current[[1]])
summary(env.grids.current[[5]])
summary(env.grids.current[[11]])





#########################################################################################################################
## 2). PROJECT MAXENT MODELS FOR MULTIPLE CLIMATE SCEANARIOS, currently 2050 AND 2070
#########################################################################################################################


## Covert the test species into directory names
test_spp = gsub(" ", "_", test.spp)


#########################################################################################################################
## For each species, use a function to create raster files and maps of all six GCMs.
## Note that some of the experimental species - e.g. Kennedia_beckxiana - still have to be modelled
env.grids.2050 = project.grids.2050(scen_2050, test_spp)
env.grids.2070 = project.grids.2070(scen_2070, test_spp)


#########################################################################################################################
## Then, make a list all the directories containing the individual GCM rasters
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
## 3). CREATE TABLE OF MAXENT RESULTS FOR CURRENT CONDITIONS
#########################################################################################################################


#########################################################################################################################
## First, read in the list of files for the current models, and specify the file path
table.list = list.files("./output/maxent/STD_VAR_ALL/")
path       = "./output/maxent/STD_VAR_ALL/"


## Could turn this into a function, and loop over a list of subfolders...
MAXENT.STD.VAR.SUMMARY <- table.list[c(1:length(table.list))] %>%
  
  ## pipe the list into lapply
  lapply(function(x) {
    
    ## create the character string
    f <- paste0(path, x, "/full/maxentResults.csv")
    
    ## load each .RData file
    d <- read.csv(f)
    
    ## now add a model column
    d = cbind(GBIF_Taxon = x, Model_run  = path, d) 
    dim(d)
    
    ## Remove path gunk, and species
    #d$GBIF_Taxon = gsub("_", " ", d$GBIF_Taxon)
    d$Model_run  = gsub("./output/maxent/", "", d$Model_run)
    d$Model_run  = gsub("/", "", d$Model_run)
    d$Species    = NULL
    d
  }) %>%
  
  ## finally, bind all the rows together
  bind_rows


## This is a summary of maxent output for current conditions
dim(MAXENT.STD.VAR.SUMMARY)
head(MAXENT.STD.VAR.SUMMARY)[1:8]


## Now check the match between the species list, and the results list. These need to match, so we can access
## the right threshold for each species.
length(intersect(test_spp, MAXENT.STD.VAR.SUMMARY$GBIF_Taxon)) ## accesssing the files from these directories... 
MAXENT.SUM.TEST  =  MAXENT.STD.VAR.SUMMARY[MAXENT.STD.VAR.SUMMARY$GBIF_Taxon %in% test_spp, ] 
head(MAXENT.SUM.TEST, 10)["Maximum.training.sensitivity.plus.specificity.Logistic.threshold"];head(MAXENT.SUM.TEST, 10)[1]
comb_spp = unique(MAXENT.SUM.TEST$GBIF_Taxon)


#########################################################################################################################
## Save results::
write.csv(MAXENT.STD.VAR.SUMMARY, "./output/maxent/MAXENT_STD_VAR_SUMMARY.csv", row.names = FALSE)





#########################################################################################################################
## 4). SUMMARIZE MAXENT RESULTS FOR EACH SPECIES ACROSS MULTIPLE GCMs
#########################################################################################################################


#########################################################################################################################
## Now for each species, use a function to combine the raster files for each climate scenario into one layer. Use the 
## species-specific thresholds that are produced by rmaxent (i.e. the columns in MAXENT.STD.VAR.SUMMARY). Linda and John 
## use these two:


## Maximum training sensitivity plus specificity Logistic threshold 
## 10 percentile training presence training omission. EG:
summary(MAXENT.SUM.TEST["Maximum.training.sensitivity.plus.specificity.Logistic.threshold"])
summary(MAXENT.SUM.TEST["X10.percentile.training.presence.training.omission"])


## Turn the maxent results into lists
thresh.max.train  = as.list(MAXENT.SUM.TEST["Maximum.training.sensitivity.plus.specificity.Logistic.threshold"]) 
thresh.max.train  = thresh.max.train$Maximum.training.sensitivity.plus.specificity.Logistic.threshold

percent.10.omiss  = as.list(MAXENT.SUM.TEST["X10.percentile.training.presence.training.omission"])
percent.10.omiss  = percent.10.omiss$X10.percentile.training.presence.training.omission
length(percent.10.omiss);length(thresh.max.train) 


#########################################################################################################################
## So we can use the values in these columns to threhsold the rasters of habitat suitability (0-1) when combining them.
## For each species, we will create a binary raster with cell values between 0-6. These cell values represent the number of GCMs 
## where that cell had a suitability value above the threshold determined by maxent (e.g. the max training sensitivity + specif).  


## These thresholded predictions of habitat suitability could be used to determine the loss or gain of species within areal units
## (e.g. significant urban areas or LGAs), between time periods (current, 2030, 2070).


## Why can't we , define the SUA in the global environment, outside the function?
# SUA = readOGR("F:/green_cities_sdm/data/base/CONTEXTUAL/IN.SUA.shp", layer = "IN.SUA")
# CRS.new  <- CRS("+init=epsg:4326") # EPSG:3577
# SUA.WGS  = spTransform(SUA, CRS.new)
# writeOGR(obj = SUA.WGS, dsn = "./data/base/CONTEXTUAL", layer = "IN_SUA_WGS", driver = "ESRI Shapefile")
# 
# areal_unit = SUA.WGS


#########################################################################################################################
## Combine output and calculate gain and loss for 2050 
DIR        = SDM.RESULTS.DIR[1] 
species    = comb_spp[1] 
thresh     = thresh.max.train[1] 
percent    = percent.10.omiss[1]
time_slice = 50


suitability.2050 = mapply(combine_gcm_threshold, 
                          DIR_list     = SDM.RESULTS.DIR[1], 
                          species_list = comb_spp[1], 
                          thresholds   = thresh.max.train[1], 
                          percentiles  = percent.10.omiss[1],
                          time_slice   = 50)


## Combine output and calculate gain and loss for 2070 
suitability.2070 = mapply(combine_gcm_threshold, 
                          DIR_list     = SDM.RESULTS.DIR[1], 
                          species_list = comb_spp[1], 
                          thresholds   = thresh.max.train[1], 
                          percentiles  = percent.10.omiss[1],
                          time_slice   = 70)


## Combine output and calculate gain and loss for 2030 
# suitability.2030 = mapply(combine_gcm_threshold, 
#                           DIR_list     = SDM.RESULTS.DIR, 
#                           species_list = comb_spp, 
#                           thresholds   = thresh.max.train, 
#                           percentiles  = percent.10.omiss,
#                           time_slice   = 30)





#########################################################################################################################
## 5). COMBINE MAXENT TABLES FOR ALL SPECIES ACROSS MULTIPLE GCMs
#########################################################################################################################


## Could turn this into a function, and loop over a list of subfolders...
MAXENT.STD.VAR.SUMMARY <- table.list[c(1:length(table.list))] %>%
  
  ## pipe the list into lapply
  lapply(function(x) {
    
    ## create the character string
    f <- paste0(path, x, "/full/maxentResults.csv")
    
    ## load each .RData file
    d <- read.csv(f)
    
    ## now add a model column
    d = cbind(GBIF_Taxon = x, Model_run  = path, d) 
    dim(d)
    
    ## Remove path gunk, and species
    #d$GBIF_Taxon = gsub("_", " ", d$GBIF_Taxon)
    d$Model_run  = gsub("./output/maxent/", "", d$Model_run)
    d$Model_run  = gsub("/", "", d$Model_run)
    d$Species    = NULL
    d
  }) %>%
  
  ## finally, bind all the rows together
  bind_rows


## This is a summary of maxent output for current conditions
dim(MAXENT.STD.VAR.SUMMARY)
head(MAXENT.STD.VAR.SUMMARY)[1:8]





#########################################################################################################################
## OUTSTANDING PREDICTION TASKS:
#########################################################################################################################


#########################################################################################################################
## The no-data areas need to be masked out as part of the calculations.

## Decide on the method for combining the GCMs

## We are missing two scenarios recommended for Australia: CanESM2 & CESM1-CAM5. These two are not on the worldclim list:

## https://www.climatechangeinaustralia.gov.au/en/support-and-guidance/faqs/eight-climate-models-data/

## http://www.worldclim.org/cmip5_30s

## Consider the final format, which files: table, plots/maps, files. Disk space important for both local and web... 
## What are the options for presenting consensus layers of all scenarios? E.G. Are there some templates that John
## has used for previous work?


#########################################################################################################################
## When we subtract the current binary layer (1, 0) from the future binary layer, we have: 

## F0 - C0 =  0 (no data in either layer)
## F0 - C1 = -1 (LOSS across all GCMs)
## F1 - C0 =  1 (GAIN across all GCMs) 
## F1 - C1 =  0 (NO CHANGE: also no data before the overlay...)


## However, this changes when we subtract the current binary layer (taking values 1, 0) from the future integer layer 
## (taking values from 0 to 6). Here, any negative value (-1) is a loss: the species was predicted to be present in that cell 
## based on current conditions, but predicted to be absent under future conditions (i.e. the future layer did not meet the suitabiity
## threshold in that location at that time.

## However, positive values (1-6) can either be a gain, or, no change. The difference needs to be accounted for, because otherwise
## we don't know if the species was predicted to occur. Can we fix this by masking the overlay to just cells with data?
## Effectively this means excluding 0 values from the overlay.

## Not all these possibilities will occur, but they are : 

## F0 - C0 =  0 no data in either layer
## F0 - C0 =  0 (NO CHANGE according to all GCMs)
## F0 - C1 =  -1 (LOSS according to all GCMs)

## F1 - C1 =  0 (NO CHANGE according to one GCM: also, no data before the overlay)
## F1 - C0 =  1 (GAIN according to one GCM)

## F2 - C1 =  1 (NO CHANGE, according to two GCMs)
## F2 - C0 =  2 (GAIN according to two GCMs)

## F3 - C1 =  2 (NO CHANGE, according to three GCMs)
## F3 - C0 =  3 (GAIN, according to three GCMs)

## F4 - C1 =  3 (NO CHANGE, according to four GCMs)
## F4 - C0 =  4 (GAIN, according to four GCMs)

## F5 - C1 =  4 (NO CHANGE, according to five GCMs)
## F5 - C0 =  5 (GAIN according to five GCMs)

## F6 - C1 =  5 (NO CHANGE, according to six GCMs?)
## F6 - 0  =  6 (GAIN, according to six GCMs?)





#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################