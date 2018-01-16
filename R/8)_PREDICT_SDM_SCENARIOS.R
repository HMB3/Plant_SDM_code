#########################################################################################################################
###################### PREDICT MAXENT TO FUTURE CLIMATES AND SUMMARISE THE RESULTS ###################################### 
#########################################################################################################################


#########################################################################################################################
## This code takes the output of the Maxent species distribution models (i.e. using current conditions), and generates 
## a prediction of habitat suitability for current and future environmental conditions. The input data table is in the 
## format of all species occurrences (rows) and environmental variables (columns).

## Finally, the predictions from all GCMs (currently using six) are combined into a single habitat suitability layer


#########################################################################################################################
## 1). READ IN THE CURRENT AND FUTURE ENVIRONMENTAL DATA FOR MULTIPLE GCMs
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


## Create a lookup table of GCMs using the WORLDCLIM website
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


## What do the grids look like? 11,177,684 NAs???!!!!!
summary(env.grids.current[[1]])
summary(env.grids.current[[5]])
summary(env.grids.current[[11]])





#########################################################################################################################
## 2). PROJECT MAXNET MODELS FOR MULTIPLE CLIMATE SCEANARIOS currently 2050 AND 2070
#########################################################################################################################


## Make the test species directory names
test_spp = gsub(" ", "_", test.spp)


#########################################################################################################################
## For each species, use a function to create raster files and maps of all six climate scenarios (GCMs).
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
## First, Read in the list of files for the current models, and specify the file path
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


## Now check the match between the species list, and the results this. These need to match, so we can access
## the right threshold for each species.
intersect(test_spp, MAXENT.STD.VAR.SUMMARY$GBIF_Taxon) ## accesssing the files from these directories... 



#########################################################################################################################
## 4). SUMMARIZE MAXENT RESULTS FOR EACH SPECIES ACROSS MULTIPLE GCMs
#########################################################################################################################


#########################################################################################################################
## Now for each species, use a function to combine the raster files for each climate scenario into one layer

## Here we are using the species-specific thresholds that are produced by rmaxent (i.e. the columns in MAXENT.STD.VAR.SUMMARY)
## Linda and John use these two:

## Maximum training sensitivity plus specificity Logistic threshold
## 10 percentile training presence training omission. EG:
summary(MAXENT.STD.VAR.SUMMARY["Maximum.training.sensitivity.plus.specificity.Logistic.threshold"])
summary(MAXENT.STD.VAR.SUMMARY["X10.percentile.training.presence.training.omission"])


## Turn these into lists
thresh.max.train = as.list(MAXENT.STD.VAR.SUMMARY["Maximum.training.sensitivity.plus.specificity.Logistic.threshold"]) 
thresh.max.train = thresh.max.train$Maximum.training.sensitivity.plus.specificity.Logistic.threshold

thresh.10.omiss  = as.list(MAXENT.STD.VAR.SUMMARY["X10.percentile.training.presence.training.omission"])
thresh.10.omiss  = thresh.10.omiss$X10.percentile.training.presence.training.omission


## So we can use the values in these columns to threhsold the rasters of habitat suitability (0-1) when combining them.
## For each species, we will create a binary raster with cell values between 0-6. These cell values represent the number of GCMs 
## where that cell had a suitability value above the threshold determined by maxent (e.g. the max training sensitivity + specif).  


## These thresholded predictions of habitat suitability could be used to determine the loss or gain of species within areal units
## (e.g. significant urban areas or LGAs), between time periods (current, 2030, 2070). So the ingredients here are the results 
## for current and future predictions. 


## Now change the function to use the species' thresholds
suitability.2050 = mapply(combine_maxent_predictions, SDM.RESULTS.DIR, test_spp, period = 50)
suitability.2070 = mapply(combine_maxent_predictions, SDM.RESULTS.DIR, test_spp, period = 70)


#########################################################################################################################
## Save results
write.csv(MAXENT.STD.VAR.SUMMARY, "./output/maxent/MAXENT_STD_VAR_SUMMARY.csv", row.names = FALSE)





#########################################################################################################################
## OUTSTANDING PREDICTION TASKS:
#########################################################################################################################


## We are missing two scenarios recommended for Australia: CanESM2 & CESM1-CAM5. These two are not on the worldclim list:

## https://www.climatechangeinaustralia.gov.au/en/support-and-guidance/faqs/eight-climate-models-data/

## http://www.worldclim.org/cmip5_30s
  

## Consider the final format, which files: table, plots/maps, files. Disk space important for both local and web... 


## What are the options for presenting consensus layers of all scenarios? EG show the occurrences, the average and then 
## the combined layers? Need to choose the format first, as the details are tricky (e.g. scale bars, etc.)





#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################