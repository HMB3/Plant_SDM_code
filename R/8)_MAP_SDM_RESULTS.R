#########################################################################################################################
###################### PREDICT MAXENT TO FUTURE CLIMATES AND SUMMARISE THE RESULTS ###################################### 
#########################################################################################################################


#########################################################################################################################
## This code takes the output of the Maxent species distribution models (i.e. using current conditions), and generates 
## a prediction of habitat suitability for current and future environmental conditions. The input data table is in the 
## format of all species occurrences (rows) and environmental variables (columns).


## The predictions from all 6 GCMs are then combined into a single habitat suitability layer.
## And the total area of habitat gained, lost or remaining stable is calculated (i.e. for AUS). 


## Using this combined layer, the % of area occupied by species within areal units (significant urban areas or SUAs), 
## under each projection (2030, 2050 and 2070) is also be calculated. 


#########################################################################################################################
## Read in all data to run the SDM code :: species lists, shapefile, rasters & tables
#source('./R/HIA_LIST_MATCHING.R')
rasterTmpFile()


#########################################################################################################################
## 1). READ IN THE CURRENT AND FUTURE ENVIRONMENTAL DATA FOR SIX GCMs
#########################################################################################################################


#########################################################################################################################
## Create a list of GCM scenarios, which are used to create maps of habitat suitability 


## Eight of the 40 CMIP5 models assessed in this project have been selected for use in provision of application-ready data. 
## This facilitates efficient exploration of climate projections for Australia.
## https://www.climatechangeinaustralia.gov.au/en/support-and-guidance/faqs/eight-climate-models-data/

## The full list is:
# scen = c("ip85bi50", "mc85bi50", "mg85bi50", "mi85bi50", "mp85bi50", 
#          "mr85bi50", "no85bi50", "ac85bi50", "bc85bi50", "cc85bi50", 
#          "cn85bi50", "gf85bi50", "gs85bi50", "hd85bi50", "he85bi50",
#          "hg85bi50", "in85bi50")


## We can only get the bioclim variables for 6 of these projections......................................................
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


## Work around because the 2030 data comes from CC in Aus, not worldclim
id.30 = gsub("50", "30", id.50)


## Create the IDs
gcms.30 <- cbind(gcms, id.30)
gcms.30$GCM = sub(" \\(#\\)", "", gcms$GCM)

gcms.50 <- cbind(gcms, id.50)
gcms.50$GCM = sub(" \\(#\\)", "", gcms$GCM)  ## sub replaces first instance in a string, gsub = global

gcms.70 <- cbind(gcms, id.70)
gcms.70$GCM = sub(" \\(#\\)", "", gcms$GCM) 


## Now create the scenario lists across which to loop
gcms.50 ; gcms.70 ; gcms.30


## Just get the 6 models picked by CSIRO for Australia, for 2030, 2050 and 2070
scen_2030 = c("mc85bi30", "no85bi30", "ac85bi30", "cc85bi30", "gf85bi30", "hg85bi30")
scen_2050 = c("mc85bi50", "no85bi50", "ac85bi50", "cc85bi50", "gf85bi50", "hg85bi50")
scen_2070 = c("mc85bi70", "no85bi70", "ac85bi70", "cc85bi70", "gf85bi70", "hg85bi70")


## Then create a stack of current environmental conditions outside the function, and an Australia shapefile for the mapping later...
# aus <- ne_states(country = 'Australia') %>% 
#   subset(!grepl('Island', name))
# 
# shapefile = aus

## Now divide the current environmental grids by 10
env.grids.current <- stack(
  file.path('./data/base/worldclim/aus/1km/bio/current',   ## ./data/base/worldclim/aus/1km/bio
            sprintf('bio_%02d.tif', 1:19)))

for(i in 1:11) {
  
  ## simple loop
  message(i)
  env.grids.current[[i]] <- env.grids.current[[i]]/10
  
}


#########################################################################################################################
## Name the environmental grids to be used in the mapping code :: this is all 19, because we are using the directory structure
grid.names = c('Annual_mean_temp',    'Mean_diurnal_range',  'Isothermality',      'Temp_seasonality', 
               'Max_temp_warm_month', 'Min_temp_cold_month', 'Temp_annual_range',  'Mean_temp_wet_qu',
               'Mean_temp_dry_qu',    'Mean_temp_warm_qu',   'Mean_temp_cold_qu',  'Annual_precip',
               'Precip_wet_month',    'Precip_dry_month',    'Precip_seasonality', 'Precip_wet_qu',
               'Precip_dry_qu',       'Precip_warm_qu',      'Precip_col_qu')





#########################################################################################################################
## 2). PROJECT MAXENT MODELS FOR MULTIPLE CLIMATE SCEANARIOS AT 2030, 2050 AND 2070
#########################################################################################################################


#########################################################################################################################
## For each species, use a function to create raster files and maps under all six GCMs at each time step
## First remove species without data from step 7
#map_spp_list = combo_spp
no_data      <- c ("Baeckea_virgata", "Kennedia_beckxiana", "Grevillea_rivularis", "Arctostaphylos_densiflora", 
                   "Cupressocyparis_leylandii", "Eucalyptus_intermedia", "Ficus_hillii", "Pentaceras_australi", 
                   "Pentaceras_australis", "Pouteria_australis", "Pouteria_chartacea", "Pouteria_eerwah", 
                   "Radermachera_gigantea", "Randia_benthamiana", "Raphiolepis_umbellata", "Tilia_mongolica", 
                   "Trema_aspera", "Xanthostemon_verticillatus")
map_spp_list     <- map_spp_list [! map_spp_list %in% no_data]
no_data %in% map_spp_list



#########################################################################################################################
## Create 2030 maps
env.grids.2030 = tryCatch(project_maxent_grids(scen_list     = scen_2030,
                                               species_list  = map_spp_list,
                                               maxent_path   = "./output/maxent/SET_VAR_KOPPEN",
                                               climate_path  = "./data/base/worldclim/aus/1km/bio",
                                               grid_names    = grid.names,
                                               time_slice    = 30,
                                               current_grids = env.grids.current),
                          
                          ## Will this work outside a loop?
                          error = function(cond) {
                            
                            message(paste('Species skipped - check', spp))
                            
                          })


#########################################################################################################################
## Create 2050 maps
env.grids.2050 = tryCatch(project_maxent_grids(scen_list     = scen_2050,
                                               species_list  = map_spp_list,
                                               time_slice    = 50,
                                               maxent_path   = "./output/maxent/SET_VAR_KOPPEN",
                                               climate_path  = "./data/base/worldclim/aus/1km/bio",
                                               grid_names    = grid.names,
                                               current_grids = env.grids.current),
                          
                          error = function(cond) {
                            
                            message(paste('Species skipped - check', spp))
                            
                          })


#########################################################################################################################
## Create 2070 maps
env.grids.2070 = tryCatch(project_maxent_grids(scen_list     = scen_2070,
                                               species_list  = map_spp_list,
                                               time_slice    = 70,
                                               maxent_path   = "./output/maxent/SET_VAR_KOPPEN",
                                               climate_path  = "./data/base/worldclim/aus/1km/bio",
                                               grid_names    = grid.names,
                                               current_grids = env.grids.current),
                          
                          error = function(cond) {
                            
                            message(paste('Species skipped - check inputs', spp))
                            
                          })





#########################################################################################################################
## 3). CREATE TABLE OF MAXENT RESULTS FOR CURRENT CONDITIONS
#########################################################################################################################


#########################################################################################################################
## First, read in the list of files for the current models, and specify the file path
path.set.var             = "./output/maxent/SET_VAR_KOPPEN/"
#map_spp_list                 = combo_spp 


## Create an object for the maxent settings :: using the same variable for every model
model.selection.settings = "Set_variables"  
records_setting          = "COORD_CLEAN"


## Create a file list for each model run
maxent.tables = list.files(path.set.var)             ## Chagne this for each variable selection strategy
maxent.tables = intersect(maxent.tables, map_spp_list)   ## Change this for new species lists
maxent_path   = path.set.var                         ## Chagne this for each variable selection strategy
length(maxent.tables)                                ## Should match the number of taxa tested
no_data %in% maxent.tables

## In linux, check if the folders are empty, then delete if empty as this will break the code:
# cd F:/green_cities_sdm/output/maxent/STD_VAR_ALL
# find . -type d -empty -print
# find . -type d -empty -delete


## Could turn this into a function, and loop over a list of subfolders...
## Then pipe the table list into lapply
## MAXENT.SUMMARY = bind_maxent_tables(table_list)
MAXENT.SUMMARY <- maxent.tables[c(1:length(maxent.tables))] %>%         ## currently x species with enough records
  
  ## pipe the list into lapply
  lapply(function(x) {
    
    ## Create the character string
    f <- paste0(maxent_path, x, "/full/maxentResults.csv")
    
    ## Load each .RData file
    d <- read.csv(f)
    
    ##
    #m <- readRDS(sprintf('%s/%s/full/model.rds', maxent_path, x))
    m <- readRDS(sprintf('%s/%s/full/maxent_fitted.rds', maxent_path, x))
    m <- m$me_full
    number_var = length(m@lambdas) - 4                                   ## (the last 4 slots of lambdas are not variables) 
    
    ## Now add a model column
    d = cbind(searchTaxon = x, 
              Settings    = model.selection.settings, 
              Number_var  = number_var, 
              Records     = records_setting, d)  ## see step 7, make a variable for multiple runs
    dim(d)
    
    ## Remove path gunk, and species
    d$Species    = NULL
    d
    
  }) %>%
  
  ## Finally, bind all the rows together
  bind_rows


## This is a summary of maxent output for current conditions
## Also which species have AUC < 0.7?
dim(MAXENT.SUMMARY)
head(MAXENT.SUMMARY, 20)[1:9]
dim(subset(MAXENT.SUMMARY, Training.AUC < 0.7))

## What are the variables we want to see?
# View(head(MAXENT.SUMMARY, 180)[, c("searchTaxon",
#                                            "Settings",
#                                            "Number_var",
#                                            "X.Training.samples",                                                                
#                                            "Iterations",                                                                        
#                                            "Training.AUC",                                                                      
#                                            "X.Background.points",  
#                                            "Maximum.training.sensitivity.plus.specificity.Logistic.threshold")])


## Now check the match between the species list, and the results list. These need to match, so we can access
## the right threshold for each species.
length(intersect(map_spp_list, MAXENT.SUMMARY$searchTaxon)) ## Accesssing the files from these directories... 
MAXENT.SUM.TEST  =  MAXENT.SUMMARY[MAXENT.SUMMARY$searchTaxon %in% map_spp_list , ] 
map_spp = unique(MAXENT.SUM.TEST$searchTaxon)
length(map_spp)


#########################################################################################################################
## Then, make a list of all the directories containing the individual GCM rasters...path.backwards.sel
SDM.RESULTS.DIR <- map_spp[c(1:length(map_spp))] %>%
  
  ## Pipe the list into lapply
  lapply(function(species) {
    
    ## Create the character string...
    m <-   sprintf('%s%s/full/', path.set.var, species)                ## path.backwards.sel
    m 
    
  }) %>%
  
  ## Bind the list together
  c()


#########################################################################################################################
## Now combine the SDM output with the niche context data 
## Get the number of aus records too ....................................................................................
NICHE.CONTEXT = COMBO.NICHE.CONTEXT[, c("searchTaxon",      "COMBO.count",       "AUS_RECORDS",       "Plant.type",        "Origin", 
                                        "Top_200",          "Total.growers",     "Number.of.States")]


## Check with John and Linda which columns will help with model selection
MAXENT.SUMM   = MAXENT.SUMMARY[, c("searchTaxon",
                                   "Settings",
                                   "Records",
                                   "Number_var",
                                   "X.Training.samples",                                                                
                                   "Iterations",                                                                        
                                   "Training.AUC",                                                                      
                                   "X.Background.points",  
                                   "Maximum.training.sensitivity.plus.specificity.Logistic.threshold")]


## Remove the underscore, and join
MAXENT.SUMM$searchTaxon = gsub("_", " ", MAXENT.SUMM$searchTaxon)
MAXENT.CHECK.TABLE      = join(NICHE.CONTEXT, MAXENT.SUMM, type = "inner")
View(MAXENT.CHECK.TABLE)


## Also could join on other tables with different settings using rbind
## rbind(MAXENT.CHECK.TABLE, MAXENT.CHECK.TABLE.SET)
## order by species and compare three rows : setvariables, backwards selection, etc


## Save - could add date as a sprintf variable to save multiple versions?
## write.csv(MAXENT.CHECK.TABLE, "./output/maxent/MAXENT_CHECK_TABLE_APRIL_2016.csv", row.names = FALSE)





#########################################################################################################################
## 4). CREATE LISTS OF HABITAT SUITABILITY THRESHOLDS
#########################################################################################################################


#########################################################################################################################
## Maxent produces a presence threshold for each species (i.e. the columns in MAXENT.SUMMARY). 
## The trouble here is that we might need to change the threshold for different species, rather than using the same one 
## for all of them. That changes the order of lists, which is a problem for looping over them.

## Read in the table of checked maps. Then subset to just the species with dodgy maps
#MAXENT.CHECK      = read.csv("./output/maxent/MAXENT_CHECK_RATING.csv", stringsAsFactors = FALSE)
MAXENT.CHECK   = read.csv("./output/maxent/MAXENT_RATING_26_2018.csv", stringsAsFactors = FALSE)
MAXENT.CHECK   = join(MAXENT.CHECK, TOT.GROW)
MAXENT.CHECK   = MAXENT.CHECK [, c(1:7, 19, 8:18)]
MAXT.CHECK.25  = subset(MAXENT.CHECK, Total.growers >= 25 & CHECK_MAP == 1 | CHECK_MAP == 2)
MAXT.CHECK.25  = completeFun(MAXT.CHECK.25, "Total.growers")

#MAXT.CHECK.25  = head(MAXT.CHECK.25, 150)
table(MAXT.CHECK.25$CHECK_MAP)
MAXT.CHECK.25 = MAXT.CHECK.25[with(MAXT.CHECK.25 , rev(order(Total.growers))), ]
View(MAXT.CHECK.25)
dim(MAXT.CHECK.25)
summary(MAXT.CHECK.25$Total.growers)


## Now write out the species list to re-process
#write.csv(MAXT.CHECK.25, "./output/maxent/MAXENT_SUA_SPP.csv", row.names = FALSE)



#########################################################################################################################
## Isolate the species which have bad maps
# MX.POOR        = subset(MAXENT.CHECK, CHECK_MAP == 2 | CHECK_MAP == 3)#[, c("searchTaxon", "CHECK_MAP", "Top_200", "Origin")]
# MX.POOR.TOP    = subset(MX.POOR,     Top_200 == "TRUE")[, c("searchTaxon", "CHECK_MAP")]
# MX.POOR.E      = subset(MX.POOR,     Origin == "Exotic")[, c("searchTaxon", "CHECK_MAP")]
# MX.POOR.TOP.E  = subset(MX.POOR,     Origin == "Exotic" & Top_200 == "TRUE")[, c("searchTaxon", "CHECK_MAP")]
# 
# 
# ## Merge with full results
# MX.POOR.TOP.E = join(MX.POOR.TOP.E, MAXENT.CHECK.TABLE, type = "left")
# MX.POOR.TOP.E = MX.POOR.TOP.E[with(MX.POOR.TOP.E , rev(order(Total.growers))), ]
# 
# 
# 
# MX.POOR        = MX.POOR[with(MX.POOR , rev(order(Total.growers))), ]
# MX.POOR.TOP    = MX.POOR.TOP[with(MX.POOR.TOP , rev(order(Total.growers))), ]
# MX.POOR.E      = MX.POOR.E[with(MX.POOR.E , rev(order(Total.growers))), ]
# MX.POOR.TOP.E  = MX.POOR.TOP.E[with(MX.POOR.TOP.E , rev(order(Total.growers))), ]
# 
# View(MX.POOR.E)
# View(MX.POOR.TOP.E)


## How do we use this to check if the poor maps are due to errors in the code?


#########################################################################################################################
## Now create a list of thresholds to loop over using the mapping functions
spp.lower.thresh  = subset(MAXENT.CHECK, CHECK_MAP == 2 | CHECK_MAP == 3)$searchTaxon
spp_lower_thresh  = gsub(" ", "_", spp.lower.thresh)

spp.best.thresh   = subset(MAXENT.CHECK, CHECK_MAP == 1 | CHECK_MAP == 2)$searchTaxon
spp_best_thresh   = gsub(" ", "_", spp.best.thresh)

MAXENT.LOWER      = MAXENT.SUM.TEST[MAXENT.SUM.TEST$searchTaxon %in% spp_lower_thresh, ] 
MAXENT.BEST       = MAXENT.SUM.TEST[MAXENT.SUM.TEST$searchTaxon %in% spp_best_thresh, ] 
identical(spp_lower_thresh, MAXENT.LOWER$searchTaxon)
identical(spp_best_thresh,  MAXENT.BEST$searchTaxon)


## John : for AUC you can report the cross-validated test AUC (if your code currently runs a cross-validated model as well), 
## and for the model threshold (for binarising) you can just use the training value (or the crossval one...there's little 
## guidance about this and you can really get away with either).


## How do the thresholds compare?
summary(MAXENT.SUM.TEST["Maximum.training.sensitivity.plus.specificity.Logistic.threshold"])   ## .training. should be .test.
summary(MAXENT.SUM.TEST["X10.percentile.training.presence.Logistic.threshold"])
summary(MAXENT.SUM.TEST["X10.percentile.training.presence.training.omission"])


## Turn the maxent results into lists :: we can use these to generate the consensus layers 
thresh.max.train       = as.list(MAXENT.SUM.TEST["Maximum.training.sensitivity.plus.specificity.Logistic.threshold"]) 
thresh.max.train       = thresh.max.train$Maximum.training.sensitivity.plus.specificity.Logistic.threshold

thresh.max.train.best  = as.list(MAXENT.BEST["Maximum.training.sensitivity.plus.specificity.Logistic.threshold"]) 
thresh.max.train.best  = thresh.max.train.best$Maximum.training.sensitivity.plus.specificity.Logistic.threshold

percent.10.log         = as.list(MAXENT.SUM.TEST["X10.percentile.training.presence.Logistic.threshold"])  ## for the twos 
percent.10.log.low     = as.list(MAXENT.LOWER["X10.percentile.training.presence.Logistic.threshold"])
percent.10.log.best    = as.list(MAXENT.BEST["X10.percentile.training.presence.Logistic.threshold"])

percent.10.log         = percent.10.log$X10.percentile.training.presence.Logistic.threshold
percent.10.log.low     = percent.10.log.low$X10.percentile.training.presence.Logistic.threshold
percent.10.log.best    = percent.10.log.best$X10.percentile.training.presence.Logistic.threshold

percent.10.om          = as.list(MAXENT.SUM.TEST["X10.percentile.training.presence.training.omission"])   ## discount
percent.10.om.low      = as.list(MAXENT.LOWER["X10.percentile.training.presence.training.omission"])
percent.10.om.best     = as.list(MAXENT.BEST["X10.percentile.training.presence.training.omission"])

percent.10.om          = percent.10.om$X10.percentile.training.presence.training.omission
percent.10.om.low      = percent.10.om.low$X10.percentile.training.presence.training.omission
percent.10.om.best     = percent.10.om.best$X10.percentile.training.presence.training.omission


#########################################################################################################################
## Then, make a list all the directories containing the individual GCM rasters...path.backwards.sel
SDM.RESULTS.DIR.LOW <- spp_lower_thresh [c(1:length(spp_lower_thresh ))] %>%
  
  ## Pipe the list into lapply
  lapply(function(species) {
    
    ## Create the character string...
    m <-   sprintf('%s%s/full/', path.set.var, species)                ## path.backwards.sel
    m 
    
  }) %>%
  
  ## Bind the list together
  c()


## The good models
SDM.RESULTS.DIR.BEST <- spp_best_thresh [c(1:length(spp_best_thresh ))] %>%
  
  ## Pipe the list into lapply
  lapply(function(species) {
    
    ## Create the character string...
    m <-   sprintf('%s%s/full/', path.set.var, species)                ## path.backwards.sel
    m 
    
  }) %>%
  
  ## Bind the list together
  c()



## Check the order of lists match, species, SUAs, areas need to match up ................................................
## It would be safer to read in the thresholds individually, so they match the species folder exactly
length(SDM.RESULTS.DIR);length(MAXENT.SUM.TEST$searchTaxon);length(thresh.max.train);
length(percent.10.log);length(percent.10.om);length(map_spp)
identical(MAXENT.SUM.TEST$searchTaxon, map_spp)


## Also, check that the lower threshold lists are ok too
length(percent.10.log.low);length(percent.10.om.low);length(spp_lower_thresh)


## The order of the directories matches
head(SDM.RESULTS.DIR, 20);head(map_spp, 20); head(MAXENT.SUM.TEST, 20)[, c("searchTaxon",
                                                                            "Maximum.training.sensitivity.plus.specificity.Logistic.threshold", 
                                                                            "X10.percentile.training.presence.Logistic.threshold")]

tail(SDM.RESULTS.DIR, 20);tail(map_spp, 20); tail(MAXENT.SUM.TEST, 20)[, c("searchTaxon",
                                                                            "Maximum.training.sensitivity.plus.specificity.Logistic.threshold", 
                                                                            "X10.percentile.training.presence.Logistic.threshold")]





#########################################################################################################################
## 5). SUMMARIZE MAXENT RESULTS FOR EACH SPECIES ACROSS MULTIPLE GCMs
#########################################################################################################################


#########################################################################################################################
## So we can use the values in these columns to threhsold the rasters of habitat suitability (0-1) when combining them.
## For each species, we will create a binary raster with cell values between 0-6. These cell values represent the number of GCMs 
## where that cell had a suitability value above the threshold determined by maxent (e.g. the max training sensitivity + specif).


#########################################################################################################################
## Loop over directories, species and one threshold for each, also taking a time_slice argument.
# DIR        = SDM.RESULTS.DIR[3] 
# species    = map_spp[3] 
# thresh     = thresh.max.train[3] 
# percent    = percent.10.log[3]
# time_slice = 30
# area_occ   = 10


# Running zonal stats for Strelitzia_reginae | 2030 combined suitability > 0.4016
# Error in names(PERECENT.AREA) <- c("SUA_NAME11", "Absent", "Present") : 
#   'names' attribute [3] must be the same length as the vector [2]
# In addition: There were 32 warnings (use warnings() to see them)


#########################################################################################################################
## Use the strictest threshold first
#########################################################################################################################


#########################################################################################################################
## Combine output and calculate gain and loss for 2030 
# suitability.2030 = tryCatch(mapply(combine_gcm_threshold,
#                                    DIR_list     = SDM.RESULTS.DIR,
#                                    species_list = map_spp,
#                                    maxent_path  = "./output/maxent/SET_VAR_KOPPEN",
#                                    thresholds   = thresh.max.train,
#                                    percentiles  = percent.10.log,
#                                    time_slice   = 30,
#                                    area_occ     = 10),
#                             
#                             error = function(cond) {
#                               
#                               message(paste('Species skipped - check inputs', spp))
#                               
#                             })
# 
# 
# ## Combine GCM output for 2050 
# suitability.2050 = tryCatch(mapply(combine_gcm_threshold, 
#                                    DIR_list     = SDM.RESULTS.DIR, 
#                                    species_list = map_spp, 
#                                    maxent_path  = "./output/maxent/SET_VAR_KOPPEN",
#                                    thresholds   = thresh.max.train,
#                                    percentiles  = percent.10.log,
#                                    time_slice   = 50,
#                                    area_occ     = 10),
#                             
#                             error = function(cond) {
#                               
#                               message(paste('Species skipped - check inputs', spp))
#                               
#                             })
# 
# 
# ## Combine GCM output for 2070 
# suitability.2070 = tryCatch(mapply(combine_gcm_threshold, 
#                                    DIR_list     = SDM.RESULTS.DIR, 
#                                    species_list = map_spp, 
#                                    maxent_path  = "./output/maxent/SET_VAR_KOPPEN",
#                                    thresholds   = thresh.max.train,
#                                    percentiles  = percent.10.log,
#                                    time_slice   = 70,
#                                    area_occ     = 10),
#                             
#                             error = function(cond) {
#                               
#                               message(paste('Species skipped - check inputs', spp))
#                               
#                             })


#########################################################################################################################
## Then use the more forgiving threshold
#########################################################################################################################


## To test it first species works 
# species    = spp_lower_thresh[13]
# thresh     = percent.10.log.low[13]
# percent    = percent.10.om.low[13]
# time_slice = 30
# area_occ   = 10

## SDM.RESULTS.DIR.LOW length doesn' match...............................................................................


#########################################################################################################################
## Combine output and calculate gain and loss for 2030 
suitability.2030 = tryCatch(mapply(combine_gcm_threshold,
                                   DIR_list     = SDM.RESULTS.DIR.LOW,
                                   species_list = spp_lower_thresh,
                                   maxent_path  = "./output/maxent/SET_VAR_KOPPEN",
                                   thresholds   = percent.10.log.low,
                                   percentiles  = percent.10.om.low,
                                   time_slice   = 30,
                                   area_occ     = 10),
                            
                            error = function(cond) {
                              
                              message(paste('Species skipped - check inputs', spp))
                              
                            })


## Combine GCM output for 2050 
suitability.2050 = tryCatch(mapply(combine_gcm_threshold, 
                                   DIR_list     = SDM.RESULTS.DIR.LOW,
                                   species_list = spp_lower_thresh,
                                   maxent_path  = "./output/maxent/SET_VAR_KOPPEN",
                                   thresholds   = percent.10.log.low,
                                   percentiles  = percent.10.om.low,
                                   time_slice   = 50,
                                   area_occ     = 10),
                            
                            error = function(cond) {
                              
                              message(paste('Species skipped - check inputs', spp))
                              
                            })


## Combine GCM output for 2070 
suitability.2070 = tryCatch(mapply(combine_gcm_threshold, 
                                   DIR_list     = SDM.RESULTS.DIR.LOW,
                                   species_list = spp_lower_thresh,
                                   maxent_path  = "./output/maxent/SET_VAR_KOPPEN",
                                   thresholds   = percent.10.log.low,
                                   percentiles  = percent.10.om.low,
                                   time_slice   = 70,
                                   area_occ     = 10),
                            
                            error = function(cond) {
                              
                              message(paste('Species skipped - check inputs', spp))
                              
                            })




#########################################################################################################################
## OUTSTANDING PREDICTION TASKS:
#########################################################################################################################


## 3). Try to thin records for ~100 spp with boundary bias: random sampling of those species records, by state and environment

## 4). Use more forgiving thresholds (10%) for all species, OR just those with bad maps:
##    "Maximum.training.sensitivity.plus.specificity.Logistic.threshold"
##    "X10.percentile.training.presence.Logistic.threshold"
##    "X10.percentile.training.presence.training.omission"

## 5). Calculate the TSS for all species (and MESS maps for a few species)

## 7). Decide on gamma diversity for this article - 200 spp - could it go to GEB? Plan figures and tables for MS



#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################