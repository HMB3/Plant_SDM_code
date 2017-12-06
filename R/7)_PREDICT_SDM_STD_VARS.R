#########################################################################################################################
######################################### PREDICT MAXENT TO FUTURE CLIMATES ############################################# 
#########################################################################################################################


## This code takes the output of the SDM models and generates a prediction of habitat suitability for current and future
## environmental conditions. The data is the format of all species occurrences (rows) and environmental variables (columns)


#########################################################################################################################
## 1). READ IN THE CURRENT AND FUTURE ENVIRONMENTAL DATA, AND CREATE LISTS FOR MODEL RUNS 
#########################################################################################################################


#########################################################################################################################
## Load packages, functions and data
source('./R/HIA_LIST_MATCHING.R')
load("./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT.RData")


## Require packages
p <- c('ff',    'things',         'raster',    'dismo',        'sp',           'latticeExtra', 'data.table', 
       'rgdal', 'rgeos',          'gdalUtils', 'rmaxent',      'readr',        'dplyr',        'tidyr',
       'readr', 'rnaturalearth',  'rasterVis', 'RColorBrewer', 'latticeExtra', 'parallel')
sapply(p, require, character.only = TRUE)


#########################################################################################################################
## create a list of GCM scenarios: 
## Eight of the 40 CMIP5 models assessed in this project have been selected for use in provision of application-ready data. 
## This facilitates efficient exploration of climate projections for Australia.
## https://www.climatechangeinaustralia.gov.au/en/support-and-guidance/faqs/eight-climate-models-data/
# scen = c("ip85bi50", "mc85bi50", "mg85bi50", "mi85bi50", "mp85bi50", 
#          "mr85bi50", "no85bi50", "ac85bi50", "bc85bi50", "cc85bi50", 
#          "cn85bi50", "gf85bi50", "gs85bi50", "hd85bi50", "he85bi50",
#          "hg85bi50", "in85bi50")


## Just get the 8 models picked by CSIRO for Australia
scen = c("mc85bi50", "no85bi50", "ac85bi50", "cn85bi50", "gf85bi50", "hg85bi50")


## Create a lookup table of GCMs
h <- read_html('http://www.worldclim.org/cmip5_30s') # Also https://www.climatechangeinaustralia.gov.au/en/support-and-guidance/faqs/eight-climate-models-data/
gcms <- h %>% 
  html_node('table') %>% 
  html_table(header=TRUE) %>% 
  filter(rcp85 != '')

id <- h %>% 
  html_nodes(xpath = "//a[text()='bi']/@href") %>% 
  as.character %>% 
  grep('85bi50', ., value=TRUE) %>% 
  basename %>% 
  sub('\\.zip.', '', .)

gcms <- cbind(gcms, id)
gcms$GCM = sub(" \\(#\\)", "", gcms$GCM)  ## sub replaces first instance in a string, gsub = global
gcms


## Create a stack of all the scenarios, then take the AV and SD of the stack. This could provide the confidence interval?
## Or, we could use the AV for the main maps instead?


#########################################################################################################################
## Create raster stacks using files in John's directories:
## sprintf has two arguments here: the main path, then the places that the bioclim number is inserted to complete the path 
env.grids.current = stack(
  file.path('./data/base/worldclim/aus/0.5/bio/current',
            sprintf('bio_%02d.tif', 1:19)))


# Future: this list is causing problems
env.grids.future = stack(
  sprintf('./data/base/worldclim/aus/0.5/bio/2050/%s/%s%s.tif',
          scen[1], scen[1], 1:19))

# env.grids.future = lapply(scen, function(x) {
#   
#   stack(
#     sprintf('./data/base/worldclim/aus/0.5/bio/2050/%s/%s%s.tif',
#             x, x, 1:19))
#   
# })


## Now rename the list of current rasters and future rasters: names look ok
class(env.grids.current);class(env.grids.future[[1]])    #[[1:19]])
names(env.grids.current);names(env.grids.future[[1:19]]) #[[1:19]])


## save(GBIF.UNRESOLVED, file = paste("./data/base/HIA_LIST/GBIF/GBIF_UNRESOLVED.RData"))


#########################################################################################################################
## Divide the temperature values by 10, because Worldclim layers are multiplied by 10 to reduced file size.
## The faster way doesn't seem to be working, because the indexing is wrong.How can we fix this?


## Current ennvironmental conditons. 
## The alternative code fails: 
## env.grids.current[[1:11]] <- env.grids.current[[1:11]]/10
for(i in 1:11) {
  
  message(i)
  env.grids.current[[i]] <- env.grids.current[[i]]/10
  env.grids.future[[i]]  <- env.grids.future[[i]]/10  
  
}


## The alternative code fails:
# for(i in seq_along(env.grids.future)) {
#   
#   env.grids.future[[i]][[1:11]] <- env.grids.future[[i]][[1:11]]/10
#   
# }


#########################################################################################################################
## Also check the environmental data
summary(env.grids.current[[1]]);summary(env.grids.future[[1]]) #[[1:19]])


## Current
# png('F:/green_cities_sdm/output/maxent/STD_VAR_ALL/env_grids_current.png',              
#     15, 15, units = 'in', res = 600)
# 
# 
# ## Use the levelplot function to make a multipanel output
# plot(env.grids.current, pch = 20, cex = 0.5) 
# 
# 
# ## finish the PNG device
# dev.off()
# 
# 
# ## Future
# png('F:/green_cities_sdm/output/maxent/STD_VAR_ALL/env_grids_future.png',              
#     15, 15, units = 'in', res = 600)
# 
# 
# ## Use the levelplot function to make a multipanel output
# plot(env.grids.future, pch = 20, cex = 0.5) 
# 
# ## finish the PNG device
# dev.off()


## Give the current and future environmental grids the same names. 
## But we can't use the same command for a raster stack vs. the list?
names(env.grids.current) <- names(env.grids.future) <- c(
  'Annual_mean_temp',    'Mean_diurnal_range',
  'Isothermality',       'Temp_seasonality',  'Max_temp_warm_month',
  'Min_temp_cold_month', 'Temp_annual_range', 'Mean_temp_wet_qu',
  'Mean_temp_dry_qu',    'Mean_temp_warm_qu', 'Mean_temp_cold_qu',  'Annual_precip',
  'Precip_wet_month',    'Precip_dry_month',  'Precip_seasonality', 'Precip_wet_qu',
  'Precip_dry_qu',       'Precip_warm_qu',    'Precip_col_qu')


## Check the rasters: list of 19 with 11 slots
names(env.grids.current);names(env.grids.future)


## Read in a shapefile for Australia, but exclude the Islands
aus <- ne_states(country = 'Australia') %>% 
  subset(!grepl('Island', name))


## Create a list of all the species folders which contain the fitted models: these were run in the previous step
species_list  <- basename(list.dirs('F:/green_cities_sdm/output/maxent/STD_VAR_ALL',   recursive = FALSE))
#save.image("STEP_7_PREDICT.RData")





#########################################################################################################################
## 2). CREATE LISTS FOR MODEL RUNS 
#########################################################################################################################
spp.all  <- unique(COMBO.NICHE.CONTEXT$searchTaxon)
str(spp.all)                 ## 6782


## The trial species
test.spp = sort(unique(c(renee.full$Species, "Betula pendula", "Fraxinus excelsior", "Quercus robur", "Fagus sylvatica")))
test.spp 


## Combine with 45 from the main list
HIA.SAMPLE = head(COMBO.NICHE.CONTEXT, 53)[, c("searchTaxon")]
test.spp   = sort(unique(c(test.spp, HIA.SAMPLE)))


## Now restrict the species_list to the test.spp
test_spp = gsub(" ", "_", test.spp)
test_spp = intersect(test_spp, species_list)





#########################################################################################################################
## 3). CREATE MAPS OF CURRENT AND FUTURE HABITAT SUITABILITY FOR ALL RECORDS AND SELECTED VARIABLES
#########################################################################################################################


#########################################################################################################################
## Use lappy to loop over a list of species
## Test on one species and scenario:
#load("STEP_7_PREDICT.RData")
#species = species_list[55]
scen_i = scen[1]


#########################################################################################################################
## Now run the code over a list of species (not enough memory for extra directories)
lapply(test_spp, function(species) {

  ## First check if the species projection has already been run...
  if(!file.exists(sprintf('F:/green_cities_sdm/output/maxent/STD_VAR_ALL/%s/full/%s_%s.tif',
                                                   species, species, scen_i))) {
    message('Doing ', species) 
  
  ## Get the scenario looping code working eventually  
  # lapply(scen, function(scen_i) {
  #   message('  Doing ', scen_i)
  
  ## Assign the scenario name to the final plot
  scen_name = gcms$GCM[gcms$id == scen_i] 
  
  ## Read in the fitted models using sprintf
  m <- readRDS(sprintf('F:/green_cities_sdm/output/maxent/STD_VAR_ALL/%s/maxent_fitted.rds', species))    
  
  # ## These numbers don't look right for precip wet month and Precip_seasonality?
  # str(m);names(m)
  # env.grids.current[[colnames(m$me_full@presence)]]
  # env.grids.future[[colnames(m$me_full@presence)]]
  
  ## Read in the occurrence files from the output directory using sprintf
  # taxa = gsub("_", " ", species)
  # occ  = subset(COMBO.RASTER.CONTEXT, searchTaxon == taxa)[, c("lon", "lat")]
  occ <- readRDS(sprintf('F:/green_cities_sdm/output/maxent/STD_VAR_ALL/%s/occ.rds', species)) %>%
    spTransform(CRS('+init=epsg:4326'))

  ## Create rasters for the current and future climate: 
  ## problems are to do with the indexing of raster vs a list of rasters...
  
  ## Also, why are there different numbers of features each time? Is this because of model selection? 
  ## So n predictors selected by the model * n features, etc?
  pred.current <- rmaxent::project(m$me_full, env.grids.current[[colnames(m$me_full@presence)]])$prediction_logistic
  pred.future  <- rmaxent::project(m$me_full, env.grids.future[[colnames(m$me_full@presence)]])$prediction_logistic
  
  ## The indexing below using scen_i is not working...
  #pred.future  <- rmaxent::project(m$me_full, env.grids.future[[scen_i]][[colnames(m$me_full@presence)]])$prediction_logistic
  
  ## Write the raster of suitability to current conditions out to file
  writeRaster(pred.current, sprintf('F:/green_cities_sdm/output/maxent/STD_VAR_ALL/%s/full/%s_current.tif',  ## add scenario to file name
                                    species, species), overwrite = TRUE)

  ## Write the raster of suitability to future conditions out to file
  writeRaster(pred.future, sprintf('F:/green_cities_sdm/output/maxent/STD_VAR_ALL/%s/full/%s_%s.tif',        ## add scenario to file name
                                   species, species, scen_i), overwrite = TRUE)
  
  ## Create an empty raster based on the future prediction
  empty <- init(pred.future, function(x) NA)
  
  #########################################################################################################################
  ## Create map of habitat suitability...the first line starts the PNG device
  ## This won't work in the loop, but it works when run individually
  png(sprintf('F:/green_cities_sdm/output/maxent/STD_VAR_ALL/%s/full/%s_%s.png', species, species, scen_i),      ## add scenario to file name
      11, 4, units = 'in', res = 300)
  
  ## Use the levelplot function to make a multipanel output
  ## Error in .Call.graphics(C_palette2, .Call(C_palette2, NULL)) : 
  ## invalid graphics state
  print(levelplot(stack(empty,
                  pred.current,
                  pred.future), margin = FALSE,
            
            ## Create a colour scheme using colbrewer: 100 is to make it continuos
            ## Also, make it a one-directional colour scheme
            scales      = list(draw = FALSE), 
            at = seq(0, 1, length = 100),
            col.regions = colorRampPalette(rev(brewer.pal(11, 'Spectral'))),
            
            ## Give each plot a name
            names.attr = c('Occurrence', 'Current', sprintf('%s, 2050, RCP8.5', scen_name)),
            colorkey   = list(height = 0.6, width = 3), xlab = '', ylab = '',
            main       = list(gsub('_', ' ', species), font = 4, cex = 2)) +
    
    ## Plot the Aus shapefile with the occurrence points for reference
    ## Why don't the points print out inside the loop?
    layer(sp.polygons(aus)) +
    layer(sp.points(occ, pch = 20, cex = 0.6, 
                    col = c('red', 'transparent', 'transparent')[panel.number()]), data = list(occ = occ)))
  
  # Why this Warning messages? 
  #   1: In min(x) : no non-missing arguments to min; returning Inf  ##
  #   2: In max(x) : no non-missing arguments to max; returning -Inf
  
  ## finish the PNG device
  dev.off()
  
  } else {
    
    message(species, ' skipped - prediction already run')         ## This condition ignores species which have no data
    
  } 
  
})

  
# })
  
  
 

  

  
  



#########################################################################################################################
## OUTSTANDING PREDICTION TASKS:
#########################################################################################################################


## Fix raster lists so all scenarios can be iterated over, indexing is causing problems


## Why do some species produce weird-looking maps. Is this because the values are curren/future environmental 
## values wrong, the scenarios are weird or there aren't enough records...or all 3?


## How can we take an consensus layer of all the scenarios to create a confidence interval?


  
  
#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################