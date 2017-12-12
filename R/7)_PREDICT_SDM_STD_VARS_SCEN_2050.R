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
source('./R/HIA_LIST_MATCHING.R')
load("./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT.RData")


## Require packages
p <- c('ff',    'things',         'raster',    'dismo',        'sp',           'latticeExtra', 'data.table', 
       'rgdal', 'rgeos',          'gdalUtils', 'rmaxent',      'readr',        'dplyr',        'tidyr',
       'readr', 'rnaturalearth',  'rasterVis', 'RColorBrewer', 'latticeExtra', 'parallel')
sapply(p, require, character.only = TRUE)


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
h <- read_html('http://www.worldclim.org/cmip5_30s') # Also https://www.climatechangeinaustralia.gov.au/en/support-and-guidance/faqs/eight-climate-models-data/
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
species_list  <- basename(list.dirs('F:/green_cities_sdm/output/maxent/STD_VAR_ALL',   recursive = FALSE))
species_rev   = sort(species_list, decreasing = TRUE)


spp.all  <- unique(COMBO.NICHE.CONTEXT$searchTaxon)
str(spp.all)                 ## 6782


## The trial species
test.spp = sort(unique(c(renee.full$Species, "Betula pendula", "Fraxinus excelsior", "Quercus robur", "Fagus sylvatica")))
test.spp 


## Combine with 45 from the main list
HIA.SAMPLE = head(COMBO.NICHE.CONTEXT, 53)[, c("searchTaxon")]
test.spp   = sort(unique(c(test.spp, HIA.SAMPLE)))


## Now make the test species directory names
test_spp = gsub(" ", "_", test.spp)
test_spp = intersect(test_spp, species_list)
test_rev = sort(test_spp, decreasing = TRUE)




#########################################################################################################################
## 3). PROJECT MAXNET MODELS FOR 2050
#########################################################################################################################


#########################################################################################################################
## First, run a loop over each scenario: options(error = recover)
env.grids.2050 = lapply(scen_2050, function(x) {
  
  ## Assign the scenario name (to use later in the plot)
  scen_name = gcms.50$GCM[gcms.50$id == x]
  
  ## Create a raster stack for each 2050 GCM - also an empty raster for the final plot
  s <- stack(
    sprintf('./data/base/worldclim/aus/0.5/bio/2050/%s/%s%s.tif',
            x, x, 1:19))
  #empty <- init(s[[1]], function(x) NA)
  
  # nm <- sub('.*bi\\d0(.*)', '\\1', names(s))
  # names(s) <- sprintf('bio%02d', as.numeric(nm))
  
  ## Rename both the current and future environmental stack...
  names(s) <- names(env.grids.current) <- c(
    'Annual_mean_temp',    'Mean_diurnal_range',
    'Isothermality',       'Temp_seasonality',  'Max_temp_warm_month',
    'Min_temp_cold_month', 'Temp_annual_range', 'Mean_temp_wet_qu',
    'Mean_temp_dry_qu',    'Mean_temp_warm_qu', 'Mean_temp_cold_qu',  'Annual_precip',
    'Precip_wet_month',    'Precip_dry_month',  'Precip_seasonality', 'Precip_wet_qu',
    'Precip_dry_qu',       'Precip_warm_qu',    'Precip_col_qu')
  
  ########################################################################################################################
  ## Divide the temperature rasters by 10: 11 million NA values?
  ## s[[1:11]] <- s[[1:11]]/10 ## that code doesn't work, this is a work-around...
  s[[1]]  = s[[1]]/10
  s[[2]]  = s[[2]]/10
  s[[3]]  = s[[3]]/10
  s[[4]]  = s[[4]]/10
  s[[5]]  = s[[5]]/10
  s[[6]]  = s[[6]]/10
  s[[7]]  = s[[7]]/10
  s[[8]]  = s[[8]]/10
  s[[9]]  = s[[9]]/10
  s[[10]] = s[[10]]/10
  s[[11]] = s[[11]]/10
  
  ########################################################################################################################
  ## Now loop over the species...   
  lapply(species_rev, function(species) {
    
    ## First check if the species projection has already been run...
    if(!file.exists(sprintf('F:/green_cities_sdm/output/maxent/STD_VAR_ALL/%s/full/%s_%s.tif',
                            species, species, x))) {
      message('Doing ', species) 
      
      ## Read in the SDM model calibrated on current conditions
      m <- readRDS(sprintf('F:/green_cities_sdm/output/maxent/STD_VAR_ALL/%s/maxent_fitted.rds', species)) 
      
      ## Read in the occurrence points used to create the SDM
      occ <- readRDS(
        sprintf('F:/green_cities_sdm/output/maxent/STD_VAR_ALL/%s/occ.rds', 
                species)) %>%
        
      ########################################################################################################################
      ## If the current raster doesn't exist, create it
      spTransform(CRS('+init=epsg:4326'))
      f_current <- sprintf('F:/green_cities_sdm/output/maxent/STD_VAR_ALL/%s/full/%s_current.tif', 
                           species, species)
      
      if(!file.exists(f_current)) {
        
        ## Report which prediction is in progress
        message('Running current prediction for', species) 
        
        pred.current <- rmaxent::project(
          m$me_full, env.grids.current[[colnames(m$me_full@presence)]])$prediction_logistic
        writeRaster(pred.current, f_current, overwrite = TRUE)
        
      } else {
        
        pred.current = raster(sprintf('F:/green_cities_sdm/output/maxent/STD_VAR_ALL/%s/full/%s_current.tif', 
                                      species, species))
      }
      
      ########################################################################################################################
      ## If the future raster doesn't exist, create it 
      f_future <- sprintf('F:/green_cities_sdm/output/maxent/STD_VAR_ALL/%s/full/%s_%s.tif', 
                          species, species, x)
      
      if(!file.exists(f_future)) {
        
        ## Report which prediction is in progress
        message('Running future prediction for', species, ' ', x) 
        
        pred.future <- rmaxent::project(
          m$me_full, s[[colnames(m$me_full@presence)]])$prediction_logistic
        writeRaster(pred.future, f_future, overwrite = TRUE)
        
        ## Now create the empty panel just before plotting
        empty <- init(pred.future, function(x) NA)
        
        ########################################################################################################################
        ## Use the levelplot function to make a multipanel output: occurrence points, current raster and future raster
        png(sprintf('F:/green_cities_sdm/output/maxent/STD_VAR_ALL/%s/full/%s_%s.png', species, species, x),      
            11, 4, units = 'in', res = 300)
        
        ## Need an empty frame
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
                        colorkey   = list(height = 0.5, width = 3), xlab = '', ylab = '',
                        main       = list(gsub('_', ' ', species), font = 4, cex = 2)) +
                
                ## Plot the Aus shapefile with the occurrence points for reference
                ## Can the points be made more legible for both poorly and well recorded species?
                layer(sp.polygons(aus)) +
                layer(sp.points(occ, pch = 20, cex = 0.4, 
                                col = c('red', 'transparent', 'transparent')[panel.number()]), data = list(occ = occ)))
        dev.off()
        
      }
      
    } else {
      
      message(species, ' ', x, ' skipped - prediction already run')   ## Ignore species which have already been run - maybe remove
      
    }
    
  })
  
})





#########################################################################################################################
## 4). PROJECT MODELS FOR 2070
#########################################################################################################################


## First, check the extents of the 2050 and 2070 grids match
ras.50  = raster("F:/green_cities_sdm/data/base/worldclim/aus/0.5/bio/2050/ac85bi50/ac85bi501.tif")
ras.70  = raster("F:/green_cities_sdm/data/base/worldclim/aus/0.5/bio/2070/ac85bi70/ac85bi701.tif")

## The extents are now identical, and the number of cells matches
identical(extent(ras.50), extent(ras.70))
ras.50;ras.70



#########################################################################################################################
## Now run a loop over each 2070 scenario
env.grids.2070 = lapply(scen_2070, function(x) {
  
  ## Assign the scenario name (to use later in the plot)
  scen_name = gcms.70$GCM[gcms.70$id == x]
  
  ## Create a raster stack for each 2050 GCM 
  s <- stack(
    sprintf('./data/base/worldclim/aus/0.5/bio/2070/%s/%s%s.tif',
            x, x, 1:19))
  
  ## Also create an empty raster for the final plot
  empty <- init(s[[1]], function(x) NA)
  
  # nm <- sub('.*bi\\d0(.*)', '\\1', names(s))
  # names(s) <- sprintf('bio%02d', as.numeric(nm))
  
  ## Rename the future environmental stack...
  names(s) <- names(env.grids.current) <- c(
    'Annual_mean_temp',    'Mean_diurnal_range',
    'Isothermality',       'Temp_seasonality',  'Max_temp_warm_month',
    'Min_temp_cold_month', 'Temp_annual_range', 'Mean_temp_wet_qu',
    'Mean_temp_dry_qu',    'Mean_temp_warm_qu', 'Mean_temp_cold_qu',  'Annual_precip',
    'Precip_wet_month',    'Precip_dry_month',  'Precip_seasonality', 'Precip_wet_qu',
    'Precip_dry_qu',       'Precip_warm_qu',    'Precip_col_qu')
  
  ## Divide the temperature rasters by 10: 11 million NA values?
  ## s[[1:11]] <- s[[1:11]]/10 ## that code doesn't work, this is a work-around...
  s[[1]]  = s[[1]]/10
  s[[2]]  = s[[2]]/10
  s[[3]]  = s[[3]]/10
  s[[4]]  = s[[4]]/10
  s[[5]]  = s[[5]]/10
  s[[6]]  = s[[6]]/10
  s[[7]]  = s[[7]]/10
  s[[8]]  = s[[8]]/10
  s[[9]]  = s[[9]]/10
  s[[10]] = s[[10]]/10
  s[[11]] = s[[11]]/10
  
  ########################################################################################################################
  ## Now loop over the species...   
  lapply(test_rev, function(species) {
    
    ## First check if the species projection has already been run...
    if(!file.exists(sprintf('F:/green_cities_sdm/output/maxent/STD_VAR_ALL/%s/full/%s_%s.tif',
                            species, species, x))) {
      message('Doing ', species) 
      
      ## Read in the SDM model calibrated on current conditions
      m <- readRDS(sprintf('F:/green_cities_sdm/output/maxent/STD_VAR_ALL/%s/maxent_fitted.rds', species)) 
      
      ## Read in the occurrence points used to create the SDM
      occ <- readRDS(
        sprintf('F:/green_cities_sdm/output/maxent/STD_VAR_ALL/%s/occ.rds', 
                species)) %>%
        
      ########################################################################################################################
      ## If the current raster doesn't exist, create it
      spTransform(CRS('+init=epsg:4326'))
      f_current <- sprintf('F:/green_cities_sdm/output/maxent/STD_VAR_ALL/%s/full/%s_current.tif', 
                           species, species)
      
      if(!file.exists(f_current)) {
        
        ## Report which prediction is in progress
        message('Running future prediction for', species, ' ', x) 
        
        pred.current <- rmaxent::project(
          m$me_full, env.grids.current[[colnames(m$me_full@presence)]])$prediction_logistic
        writeRaster(pred.current, f_current, overwrite = TRUE)
        
      } else {
        
        ## Otherwise just read it in
        pred.current = raster(sprintf('F:/green_cities_sdm/output/maxent/STD_VAR_ALL/%s/full/%s_current.tif', 
                                      species, species))
      }
      
      ########################################################################################################################
      ## If the future raster doesn't exist, create it 
      f_future <- sprintf('F:/green_cities_sdm/output/maxent/STD_VAR_ALL/%s/full/%s_%s.tif', 
                          species, species, x)
      
      if(!file.exists(f_future)) {
        
        ## Report which prediction is in progress
        message('Doing future prediction for', species, ' ', x) 
        
        pred.future <- rmaxent::project(
          m$me_full, s[[colnames(m$me_full@presence)]])$prediction_logistic
        writeRaster(pred.future, f_future, overwrite = TRUE)
        
        ## Now create the empty panel just before plotting...this may have been causing problems!
        empty <- init(pred.future, function(x) NA)
        
        ########################################################################################################################
        ## Use the levelplot function to make a multipanel output: occurrence points, current raster and future raster
        png(sprintf('F:/green_cities_sdm/output/maxent/STD_VAR_ALL/%s/full/%s_%s.png', species, species, x),      
            11, 4, units = 'in', res = 300)
        
        ## Need an empty frame
        print(levelplot(stack(empty,
                              pred.current,
                              pred.future), margin = FALSE,
                        
                        ## Create a colour scheme using colbrewer: 100 is to make it continuos
                        ## Also, make it a one-directional colour scheme
                        scales      = list(draw = FALSE), 
                        at = seq(0, 1, length = 100),
                        col.regions = colorRampPalette(rev(brewer.pal(11, 'Spectral'))),
                        
                        ## Give each plot a name
                        names.attr = c('Occurrence', 'Current', sprintf('%s, 2070, RCP8.5', scen_name)),
                        colorkey   = list(height = 0.5, width = 3), xlab = '', ylab = '',
                        main       = list(gsub('_', ' ', species), font = 4, cex = 2)) +
                
                ## Plot the Aus shapefile with the occurrence points for reference
                ## Can the points be made more legible for both poorly and well recorded species?
                layer(sp.polygons(aus)) +
                layer(sp.points(occ, pch = 20, cex = 0.4, 
                                col = c('red', 'transparent', 'transparent')[panel.number()]), data = list(occ = occ)))
        dev.off()
        
      }
      
    } else {
      
      message(species, ' ', x, ' skipped - prediction already run')   ## Ignore species which have already been run - maybe remove
      
    }
    
  })
  
})





#########################################################################################################################
## OUTSTANDING PREDICTION TASKS:
#########################################################################################################################


## 2070 raster crop is failing: need to set the extent and resolution to be exactly the same...


## We are missing two scenarios recommended for Australia: CanESM2 & CESM1-CAM5. These two are not on the worldclim list:

## https://www.climatechangeinaustralia.gov.au/en/support-and-guidance/faqs/eight-climate-models-data/

## http://www.worldclim.org/cmip5_30s
  


## Consider the final format needed: which files: table, plots/maps, files. Space important for both local and web 


## How can we take an consensus layer of all the scenarios, to create a confidence interval? See ensemble.raster





#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################