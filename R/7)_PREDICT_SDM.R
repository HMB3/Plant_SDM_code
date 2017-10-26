#########################################################################################################################
######################################### PREDICT MAXENT TO FUTURE CLIMATES ############################################# 
#########################################################################################################################


## This code takes a table of all species occurrences (rows) and environmental values "columns", and runs a maxent model 
## for each species. 


#########################################################################################################################
## 1). READ IN THE CURRENT AND FUTURE ENVIRONMENTAL DATA, AND CREATE LISTS FOR MODEL RUNS 
#########################################################################################################################


#########################################################################################################################
## Create raster stacks
env.grids.current = stack(
  file.path('//sci-7910/F/data/worldclim/aus/0.5/bio/current',
            sprintf('bio_%02d.tif', 1:19)))

env.grids.future = stack(
  sprintf('//sci-7910/F/data/worldclim/aus/0.5/bio/2050/%s/%s%s.tif',
          scen, scen, 1:19))


env.grids.future  = stack (c("//SCI-7910/f/data/worldclim/world/0.5/bio/2050/ac85bi50/ac85bi501.tif",
                             "//SCI-7910/f/data/worldclim/world/0.5/bio/2050/ac85bi50/ac85bi502.tif",
                             "//SCI-7910/f/data/worldclim/world/0.5/bio/2050/ac85bi50/ac85bi503.tif",
                             "//SCI-7910/f/data/worldclim/world/0.5/bio/2050/ac85bi50/ac85bi504.tif",
                             "//SCI-7910/f/data/worldclim/world/0.5/bio/2050/ac85bi50/ac85bi505.tif",
                             "//SCI-7910/f/data/worldclim/world/0.5/bio/2050/ac85bi50/ac85bi506.tif",
                             "//SCI-7910/f/data/worldclim/world/0.5/bio/2050/ac85bi50/ac85bi507.tif",
                             "//SCI-7910/f/data/worldclim/world/0.5/bio/2050/ac85bi50/ac85bi508.tif",
                             "//SCI-7910/f/data/worldclim/world/0.5/bio/2050/ac85bi50/ac85bi509.tif",
                             "//SCI-7910/f/data/worldclim/world/0.5/bio/2050/ac85bi50/ac85bi5010.tif",
                             "//SCI-7910/f/data/worldclim/world/0.5/bio/2050/ac85bi50/ac85bi5011.tif",
                             "//SCI-7910/f/data/worldclim/world/0.5/bio/2050/ac85bi50/ac85bi5012.tif",
                             "//SCI-7910/f/data/worldclim/world/0.5/bio/2050/ac85bi50/ac85bi5013.tif",
                             "//SCI-7910/f/data/worldclim/world/0.5/bio/2050/ac85bi50/ac85bi5014.tif",
                             "//SCI-7910/f/data/worldclim/world/0.5/bio/2050/ac85bi50/ac85bi5015.tif",
                             "//SCI-7910/f/data/worldclim/world/0.5/bio/2050/ac85bi50/ac85bi5016.tif",
                             "//SCI-7910/f/data/worldclim/world/0.5/bio/2050/ac85bi50/ac85bi5017.tif",
                             "//SCI-7910/f/data/worldclim/world/0.5/bio/2050/ac85bi50/ac85bi5018.tif",
                             "//SCI-7910/f/data/worldclim/world/0.5/bio/2050/ac85bi50/ac85bi5019.tif"))


#########################################################################################################################
## Divide the temperature values by 10, because Worldclim are multiplied by 10 to reduced file size
for(i in 1:11) {
  
  message(i)
  env.grids.current[[i]] <- env.grids.current[[i]]/10
  env.grids.future[[i]]  <- env.grids.future[[i]]/10  
  
}


## give the current and future environmental grids the same names
names(env.grids.current) <- names(env.grids.future) <- c(
  'Annual_mean_temp',    'Mean_diurnal_range',
  'Isothermality',       'Temp_seasonality',  'Max_temp_warm_month',
  'Min_temp_cold_month', 'Temp_annual_range', 'Mean_temp_wet_qu',
  'Mean_temp_dry_qu',    'Mean_temp_warm_qu', 'Mean_temp_cold_qu',  'Annual_precip',
  'Precip_wet_month',    'Precip_dry_month',  'Precip_seasonality', 'Precip_wet_qu',
  'Precip_dry_qu',       'Precip_warm_qu',    'Precip_col_qu')


## Read in a shapefile for Australia, but exclude the Islands
aus <- ne_states(country = 'Australia') %>% 
  subset(!grepl('Island', name))


## Create a list of all the species folders which contain the fitted models
species_list  <- basename(list.dirs('F:/green_cities_sdm/output/maxent/baseline',   recursive = FALSE))
scenario_list <- basename(list.dirs('//sci-7910/F/data/worldclim/aus/0.5/bio/2050', recursive = FALSE))


## Now save/load .RData file for the next session
save.image("STEP_7_PREDICT_SDM.RData")
load("STEP_7_PREDICT_SDM.RData")





#########################################################################################################################
## 2). CREATE MAPS OF CURRENT AND FUTURE HABITAT SUITABILITY
#########################################################################################################################


# Error in m[, i] <- getValues(x@layers[[i]]) : 
#   number of items to replace is not a multiple of replacement length

#########################################################################################################################
## Use lappy to loop over a list of species
## Test on one species:
species = species_list[320] # [1] "Lomandra_longifolia"


lapply(species_list, function(species)) {
  message('Doing ', species)
  
  lapply(scenario_list, function(scen) {
    message('  Doing ', scen)
    
    ## Read in the fitted models using sprintf
    m <- readRDS(sprintf('F:/green_cities_sdm/output/maxent/baseline/%s/maxent_fitted.rds', species))
    
    ## Read in the occurrence files using sprintf
    occ <- readRDS(sprintf('F:/green_cities_sdm/output/maxent/baseline/%s/occ.rds', species)) %>% 
      spTransform(CRS('+init=epsg:4326'))
    
    ## Create rasters for the current and future climate
    pred.current <- rmaxent::project(m$me_full, env.grids.current[[colnames(m$me_full@presence)]])
    pred.future  <- rmaxent::project(m$me_full, env.grids.future[[colnames(m$me_full@presence)]])
    
    
    ## Make a list of current: Sprintf has three function arguments here...
    ## The folders will need to change
    writeRaster(pred.current, sprintf('F:/green_cities_sdm/output/maxent/baseline/%s/full/%s_current.tif', 
                                      species, species))
    
    ## Make a list of future
    writeRaster(pred.future, sprintf('F:/green_cities_sdm/output/maxent/baseline/%s/full/%s_%s.tif', 
                                     species, species, scen))
    
    ## Create an empty raster
    empty <- init(pred.future$prediction_logistic, function(x) NA)
    
    #########################################################################################################################
    ## Create map of habitat suitability...the first line starts the PNG device
    png(sprintf('F:/green_cities_sdm/output/maxent/baseline/%s/full/%s.png', species, species), 
        11, 4, units = 'in', res = 300)
    
    ## Use the levelplot function to make a multipanel output
    levelplot(stack(empty,
                    pred.current$prediction_logistic,
                    pred.future$prediction_logistic), margin = FALSE, 
              
              ## Create a colour scheme using colbrewer
              scales      = list(draw = FALSE), at = seq(0, 1, length = 100),
              col.regions =colorRampPalette(rev(brewer.pal(11, 'Spectral'))),
              
              ## Give each plot a name
              names.attr = c('Occurrence', 'Current', 'CSIRO Access 1.0, 2050, RCP8.5'),
              colorkey   = list(height = 0.6, width = 3), xlab = '', ylab = '',
              main       = list(gsub('_', ' ', species), font = 4, cex = 2)) +
      
      ## Plot the Aus shapefile with the occurrence points for reference
      layer(sp.polygons(aus)) +
      layer(sp.points(occ, pch = 20, cex = 0.5, 
                      col = c('red', 'transparent', 'transparent')[panel.number()]))
    
    ## finish the PNG device
    dev.off()
    
  }
  
)}



## Now save .RData file for the next session
save.image("STEP_6_SDM.RData")






#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################
