#########################################################################################################################
######################################### PREDICT MAXENT TO FUTURE CLIMATES ############################################# 
#########################################################################################################################


## This code takes the output of the SDM models and generates a prediction of habitat suitability for current and future
## conditions. 


#########################################################################################################################
## 1). READ IN THE CURRENT AND FUTURE ENVIRONMENTAL DATA, AND CREATE LISTS FOR MODEL RUNS 
#########################################################################################################################


#########################################################################################################################
## create scenario list first



## Create raster stacks: 
env.grids.current = stack(
  file.path('//sci-7910/F/data/worldclim/aus/0.5/bio/current',
            sprintf('bio_%02d.tif', 1:19)))

## Future: the problem is occurring in here...
env.grids.future = stack(
  sprintf('//sci-7910/F/data/worldclim/aus/0.5/bio/2050/%s/%s%s.tif',
          scen, scen, 1:19))


env.grids.future  = c("//SCI-7910/f/data/worldclim/world/0.5/bio/2050/ac85bi50/ac85bi501.tif",
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
                      "//SCI-7910/f/data/worldclim/world/0.5/bio/2050/ac85bi50/ac85bi5019.tif")


## Convert all the rasters to a stack
env.grids.future <- stack(env.grids.future)


## str(env.grids.current)
## str(env.grids.future)


#########################################################################################################################
## Divide the temperature values by 10, because Worldclim layers are multiplied by 10 to reduced file size.


## R is writing a version of these files to memory for some reason...in this directory:
## C:\Users\user\AppData\Local\Temp\Rtmpiwken9\raster
## tempdir() ## set.tempdir("F:/RTEMP") does not work

## Create a file called .Renviron in the directory given by Sys.getenv('R_USER') and save it with the line TMP = '<your-desired-tempdir>'
## write("TMP = '<your-desired-tempdir>'", file=file.path(Sys.getenv('R_USER'), '.Renviron'))


## Is there a faster way to divide the rasters?
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


## Check the rasters
str(env.grids.current)
str(env.grids.future)


## Read in a shapefile for Australia, but exclude the Islands
aus <- ne_states(country = 'Australia') %>% 
  subset(!grepl('Island', name))


## Create a list of all the species folders which contain the fitted models: these were run in the previous step
species_list  <- basename(list.dirs('F:/green_cities_sdm/output/maxent/baseline',   recursive = FALSE))
scenario_list <- basename(list.dirs('//sci-7910/F/data/worldclim/aus/0.5/bio/2050', recursive = FALSE))


#########################################################################################################################
## Now save/load .RData file for the next session
save.image("STEP_7_PREDICT_SDM.RData")
load("STEP_7_PREDICT_SDM.RData")





#########################################################################################################################
## 2). CREATE MAPS OF CURRENT AND FUTURE HABITAT SUITABILITY
#########################################################################################################################


# Error in m[, i] <- getValues(x@layers[[i]]) : 
#   number of items to replace is not a multiple of replacement length
## debugonce(project)
load("STEP_7_PREDICT_SDM.RData")


##
env.grids.current[[colnames(m$me_full@presence)]]
env.grids.future[[colnames(m$me_full@presence)]]


## Ok so why does R need to save a version of each raster to a temporary folder?
# Error in file(fn, "rb") : cannot open the connection
# In addition: Warning message:
#   In file(fn, "rb") :
#   cannot open file 'C:\Users\user\AppData\Local\Temp\Rtmpiwken9\raster\r_tmp_2017-10-26_143723_11284_83052.gri': No such file or directory


#########################################################################################################################
## Use lappy to loop over a list of species
## Test on one species:
species = species_list[320] # [1] "Lomandra_longifolia"


##
lapply(species_list, function(species) {
  message('Doing ', species)
  
  lapply(scenario_list, function(scen) {
    message('  Doing ', scen)
    
    ## Read in the fitted models using sprintf
    m <- readRDS(sprintf('F:/green_cities_sdm/output/maxent/baseline/%s/maxent_fitted.rds', species))
    
    # ## These numbers don't look right 
    # env.grids.current[[colnames(m$me_full@presence)]]
    # env.grids.future[[colnames(m$me_full@presence)]]
    
    ## Read in the occurrence files using sprintf
    occ <- readRDS(sprintf('F:/green_cities_sdm/output/maxent/baseline/%s/occ.rds', species)) %>% 
      spTransform(CRS('+init=epsg:4326'))
    
    ## Create rasters for the current and future climate
    ## Calculating contribution of feature 11 of 11..............why 11 and not 19?
    ## Or are features the maxent setting
    ## Also this takes a lot of time...why is the future prediction so slow?
    pred.current <- rmaxent::project(m$me_full, env.grids.current[[colnames(m$me_full@presence)]])
    pred.future  <- rmaxent::project(m$me_full, env.grids.future[[colnames(m$me_full@presence)]])
    
    ## What is going wrong here?
    # Error in m[, i] <- getValues(x@layers[[i]]) : 
    #   number of items to replace is not a multiple of replacement length
    ## debugonce(project)
    
    
    ## Write the current raster out
    ## printf has three function arguments here. The folders will need to change
    writeRaster(pred.current, sprintf('F:/green_cities_sdm/output/maxent/baseline/%s/full/%s_current.tif', 
                                      species, species))
    
    ## Write the future raster out
    writeRaster(pred.future, sprintf('F:/green_cities_sdm/output/maxent/baseline/%s/full/%s_%s.tif', 
                                     species, species, scen))
    
    ## Create an empty raster based on the future prediction
    empty <- init(pred.future$prediction_logistic, function(x) NA)
    
    
    # Warning message:
    #   In .rasterFromRasterFile(grdfile, band = band, objecttype, ...) :
    #   size of values file does not match the number of cells (given the data type)
    
    
    #########################################################################################################################
    ## Create map of habitat suitability...the first line starts the PNG device
    
    # Error in compareRaster(x) : different extent
    
    png(sprintf('F:/green_cities_sdm/output/maxent/baseline/%s/full/%s.png', species, species), 
        11, 4, units = 'in', res = 300)
    
    ## Use the levelplot function to make a multipanel output
    levelplot(stack(empty,
                    pred.current$prediction_logistic,
                    pred.future$prediction_logistic), margin = FALSE, 
              
              ## Create a colour scheme using colbrewer
              scales      = list(draw = FALSE), at = seq(0, 1, length = 100),
              col.regions = colorRampPalette(rev(brewer.pal(11, 'Spectral'))),
              
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
#save.image('STEP_7_SDM.RData')
  
  
## What is the output?
  
  
  
  
#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################
  