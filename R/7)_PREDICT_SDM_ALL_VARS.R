#########################################################################################################################
######################################### PREDICT MAXENT TO FUTURE CLIMATES ############################################# 
#########################################################################################################################


## This code takes the output of the SDM models and generates a prediction of habitat suitability for current and future
## environmental conditions. 


#########################################################################################################################
## 1). READ IN THE CURRENT AND FUTURE ENVIRONMENTAL DATA, AND CREATE LISTS FOR MODEL RUNS 
#########################################################################################################################


#########################################################################################################################
## create scenario list first: or just one at a time
#load("STEP_6_SDM.RData")
scen = c("ip85bi50", "mc85bi50", "mg85bi50", "mi85bi50", "mp85bi50", 
         "mr85bi50", "no85bi50", "ac85bi50", "bc85bi50", "cc85bi50", 
         "cn85bi50", "gf85bi50", "gs85bi50", "hd85bi50", "he85bi50",
         "hg85bi50", "in85bi50")

## create lookup table
library(rvest)
h <- read_html('http://www.worldclim.org/cmip5_30s')
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


#########################################################################################################################
## Create raster stacks using files in John's directories:
## sprintf has two arguments here: the main path, then the places that the bioclim number is inserted to complete the path 
env.grids.current = stack(
  file.path('//SCI-7910/F/data/worldclim/aus/0.5/bio/current',
            sprintf('bio_%02d.tif', 1:19)))


## Future: just running code for one scenario at a time
env.grids.future = lapply(scen, function(x) {
  stack(
    sprintf('//SCI-7910/F/data/worldclim/aus/0.5/bio/2050/%s/%s%s.tif',
            x, x, 1:19))
})
names(env.grids.future) <- scen


##
str(env.grids.current)
str(env.grids.future)


#########################################################################################################################
## Divide the temperature values by 10, because Worldclim layers are multiplied by 10 to reduced file size.


## R is writing a version of these files to memory for some reason...in this directory:
## C:\Users\user\AppData\Local\Temp\Rtmpiwken9\raster
## tempdir() ## set.tempdir("F:/RTEMP") does not work

## Create a file called .Renviron in the directory given by Sys.getenv('R_USER') and save it with the line TMP = '<your-desired-tempdir>'
## write("TMP = '<your-desired-tempdir>'", file=file.path(Sys.getenv('R_USER'), '.Renviron'))


## Is there a faster way to divide the rasters?
env.grids.current[[1:11]] <- env.grids.current[[1:11]]/10
for(i in seq_along(env.grids.future)) {
  env.grids.future[[i]][[1:11]] <- env.grids.future[[i]][[1:11]]/10
}


## give the current and future environmental grids the same names
names(env.grids.current) <- names(env.grids.future) <- c(
  'Annual_mean_temp',    'Mean_diurnal_range',
  'Isothermality',       'Temp_seasonality',  'Max_temp_warm_month',
  'Min_temp_cold_month', 'Temp_annual_range', 'Mean_temp_wet_qu',
  'Mean_temp_dry_qu',    'Mean_temp_warm_qu', 'Mean_temp_cold_qu',  'Annual_precip',
  'Precip_wet_month',    'Precip_dry_month',  'Precip_seasonality', 'Precip_wet_qu',
  'Precip_dry_qu',       'Precip_warm_qu',    'Precip_col_qu')


## Check the rasters: list of 19 with 11 slots
str(env.grids.current)
str(env.grids.future)


## Read in a shapefile for Australia, but exclude the Islands
aus <- ne_states(country = 'Australia') %>% 
  subset(!grepl('Island', name))


## Create a list of all the species folders which contain the fitted models: these were run in the previous step
species_list  <- basename(list.dirs('F:/green_cities_sdm/output/maxent/baseline',   recursive = FALSE))
# scenario_list <- basename(list.dirs('//sci-7910/F/data/worldclim/aus/0.5/bio/2050', recursive = FALSE))


## The trial species based on Renee's species
test.spp = sort(unique(c(renee.full$Species, "Betula pendula", "Fraxinus excelsior", "Quercus robur", "Fagus sylvatica")))
test.spp[35]


#########################################################################################################################
## Now save/load .RData file for the next session
save.image("STEP_7_PREDICT_SDM.RData")
load("STEP_7_PREDICT_SDM.RData")





#########################################################################################################################
## 2). CREATE MAPS OF CURRENT AND FUTURE HABITAT SUITABILITY FOR ALL RECORDS AND SELECTED VARIABLES
#########################################################################################################################


## Currenlty, this code is not working on a single value,  


#########################################################################################################################
## Use lappy to loop over a list of species
## Test on one species:
species = species_list[458] # [1] "Lomandra_longifolia"


#########################################################################################################################
## Also, create a list of directories to loop over
## Now run the code over a list of species...
lapply(species_list, function(species) {
  message('Doing ', species)
  
  lapply(scen, function(scen_i) {
    message('  Doing ', scen_i)
    
    scen_name = gcms$GCM[gcms$id == scen_i] 
    
    ## Read in the fitted models using sprintf
    m <- readRDS(sprintf('F:/green_cities_sdm/output/maxent/SEL_VAR_ALL/%s/maxent_fitted.rds', species))    ## Change dir
    
    # ## These numbers don't look right 
    # str(m)
    # env.grids.current[[colnames(m$me_full@presence)]]
    # env.grids.future[[colnames(m$me_full@presence)]]
    
    ## Read in the occurrence files from the output directory using sprintf
    occ <- readRDS(sprintf('F:/green_cities_sdm/output/maxent/SEL_VAR_ALL/%s/occ.rds', species)) %>%        ## Change dir
      spTransform(CRS('+init=epsg:4326'))
    
    ## str(occ)
    ## Create rasters for the current and future climate
    ## Calculating contribution of feature 11 of 11..............why 11 and not 19?
    ## Or are features the maxent setting
    ## Also this takes a lot of time...why is the future prediction so slow?
    pred.current <- rmaxent::project(m$me_full, env.grids.current[[colnames(m$me_full@presence)]])$prediction_logistic
    pred.future  <- rmaxent::project(m$me_full, env.grids.future[[scen_i]][[colnames(m$me_full@presence)]])$prediction_logistic
    
    ## str(pred.current)
    ## str(pred.future)
    
    ## Write the current raster out
    ## The folders will need to change as I add model runs for all variables vs select, 
    ## all records vs. cultivated and non-cultivated
    ## does this need to be indexed? pred.current[[1]]
    writeRaster(pred.current, sprintf('F:/green_cities_sdm/output/maxent/SEL_VAR_ALL/%s/full/%s_current.tif',  ## Change dir
                                      species, species))
    # 
    # # Warning message:
    # #   In unlist(lapply(elist, findLocals1, shadowed, cntxt)) :
    # #   closing unused connection 3 (F:/green_cities_sdm/RTEMP/RtmpQxrX5c/raster/r_tmp_2017-11-01_165316_9480_78386.gri)
    # 
    # ## Write the future raster out: does this need to be indexed? pred.future[[1]]
    writeRaster(pred.future, sprintf('F:/green_cities_sdm/output/maxent/SEL_VAR_ALL/%s/full/%s_%s.tif',        ## Change dir
                                     species, species, scen_i))
    
    ## Create an empty raster based on the future prediction
    empty <- init(pred.future, function(x) NA)
    
    # Warning message:
    #   In .rasterFromRasterFile(grdfile, band = band, objecttype, ...) :
    #   size of values file does not match the number of cells (given the data type)
    
    
    #########################################################################################################################
    ## Create map of habitat suitability...the first line starts the PNG device
    # Error in compareRaster(x) : different extent
    png(sprintf('F:/green_cities_sdm/output/maxent/SEL_VAR_ALL/%s/full/%s.png', species, species),              ## Change dir
        11, 4, units = 'in', res = 300)
    
    ## Use the levelplot function to make a multipanel output
    levelplot(stack(empty,
                    pred.current,
                    pred.future), margin = FALSE, 
              
              ## Create a colour scheme using colbrewer: 100 is to make it continuos
              scales      = list(draw = FALSE), at = seq(0, 1, length = 100),
              col.regions = colorRampPalette(rev(brewer.pal(11, 'Spectral'))),
              
              ## Give each plot a name
              names.attr = c('Occurrence', 'Current', sprintf('%s, 2050, RCP8.5', scen_name)),
              colorkey   = list(height = 0.6, width = 3), xlab = '', ylab = '',
              main       = list(gsub('_', ' ', species), font = 4, cex = 2)) +
      
      ## Plot the Aus shapefile with the occurrence points for reference
      layer(sp.polygons(aus)) +
      layer(sp.points(occ, pch = 20, cex = 0.5, 
                      col = c('red', 'transparent', 'transparent')[panel.number()]))
    
    ## finish the PNG device
    dev.off()
    
  })
})
  
  
 

  

  
  



#########################################################################################################################
## OUTSTANDING PREDICTION TASKS:
#########################################################################################################################


## Create a switch to skip files that exist

## Create folder structure to hold the output

## 


  
  
  
  
#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################