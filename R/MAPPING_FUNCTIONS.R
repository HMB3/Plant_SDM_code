#########################################################################################################################
############################################# FUNCTIONS FOR MAPPING SDMS ################################################ 
#########################################################################################################################


#########################################################################################################################
## PROJECT MAXENT MODELS
#########################################################################################################################

## flag isses with.......................................................................................................

#########################################################################################################################
## Create maxent maps for a given time period 
# x             = scen_2030[1]
# species       = kop_spp[1]
# time_slice    = 30
# maxent_path   = "./output/maxent/SET_VAR_COORDCLEAN"
# climate_path  = "./data/base/worldclim/aus/1km/bio"
# grid_names    = grid.names
# current_grids = env.grids.current


## One function for all time periods
project_maxent_grids = function(scen_list, species_list, maxent_path, climate_path, grid_names, time_slice, current_grids) {
  
  ## Read in Australia
  aus = readRDS("F:/green_cities_sdm/data/base/CONTEXTUAL/aus_states.rds") %>%
    spTransform(CRS('+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'))
  
  ## First, run a loop over each scenario:    
  lapply(scen_list, function(x) {
    
    ## Create a raster stack for each of the 6 GCMs, not for each species
    s <- stack(sprintf('%s/20%s/%s/%s%s.tif', climate_path, time_slice, x, x, 1:19))
    #s <- stack(sprintf('%s/20%s/%s/%s%s.tif', climate_path, time_slice, x, x, 1:19))
    projection(s)
    
    ## Rename both the current and future environmental stack...
    names(s) <- names(current_grids) <- grid_names 
    
    ########################################################################################################################
    ## Divide the 11 temperature rasters by 10: NA values are the ocean
    ## s[[1:11]] <- s[[1:11]]/10 ## That code doesn't work
    message('20', time_slice, ' rasters / 10 ', x)
    for(i in 1:11) {
      
      ## Simple loop
      message(i)
      s[[i]] <- s[[ i]]/10
      
    }
    
    ## Then apply each GCM to each species
    lapply(species_list, function(species) {
      
      ## First, check if the maxent model exists
      save_name = gsub(' ', '_', species)
      if(file.exists(sprintf('%s/%s/full/maxent_fitted.rds', maxent_path, species))) {
        message('Doing ', species)
        
        ## Then, check if the species projection has already been run...
        if(!file.exists(sprintf('%s/%s/full/%s_%s.tif', maxent_path, species, species, x))) {
          
          ## Assign the scenario name (to use later in the plot)
          scen_name = eval(parse(text = sprintf('gcms.%s$GCM[gcms.%s$id == x]', time_slice, time_slice)))           
          
          ########################################################################################################################
          ## Now read in the SDM model calibrated on current conditions  ## maxent_fitted.rds
          #m <- readRDS(sprintf('%s/%s/full/model.rds', maxent_path, species)) 
          m <- readRDS(sprintf('%s/%s/full/maxent_fitted.rds', maxent_path, species)) 
          m <- m$me_full  ## class(m);View(m)
          
          ## Read in the occurrence points used to create the SDM :: need the transform to plot later
          occ <- readRDS(sprintf('%s/%s/%s_occ.rds', maxent_path, species, save_name)) %>%
            spTransform(CRS('+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'))  
          ## '+init=epsg:4326' ## '+init=ESRI:3577' ## '+init=ESRI:54009'
          
          ## Create file path fpr current raster doesn't exist, create it
          f_current <- sprintf('%s/%s/full/%s_current.tif', maxent_path, species, species)
          
          if(!file.exists(f_current)) {
            
            ## Report which prediction is in progress :: m$me_full, m$me_full@presence
            message('Running current prediction for ', species) 
            
            pred.current <- rmaxent::project(
              m, current_grids[[colnames(m@presence)]])$prediction_logistic
            writeRaster(pred.current, f_current, overwrite = TRUE)
            
          } else {
            
            pred.current = raster(sprintf('%s/%s/full/%s_current.tif',
                                          maxent_path, species, species))
          }
          
          ########################################################################################################################
          ## Create file path for future raster doesn't exist, create it 
          f_future <- sprintf('%s/%s/full/%s_%s.tif', 
                              maxent_path, species, species, x)
          
          if(!file.exists(f_future)) {
            
            ## Report which prediction is in progress
            message('Running future prediction for ', species, ' ', x) 
            
            ## Create the future raster
            pred.future <- rmaxent::project(
              m, s[[colnames(m@presence)]])$prediction_logistic
            writeRaster(pred.future, f_future, overwrite = TRUE)
            
            ## Now create the empty panel just before plotting
            empty_ras <- init(pred.current, function(x) NA) 
            
            ## Check projections ..................................................................................................
            projection(aus);projection(occ);projection(empty_ras)
            projection(pred.current);projection(pred.future)
            
            ########################################################################################################################
            ## Use the levelplot function to make a multipanel output: occurrence points, current raster and future raster
            png(sprintf('%s/%s/full/%s_%s.png', maxent_path, species, species, x),      
                11, 4, units = 'in', res = 300)
            
            ## Need an empty frame
            print(levelplot(stack(empty_ras,
                                  pred.current, 
                                  pred.future, quick = TRUE), margin = FALSE,
                            
                            ## Create a colour scheme using colbrewer: 100 is to make it continuos
                            ## Also, make it a one-directional colour scheme
                            scales      = list(draw = FALSE), 
                            at = seq(0, 1, length = 100),
                            col.regions = colorRampPalette(rev(brewer.pal(11, 'Spectral'))),
                            
                            ## Give each plot a name: the third panel is the GCM
                            names.attr = c('Australian records', 'Current', sprintf('%s, 20%s, RCP8.5', scen_name, time_slice)),
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
          
          message(species, ' ', x, ' skipped - prediction already run') ## Skip species which have already been projected
          
        }
        
      } else {
        
        message(species, ' ', x, ' skipped - SDM not yet run')          ## Skip species with no existing maxent model
        
      }
      
    })
    
  })
  
}





#########################################################################################################################
## COMBINE GCM FUNCTIONS
#########################################################################################################################


#########################################################################################################################
## Loop over directories, species and one threshold for each, also taking a time_slice argument. Next, make the lists generic too
## pervious version in R/old/model_combine.R
combine_gcm_threshold = function(DIR_list, species_list, maxent_path, thresholds, percentiles, time_slice, area_occ) {
  
  ## How can the shapefiles be read in once, not for each species?...................................................
  
  ###################################################################################################################
  ## Read in shapefiles :: this should be done outside the loop
  aus        = readRDS("F:/green_cities_sdm/data/base/CONTEXTUAL/aus_states.rds") %>%
    spTransform(CRS('+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'))
  
  LAND       = readRDS("F:/green_cities_sdm/data/base/CONTEXTUAL/LAND_world.rds")
  
  areal_unit = readRDS("F:/green_cities_sdm/data/base/CONTEXTUAL/SUA.rds") %>%
    spTransform(CRS('+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'))
  
  areal_unit = areal_unit[order(areal_unit$SUA_NAME11),]
  
  ## Loop over each directory
  lapply(DIR_list, function(DIR) { 
    
    ## And each species - although we don't want all possible combinations. use mapply in the function call
    lapply(species_list, function(species) {
      
      # for (slice in time_slice) {
      
      ###################################################################################################################
      ## Create a list of the rasters in each directory, then take the mean. Tidy this up with %, etc
      message('Calcualting mean of GCMs for ', species)
      
      ## Check if the mean GCM raster exists
      f_mean = sprintf('%s/%s/full/%s_20%s_suitability_mean.tif', maxent_path, species, species, time_slice)
      
      ## The mean of the GCMs doesn't exist, create it
      if(!file.exists(f_mean)) {
        
        raster.list       = list.files(as.character(DIR), pattern = sprintf('bi%s.tif', time_slice), full.names = TRUE)  
        suit              = stack(raster.list)
        suit.list         = unstack(suit)
        combo_suit_mean   = mean(suit)                            ## plot(mean.suit)
        
        writeRaster(combo_suit_mean , sprintf('%s/%s/full/%s_20%s_suitability_mean.tif', 
                                              maxent_path, species, species, time_slice), overwrite = TRUE)
        
      } else {
        
        ## Create another level without the mean calculation
        raster.list = list.files(as.character(DIR), pattern = sprintf('bi%s.tif', time_slice), full.names = TRUE)  
        suit        = stack(raster.list)
        suit.list   = unstack(suit)
        
      }
      
      #########################################################################################################################
      ## Then, create rasters that meet habitat suitability criteria thresholds, determined by the rmaxent function
      for (thresh in thresholds) {
        
        for (percent in percentiles) { 
          
          ## Can we set 0 to NA in the rasters before running the calculations?
          ## Check if the combined suitability raster exists
          f_max_train_suit <- sprintf('%s/%s/full/%s_20%s%s%s.tif', maxent_path,
                                      species, species, time_slice, "_Max_train_sensit_above_", thresh)
          
          ## If it doesn't exist, create the suitability raster
          if(!file.exists(f_max_train_suit)) {
            
            ## Print the species being analysed
            message('doing ', species, ' | Max train sensit > ', thresh, ' for 20', time_slice)
            
            ## Read in the current suitability raster
            f_current <- raster(sprintf('%s/%s/full/%s_current.tif', 
                                        maxent_path, species, species))
            
            ## First, create a simple function to threshold each of the rasters in raster.list
            thresh_greater  = function (x) {x > thresh}
            percent_greater = function (x) {x > percent}
            
            ## Then apply this to just the current suitability raster. These functions use the : 
            ## Maximum training sensitivity plus specificity Logistic threshold
            ## 10th percentile training presence training omission
            current_suit_thresh  = thresh_greater(f_current)
            current_suit_percent = percent_greater(f_current) 
            
            ## Now, apply these functions to the list of rasters (6 GCMs) for each species
            ## Also, the 'init' function initializes a raster object with values
            ## suit_ras_greater    = reduce(suit.list, thresh_above_fun, .init = suit.list[[1]] > thresh)
            
            ## Check the logic of removing the zeros.....................................................
            ## If we are just using the combined rasters without calculating the difference, don't worry about the zeros   
            
            #########################################################################################################################
            ## First, calculate the cells which are greater that the: 
            ## Maximum training sensitivity plus specificity Logistic threshold
            message('Running thresholds for ', species, ' | 20', time_slice, ' combined suitability > ', thresh)
            
            suit_ras1_thresh   = thresh_greater(suit.list[[1]])   ## Abbreviate these...
            suit_ras2_thresh   = thresh_greater(suit.list[[2]])
            suit_ras3_thresh   = thresh_greater(suit.list[[3]])
            suit_ras4_thresh   = thresh_greater(suit.list[[4]])   
            suit_ras5_thresh   = thresh_greater(suit.list[[5]])
            suit_ras6_thresh   = thresh_greater(suit.list[[6]])
            
            ## Then calculate the cells which are greater than the 10th percentile training presence training omission
            suit_ras1_percent  = percent_greater(suit.list[[1]])
            suit_ras2_percent  = percent_greater(suit.list[[2]])
            suit_ras3_percent  = percent_greater(suit.list[[3]])
            suit_ras4_percent  = percent_greater(suit.list[[4]])
            suit_ras5_percent  = percent_greater(suit.list[[5]])
            suit_ras6_percent  = percent_greater(suit.list[[6]])
            
            #########################################################################################################################
            ## Then sum them up: All the threshholds
            combo_suit_thresh   =  Reduce("+", list(suit_ras1_thresh, suit_ras2_thresh, suit_ras3_thresh,
                                                    suit_ras4_thresh, suit_ras5_thresh, suit_ras6_thresh))
            
            ## All the percentiles
            combo_suit_percent  =  Reduce("+", list(suit_ras1_percent, suit_ras2_percent, suit_ras3_percent,
                                                    suit_ras4_percent, suit_ras5_percent, suit_ras6_percent))

            #########################################################################################################################
            ## For each species, create a binary raster with cells > 4 GCMs above the maxent threshold = 1, and cells with < 4 GCMs = 0. 
            message('Calculating change for ', species, ' | 20', time_slice, ' combined suitability > ', thresh)
            
            ## Functions for thresholding rasters
            #band_ras     <- function(x) {ifelse(x >=  4, 4, ifelse(x > 0 & x < 4, 3, x)) }
            band_4           <- function(x) {ifelse(x >=  4, 1, 0) }
            combo_suit_4GCM  <- calc(combo_suit_thresh, fun = band_4)
            
            #########################################################################################################################
            ## Subtract the current from the future
            Combo_future_minus_current = combo_suit_4GCM - current_suit_thresh
            
            ## Plot to check
            plot(current_suit_thresh, main = gsub('_', ' ', (sprintf('%s current max_train_sensit > %s', species, thresh))))
            plot(combo_suit_percent,  main = gsub('_', ' ', (sprintf('%s future 10th percentile > %s',   species, percent))))
            plot(combo_suit_thresh,   main = gsub('_', ' ', (sprintf('%s future max_train_sensit > %s',  species, thresh))))
            plot(combo_suit_4GCM,     main = gsub('_', ' ', (sprintf('%s 4+ GCMs > %s',  species, thresh))))
            plot(Combo_future_minus_current,     main = gsub('_', ' ', (sprintf('%s 4+ GCMs > %s - current',  species, thresh))))
            
            #########################################################################################################################
            ## For each species, calculate the projected rainfall and temperature increase and decreas for each GCM? Could plot this
            ## as part of the final figure.
            
            #########################################################################################################################
            ## Then using this GCM consensus, calculate whether the species is likely to be present in each SUA.
            ## Decide on a threshold of % area (10?) of the SUA that needs to be occupied, for each species to be considered present. 
            message('Running SUA area calculations for ', species, ' | 20', time_slice, ' combined suitability > ', thresh)
            
            ## Check the order of lists match, species, SUAs, areas need to match up ................................................
            
            ## Should we also create a simple count of the no. of cells per SUA with where > 4 GCMs met the suitability threshold?
            ## Need to fix this so that it has the same order as the shapefile
            # sp.count <- spatialEco::zonal.stats(x = areal_unit, y = combo_suit_4GCM, stat = sum, trace = TRUE, plot = TRUE) 
            # z.med   <- spatialEco::zonal.stats(x = areal_unit, y = combo_suit_4GCM, stat = median, trace = TRUE, plot = TRUE)
            
            #########################################################################################################################
            ## Then, extract the values of the presence raster for each areal unit: generates a list
            ext.current  <- extract(current_suit_thresh, areal_unit, method = 'simple')
            ext.future   <- extract(combo_suit_4GCM,     areal_unit, method = 'simple')
            
            ## A function to tabulate the raster values by aerial unit, returning a data frame
            tabFunc <- function(indx, extracted, region, regname) {
              
              dat<-as.data.frame(table(extracted[[indx]]))
              dat$name<-region[[regname]][[indx]]
              return(dat)
              
            }
            
            ## Run through each areal unit and calculate a table of the count of raster cells
            tabs.current <- lapply(seq(ext.current), tabFunc, ext.current, areal_unit, "SUA_NAME11")
            tabs.current <- do.call("rbind", tabs.current )
            
            tabs.future  <- lapply(seq(ext.future), tabFunc, ext.future, areal_unit, "SUA_NAME11")
            tabs.future  <- do.call("rbind", tabs.future )
            
            ## Can we get the count here, so we don't need to run zonal stats as well
            # tabs.count      = tabs
            # tabs.count$Freq = ifelse(tabs.count$Var1 == 0, 0, tabs.count$Freq) 
            # sp.count <- aggregate(tabs.count$Freq, by = list(Category = tabs$name), FUN = sum)
            # names(sp.count) = c('SUA_NAME11', 'COUNT')
            # head(sp.count)
            
            #########################################################################################################################
            ## Now mutate the current table
            PERECENT.AREA.CURRENT <- tabs.current %>%
              
              group_by(name) %>%                                          ## Group by region
              mutate(totcells = sum(Freq),                                ## How many cells overall?
                     percent.area = round(100 * Freq / totcells, 2)) %>%  ## Cells divided by total cells
              
              dplyr::select(-c(Freq, totcells)) %>%                       ## There is a select func in raster so need to specify
              spread(key = Var1, value = percent.area, fill = 0)     %>%  ## Make wide format
              as.data.frame()
            
            #########################################################################################################################
            ## Now mutate the future table
            PERECENT.AREA.FUTURE <- tabs.future %>%
              
              group_by(name) %>%                                          ## Group by region
              mutate(totcells = sum(Freq),                                ## How many cells overall?
                     percent.area = round(100 * Freq / totcells, 2)) %>%  ## Cells divided by total cells
              
              dplyr::select(-c(Freq, totcells)) %>%                       ## There is a select func in raster so need to specify
              spread(key = Var1, value = percent.area, fill = 0)     %>%  ## Make wide format
              as.data.frame()
            
            ## Rename and create a column for whether or not the species occupies that area 
            names(PERECENT.AREA.CURRENT) =  c('SUA_NAME11', 'Absent', 'Present') 
            names(PERECENT.AREA.FUTURE)  =  c('SUA_NAME11', 'Absent', 'Present') 
            
            PERECENT.AREA.CURRENT$species_present = ifelse(PERECENT.AREA.CURRENT$Present >= area_occ, 1, 0)
            PERECENT.AREA.FUTURE$species_present  = ifelse(PERECENT.AREA.FUTURE$Present  >= area_occ, 1, 0)
            
            head(PERECENT.AREA.CURRENT)
            head(PERECENT.AREA.FUTURE)
            
            ## Create a table with the columns: REGION, AREA SPECIES, STATS (for each time slice)
            ## ifelse(<condition>, <yes>, ifelse(<condition>, <yes>, <no>))
            GCM.AREA.SUMMARY <- data.frame(SUA          = areal_unit$SUA_NAME11, 
                                           AREA_SQKM    = areal_unit$AREA_SQKM,
                                           SPECIES      = species,
                                           PERIOD       = time_slice,
                                           AREA_THRESH  = area_occ,
                                           MAX_TRAIN    = thresh,
                                           PERCENT_AREA = PERECENT.AREA.FUTURE$Present,
                                           AREA_CHANGE  = PERECENT.AREA.FUTURE$Present - PERECENT.AREA.CURRENT$Present,
                                           PRESENT      = PERECENT.AREA.FUTURE$species_present)
            
            
            ## Can we calculate these in the table?
            ## Gain (ie not suitable now but suitable in future)
            ## Loss (ie suitable now but not in future)
            ## Stable (suitable in both)
            GCM.AREA.SUMMARY$GAIN_LOSS  <- with(GCM.AREA.SUMMARY, ifelse(
              GCM.AREA.SUMMARY$AREA_CHANGE <= -1 & GCM.AREA.SUMMARY$AREA_CHANGE < 0, 'LOSS', ifelse(
                GCM.AREA.SUMMARY$AREA_CHANGE >= 1, 'GAIN', ifelse(
                  GCM.AREA.SUMMARY$AREA_CHANGE == 0, 'STABLE', 'GAIN'))))
            
            
            ##
            GCM.AREA.SUMMARY$GAIN_LOSS 
            View(GCM.AREA.SUMMARY) 
            
            #########################################################################################################################
            ## Then save the table of SUA results for all species to a datafile...
            ## This would be the file to loop over to create a summary of species per SUA
            write.csv(GCM.AREA.SUMMARY, sprintf('%s/%s/full/%s_20%s_%s%s.csv', maxent_path,
                                                species, species, time_slice, area_occ, "pc_area_SUA_summary"), row.names = FALSE)
            
            #########################################################################################################################
            #########################################################################################################################
            ## If they don't exist, write the rasters for each species/threshold
            if(!file.exists(f_max_train_suit)) {
              
              ## Write the current suitability raster, thresholded using the Maximum training sensitivity plus specificity Logistic threshold
              message('Writing ', species, ' current', ' max train > ', thresh) 
              writeRaster(current_suit_thresh, sprintf('%s/%s/full/%s_%s%s.tif', maxent_path,
                                                       species, species, "current_suit_above_", thresh), overwrite = TRUE) 
              
              ## Write the combined suitability raster, thresholded using the maximum training value
              message('Writing ', species, ' | 20', time_slice, ' max train > ', thresh) 
              writeRaster(combo_suit_thresh, sprintf('%s/%s/full/%s_20%s%s%s.tif', maxent_path,
                                                     species, species, time_slice, "_Max_train_sensit_above_", thresh), overwrite = TRUE)
              
              ## Write the combined suitability raster, thresholded using the percentile value
              message('Writing ', species, ' | 20', time_slice, ' 10th percentile > ', percent)
              writeRaster(combo_suit_percent, sprintf('%s/%s/full/%s_20%s%s%s.tif', maxent_path,
                                                      species, species, time_slice, "_10_percentile_omiss_above_", percent), overwrite = TRUE)
              
              ## Write the combined future raster with > 4 GCMs above the maximum training value
              message('Writing ', species, ' | 20', time_slice, ' 10th percentile > ', percent) 
              writeRaster(combo_suit_4GCM, sprintf('%s/%s/full/%s_20%s%s%s.tif', maxent_path,
                                                   species, species, time_slice, "_4GCMs_above_", thresh), overwrite = TRUE)
              
              ## Write the future - current thresholded rasters to file
              message('Writing ', species, ' | 20', time_slice, ' Future - current', thresh) 
              writeRaster(Combo_future_minus_current, sprintf('%s/%s/full/%s_20%s%s%s.tif', maxent_path,
                                                              species, species, time_slice, "_minus_current_", thresh), overwrite = TRUE)
              
            } else {
              
              message(species, ' 20', time_slice, ' combined suitability > ', thresh, ' skipped - already exists')   ## 
              
            }
            
            ########################################################################################################################
            ## Now create the empty panel just before plotting, and read in the occurrence and background points and original model
            
            
            # We need to be able to flick through maps. So we want a three panelled page with 
            # a) occ points overlaid on Koppen zones, 
            # b) current continuous suitability, 
            # c) thresholded suitability.
            
            ## And change the projection for the operators...........................................................................
            ## '+init=esri:54009'
            save_name = gsub(' ', '_', species)
            empty <- init(combo_suit_thresh, function(x) NA)
            
            occ <- readRDS(sprintf('%s/%s/%s_occ_swd.rds', maxent_path, species, save_name)) %>%
              spTransform(CRS('+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'))  
            
            bg <- readRDS(sprintf('%s/%s/%s_bg_swd.rds', maxent_path, species, save_name)) %>%
              spTransform(CRS('+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'))
            
            #m <- readRDS(sprintf('%s/%s/full/model.rds', maxent_path, species)) 
            m <- readRDS(sprintf('%s/%s/full/maxent_fitted.rds', maxent_path, species)) 
            m <- m$me_full  ## class(m);View(m)
            
            ## Use the 'levelplot' function to make a multipanel output: occurences, percentiles and thresholds
            message('Writing figure for ', species, ' | 20', time_slice, ' > ', thresh) 
            
            png(sprintf('%s/%s/full/%s_20%s_suitability_above_%s.png', maxent_path,
                        species, species, time_slice, thresh),
                11, 4, units = 'in', res = 300)
            
            ## Need an empty frame
            print(levelplot(stack(empty,                ## needs to have a different colour scale,
                                  combo_suit_percent,
                                  combo_suit_thresh), margin = FALSE,
                            
                            ## Create a colour scheme using colbrewer: 100 is to make it continuos
                            ## Also, make it a one-directional colour scheme
                            scales      = list(draw = FALSE), 
                            at = seq(0, 6, length = 100),
                            col.regions = colorRampPalette(rev(brewer.pal(11, 'Spectral'))),
                            
                            ## Give each plot a name
                            names.attr = c('Aus occurrences', 
                                           sprintf('20%s GCM 10thp train omission > %s', time_slice, percent), 
                                           sprintf('20%s GCM Max train logis > %s',      time_slice, thresh)),
                            
                            colorkey   = list(height = 0.5, width = 3), xlab = '', ylab = '',
                            main       = list(gsub('_', ' ', species), font = 4, cex = 2)) +
                    
                    ## Plot the Aus shapefile with the occurrence points for reference
                    ## Can we assign different shapefiles to different panels, rather than to them all?
                    
                    layer(sp.polygons(aus)) + ## sp.polygons(aus))
                    layer(sp.points(occ, pch = 20, cex = 0.4, 
                                    col = c('red', 'transparent', 'transparent')[panel.number()]), data = list(occ = occ)))
            
            ## Finish the device
            dev.off()
            
            ########################################################################################################################
            ## Another .png for the global records: str(LAND$long)
            # if(!file.exists(sprintf('%s/%s/full/%s_%s.png', maxent_path, species, species, "global_records"))) {
            #   
            #   png(sprintf('%s/%s/full/%s_%s.png', maxent_path,
            #               species, species, "global_records"),
            #       16180, 10000, units = 'px', res = 600)
            #   
            #   ## How do we locate bad records in the dataset after spotting them?
            #   plot(LAND, 
            #        lwd = 0.5, asp = 1, axes = TRUE, cex.axis = 3.5,
            #        col = 'darkolivegreen3', bg = 'lightblue', cex.lab = 3)
            #   
            #   points(occ, pch = ".", cex = 3.5, col = "red", cex.lab = 3, cex.main = 4, cex.axis = 2, 
            #          main = paste0("Global occurrences for ", species), 
            #          xlab = "", ylab = "", asp = 1)
            #   
            #   ## Title 
            #   title(paste0("Global points for ", species),
            #         cex.main = 4,   font.main = 4, col.main = "blue")
            #   
            #   ## Finsh the device
            #   dev.off()
            #   
            # } else {
            #   
            #   message("Global records maps exists for ", species)
            #   
            # }
            
            ########################################################################################################################
            ## Another PNG for the background points....
            # if(!file.exists(sprintf('%s/%s/full/%s_%s.png', maxent_path, species, species, "background_records"))) {
            #   
            #   png(sprintf('%s/%s/full/%s_%s.png', maxent_path,
            #               species, species, "background_records"),
            #       16180, 10000, units = 'px', res = 600)
            #   
            #   ## How do we locate bad records in the dataset after spotting them?
            #   plot(LAND,  
            #        lwd = 0.5, asp = 1, axes = TRUE, cex.axis = 3.5,
            #        col = 'darkolivegreen3', bg = 'lightblue', cex.lab = 3)
            #   
            #   points(bg, pch = ".", cex = 1.6, col = "blue", cex.lab = 3, cex.main = 4, cex.axis = 2, 
            #          main = paste0("Global occurrences for ", species), 
            #          xlab = "", ylab = "", asp = 1)
            #   
            #   ## Title 
            #   title(paste0("Bacground points for ", species),
            #         cex.main = 4,   font.main = 4, col.main = "blue")
            #   
            #   ## Finish the device
            #   dev.off() 
            #   
            # } else {
            #   
            #   message("Background records maps exists for ", species)
            #   
            # }
            
            ########################################################################################################################
            ## Plot the models: can two plots be combined into one?
            ## Make these unique names, and they can be searched in windows. Otherwise, we can just click into each subfolder. 
            ## To sort, names would need to be: spp + unique_extension
            if(!file.exists(sprintf('%s/%s/full/%s_current.png', maxent_path, species, species, "variable_contribution"))) {
              
              png(sprintf('%s/%s/full/%s_current.png', maxent_path,
                          species, species, "variable_contribution"),
                  3236, 2000, units = 'px', res = 300)
              
              ## Set the margins
              # par(mgp      = c(10, 4, 0), 
              #     oma      = c(1.5, 1.5, 1.5, 1.5),
              #     font.lab = 2)
              
              plot(m, col = "blue", pch = 19, cex.lab = 2, cex.axis = 5, cex.main = 2, 
                   main   = paste0("Variables for ", species), 
                   xlab   = "Maxent contribution (%)")
              
              ## Finish the device
              dev.off()
              
              ## Plot the response curves too
              png(sprintf('%s/%s/full/%s_%s.png', maxent_path,
                          species, species, "response_curves"),
                  3236, 2000, units = 'px', res = 300)
              
              ## Add detail to the response plot
              response(m, pch = 19, cex.lab = 2, cex.axis = 1.5, lwd = 2) 
              #ylab   = "",
              #main   = paste0(species, " responses"))
              
              ## Finish the device
              dev.off()
              
            } else {
              
              message("Variable contribution plot exists for ", species)
              
            }
            
            ## current.list = as.data.frame(env.grids.current)
            ## test = similarity(occ, env.grids.current, full = FALSE)
            
          } else {
            
            message(species, ' 20', time_slice, ' combined suitability > ', thresh, ' skipped - already exists')   ## 
            
          }
          
        }
        
      }
      
    })
    
  })
  
}





#########################################################################################################################
## ANOMALY FUNCTIONS
#########################################################################################################################


#########################################################################################################################
## For each scenario, calcualte the mean annual temperature and annual rainfall anomaly
calculate.anomaly = function(scen_list, time_slice, climate_path) {
  
  #########################################################################################################################
  ## Create current rasters :: can the raster creation happen just once?
  ## And read in the current data. Can the specific "current" be removed withouth makin a seconf path argument?
  message('current rasters / 10')
  
  env.grids.current       = stack(file.path(sprintf('%s/current/bio_%02d.tif', climate_path, 1:19)))
  env.grids.current[[1]]  = env.grids.current[[1]]/10
  current.bio1            = env.grids.current[[1]]
  current.bio12           = env.grids.current[[12]]
  
  ## First, run a loop over each scenario :: 
  lapply(scen_list, function(x) {
    
    ## Assign the scenario name (to use later in the plot)
    scen_name = eval(parse(text = sprintf('gcms.%s$GCM[gcms.%s$id == x]', time_slice, time_slice)))
    
    #########################################################################################################################
    ## Create a raster stack for each 2050 GCM - also an empty raster for the final plot
    s <- stack(sprintf('%s/20%s/%s/%s%s.tif', climate_path, time_slice, x, x, 1:19))
    
    ########################################################################################################################
    ## Create future rasters
    message('20', time_slice, ' rasters / 10')
    
    s[[1]]       = s[[1]]/10
    future.bio1  = s[[1]]
    future.bio12 = s[[12]]
    
    ## Now subtract the current from the future...
    temp.anomaly = future.bio1  - current.bio1
    rain.anomaly = future.bio12 - current.bio12
    
    plot(temp.anomaly, main = paste0("BIO1  anomaly ", x))
    plot(rain.anomaly, main = paste0("BIO12 anomaly ", x))
    
    ########################################################################################################################
    ## Write the rasters for each species/threshold
    message('Writing BIO1/12 anomaly rasters for 20', time_slice, ' ', x)
    
    if(!file.exists(sprintf('%s/anomalies/%s_%s.tif', climate_path, "BIO1_anomaly",  x))) {
      
      ## Temp anomalies
      writeRaster(temp.anomaly, sprintf('%s/anomalies/%s_%s.tif', climate_path, "BIO1_anomaly",  x), overwrite = TRUE)
      
    } else {
      
      message('Skipped ', "BIO1 anomaly for ", x, ', already exists')
      
    }
    
    if(!file.exists(sprintf('%s/anomalies/%s_%s.tif', climate_path, "BIO12_anomaly",  x))) {
      
      ## Rain anomalies
      writeRaster(rain.anomaly, sprintf('%s/anomalies/%s_%s.tif', climate_path, "BIO12_anomaly", x), overwrite = TRUE)
      
    } else {
      
      message('Skipped ', "BIO12 anomaly for ", x, ', already exists')
      
    }
    
  })
  
} 





#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################