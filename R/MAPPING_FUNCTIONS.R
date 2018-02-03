#########################################################################################################################
############################################# FUNCTIONS FOR MAPPING SDMS ################################################ 
#########################################################################################################################


#########################################################################################################################
## PROJECT MAXENT MODELS
#########################################################################################################################

## flag isses with.......................................................................................................

#########################################################################################################################
## 2050
project.grids.2050 = function(scen_2050, test_spp) {
  
  ##  First, run a loop over each scenario: options(error = recover)    
  lapply(scen_2050, function(x) {
    
    ## Assign the scenario name (to use later in the plot)
    scen_name = gcms.50$GCM[gcms.50$id == x]                       
    
    ## Create a raster stack for each 2050 GCM - also an empty raster for the final plot
    s <- stack(
      sprintf('./data/base/worldclim/aus/0.5/bio/2050/%s/%s%s.tif',
              x, x, 1:19))
    
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
    
    # for(i in 1:11) {
    #   
    #   message(i)
    #   s[[i]] <- s[[ i]]/10
    #   
    # }
    
    message('Doing raster calculation 2050')
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
    lapply(test_spp, function(species) {
      
      ## First check if the species projection has already been run...
      ## Check model exists: needs this condition..........................................................................
      if(file.exists(sprintf('./output/maxent/STD_VAR_ALL/%s/full/model.rds', species))) {
        message('Doing ', species)
        
        ## immediately after?
        if(!file.exists(sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_%s.tif',
                                species, species, x))) {
          message('Projecting ', species) 
          
          ## Read in the SDM model calibrated on current conditions
          #m <- readRDS(sprintf('./output/maxent/STD_VAR_ALL/%s/maxent_fitted.rds', species))
          m <- readRDS(sprintf('./output/maxent/STD_VAR_ALL/%s/full/model.rds', species)) 
          
          ## Read in the occurrence points used to create the SDM
          ## And if the current raster doesn't exist, create it
          occ <- readRDS(
            sprintf('./output/maxent/STD_VAR_ALL/%s/occ.rds', 
                    species)) %>%
            
            spTransform(CRS('+init=epsg:4326'))
          f_current <- sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_current.tif', 
                               species, species)
          
          if(!file.exists(f_current)) {
            
            ## Report which prediction is in progress m$me_full, m$me_full@presence
            message('Running current prediction for ', species) 
            
            pred.current <- rmaxent::project(
              m, env.grids.current[[colnames(m@presence)]])$prediction_logistic
            writeRaster(pred.current, f_current, overwrite = TRUE)
            
          } else {
            
            pred.current = raster(sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_current.tif', 
                                          species, species))
          }
          
          ########################################################################################################################
          ## If the future raster doesn't exist, create it 
          f_future <- sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_%s.tif', 
                              species, species, x)
          
          if(!file.exists(f_future)) {
            
            ## Report which prediction is in progress
            message('Running future prediction for ', species, ' ', x) 
            
            pred.future <- rmaxent::project(
              m, s[[colnames(m@presence)]])$prediction_logistic
            writeRaster(pred.future, f_future, overwrite = TRUE)
            
            ## Now create the empty panel just before plotting
            empty <- init(pred.future, function(x) NA)
            
            ########################################################################################################################
            ## Use the levelplot function to make a multipanel output: occurrence points, current raster and future raster
            png(sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_%s.png', species, species, x),      
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
          
          message(species, ' ', x, ' skipped - prediction already run') ## Skip species which have already been projected
          
        }
        
      } else {
        
        message(species, ' ', x, ' skipped - SDM not yet run')   ## Skip species with no existing SDM
        
      }
      
    })
    
  })
  
}





#########################################################################################################################
## 2070
project.grids.2070 = function(scen_2070, test_spp, time_slice) {
  
  ##  First, run a loop over each scenario: options(error = recover)    
  lapply(scen_2070, function(x) {
    
    ## Assign the scenario name (to use later in the plot)
    scen_name = gcms.70$GCM[gcms.70$id == x]
    
    ## Create a raster stack for each 2070 GCM - also an empty raster for the final plot
    s <- stack(
      sprintf('./data/base/worldclim/aus/0.5/bio/2070/%s/%s%s.tif',
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
    message('Doing raster calculation 2070')
    
    # for(i in 1:11) {
    #   
    #   message(i)
    #   s[[i]] <- s[[ i]]/10
    #   
    # }
    
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
    lapply(test_spp, function(species) {
      
      ## First check if the species projection has already been run...
      if(file.exists(sprintf('./output/maxent/STD_VAR_ALL/%s/full/model.rds', species))) {
        message('Doing ', species)
        
        ## The, run the projections if each projections doesn't exist
        if(!file.exists(sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_%s.tif',
                                species, species, x))) {
          message('Projecting ', species) 
          
          ## Read in the SDM model calibrated on current conditions
          #m <- readRDS(sprintf('./output/maxent/STD_VAR_ALL/%s/maxent_fitted.rds', species))
          m <- readRDS(sprintf('./output/maxent/STD_VAR_ALL/%s/full/model.rds', species)) 
          
          ## Read in the occurrence points used to create the SDM
          ## And if the current raster doesn't exist, create it
          occ <- readRDS(
            sprintf('./output/maxent/STD_VAR_ALL/%s/occ.rds', 
                    species)) %>%
            
            spTransform(CRS('+init=epsg:4326'))
          f_current <- sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_current.tif', 
                               species, species)
          
          if(!file.exists(f_current)) {
            
            ## Report which prediction is in progress m$me_full, m$me_full@presence
            message('Running current prediction for ', species) 
            
            pred.current <- rmaxent::project(
              m, env.grids.current[[colnames(m@presence)]])$prediction_logistic
            writeRaster(pred.current, f_current, overwrite = TRUE)
            
          } else {
            
            pred.current = raster(sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_current.tif', 
                                          species, species))
          }
          
          ########################################################################################################################
          ## If the future raster doesn't exist, create it 
          f_future <- sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_%s.tif', 
                              species, species, x)
          
          if(!file.exists(f_future)) {
            
            ## Report which prediction is in progress
            message('Running future prediction for ', species, ' ', x) 
            
            pred.future <- rmaxent::project(
              m, s[[colnames(m@presence)]])$prediction_logistic
            writeRaster(pred.future, f_future, overwrite = TRUE)
            
            ## Now create the empty panel just before plotting
            empty <- init(pred.future, function(x) NA)
            
            ########################################################################################################################
            ## Use the levelplot function to make a multipanel output: occurrence points, current raster and future raster
            png(sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_%s.png', species, species, x),      
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
          
          message(species, ' ', x, ' skipped - prediction already run') ## Skip species which have already been projected
          
        }
        
      } else {
        
        message(species, ' ', x, ' skipped - SDM not yet run')   ## Skip species with no existing SDM
        
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



#########################################################################################################################
## Loop over directories, species and one threshold for each, also taking a time_slice argument. Next, make the lists generic too
## pervious version in R/old/model_combine.R
combine_gcm_threshold = function(DIR_list, species_list, thresholds, percentiles, time_slice, area_occ) {
  
  lapply(DIR_list, function(DIR) { 
    
    ## And each species - although we don't want all possible combinations. How can this loop be improved?
    lapply(species_list, function(species) {
      
      # for (slice in time_slice) {
      
      ###################################################################################################################
      ## Create a list of the rasters in each directory, then take the mean. Tidy this up with %, etc
      message('Calcualting mean of GCMs for ', species)
      
      raster.list = list.files(as.character(DIR), pattern = sprintf('bi%s.tif', time_slice), full.names = TRUE)  
      suit        = stack(raster.list)
      suit.list   = unstack(suit)
      mean.suit   = mean(suit)                            ## plot(mean.suit)
      #median.suit = median(suit)
      
      ## Write the mean to file
      f_mean = sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_20%s_suitability_mean.tif', 
                       species, species, time_slice)
      
      if(!file.exists(f_mean)) {
        
        writeRaster(mean.suit, sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_20%s_suitability_mean.tif', 
                                       species, species, time_slice), overwrite = TRUE)
        
      } else {
        
        message(species, ' 20', time_slice, ' mean suitability skipped - already exists')   ## 
        
      }
      
      #########################################################################################################################
      #########################################################################################################################
      ## Then create rasters that meet habitat suitability criteria thresholds, determined by the rmaxent function
      for (thresh in thresholds) {
        
        for (percent in percentiles) { 
          
          ## Can we set 0 to NA in the rasters before running the calculations?
          ## Check if the combined suitability raster exists
          f_max_train_suit <- sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_20%s%s%s.tif',
                                      species, species, time_slice, "_Max_train_sensit_above_", thresh)
          
          ## If it doesn't exist, create the suitability raster
          if(!file.exists(f_max_train_suit)) {
            
            ## Print the species being analysed
            message('doing ', species, ' | Max train sensit > ', thresh, ' for 20', time_slice)
            
            ## Read in the current suitability raster
            f_current <- raster(sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_current.tif', 
                                        species, species))
            
            ## First, create a simple function to threshold each of the rasters in raster.list
            thresh_greater  = function (x) {x > thresh}
            percent_greater = function (x) {x > percent}
            
            ## Then apply this to just the current suitability raster. These functions use the : 
            ## Maximum training sensitivity plus specificity Logistic threshold and the :
            ## 10th percentile training presence training omission
            current_suit_thresh  = thresh_greater(f_current)
            current_suit_percent = percent_greater(f_current) 
            
            ## Now, apply these functions to the list of rasters (6 GCMs) for each species
            ## Also, the 'init' function initializes a raster object with values
            ## suit_ras_greater    = reduce(suit.list, thresh_above_fun, .init = suit.list[[1]] > thresh)
            
            ## Check the logic of removing the zeros.....................................................
            ## If we are just using the combined rasters without calculating the difference, do we need to worry about the zeros?   
            
            #########################################################################################################################
            ## First, calculate the cells which are greater that the: 
            ## Maximum training sensitivity plus specificity Logistic threshold
            message('Running thresholds for ', species, ' | 20', time_slice, ' combined suitability > ', thresh)
            
            suit_ras1_thresh   = thresh_greater(suit.list[[1]])  
            suit_ras2_thresh   = thresh_greater(suit.list[[2]])
            suit_ras3_thresh   = thresh_greater(suit.list[[3]])
            suit_ras4_thresh   = thresh_greater(suit.list[[4]])   ## Abbreviate this...
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
            
            ## Plot to check
            plot(current_suit_thresh, main = gsub('_', ' ', (sprintf('%s current max_train_sensit > %s', species, thresh))))
            plot(combo_suit_percent,  main = gsub('_', ' ', (sprintf('%s future 10th percentile > %s',   species, percent))))
            plot(combo_suit_thresh,   main = gsub('_', ' ', (sprintf('%s future max_train_sensit > %s',  species, thresh))))
            plot(combo_suit_4GCM,     main = gsub('_', ' ', (sprintf('%s 4+ GCMs > %s',  species, thresh))))
            
            #########################################################################################################################
            ## For each species, calculate the projected rainfall and temperature increase and decreas for each GCM? 
            ## No, calculate this once for each GCM.
            
            #########################################################################################################################
            ## Then using this GCM consensus, calculate whether the species is likely to be present in each SUA.
            ## Decide on a threshold of % area (10?) of the SUA that needs to be occupied, for each species to be considered present. 
            message('Running zonal stats for ', species, ' | 20', time_slice, ' combined suitability > ', thresh)
            
            ## Check the order of lists match, species, SUAs, areas need to match up ................................................
            
            ## First read in the shapefile - can we do this outside the function?
            ## Also sort the attribute table so that it matches the list 
            areal_unit = readOGR("F:/green_cities_sdm/data/base/CONTEXTUAL/IN_SUA_WGS.shp", layer = "IN_SUA_WGS")
            areal_unit = areal_unit[order(areal_unit$SUA_NAME11),] 
            
            ## Should we also create a simple count of the no. of cells per SUA with where > 4 GCMs met the suitability threshold?
            ## Need to fix this so that it has the same order as the shapefile
            #sp.count <- spatialEco::zonal.stats(x = areal_unit, y = combo_suit_4GCM, stat = sum, trace = TRUE, plot = TRUE) 
            # z.med   <- spatialEco::zonal.stats(x = areal_unit, y = combo_suit_4GCM, stat = median, trace = TRUE, plot = TRUE)
            
            #########################################################################################################################
            ## Then, extract the values of the presence raster for each areal unit: generates a list: 32 seconds
            ext  <- extract(combo_suit_4GCM, areal_unit, method = 'simple')
            
            ## A function to tabulate the raster values by aerial unit, returning a data frame
            tabFunc<-function(indx, extracted, region, regname) {
              
              dat<-as.data.frame(table(extracted[[indx]]))
              dat$name<-region[[regname]][[indx]]
              return(dat)
              
            }
            
            ## Run through each areal unit and calculate a table of the count of raster cells
            tabs <- lapply(seq(ext), tabFunc, ext, areal_unit, "SUA_NAME11")
            tabs <- do.call("rbind", tabs )
            
            ## Can we get the count here, so we don't need to run zonal stats as well
            # tabs.count      = tabs
            # tabs.count$Freq = ifelse(tabs.count$Var1 == 0, 0, tabs.count$Freq) 
            # sp.count <- aggregate(tabs.count$Freq, by = list(Category = tabs$name), FUN = sum)
            # names(sp.count) = c('SUA_NAME11', 'COUNT')
            # head(sp.count)
            
            #########################################################################################################################
            ## Now mutate the table
            PERECENT.AREA <- tabs %>%
              
              group_by(name) %>%                                          ## group by region
              mutate(totcells = sum(Freq),                                ## how many cells overall
                     percent.area = round(100 * Freq / totcells, 2)) %>%  ## cells /total cells
              
              dplyr::select(-c(Freq, totcells)) %>%                       ## there is a select func in raster so need to specify
              spread(key = Var1, value = percent.area, fill = 0)  %>%     ## make wide format
              as.data.frame()
            
            ## Rename and create a column for whether or not the species occupies that area 
            names(PERECENT.AREA) =  c('SUA_NAME11', 'Absent', 'Present') 
            PERECENT.AREA$species_present = ifelse(PERECENT.AREA$Present >= area_occ, 1, 0)
            head(PERECENT.AREA)
            
            ## Create a table with the columns: REGION, AREA SPECIES, STATS (for each time slice)
            ## Change this so that we can summarise across all species as suggested by Linda 
            GCM.AREA.SUMMARY <- data.frame(SUA        = areal_unit$SUA_NAME11, 
                                           AREA_SQKM  = areal_unit$AREA_SQKM,
                                           SPECIES    = species,
                                           #CELL_COUNT = sp.count$COUNT,
                                           PERCENT_AREA = PERECENT.AREA$Present,
                                           PRESENT      = PERECENT.AREA$species_present)
            
            ## Rename columns using sprintf, so we can include the suitability threshold and the time slice
            names(GCM.AREA.SUMMARY) <-  c('SUA',
                                          'AREA_SQKM',
                                          'SPECIES',
                                          #'CELL_COUNT',
                                          sprintf("Percent area where 4GCMs > %s in 20%s",   thresh, time_slice),
                                          sprintf('Species present 4GCMs > %s in 20%s',   thresh, time_slice))
            View(GCM.AREA.SUMMARY)
            
            ##
            dim(GCM.AREA.SUMMARY)
            unique(GCM.AREA.SUMMARY$SPECIES)
            
            #########################################################################################################################
            ## Then save the table of SUA results for all species to a datafile...
            ## This would be the file to loop over to create a summary of species per SUA
            write.csv(GCM.AREA.SUMMARY, sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_%s.csv',
                                                species, species, "SUA_summary"), row.names = FALSE)
            
            #########################################################################################################################
            #########################################################################################################################
            ## Write the rasters for each species/threshold
            message('Writing ', species, ' current', ' max train > ', thresh) 
            writeRaster(current_suit_thresh, sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_%s%s.tif',
                                                     species, species, "current_suit_above_", thresh), overwrite = TRUE) 
            
            message('Writing ', species, ' | 20', time_slice, ' max train > ', thresh) 
            writeRaster(combo_suit_thresh, sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_20%s%s%s.tif',
                                                   species, species, time_slice, "_Max_train_sensit_above_", thresh), overwrite = TRUE)
            
            ## Write the percentile raster
            message('Writing ', species, ' | 20', time_slice, ' 10th percentile > ', percent) 
            writeRaster(combo_suit_percent, sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_20%s%s%s.tif',
                                                    species, species, time_slice, "_10_percentile_omiss_above_", percent), overwrite = TRUE)
            
            ########################################################################################################################
            ## Now create the empty panel just before plotting, and read in the occurrence data
            empty <- init(combo_suit_thresh, function(x) NA)
            occ   <- readRDS(sprintf('./output/maxent/STD_VAR_ALL/%s/occ.rds', species)) 
            
            ## Use the 'levelplot' function to make a multipanel output: occurences, percentiles and thresholds
            message('Writing figure for ', species, ' | 20', time_slice, ' > ', thresh) 
            
            png(sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_20%s_suitability_above_%s.png',
                        species, species, time_slice, thresh),
                11, 4, units = 'in', res = 300)
            
            ## Need an empty frame
            print(levelplot(stack(empty,
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
                    ## Why are the points not working using "occ" in this way?
                    layer(sp.polygons(aus)) +
                    layer(sp.points(occ, pch = 20, cex = 0.4, 
                                    col = c('red', 'transparent', 'transparent')[panel.number()]), data = list(occ = occ)))
            
            ## Finish the device
            dev.off()
            
          } else {
            
            message(species, ' 20', time_slice, ' combined suitability > ', thresh, ' skipped - already exists')   ## 
            
          }
          
        }
        
      }
      
    })
    
  })
  
}
    




#########################################################################################################################
## TABLE FUNCTIONS
#########################################################################################################################


## Could turn this into a function, and loop over a list of subfolders...
bind_maxent_tables = function(x) {
  
  ## First check the table exists
  maxent.tables = list.files("./output/maxent/STD_VAR_ALL/")
  path       = "./output/maxent/STD_VAR_ALL/"
  
  f <- paste0(path, x, "/full/maxentResults.csv")
  if(file.exists(f)){
    
    ## Then pipe the table list into lapply
    MAXENT.STD.VAR.SUMMARY <- table_list[c(1:length(table_list))] %>%           ## Change to 1:20 if SDMs not complete 
      
      ## pipe the list into lapply
      lapply(function(x) {
        
        ## Create the character string
        f <- paste0(path, x, "/full/maxentResults.csv")
        
        ## load each .RData file
        d <- read.csv(f)
        
        ## now add a model column
        d = cbind(GBIF_Taxon = x, Model_run  = path, d) 
        dim(d)
        
        ## Remove path gunk, and species
        #d$GBIF_Taxon = gsub("_", " ", d$GBIF_Taxon)
        d$Model_run  = gsub("./output/maxent/", "", d$Model_run)  ## Model_run needed, test effect of different settings
        d$Model_run  = gsub("/", "", d$Model_run)
        d$Species    = NULL
        d
        
      }) %>%
      
      ## Finally, bind all the rows together
      bind_rows
    
  } else {
    
    message(x, ' table does not exist')   ## Skip species with no existing SDM
    
  }
  
} 





#########################################################################################################################
## ANOMALY FUNCTIONS
#########################################################################################################################
    

#########################################################################################################################
## For each scenario, calcualte the mean annual temperature and annual rainfall anomaly
calculate.anomaly.2050 = function(scen_2050) {
  
  ##  First, run a loop over each scenario: options(error = recover)    
  lapply(scen_2050, function(x) {
    
    ## Assign the scenario name (to use later in the plot)
    scen_name = gcms.50$GCM[gcms.50$id == x]
    
    ## create the current raster
    aus <- ne_states(country = 'Australia') %>% 
      subset(!grepl('Island', name))
    
    env.grids.current <- stack(
      file.path('./data/base/worldclim/aus/0.5/bio/current',
                sprintf('bio_%02d.tif', 1:19)))
    
    #########################################################################################################################
    ## Create current rasters 
    env.grids.current[[1]]  = env.grids.current[[1]]/10
    current.bio1            = env.grids.current[[1]]
    current.bio12           = env.grids.current[[12]]
    
    ## Create a raster stack for each 2050 GCM - also an empty raster for the final plot
    s <- stack(
      sprintf('./data/base/worldclim/aus/0.5/bio/2050/%s/%s%s.tif',
              x, x, 1:19))

    ########################################################################################################################
    ## Create future rasters
    s[[1]]  = s[[1]]/10
    future.bio1  = s[[1]]
    future.bio12 = s[[12]]
    
    ## Now subtract the current from the future...
    temp.anomaly = future.bio1  - current.bio1
    rain.anomaly = future.bio12 - current.bio12
    
    plot(temp.anomaly, main = paste0("BIO1  anomaly ", x))
    plot(rain.anomaly, main = paste0("BIO12 anomaly ", x))
    
    ########################################################################################################################
    ## Write the rasters for each species/threshold
    message('Writing BIO1/12 anomaly rasters for 2050 ', x) 
    writeRaster(temp.anomaly, sprintf('./data/base/worldclim/aus/0.5/bio/anomalies/%s_%s.tif',
                                      "BIO1_anomaly", x), overwrite = TRUE) 
    
    writeRaster(rain.anomaly, sprintf('./data/base/worldclim/aus/0.5/bio/anomalies/%s_%s.tif',
                                      "BIO12_anomaly", x), overwrite = TRUE)
    
  })
  
} 




#########################################################################################################################
## For each scenario, calcualte the mean annual temperature and annual rainfall anomaly
calculate.anomaly.2070 = function(scen_2070) {
  
  ##  First, run a loop over each scenario: options(error = recover)    
  lapply(scen_2070, function(x) {
    
    ## Assign the scenario name (to use later in the plot)
    scen_name = gcms.50$GCM[gcms.50$id == x]
    
    ## create the current raster
    aus <- ne_states(country = 'Australia') %>% 
      subset(!grepl('Island', name))
    
    env.grids.current <- stack(
      file.path('./data/base/worldclim/aus/0.5/bio/current',
                sprintf('bio_%02d.tif', 1:19)))
    
    #########################################################################################################################
    ## Create current rasters 
    env.grids.current[[1]]  = env.grids.current[[1]]/10
    current.bio1            = env.grids.current[[1]]
    current.bio12           = env.grids.current[[12]]
    
    ## Create a raster stack for each 2050 GCM - also an empty raster for the final plot
    s <- stack(
      sprintf('./data/base/worldclim/aus/0.5/bio/2070/%s/%s%s.tif',
              x, x, 1:19))
    
    ########################################################################################################################
    ## Create future rasters
    s[[1]]  = s[[1]]/10
    future.bio1  = s[[1]]
    future.bio12 = s[[12]]
    
    ## Now subtract the current from the future...
    temp.anomaly = future.bio1  - current.bio1
    rain.anomaly = future.bio12 - current.bio12
    
    plot(temp.anomaly, main = paste0("BIO1  anomaly ", x))
    plot(rain.anomaly, main = paste0("BIO12 anomaly ", x))
    
    ########################################################################################################################
    ## Write the rasters for each species/threshold
    message('Writing BIO1/12 anomaly rasters for 2070 ', x) 
    writeRaster(temp.anomaly, sprintf('./data/base/worldclim/aus/0.5/bio/anomalies/%s_%s.tif',
                                      "BIO1_anomaly", x), overwrite = TRUE) 
    
    writeRaster(rain.anomaly, sprintf('./data/base/worldclim/aus/0.5/bio/anomalies/%s_%s.tif',
                                      "BIO12_anomaly", x), overwrite = TRUE)
    
  })
  
}  





#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################