#########################################################################################################################
############################################# FUNCTIONS FOR MAPPING SDMS ################################################ 
#########################################################################################################################


#########################################################################################################################
## PROJECT MAXENT MODELS
#########################################################################################################################


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
      if(!file.exists(sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_%s.tif',
                              species, species, x))) {
        message('Doing ', species) 
        
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
          
          ## Report which prediction is in progress m$me_full
          message('Running current prediction for ', species) 
          
          pred.current <- rmaxent::project(
            m, env.grids.current[[colnames(m$me_full@presence)]])$prediction_logistic
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
            m$me_full, s[[colnames(m$me_full@presence)]])$prediction_logistic
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
        
        message(species, ' ', x, ' skipped - prediction already run')   ## Ignore species which have already been run - maybe remove
        
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
      if(!file.exists(sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_%s.tif',
                              species, species, x))) {
        message('Doing ', species) 
        
        ## Read in the SDM model calibrated on current conditions
        #m <- readRDS(sprintf('./output/maxent/STD_VAR_ALL/%s/maxent_fitted.rds', species))
        m <- readRDS(sprintf('./output/maxent/STD_VAR_ALL/%s/full/model.rds', species)) 
        
        ## Read in the occurrence points used to create the SDM
        occ <- readRDS(
          sprintf('./output/maxent/STD_VAR_ALL/%s/occ.rds', 
                  species)) %>%
          
        ########################################################################################################################
        ## If the current raster doesn't exist, create it
        spTransform(CRS('+init=epsg:4326'))
        f_current <- sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_current.tif', 
                             species, species)
        
        if(!file.exists(f_current)) {
          
          ## Report which prediction is in progress
          message('Running current prediction for ', species) 
          
          pred.current <- rmaxent::project(
            m$me_full, env.grids.current[[colnames(m$me_full@presence)]])$prediction_logistic
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
            m$me_full, s[[colnames(m$me_full@presence)]])$prediction_logistic
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
  
}
 




#########################################################################################################################
## COMBINE GCM FUNCTIONS
#########################################################################################################################


#########################################################################################################################
## Loop over directories, species and one threshold for each, also taking a time_slice argument. Next, make the lists generic too
combine_gcm_threshold = function(DIR_list, species_list, thresholds, percentiles, time_slice, area_occ) {
  
  lapply(DIR_list, function(DIR) { 
    
    ## And each species - although we don't want all possible combinations. How can this loop be improved?
    lapply(species_list, function(species) {
      
      # for (slice in time_slice) {
        
        ###################################################################################################################
        ## Create a list of the rasters in each directory, then take the mean. Tidy this up with %, etc
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
            f_suit <- sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_20%s%s%s.tif',
                              species, species, time_slice, "_combined_suitability_above_", thresh)
            
            ## If it doesn't exist, create the suitability raster
            if(!file.exists(f_suit)) {
              
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
              
              ## Then greater than the 10 percentile training presence training omission
              suit_ras1_percent  = percent_greater(suit.list[[1]])
              suit_ras2_percent  = percent_greater(suit.list[[2]])
              suit_ras3_percent  = percent_greater(suit.list[[3]])
              suit_ras4_percent  = percent_greater(suit.list[[4]])
              suit_ras5_percent  = percent_greater(suit.list[[5]])
              suit_ras6_percent  = percent_greater(suit.list[[6]])
              
              #########################################################################################################################
              ## Then sum them up: could use % to compress this..
              ## All the threshholds
              combo_suit_thresh   =  Reduce("+", list(suit_ras1_thresh, suit_ras2_thresh, suit_ras3_thresh,
                                                      suit_ras4_thresh, suit_ras5_thresh, suit_ras6_thresh))
              
              ## All the percentiles
              combo_suit_percent  =  Reduce("+", list(suit_ras1_percent, suit_ras2_percent, suit_ras3_percent,
                                                      suit_ras4_percent, suit_ras5_percent, suit_ras6_percent))
              
              #########################################################################################################################
              #########################################################################################################################
              ## Next, calcualte the loss or gain between the two time periods. Create a binary raster for GCM layer
              message('Calculating change for ', species, ' | 20', time_slice, ' combined suitability > ', thresh)
              
              ## Functions for thresholding rasters
              binary_ras   <- function(x) {ifelse(x >=  1, 1, 0) }
              binary_4_ras <- function(x) {ifelse(x >=  4, 1, 0) }
              band_ras     <- function(x) {ifelse(x >=  4, 4, ifelse(x > 0 & x < 4, 3, x)) }
              
              ## Function to tabulate raster values by aerial unit (e.g. SUA) and return a data.frame
              tabFunc <- function(indx, extracted, region, regname) {
                dat<-as.data.frame(table(extracted[[indx]]))
                dat$name<-region[[regname]][[indx]]
                return(dat)
              }
              
              combo_suit_band  <- calc(combo_suit_thresh, fun = band_ras)
              combo_suit4_band <- calc(combo_suit_thresh, fun = binary_4_ras)
              #combo_suit_binary <- calc(combo_suit_thresh, fun = binary)
              
              ## Then subtract the binary current layer from the binary future layer. So need to mask the no-data from this overlay calc.  
              # binary_future_minus_current     = overlay(combo_suit_binary,
              #                                           current_suit_thresh,
              #                                           fun = function(r1, r2) {return (r1 - r2)})
              
              #########################################################################################################################
              ## Subtract the binary current layer from the discrete future layer 
              discrete_future_minus_current = overlay(combo_suit_band,
                                                      current_suit_thresh,
                                                      fun = function(r1, r2) {return (r1 - r2)})                                           
              
              #########################################################################################################################
              ## Subtract the binary current layer from the integer future layer 
              integer_future_minus_current  = overlay(combo_suit_thresh,
                                                      current_suit_thresh,
                                                      fun = function(r1, r2) {return (r1 - r2)})
              
              ## Plot the difference between future and current layers...
              plot(current_suit_thresh, main = gsub('_', ' ', (sprintf('%s current Max_train_sensit > %s', species, thresh))))
              plot(combo_suit_thresh,   main = gsub('_', ' ', (sprintf('%s future Max_train_sensit > %s',  species, thresh))))
              
              ## Plot the integer difference
              plot(integer_future_minus_current,
                   main = gsub('_', ' ', (sprintf('%s future - current  Max_train_sensit > %s', species, thresh))))
              
              ## Plot the band difference
              plot(discrete_future_minus_current,
                   main = gsub('_', ' ', (sprintf('%s future - discrete  Max_train_sensit > %s', species, thresh))))
              
              #########################################################################################################################
              ## Now, mask out the zero values?
              # plus    = overlay(combo_suit_thresh,
              #                   current_suit_thresh,
              #                   fun = function(r1, r2) {return (r1 + r2)})
              # 
              # mask <- calc(plus, fun = rc)
              # mask[mask < 0] <- NA
              # 
              # pol <- rasterToPolygons(mask, fun = function(x){x>0})
              # 
              # 
              # 
              # test   =  polygonizer(plus, outshape = NULL, pypath = "F:/green_cities_sdm/R")
              # crs(b) <- crs(r)
              # crop   <- crop(r, b)
              
              
              #########################################################################################################################
              ## Then calculate the loss or gain within a given areal unit. Use the SUA's, but could be anything!
              ## How can this be loaded outside the function?
              message('Running zonal stats for ', species, ' | 20', time_slice, ' combined suitability > ', thresh)
              areal_unit = readOGR("F:/green_cities_sdm/data/base/CONTEXTUAL/IN_SUA_WGS.shp", layer = "IN_SUA_WGS")
              areal_unit = areal_unit[order(areal_unit$SUA_NAME11),] 
              
              ## For each species, create a binary raster with cells > 4 GCMs above the maxent threshold = 1, and cells with < 4 GCMs = 0. 
              ## Decide on a threshold of % area (10?) of the SUA that needs to be occupied, for each species to be considered present. 
 
              ## First, create a simple count of the no. of cells per SUA with where > 4 GCMs met the suitability threshold
              ## Need to fix this so that it has the same order as 
              #sp.count <- spatialEco::zonal.stats(x = areal_unit, y = combo_suit4_band, stat = sum, trace = TRUE, plot = TRUE) 
              # z.max  <- spatialEco::zonal.stats(x = areal_unit, y = intergter_future_minus_current, stat = max,    trace = TRUE, plot = TRUE)
              # z.min  <- spatialEco::zonal.stats(x = areal_unit, y = intergter_future_minus_current, stat = min,    trace = TRUE, plot = TRUE)
              # z.med  <- spatialEco::zonal.stats(x = areal_unit, y = intergter_future_minus_current, stat = median, trace = TRUE, plot = TRUE)
              
              ## Then, extract the values of the presence raster for each areal unit: generates a list: 32 seconds
              ext  <- extract(combo_suit4_band, areal_unit, method = 'simple')
              
              ## Run through each areal unit and calculate a table of the count of raster cells
              tabs <- lapply(seq(ext), tabFunc, ext, areal_unit, "SUA_NAME11")
              tabs <- do.call("rbind", tabs )
              
              ## Get the count here so we don't need to run zonal stats as well
              # tabs.count      = tabs
              # tabs.count$Freq = ifelse(tabs.count$Var1 == 0, 0, tabs.count$Freq) 
              # sp.count <- aggregate(tabs.count$Freq, by = list(Category = tabs$name), FUN = sum)
              # names(sp.count) = c('SUA_NAME11', 'COUNT')
              # head(sp.count)
              
              ## Now mutate the table
              PERECENT.AREA <- tabs %>%
                
                group_by(name) %>%                                          ## group by region
                mutate(totcells = sum(Freq),                                ## how many cells overall
                       percent.area = round(100 * Freq / totcells, 2)) %>%  ## cells by landuse/total cells
   
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
                                             PRESENT    = PERECENT.AREA$species_present)
  
              ## Rename columns using sprintf, so we can include the suitability threshold and the time slice
              names(GCM.AREA.SUMMARY) <-  c('SUA',
                                            'AREA_SQKM',
                                            'SPECIES',
                                            #'CELL_COUNT',
                                            sprintf('PRESENT_4_GCMs > %s in 20%s',   thresh, time_slice))
              View(GCM.AREA.SUMMARY)
 
              ##
              dim(GCM.AREA.SUMMARY)
              unique(GCM.AREA.SUMMARY$SPECIES)
              
              ## Then save the table of SUA results for all species to a datafile... 
              write.csv(GCM.AREA.SUMMARY, sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_%s%s.csv',
                                                  species, species, "Max_train_sensit_above_", thresh), row.names = FALSE)

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
              
              #########################################################################################################################
              ## Write the raster for discrete loss/gain
              message('Writing ', species, ' 20', time_slice, ' loss/gain max train > ', thresh) 
              writeRaster(combo_suit4_band, sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_20%s%s%s.tif',
                                                               species, species, time_slice, "_combo_suit4_band_", thresh), overwrite = TRUE)
              
              ## Write the rasters for binary loss/gain for each species/threshold
              message('Writing ', species, ' 20', time_slice, ' loss/gain max train > ', thresh) 
              writeRaster(binary_future_minus_current, sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_20%s%s%s.tif',
                                                               species, species, time_slice, "_binary_future_minus_current_", thresh), overwrite = TRUE)
              
              ## Write the rasters for interger loss/gain for each species/threshold
              message('Writing ', species, ' 20', time_slice, ' loss/gain 10th percentile > ', percent) 
              writeRaster(integer_future_minus_current, sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_20%s%s%s.tif',
                                                                 species, species, time_slice, "_integer_future_minus_current_", thresh), overwrite = TRUE)
              
              ########################################################################################################################
              ## Now create the empty panel just before plotting, and read in the occurrence dat
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
## Take a period argument. Next, make the lists generic too
# combine_maxent_predictions = function(SDM.RESULTS.DIR, test_spp, period){
#   
#   lapply(SDM.RESULTS.DIR, function(DIR) { 
#     
#     ## And each species - although we don't want all possible combinations. How can this loop be improved?
#     lapply(test_spp, function(species) {
#       
#       ###################################################################################################################
#       ## Create a list of the rasters in each directory, then take the mean. How long does the mean calculation take?
#       ## Tidy this up with %, etc
#       raster.list = list.files(as.character(DIR), pattern = sprintf('bi%s.tif', period), full.names = TRUE)  ##  sprintf('bi%s.tif', period)
#       suit        = stack(raster.list)
#       suit.list   = unstack(suit)
#       mean.suit   = mean(suit)   ## plot(mean.suit)
#       
#       ## Write the mean to file
#       f_mean = sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_20%s_suitability_mean.tif', 
#                        species, species, period)
#       
#       
#       if(!file.exists(f_mean)) {
#         
#         writeRaster(mean.suit, sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_20%s_suitability_mean.tif', 
#                                        species, species, period), overwrite = TRUE)
#         
#       } else {
#         
#         message(species, ' 20', period, ' mean suitability skipped - already exists')   ## 
#         
#       }
#       
#       ###################################################################################################################
#       ## Then create rasters that meet habitat suitability criteria thresholds
#       for (thresh in c(0.7)) {  ## for (thresh in c(0.5, 0.7, 0.8)) {
#         
#         ## Check if the combined suitability raster exists
#         f_suit <- sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_20%s%s%s.tif',
#                           species, species, period, "_combined_suitability_above_", thresh)
#         
#         
#         ## If it doesn't exist, create the suitability raster
#         if(!file.exists(f_suit)) {
#           
#           ## First create a simple function to threshold each of the rasters in raster.list
#           ## So how does this change to add them up? First
#           scen_greater = function (x) {x > thresh}
#           
#           ## Then apply the function to the GCM list for each species
#           ## The Init function initializes a raster object with values
#           #suit_ras_greater    = reduce(suit.list, thresh_above_fun, .init = suit.list[[1]] > thresh)
#           
#           suit_ras1_greater  = scen_greater(suit.list[[1]])   ## do this better...
#           suit_ras2_greater  = scen_greater(suit.list[[2]])
#           suit_ras3_greater  = scen_greater(suit.list[[3]])
#           suit_ras4_greater  = scen_greater(suit.list[[4]])
#           suit_ras5_greater  = scen_greater(suit.list[[5]])
#           suit_ras6_greater  = scen_greater(suit.list[[6]])
#           
#           ## Then sum them up: could use a function or magrittr to compress this..
#           suit_ras_combined_sum   =  Reduce("+", list(suit_ras1_greater, suit_ras2_greater, suit_ras3_greater,
#                                                        suit_ras4_greater, suit_ras5_greater, suit_ras6_greater))
#           ## plot(suit_ras_combined_sum )
#           
#           ## Write the raster for each species and threshold inside the loop. But how to access the rasters for plotting?
#           message('Writing ', species, ' 20', period, ' combined suitability > ', thresh) 
#           writeRaster(suit_ras_combined_sum, sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_20%s%s%s.tif',
#                                                       species, species, period, "_combined_suitability_above_", thresh), overwrite = TRUE)
#           
#           
#         } else {
#           
#           message(species, ' 20', period, ' combined suitability > ', thresh, ' skipped - already exists')   ## 
#           
#         }
#         
#       }
#       
#       ## Now create the empty panel just before plotting, and read in the occurrence data
#       f_current <- sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_current.tif', 
#                            species, species)
#       
#       empty <- init(f_current, function(x) NA)
#       occ <- readRDS(sprintf('./output/maxent/STD_VAR_ALL/%s/occ.rds', species)) 
#       
#       ########################################################################################################################
#       ## Use the levelplot function to make a multipanel output: average, threshold 1, threshold 2
#       png(sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_20%s_suitability_above_%s.png',
#                   species, species, period, thresh),
#           11, 4, units = 'in', res = 300)
#       
#       ## Need an empty frame
#       print(levelplot(stack(empty,
#                             mean.suit,
#                             suit_ras_combined_sum), margin = FALSE,
#                       
#                       ## Create a colour scheme using colbrewer: 100 is to make it continuos
#                       ## Also, make it a one-directional colour scheme
#                       scales      = list(draw = FALSE), 
#                       at = seq(0, 6, length = 100),
#                       col.regions = colorRampPalette(rev(brewer.pal(11, 'Spectral'))),
#                       
#                       ## Give each plot a name
#                       names.attr = c('Aus occurrences', sprintf('20%s GCM average', period), sprintf('20%s GCM > %s', period, thresh)),
#                       colorkey   = list(height = 0.5, width = 3), xlab = '', ylab = '',
#                       main       = list(gsub('_', ' ', species), font = 4, cex = 2)) +
#               
#               ## Plot the Aus shapefile with the occurrence points for reference
#               ## Why are the points not working using "occ" in this way?
#               layer(sp.polygons(aus)) +
#               layer(sp.points(occ, pch = 20, cex = 0.4, 
#                               col = c('red', 'transparent', 'transparent')[panel.number()]), data = list(occ = occ)))
#       
#       ## Finish the device
#       dev.off()
#       
#     })
#     
#   })
#   
# }
    

    
  

#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################