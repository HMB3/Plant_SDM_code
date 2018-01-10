#########################################################################################################################
############################################# FUNCTIONS FOR MAPPING SDMS ################################################ 
#########################################################################################################################


#########################################################################################################################
## PROJECT MAXENT MODELS
#########################################################################################################################


#########################################################################################################################
## 2050
project.grids.2050 = function(scen_2050, test_spp){
  
  ##  First, run a loop over each scenario: options(error = recover)    
  lapply(scen_2050, function(x) {
    
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
    lapply(test_spp, function(species) {
      
      ## First check if the species projection has already been run...
      if(!file.exists(sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_%s.tif',
                              species, species, x))) {
        message('Doing ', species) 
        
        ## Read in the SDM model calibrated on current conditions
        m <- readRDS(sprintf('./output/maxent/STD_VAR_ALL/%s/maxent_fitted.rds', species)) 
        
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
project.grids.2070 = function(scen_2070, test_spp){
  
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
        m <- readRDS(sprintf('./output/maxent/STD_VAR_ALL/%s/maxent_fitted.rds', species)) 
        
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
## 2050 combine
ensemble.2050 = function(SDM.RESULTS.DIR, test_spp){
  
  lapply(SDM.RESULTS.DIR, function(DIR) { 
    
    ## And each species - although we don't want all possible combinations. How can this loop be improved?
    lapply(test_spp, function(species) {
      
      ###################################################################################################################
      ## Create a list of the rasters in each directory, then take the mean. How long does the mean calculation take?
      ## Tidy this up with %, etc
      raster.list = list.files(as.character(DIR), pattern = "bi50.tif", full.names = TRUE)
      suit        = stack(raster.list)
      suit.list   = unstack(suit)
      mean.suit   = mean(suit)   ## plot(mean.suit)
      
      ## Write the mean to file
      f_mean = sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_2050_suitability_mean.tif', 
                       species, species)
      
      
      if(!file.exists(f_mean)) {
        
        writeRaster(mean.suit, sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_2050_suitability_mean.tif', 
                                       species, species), overwrite = TRUE)
        
      } else {
        
        message(species, ' 2050 mean suitability skipped - already exists')   ## 
        
      }
      
      ###################################################################################################################
      ## Then create rasters that meet habitat suitability criteria thresholds
      for (thresh in c(0.7)) {  ## for (thresh in c(0.5, 0.7, 0.8)) {
        
        ## Check if the combined suitability raster exists
        f_suit <- sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_%s%s.tif',
                          species, species, "2050_suitability_consensus_greater_", thresh)
        
        
        ## If it doesn't exist, create the suitability raster
        if(!file.exists(f_suit)) {
          
          ## First create a simple function to threshold each of the rasters in raster.list
          ## So how does this change to add them up? First
          scen_greater = function (x) {x > thresh}
          
          ## Then apply the function to the GCM list for each species
          ## The Init function initializes a raster object with values
          #suit_ras_greater    = reduce(suit.list, thresh_greater_fun, .init = suit.list[[1]] > thresh)
          
          suit_ras1_greater  = scen_greater(suit.list[[1]])   ## do this better...
          suit_ras2_greater  = scen_greater(suit.list[[2]])
          suit_ras3_greater  = scen_greater(suit.list[[3]])
          suit_ras4_greater  = scen_greater(suit.list[[4]])
          suit_ras5_greater  = scen_greater(suit.list[[5]])
          suit_ras6_greater  = scen_greater(suit.list[[6]])
          
          ## Then sum them up: could use a function or magrittr to compress this..
          suit_ras_consensus_sum   =  Reduce("+", list(suit_ras1_greater, suit_ras2_greater, suit_ras3_greater,
                                                       suit_ras4_greater, suit_ras5_greater, suit_ras6_greater))
          ## plot(suit_ras_consensus_sum )
          
          ## Write the raster for each species and threshold inside the loop. But how to access the rasters for plotting?
          message('Writing ', species, ' 2050 suitability > ', thresh) 
          writeRaster(suit_ras_consensus_sum, sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_%s%s.tif',
                                                      species, species, "2050_suitability_consensus_greater_", thresh), overwrite = TRUE)
          
          
        } else {
          
          message(species, ' 2050 suitability consensus > ', thresh, ' skipped - already exists')   ## 
          
        }
        
      }
      
      ## Now create the empty panel just before plotting, and read in the occurrence data
      empty <- init(suit_ras_consensus_sum, function(x) NA)
      occ <- readRDS(sprintf('./output/maxent/STD_VAR_ALL/%s/occ.rds', species)) 
      
      ########################################################################################################################
      ## Use the levelplot function to make a multipanel output: average, threshold 1, threshold 2
      png(sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_suitability_above_%s.tif',
                  species, species, thresh),
          11, 4, units = 'in', res = 300)
      
      ## Need an empty frame
      print(levelplot(stack(empty,
                            mean.suit,
                            suit_ras_consensus_sum), margin = FALSE,
                      
                      ## Create a colour scheme using colbrewer: 100 is to make it continuos
                      ## Also, make it a one-directional colour scheme
                      scales      = list(draw = FALSE), 
                      at = seq(0, 6, length = 100),
                      col.regions = colorRampPalette(rev(brewer.pal(11, 'Spectral'))),
                      
                      ## Give each plot a name
                      names.attr = c('Occurrence', 'GCM average', 'GCM > 0.7'),
                      colorkey   = list(height = 0.5, width = 3), xlab = '', ylab = '',
                      main       = list(gsub('_', ' ', species), font = 4, cex = 2)) +
              
              ## Plot the Aus shapefile with the occurrence points for reference
              ## Why are the points not working?
              layer(sp.polygons(aus)) +
              layer(sp.points(occ, pch = 20, cex = 0.4, 
                              col = c('red', 'transparent', 'transparent')[panel.number()]), data = list(occ = occ)))
      
      ## Finish the device
      dev.off()
      
    })
    
  })
  
}




#########################################################################################################################
## 2070 combine
ensemble.2070 = function(SDM.RESULTS.DIR, test_spp){
  
  lapply(SDM.RESULTS.DIR, function(DIR) { 
    
    ## And each species - although we don't want all possible combinations. How can this loop be improved?
    lapply(test_spp, function(species) {
      
      ###################################################################################################################
      ## Create a list of the rasters in each directory, then take the mean. How long does the mean calculation take?
      ## Tidy this up with %, etc
      raster.list = list.files(as.character(DIR), pattern = "bi70.tif", full.names = TRUE)
      suit        = stack(raster.list)
      suit.list   = unstack(suit)
      mean.suit   = mean(suit)   ## plot(mean.suit)
      
      ## Write the mean to file
      f_mean = sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_2070_suitability_mean.tif', 
                       species, species)
      
      
      if(!file.exists(f_mean)) {
        
        writeRaster(mean.suit, sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_2070_suitability_mean.tif', 
                                       species, species), overwrite = TRUE)
        
      } else {
        
        message(species, ' 2070 mean suitability skipped - already exists')   ## 
        
      }
      
      ###################################################################################################################
      ## Then create rasters that meet habitat suitability criteria thresholds
      for (thresh in c(0.7)) {  ## for (thresh in c(0.5, 0.7, 0.8)) {
        
        ## Check if the combined suitability raster exists
        f_suit <- sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_%s%s.tif',
                          species, species, "2070_suitability_consensus_greater_", thresh)
        
        
        ## If it doesn't exist, create the suitability raster
        if(!file.exists(f_suit)) {
          
          ## First create a simple function to threshold each of the rasters in raster.list
          ## So how does this change to add them up? First
          scen_greater = function (x) {x > thresh}
          
          ## Then apply the function to the GCM list for each species
          ## The Init function initializes a raster object with values
          #suit_ras_greater    = reduce(suit.list, thresh_greater_fun, .init = suit.list[[1]] > thresh)
          
          suit_ras1_greater  = scen_greater(suit.list[[1]])   ## do this better...
          suit_ras2_greater  = scen_greater(suit.list[[2]])
          suit_ras3_greater  = scen_greater(suit.list[[3]])
          suit_ras4_greater  = scen_greater(suit.list[[4]])
          suit_ras5_greater  = scen_greater(suit.list[[5]])
          suit_ras6_greater  = scen_greater(suit.list[[6]])
          
          ## Then sum them up: could use a function or magrittr to compress this..
          suit_ras_consensus_sum   =  Reduce("+", list(suit_ras1_greater, suit_ras2_greater, suit_ras3_greater,
                                                       suit_ras4_greater, suit_ras5_greater, suit_ras6_greater))
          ## plot(suit_ras_consensus_sum )
          
          ## Write the raster for each species and threshold inside the loop. But how to access the rasters for plotting?
          message('Writing ', species, ' 2070 suitability > ', thresh) 
          writeRaster(suit_ras_consensus_sum, sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_%s%s.tif',
                                                      species, species, "2070_suitability_consensus_greater_", thresh), overwrite = TRUE)
          
          
        } else {
          
          message(species, ' 2070 suitability consensus > ', thresh, ' skipped - already exists')   ## 
          
        }
        
      }
      
      ## Now create the empty panel just before plotting, and read in the occurrence data
      empty <- init(suit_ras_consensus_sum, function(x) NA)
      occ <- readRDS(sprintf('./output/maxent/STD_VAR_ALL/%s/occ.rds', species)) 
      
      ########################################################################################################################
      ## Use the levelplot function to make a multipanel output: average, threshold 1, threshold 2
      png(sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_suitability_above_%s.tif',
                  species, species, thresh),
          11, 4, units = 'in', res = 300)
      
      ## Need an empty frame
      print(levelplot(stack(empty,
                            mean.suit,
                            suit_ras_consensus_sum), margin = FALSE,
                      
                      ## Create a colour scheme using colbrewer: 100 is to make it continuos
                      ## Also, make it a one-directional colour scheme
                      scales      = list(draw = FALSE), 
                      at = seq(0, 6, length = 100),
                      col.regions = colorRampPalette(rev(brewer.pal(11, 'Spectral'))),
                      
                      ## Give each plot a name
                      names.attr = c('Occurrence', 'GCM average', 'GCM > 0.7'),
                      colorkey   = list(height = 0.5, width = 3), xlab = '', ylab = '',
                      main       = list(gsub('_', ' ', species), font = 4, cex = 2)) +
              
              ## Plot the Aus shapefile with the occurrence points for reference
              ## Why are the points not working?
              layer(sp.polygons(aus)) +
              layer(sp.points(occ, pch = 20, cex = 0.4, 
                              col = c('red', 'transparent', 'transparent')[panel.number()]), data = list(occ = occ)))
      
      ## Finish the device
      dev.off()
      
    })
    
  })
  
}
    

    
  

#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################