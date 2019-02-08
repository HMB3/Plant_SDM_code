## One function for all time periods
## Maybe best to 
project_maxent_grids_mess = function(shp, scen_list, species_list, threshold_list, maxent_path, 
                                     climate_path, grid_names, time_slice, current_grids) {
  
  ## Read in Australia
  aus = shp %>%
    spTransform(ALB.CONICAL)
  
  ## First, run a loop over each scenario:    
  lapply(scen_list, function(x) {
    
    ## Create a raster stack for each of the 6 GCMs, not for each species
    s <- stack(sprintf('%s/20%s/%s/%s%s.tif', climate_path, time_slice, x, x, 1:19))
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
    mapply(function(species, threshold) {
      
      ## First, check if the maxent model exists
      ## Can we skip the species before dividing the rasters?................................................................
      save_name = gsub(' ', '_', species)
      if(file.exists(sprintf('%s/%s/full/maxent_fitted.rds', maxent_path, species))) {
        message('Doing ', species)
        
        ## Then, check if the species projection has already been run...
        if(!file.exists(sprintf('%s/%s/full/%s_%s.tif', maxent_path, species, species, x))) {
          
          ## Assign the scenario name (to use later in the plot)
          scen_name = eval(parse(text = sprintf('gcms.%s$GCM[gcms.%s$id == x]', time_slice, time_slice)))           
          
          ########################################################################
          ## Now read in the SDM model, calibrated on current conditions
          m   <- readRDS(sprintf('%s/%s/full/maxent_fitted.rds', maxent_path, species))
          m   <- m$me_full  ## 
          swd <- readRDS(sprintf('%s%s/swd.rds',                maxent_path, species, species))
          occ <- readRDS(sprintf('%s/%s/%s_occ.rds', maxent_path, species, save_name)) %>%
            spTransform(ALB.CONICAL)  
          
          ## Read in the current and future raster
          f_current  <- sprintf('%s/%s/full/%s_current.tif', maxent_path, species, species)
          hs_current <- raster(sprintf('%s%s/full/%s_%s%s.tif',    maxent_path,
                                       species, species, "current_suit_above_", threshold))
          
          if(!file.exists(f_current)) {
            
            ## Report which prediction is in progress :: m$me_full, m$me_full@presence
            message('Running current prediction for ', species) 
            
            pred.current <- rmaxent::project(
              m, current_grids[[colnames(m@presence)]])$prediction_logistic
            writeRaster(pred.current, f_current, overwrite = TRUE)
            
            
            #####################################################################
            ## Report current mess map in progress
            message('Running current mess map for ', species) 
            
            grid_names           = sdm.predictors
            current_grids        = aus.grids.current
            names(current_grids) = grid_names 
            current_grids        = subset(current_grids, intersect(names(current_grids), sdm.select))
            
            ## Create the current mess map :: check the variable names are the same
            mess_current  <- similarity(current_grids, swd[, names(current_grids)], full = TRUE)
            novel_current <- mess_current$similarity_min   < 0  ##   All novel environments are < 0
            novel_current[novel_current==0] <- NA               ##   0 values are NA
            
            ##################################################################
            ## Write out the current mess maps 
            writeRaster(mess_current$similarity_min, sprintf('%s%s/full/%s_%s.tif', 
                                                             maxent_path, species, species, "current_mess_map"), 
                        overwrite = TRUE, datatype = 'INT2S')
            
            ##################################################################
            ## Create a PNG file of all the CURRENT MESS output
            message('Creating mess maps for each environmental predictor for', species)
            mapply(function(r, name) {
              
              p <- levelplot(r, margin = FALSE, scales = list(draw = FALSE),
                             at = seq(minValue(r), maxValue(r), len = 100),
                             colokey = list(height = 0.6), main = gsub('_', ' ', sprintf('%s (%s)', name, species))) +
                
                latticeExtra::layer(sp.polygons(poly), data = list(poly = poly))   ## Use this in previous functions
              
              p <- diverge0(p, 'RdBu')
              f <- sprintf('%s%s/full/%s_messCurrent__%s.png', maxent_path, species, species, name)
              
              png(f, 8, 8, units = 'in', res = 300, type = 'cairo')
              print(p)
              dev.off()
              
              
            }, unstack(mess_current$similarity), names(mess_current$similarity))
            
            ## Write the raster of novel environments to the maxent directory 
            ## The "full" directory is getting full, could create a sub dir for MESS maps
            message('Writing maps of novel environments to file for', species) 
            writeRaster(novel_current, sprintf('%s%s/full/%s_%s.tif', 
                                               maxent_path, species, species, "current_novel_map"), 
                        overwrite = TRUE, datatype = 'INT2S')
            
            ##################################################################
            # mask out novel environments 
            # is.na(novel_current) is a binary layer showing 
            # not novel [=1] vs novel [=0], 
            # so multiplying this with hs_current will mask out novel
            hs_current_notNovel <- hs_current * is.na(novel_current) 
            
            ## Write out not-novel raster
            message('Writing maps of un- novel environments to file for', species) 
            
            writeRaster(hs_current_notNovel, sub('\\.tif', '_notNovel.tif', hs_current@file@name), 
                        overwrite = TRUE, datatype = 'INT2S')
            
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
            
            #####################################################################
            ## Report future mess map in progress
            message('Running future mess map for ', species) 
            
            grid_names          = sdm.predictors
            future_grids        = s
            names(future_grids) = grid_names 
            future_grids        = subset(future_grids, intersect(names(future_grids), sdm.select))
            
            ## Create the future mess map :: check the variable names are the same
            mess_future  <- similarity(future_grids, swd[, names(future_grids)], full = TRUE)
            novel_future <- mess_future$similarity_min   < 0  ##   All novel environments are < 0
            novel_future[novel_future==0] <- NA               ##   0 values are NA
            
            ##################################################################
            ## Write out the future mess maps 
            writeRaster(mess_future$similarity_min, sprintf('%s%s/full/%s_%s.tif', 
                                                             maxent_path, species, species, "future_mess_map"), 
                        overwrite = TRUE, datatype = 'INT2S')
            
            ##################################################################
            ## Create a PNG file of all the future MESS output: Check what r is here again?
            ## This should include the scenario name in the final output :: x
            message('Creating mess maps for each future environmental predictor for ', x, ' scenario for ', species)
            mapply(function(r, name) {
              
              p <- levelplot(r, margin = FALSE, scales = list(draw = FALSE),
                             at = seq(minValue(r), maxValue(r), len = 100),
                             colokey = list(height = 0.6), main = gsub('_', ' ', sprintf('%s (%s)', name, species))) +
                
                latticeExtra::layer(sp.polygons(poly), data = list(poly = poly))   ## Use this in previous functions
              
              p <- diverge0(p, 'RdBu')
              f <- sprintf('%s%s/full/%s_messFuture__%s.png', maxent_path, species, species, name)
              
              png(f, 8, 8, units = 'in', res = 300, type = 'cairo')
              print(p)
              dev.off()
              
              
            }, unstack(mess_future$similarity), names(mess_future$similarity))
            
            ## Write the raster of novel environments to the maxent directory 
            ## The "full" directory is getting full, could create a sub dir for MESS maps
            message('Writing maps of novel environments to file for', species) 
            writeRaster(novel_future, sprintf('%s%s/full/%s_%s.tif', 
                                               maxent_path, species, species, "future_novel_map"), 
                        overwrite = TRUE, datatype = 'INT2S')
            
            ##################################################################
            # mask out novel environments 
            # is.na(novel_future) is a binary layer showing 
            # not novel [=1] vs novel [=0], 
            # so multiplying this with hs_future will mask out novel
            hs_future_notNovel <- hs_future * is.na(novel_future) 
            
            ## Write out not-novel raster
            message('Writing maps of un- novel environments to file for', species) 
            
            writeRaster(hs_future_notNovel, sub('\\.tif', '_notNovel.tif', hs_future@file@name), 
                        overwrite = TRUE, datatype = 'INT2S')
            
            
            ########################################################################################################################
            ## Now create the empty panel just before plotting
            empty_ras <- init(pred.current, function(x) NA) 
            
            ## Check exents
            projection(aus);projection(occ);projection(empty_ras)
            projection(pred.current);projection(pred.future)
            
            identical(extent(pred.current), extent(pred.future))
            
            ########################################################################################################################
            ## Use the levelplot function to make a multipanel output: occurrence points, current raster and future raster
            png(sprintf('%s/%s/full/%s_%s.png', maxent_path, species, species, x),      
                11, 4, units = 'in', res = 300)
            
            ## If possible, add a hatched section to each current and future panel, which shows the novel environments
            ## .....................................................................................................................
            
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
                    ## layer(sp.polygons(aus_albers), data = list(aus_albers = aus_albers))
                    latticeExtra::layer(sp.polygons(aus), data = list(aus = aus)) +
                    latticeExtra::layer(sp.points(occ, pch = 19, cex = 0.15, 
                                                  col = c('red', 'transparent', 'transparent')[panel.number()]), data = list(occ = occ)))
            dev.off()
            
          }
          
        } else {
          
          message(species, ' ', x, ' skipped - prediction already run') ## Skip species which have already been projected
          
        }
        
      } else {
        
        message(species, ' ', x, ' skipped - SDM not yet run')          ## Skip species with no existing maxent model
        
      }
      
    }, species_list, threshold_list, SIMPLIFY = FALSE)
    
  })
  
}