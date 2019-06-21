# env.grids.2030 = tryCatch(project_maxent_grids_mess(poly          = AUS,          ## A shapefile, e.g. Australia
#                                                     scen_list     = scen_2030,    ## A list of climate scenarios
#                                                     species_list  = map_spp_list, ## A list of species folders with maxent models
#                                                     maxent_path   = maxent_path,  ## the output folder
#                                                     climate_path  = "./data/base/worldclim/aus/1km/bio", ## climate data
#                                                     grid_names    = grid.names,   ## names of the predictor grids
#                                                     time_slice    = 30,           ## Time period
#                                                     current_grids = aus.grids.current,
#                                                     create_mess   = "TRUE"),  ## predictor grids
#                           
#                           ## Skip species
#                           error = function(cond) {
#                             
#                             message(paste('Species skipped - check', spp))
#                             
#                           })
                                     

## Try to run the mess maps at the same time as the map creation?
project_maxent_grids_mess = function(poly, scen_list, species_list, maxent_path, climate_path, 
                                     grid_names, time_slice, current_grids, create_mess) {
  
  ## Read in Australia
  poly = poly %>%
    spTransform(ALB.CONICAL)
  
  # Parallelising stuff #####################
  cl <- makeCluster(6)
  clusterExport(cl, c(
    'poly', 'scen_list', 'maxent_path', 'climate_path', 
    'grid_names', 'time_slice', 'current_grids', 'create_mess',
    'hatch', 'polygonizer', 'sdm.predictors', 'aus.grids.current', 
    'sdm.select', 'diverge0'))
  clusterEvalQ(cl, {
    library(rmaxent)
    library(sp)
    library(raster)
    library(rasterVis)
    library(latticeExtra)
  })
  ###########################################
  
  
  ## First, run a loop over each scenario:    
  lapply(scen_list, function(x) {
    
    ## Create a raster stack for each of the 6 GCMs, not for each species
    s <- stack(sprintf('%s/20%s/%s/%s%s.tif', climate_path, time_slice, x, x, 1:19))
    identical(projection(s), projection(poly))
    
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
    
    ## Then apply each GCM to each species.
    ## First, check if the maxent model exists
    ## Can we skip the species before dividing the rasters?................................................................
    ## Then apply each GCM to each species
    
    #lapply(species_list, function(species) {
    
    # Parallelising stuff #####################
    clusterExport(cl, c('x', 's'))
    parLapply(cl, species_list, function(species) {
    ###########################################
      
      save_name = gsub(' ', '_', species)
      if(file.exists(sprintf('%s/%s/full/maxent_fitted.rds', maxent_path, species))) {
        message('Doing ', species)
        
        ## Then, check if the species projection has already been run...
        if(!file.exists(sprintf('%s/%s/full/%s_%s.tif', maxent_path, species, species, x))) {
          
          ## Assign the scenario name (to use later in the plot)
          scen_name = eval(parse(text = sprintf('gcms.%s$GCM[gcms.%s$id == x]', time_slice, time_slice)))           
          
          ########################################################################
          ## Now read in the SDM model, calibrated on current conditions
          m   <- readRDS(sprintf('%s/%s/full/maxent_fitted.rds', maxent_path, species))$me_full 
          swd <- readRDS(sprintf('%s%s/swd.rds',                 maxent_path, species, species))
          occ <- readRDS(sprintf('%s/%s/%s_occ.rds',             maxent_path, species, save_name)) %>%
            spTransform(ALB.CONICAL)  
          
          ## Create a file path for the current raster prediction
          f_current  <- sprintf('%s/%s/full/%s_current.tif', maxent_path, species, species)
          
          ## If the current raster prediction has not been run, run it
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
          
          #####################################################################
          ## Report current mess map in progress
          MESS_dir = sprintf('%s%s/full/%s/', 
                             maxent_path, species, 'MESS_output')
          f_mess_current = sprintf('%s%s%s.tif', MESS_dir, species, "_current_mess_map")
          
          if(create_mess == "TRUE" & !file.exists(f_mess_current)) {
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
            ## Write out the current mess maps - 
            ## create a new folder for the mess output - we are going to print it to the maps
            if(!dir.exists(MESS_dir)) {
              dir.create(MESS_dir)
              
            } else {
              
              message(species, ' MESS directory already created') 
              
            }
            
            ## Then write the mess output to a directory inside the 'full' maxent folder
            writeRaster(mess_current$similarity_min, sprintf('%s%s%s.tif', MESS_dir, species, "_current_mess_map"), 
                        overwrite = TRUE) ## datatype = 'INT2S' for integer 
            
            ##################################################################
            ## Create a PNG file of all the CURRENT MESS output
            ## r    = unstack(mess_current$similarity) :: list of environmental rasters
            ## name = names(mess_current$similarity)   :: names of the rasters
            message('Creating mess maps for each current environmental predictor for', species)
            mapply(function(r, name) {
              
              ## Create a level plot for each species
              p <- levelplot(r, margin = FALSE, scales = list(draw = FALSE),
                             at = seq(minValue(r), maxValue(r), len = 100),
                             colorkey = list(height = 0.6), 
                             main = gsub('_', ' ', sprintf('Current_mess_for_%s (%s)', name, species))) +
                
                latticeExtra::layer(sp.polygons(poly), data = list(poly = poly))  ## need list() for polygon
              
              p <- diverge0(p, 'RdBu')
              f <- sprintf('%s%s%s%s.png', MESS_dir, species, "_current_mess_", name)
              
              png(f, 8, 8, units = 'in', res = 300, type = 'cairo')
              print(p)
              dev.off()
              
            }, unstack(mess_current$similarity), names(mess_current$similarity))
            
            ## Write the raster of novel environments to the maxent directory 
            ## The "full" directory is getting full, could create a sub dir for MESS maps
            message('Writing maps of novel environments to file for', species) 
            writeRaster(novel_current, sprintf('%s%s%s.tif', MESS_dir, species, "_current_novel_map"), 
                        overwrite = TRUE, )
            
            ##################################################################
            # Now mask out novel environments.................................
            # John suggested we might not use the MESS in this way
            
            # is.na(novel_current) is a binary layer showing 
            # not novel [=1] vs novel [=0], 
            # so multiplying this with hs_current will mask out novel
            hs_current_notNovel <- pred.current * is.na(novel_current) 
            
            ## Write out not-novel raster :: this can go to the main directory
            message('Writing maps of un - novel environments to file for', species) 
            writeRaster(hs_current_notNovel, sub('\\.tif', '_notNovel.tif', hs_current@file@name), 
                        overwrite = TRUE)
            
          } else {
            
            message('Dont run MESS maps for ', species) 
            
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
            f_mess_future = sprintf('%s%s%s%s.tif', MESS_dir, species, "_future_mess_", x)
            
            if(create_mess == "TRUE" & !file.exists(f_mess_future)) {
              message('Running current mess map for ', species)
              
              grid_names          = sdm.predictors   ## same grid names
              future_grids        = s                ## the stack of 8 rasters for scenario x
              names(future_grids) = grid_names 
              future_grids        = subset(future_grids, intersect(names(future_grids), sdm.select))
              
              ## Create the future mess map :: check the variable names are the same
              mess_future  <- similarity(future_grids, swd[, names(future_grids)], full = TRUE)
              novel_future <- mess_future$similarity_min   < 0  ##   All novel environments are < 0
              novel_future[novel_future==0] <- NA               ##   0 values are NA
              
              ##################################################################
              ## Write out the future mess maps, for all variables
              writeRaster(mess_future$similarity_min, sprintf('%s%s%s%s.tif', MESS_dir, species, "_future_mess_", x), 
                          overwrite = TRUE)
              
              ##################################################################
              ## Create a PNG file of all the future MESS output: 
              ## r    = unstack(mess_current$similarity) :: list of environmental rasters
              ## name = names(mess_current$similarity)   :: names of the rasters
              message('Creating mess maps of each future environmental predictor for ', x, ' scenario for ', species)
              mapply(function(r, name) {
                
                p <- levelplot(r, margin = FALSE, scales = list(draw = FALSE),
                               at = seq(minValue(r), maxValue(r), len = 100),
                               colorkey = list(height = 0.6), 
                               main = gsub('_', ' ', sprintf('Future_mess_for_%s_%s (%s)', name, x, species, x))) +
                  
                  latticeExtra::layer(sp.polygons(poly), data = list(poly = poly))   ## Use this in previous functions
                
                p <- diverge0(p, 'RdBu')
                f <- sprintf('%s%s%s%s%s%s.png', MESS_dir, species, "_future_mess_", name, "_", x)
                
                png(f, 8, 8, units = 'in', res = 300, type = 'cairo')
                print(p)
                dev.off()
                
              }, unstack(mess_future$similarity), names(mess_future$similarity))
              
              ## Write the raster of novel environments to the maxent directory 
              ## The "full" directory is getting full, could create a sub dir for MESS maps
              message('Writing maps of novel environments to file for', species) 
              writeRaster(novel_future, sprintf('%s%s/full/%s_%s.tif', 
                                                maxent_path, species, species, "future_novel_map"), 
                          overwrite = TRUE)
              
              ##################################################################
              # mask out novel environments 
              # is.na(novel_future) is a binary layer showing 
              # not novel [=1] vs novel [=0], 
              # so multiplying this with hs_future will mask out novel
              hs_future_notNovel <- pred.future * is.na(novel_future) 
              
              ## Write out not-novel raster
              message('Writing maps of un- novel environments to file for', species) 
              writeRaster(hs_future_notNovel, sub('\\.tif', '_notNovel.tif', hs_future@file@name), 
                          overwrite = TRUE)
              
            } else {
              
              message('Dont run MESS maps for ', species) 
              
            }
            
            #####
            ## Convert binary rasters of novel climate to polygons
            novel_current_poly <- polygonizer('output/maxent/SUA_TREES_ANALYSIS/Abelia_grandiflora/full/Abelia_grandiflora_current_novel_map.tif')
            novel_future_poly <- polygonizer('output/maxent/SUA_TREES_ANALYSIS/Abelia_grandiflora/full/Abelia_grandiflora_future_novel_map.tif')
            
            ## Create a SpatialLines object that indicates novel areas (this will be overlaid)            
            # (Below we create a dummy polygon as the first list element (which is the extent
            # of the raster, expanded by 10%), to plot on panel 1)
            novel_hatch <- list(as(extent(pred.current)*1.1, 'SpatialPolygons'),  
                                hatch(novel_current_poly, 50), hatch(novel_future_poly, 50)) # 50 = approx 50 lines across the extent of the poly
            
            
            ########################################################################################################################
            ## Now create the empty panel just before plotting
            empty_ras <- init(pred.current, function(x) NA) 
            
            ## Check exents :: one of the species combinations failed this test
            projection(poly);projection(occ);projection(empty_ras)
            projection(pred.current);projection(pred.future)
            identical(extent(pred.current), extent(pred.future))
            
            ########################################################################################################################
            ## Use the levelplot function to make a multipanel output: occurrence points, current raster and future raster
            if(create_mess == "TRUE") {
              message('Running current mess map for ', species)
              
              ############################################################
              ## Create level plot including MESS maps                        
              png(sprintf('%s/%s/full/%s_%s.png', maxent_path, species, species, x),      
                  11, 4, units = 'in', res = 300)
              
              ## If possible, add a hatched section to each current and future panel, which shows the novel environments
              ## That would be the objects :: 
              ## novel_current & novel_future 
              ## .....................................................................................................................
              
              ## Need an empty frame
              ## Here can we add the novel_current & novel_future layers as a hatching?
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
                      
                                     
                      ## Try adding novel maps as vectors.....................................               
                      latticeExtra::layer(sp.polygons(poly), data = list(poly = poly)) +
                      latticeExtra::layer(sp.points(occ, pch = 19, cex = 0.15, 
                                                    col = c('red', 'transparent', 'transparent')[panel.number()]), data = list(occ = occ)) +
                      latticeExtra::layer(sp.polygons(h[[panel.number()]]), data = list(h = novel_hatch)))
              dev.off()
              
            } else {
              
              ############################################################
              ## OR Create level plot without MESS maps                        
              png(sprintf('%s/%s/full/%s_%s.png', maxent_path, species, species, x),      
                  11, 4, units = 'in', res = 300)
              
              ## Create levelplot without MESS maps
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
                      latticeExtra::layer(sp.polygons(poly), data = list(poly = poly)) +
                      latticeExtra::layer(sp.points(occ, pch = 19, cex = 0.15, 
                                                    col = c('red', 'transparent', 'transparent')[panel.number()]), data = list(occ = occ)))
              dev.off()                     
              
              message('Dont run MESS maps for ', species) 
              
            }
            
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