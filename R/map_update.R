## One function for all time periods
project_maxent_grids = function(scen_list, species_list, maxent_path, climate_path, grid_names, time_slice, current_grids) {
  
  ## Read in Australia
  aus = AUS %>%
    spTransform(ALB.CONICAL)
  
  ## First, run a loop over each GCM:    
  mapply(function(GCM, species) {
    
    ## Create a raster stack for each of the 6 GCMs, not for each species
    s <- stack(sprintf('%s/20%s/%s/%s%s.tif', climate_path, time_slice, GCM, GCM, 1:19))
    projection(s)
    
    ## Rename both the current and future environmental stack...
    names(s) <- names(current_grids) <- grid_names 
    
    ########################################################################################################################
    ## Divide the 11 temperature rasters by 10: NA values are the ocean
    ## s[[1:11]] <- s[[1:11]]/10 ## That code doesn't work
    message('20', time_slice, ' rasters / 10 ', GCM)
    for(i in 1:11) {
      
      ## Simple loop
      message(i)
      s[[i]] <- s[[ i]]/10
      
    }
    
    ## Then apply each GCM to each species     
      ## First, check if the maxent model exists
      ## Can we skip the species before dividing the rasters?................................................................
      save_name = gsub(' ', '_', species)
      if(file.exists(sprintf('%s/%s/full/maxent_fitted.rds', maxent_path, species))) {
        message('Doing ', species)
        
        ## Then, check if the species projection has already been run...
        if(!file.exists(sprintf('%s/%s/full/%s_%s.tif', maxent_path, species, species, GCM))) {
          
          ## Assign the GCM name (to use later in the plot)
          scen_name = eval(parse(text = sprintf('gcms.%s$GCM[gcms.%s$id == GCM]', time_slice, time_slice)))           
          
          ########################################################################################################################
          ## Now read in the SDM model calibrated on current conditions  ## maxent_fitted.rds
          #m <- readRDS(sprintf('%s/%s/full/model.rds', maxent_path, species)) 
          m <- readRDS(sprintf('%s/%s/full/maxent_fitted.rds', maxent_path, species)) 
          m <- m$me_full  ## class(m);View(m)
          
          ## Read in the occurrence points used to create the SDM :: need the transform to plot later
          occ <- readRDS(sprintf('%s/%s/%s_occ.rds', maxent_path, species, save_name)) %>%
            spTransform(ALB.CONICAL)  
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
                              maxent_path, species, species, GCM)
          
          if(!file.exists(f_future)) {
            
            ## Report which prediction is in progress
            message('Running future prediction for ', species, ' ', GCM) 
            
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
            png(sprintf('%s/%s/full/%s_%s.png', maxent_path, species, species, GCM),      
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
                    ## layer(sp.polygons(aus_albers), data = list(aus_albers = aus_albers))
                    latticeExtra::layer(sp.polygons(aus), data = list(aus_albers = aus_albers)) +
                    latticeExtra::layer(sp.points(occ, pch = 19, cex = 0.15, 
                                    col = c('red', 'transparent', 'transparent')[panel.number()]), data = list(occ = occ)))
            dev.off()
            
          }
          
        } else {
          
          message(species, ' ', GCM, ' skipped - prediction already run') ## Skip species which have already been projected
          
        }
        
      } else {
        
        message(species, ' ', GCM, ' skipped - SDM not run')             ## Skip species with no maxent model
        
      }
      
    
  })
  
}