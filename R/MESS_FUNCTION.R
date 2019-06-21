# mapply(function(r, name) {
#   
#   p <- levelplot(r, margin = FALSE, scales = list(draw = FALSE), 
#                  at = seq(minValue(r), maxValue(r), len = 100), 
#                  colokey = list(height = 0.6), main = gsub('_', ' ', sprintf('%s (%s)', name, species))) + 
#     
#     layer(sp.polygons(aus), data = list(aus = aus))  ## Use this in previous functions
#   
#   p <- diverge0(p, 'RdBu')
#   f <- sprintf('%s%s/full/%s_messCurrent__%s.png', maxent_path, species, species, name)
#   
#   png(f, 8, 8, units = 'in', res = 300, type = 'cairo')
#   print(p)
#   dev.off()
#   
#   ## unstack(mess_current$similarity), names(mess_current$similarity)
#   
# }, r, name)
  
## One function for all time periods
create_maxent_mess = function(scen_list, species_list, maxent_path, climate_path, grid_names, time_slice, current_grids) {
  
  ## Read in Australia
  aus = AUS %>%
    spTransform(ALB.CONICAL)
  
  ## First, run a loop over each scenario:  
  ## x = scen_list[1]  
  lapply(scen_list, function(x) {
    
    ## Create a raster stack for each of the 6 GCMs, not for each species
    s_current <- current_grids
    s_future  <- stack(sprintf('%s/20%s/%s/%s%s.tif', climate_path, time_slice, x, x, 1:19))
    projection(s_current);projection(s_future)
    
    ## Rename both the current and future environmental stack...
    names(s_future) <- names(s_future) <- grid_names 
    
    ########################################################################################################################
    ## Divide the 11 temperature rasters by 10: NA values are the ocean
    ## s[[1:11]] <- s[[1:11]]/10 ## That code doesn't work
    message('20', time_slice, ' rasters / 10 ', x)
    for(i in 1:11) {
      
      ## Simple loop
      message(i)
      s_future[[i]] <- s_future[[ i]]/10
      
    }
    
    ## Then apply each GCM to each species
    lapply(species_list, function(species) {
      
      ## Might nee to add another loop over the "threshold" of habitat suitability
      
      ## First, check if the maxent model exists
      ## Can we skip the species before dividing the rasters?................................................................
      save_name = gsub(' ', '_', species)
      if(file.exists(sprintf('%s/%s/full/maxent_fitted.rds', maxent_path, species))) {
        message('Doing ', species)
        
        ## Then, check if the mess map has already been run
        if(!file.exists(sprintf('%s/%s/full/%s_%s%s.tif', maxent_path, species, species, x, "_mess_map"))) {
          
          ## Assign the scenario name (to use later in the plot)
          scen_name = eval(parse(text = sprintf('gcms.%s$GCM[gcms.%s$id == x]', time_slice, time_slice)))           
          
          ########################################################################################################################
          ## Now read in the model, swd file and extract the varaible names 
          m    <- readRDS(sprintf('%s/%s/full/maxent_fitted.rds', maxent_path, species))[["me_full"]]
          swd  <- readRDS(sprintf('%s/%s/swd.rds', maxent_path, species, save_name)) 
          vars <- names(m@presence)
          
          ## Select only the eight bioclim variables used in the maxent models
          current_vars <- dropLayer(s_current, setdiff(names(s_current), vars))
          future_vars  <- dropLayer(s_future,  setdiff(names(s_future),  vars))
          
          ## Read in current continuous raster, and the binary raster thresholded raster (0-1)
          f_current  <- raster(sprintf('%s/%s/full/%s_current.tif', maxent_path, species, species))
          hs_current <- raster(sprintf('%s/%s/full/%s_%s%s.tif',    maxent_path,
                                       species, species, "current_suit_above_", thresh))

          ########################################################################################################################
          ## Create file path for future raster doesn't exist, create it 
          f_future <- sprintf('%s/%s/full/%s_%s%s.tif', maxent_path, species, species, x, "_mess_map")
          
          if(!file.exists(f_future)) {
            
            ## Report current mess map in progress
            message('Running current mess map for ', species) 
            
            ## Create the current mess map :: check the variable names are the same
            setdiff(names(current_vars), names(swd))
            mess_current  <- dismo::mess(current_vars, swd)   ##  time this
            novel_current <- mess_current   < 0               ## "0"-not clear
            novel_current[novel_current==0] <- NA             ## "NA"- not clear
            novel_current <- mask(novel_current, hs_current)  ##  Mask the current raster to the 
            
            ##################################################################
            ## Write out the current mess maps 
            writeRaster(mess_current, sprintf('%s/%s/full/%s_%s.tif', maxent_path,
                                               species, species, "current_mess"), 
                        overwrite = TRUE, datatype = 'INT2S')
  
            writeRaster(novel_current, sprintf('%s/%s/full/%s_%s.tif', maxent_path,
                                               species, species, "novel_current_mess"), 
                        overwrite = TRUE, datatype = 'INT2S')

            ## Report future mess map in progress
            message('Running future mess map for ', species, ' ', x) 
            
            ## Create the future mess map
            setdiff(names(future_vars), names(swd))
            mess_future  <- mess(future_vars, swd)
            novel_future <- mess_future < 0
            novel_future[novel_future==0] <- NA
            novel_future <- mask(novel_future, hs_current)    ## Mask grid to exclude PNG etc.
            
            ##################################################################
            ## Write out the current mess maps 
            writeRaster(mess_future, sprintf('%s/%s/full/%s_%s.tif', maxent_path,
                                              species, species, "future_mess"), 
                        overwrite = TRUE, datatype = 'INT2S')
            
            writeRaster(novel_future, sprintf('%s/%s/full/%s_%s.tif', maxent_path,
                                               species, species, "novel_future_mess"), 
                        overwrite = TRUE, datatype = 'INT2S')
            
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
                    layer(sp.points(occ, pch = 19, cex = 0.15, 
                                    col = c('red', 'transparent', 'transparent')[panel.number()]), data = list(occ = occ)))
            dev.off()
            
            
            ########################################################################################################################
            ## Future mess map
            message('Running current mess map for ', species, ' ', x) 

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