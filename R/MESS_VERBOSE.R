## One function for all time periods
create_maxent_mess = function(species_list, threshold_list, maxent_path, current_grids) {
  
  ## Read in Australia
  aus = AUS %>%
    spTransform(ALB.CONICAL)
  
  # import gist that plots maps with diverging colour ramps
  devtools::source_gist('306e4b7e69c87b1826db', filename='diverge0.R') # function that plots a map with a diverging colour ramp
  
  ## First, run a loop over each species  
  mapply(function(species, threshold) {
    
    ## Create a raster stack for each of the 6 GCMs, not for each species
    
    #####################################################################
    ## Now read in the model, swd file, occ data and extract the varaible names 
    m    <- readRDS(sprintf('%s%s/full/maxent_fitted.rds', maxent_path, species))[["me_full"]]
    swd  <- readRDS(sprintf('%s%s/swd.rds',                maxent_path, species, species))
    occ  <- readRDS(sprintf('%s%s/swd.rds',                maxent_path, species, species))
    occ  <- readRDS(sprintf('%s%s/%s_occ.rds',             maxent_path, species, species)) %>%
      spTransform(ALB.CONICAL) 
    
    ## First, check if the mess map has already been run
    if(!file.exists(sprintf('%s%s/full/%s_%s%s%s.tif', maxent_path, species, species, 
                            "current_suit_above_", threshold, "_notNovel"))) {
      
      ## Then read in the current continuous raster, and the binary raster thresholded raster (0-1)
      f_current  <- raster(sprintf('%s%s/full/%s_current.tif', maxent_path, species, species))
      hs_current <- raster(sprintf('%s%s/full/%s_%s%s.tif',    maxent_path,
                                   species, species, "current_suit_above_", threshold))
      
      #####################################################################
      ## Report current mess map in progress
      message('Running current mess map for ', species) 
      
      ## Create the current mess map :: check the variable names are the same
      #length(intersect(names(current_grids), names(swd)))
      mess_current <- similarity(current_grids, swd[, names(current_grids)], full = TRUE)
      novel_current <- mess_current$similarity_min   < 0  ##   All novel environments are < 0
      novel_current[novel_current==0] <- NA               ##   0 values are NA
      
      ##################################################################
      ## Write out the current mess maps 
      writeRaster(mess_current$similarity_min, sprintf('%s%s/full/%s_%s.tif', 
                                                       maxent_path, species, species, "current_mess_map"), 
                  overwrite = TRUE, datatype = 'INT2S')
      
      ##################################################################
      ## Create a PNG file of all the MESS output
      message('Creating mess maps for each environmental predictor for', species)
      
      ## First explictily create the lists of grids and names
      grids      = unstack(mess_current$similarity)
      grid_names = names(mess_current$similarity)
      
      ## Loop over each predictor and name
      for (i in 1:length(grids)) {
        
        for (name in grid_names) {
          
          ## Create a level plot
          p <- levelplot(grids[[i]], margin = FALSE, scales = list(draw = FALSE),
                         at = seq(minValue(grids[[i]]), maxValue(grids[[i]]), len = 100),
                         colokey = list(height = 0.6), main = gsub('_', ' ', sprintf('%s (%s)', name, species))) +
            
            layer(sp.polygons(aus), data = list(aus = aus))  ## Use this in previous functions
          
          p <- diverge0(p, 'RdBu')
          f <- sprintf('%s%s/full/%s_messCurrent_%s.png', maxent_path, species, species, name)
          
          png(f, 8, 8, units = 'in', res = 300, type = 'cairo')
          print(p)
          
          ## finish the device
          dev.off()
          
          ## Write the raster of novel environments to the maxent directory 
          ## The "full" directory is getting full, could create a sub dir for MESS maps
          message('Writing maps of novel environments to file for', species) 
          
          writeRaster(novel_current, sprintf('%s%s/full/%s_%s.tif', 
                                             maxent_path, species, species, "current_novel_map"), 
                      overwrite = TRUE, datatype = 'INT2S')
          
          hs_current_notNovel <- hs_current * is.na(novel_current) 
          
          
          ##################################################################
          # mask out novel environments 
          # is.na(novel_current) is a binary layer showing 
          # not novel [=1] vs novel [=0], 
          # so multiplying with hs_current will mask out novel
          message('Writing maps of un- novel environments to file for', species) 
          
          writeRaster(hs_current_notNovel, sub('\\.tif', '_notNovel.tif', hs_current@file@name), 
                      overwrite = TRUE, datatype = 'INT2S')
          
        }
        
      }
      
    } else {
      
      message(species, ' skipped - MESS map already run') 
      
    }
    
  }, species_list, threshold_list, SIMPLIFY = FALSE)
  
}