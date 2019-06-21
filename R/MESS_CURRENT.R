# species       = mess_spp[1]
# threshold     = MESS.THRESH[1]
# maxent_path   = maxent_path 
# climate_path  = "./data/base/worldclim/aus/1km/bio"
# grid_names    = grid.names
# current_grids = aus.grids.current


## One function for all time periods
create_maxent_mess = function(species_list, threshold_list, maxent_path, climate_path, grid_names, current_grids) {
  
  ## Read in Australia
  aus = AUS %>%
    spTransform(ALB.CONICAL)
  
  ## First, run a loop over each species  
  lapply(species_list, function(species) {
    
    ## Create a raster stack for each of the 6 GCMs, not for each species
    s_current <- current_grids
    projection(s_current)
    
    ## Rename both the current and future environmental stack...
    names(s_current) <- grid_names 
    
    ## Then loop over each threshold
    lapply(threshold_list, function(threshold) {
      
      ## First, check if the maxent model exists
      save_name = gsub(' ', '_', species)
      message('Doing ', species)
      
      ## Then, check if the mess map has already been run
      if(!file.exists(sprintf('%s%s/full/%s_%s.tif', maxent_path, species, species, "current_mess_map"))) {
        
        #####################################################################
        ## Now read in the model, swd file, occ data and extract the varaible names 
        m    <- readRDS(sprintf('%s%s/full/maxent_fitted.rds', maxent_path, species))[["me_full"]]
        swd  <- readRDS(sprintf('%s%s/swd.rds', maxent_path, species, save_name))
        occ  <- readRDS(sprintf('%s%s/swd.rds', maxent_path, species, save_name))
        occ  <- readRDS(sprintf('%s%s/%s_occ.rds', maxent_path, species, save_name)) %>%
          spTransform(ALB.CONICAL) 
        vars <- names(m@presence)
        
        ## Select only the eight bioclim variables used in the maxent models
        #current_vars <- dropLayer(s_current, setdiff(names(s_current), vars))
        current_vars <- subset(s_current, intersect(names(s_current), vars))


        ## Read in current continuous raster, and the binary raster thresholded raster (0-1)
        f_current  <- raster(sprintf('%s%s/full/%s_current.tif', maxent_path, species, species))
        hs_current <- raster(sprintf('%s%s/full/%s_%s%s.tif',    maxent_path,
                                     species, species, "current_suit_above_", threshold))
                                     
        ## Check the plot
        ## plot(hs_current) 
        ## points(occ, pch = 19, col = "red", cex = 0.4)
        
        #####################################################################
        ## Report current mess map in progress
        message('Running current mess map for ', species) 
        
        ## Create the current mess map :: check the variable names are the same
        setdiff(names(current_vars), names(swd))
        mess_current  <- dismo::mess(current_vars, swd)   ##  time this
        novel_current <- mess_current   < 0               ## "0"-not clear
        novel_current[novel_current==0] <- NA             ## "NA"- not clear
        novel_current <- mask(novel_current, hs_current)  ##  Mask the current raster to the current suitability
        
        ##################################################################
        ## Write out the current mess maps 
        writeRaster(mess_current, sprintf('%s%s/full/%s_%s.tif', 
                                          maxent_path, species, species, "current_mess_map"), 
                    overwrite = TRUE, datatype = 'INT2S')
        
        writeRaster(novel_current, sprintf('%s%s/full/%s_%s.tif', 
                                           maxent_path, species, species, "current_novel_map"), 
                    overwrite = TRUE, datatype = 'INT2S')
        
        ## Now create the empty panel just before plotting
        empty_ras <- init(mess_current, function(x) NA) 
        
        ## Check projections ..................................................................................................
        projection(aus);projection(occ);projection(empty_ras);projection(mess_current);projection(novel_current)
        
        ##################################################################
        ## Convert to polygon ...this may not work on your computer... 
        ## it requires GDAL and Python to be set up in a specific way
        # poly       <- polygonizer(novel_current)
        # poly_hatch <- hatch(poly, 60)
        # 
        # levelplot(hs_current, at = seq(0, 1, len = 101), 
        #           
        #           col.regions = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(100),
        #           margin      = FALSE, scales = list(draw = FALSE), 
        #           colorkey    = list(height = 0.5)) +
        #   
        #   layer(sp.polygons(aus)) +
        #   layer(sp.polygons(poly_hatch))
        
        ##################################################################
        ## Print to screen
        print(levelplot(stack(empty_ras, 
                              mess_current, 
                              novel_current,
                              hs_current, 
                              quick = TRUE), margin = FALSE,
                        
                        ## Create a colour scheme using colbrewer: 100 is to make it continuos
                        ## Also, make it a one-directional colour scheme
                        scales      = list(draw = FALSE), 
                        at = seq(0, cellStats(mess_current, stat = 'min'), length = 100),
                        col.regions = colorRampPalette(rev(brewer.pal(11, 'Spectral'))),
                        
                        ## Give each plot a name: the third panel is the GCM
                        names.attr = c('Australian records', 'Current Mess', 'Novel env', 'Current threshold'),
                        colorkey   = list(height = 0.5, width = 3), xlab = '', ylab = '',
                        main       = list(gsub('_', ' ', species), font = 4, cex = 2)) +
                
                ## Plot the Aus shapefile with the occurrence points for reference
                ## Can the points be made more legible for both poorly and well recorded species?
                layer(sp.polygons(aus)) +
                layer(sp.points(occ, pch = 19, cex = 0.15, 
                                col = c('red', 'transparent', 'transparent')[panel.number()]), data = list(occ = occ)))
        
        ##################################################################
        ## Use the levelplot function to make a multipanel output: 
        ## occurrence points, current raster and future raster
        png(sprintf('%s%s/full/%s_%s.png', maxent_path, species, species, "current_mess_map"),      
            11, 4, units = 'in', res = 300)
        
        ## Need an empty frame
        print(levelplot(stack(empty_ras, 
                              mess_current, 
                              novel_current,
                              hs_current, 
                              quick = TRUE), margin = FALSE,
                        
                        ## Create a colour scheme using colbrewer: 100 is to make it continuos
                        ## Also, make it a one-directional colour scheme
                        scales      = list(draw = FALSE), 
                        at = seq(0, cellStats(mess_current, stat = 'min'), length = 100),
                        col.regions = colorRampPalette(rev(brewer.pal(11, 'Spectral'))),
                        
                        ## Give each plot a name: the third panel is the GCM
                        names.attr = c('Australian records', 'Current Mess', 'Novel env', 'Current threshold'),
                        colorkey   = list(height = 0.5, width = 3), xlab = '', ylab = '',
                        main       = list(gsub('_', ' ', species), font = 4, cex = 2)) +
                
                ## Plot the Aus shapefile with the occurrence points for reference
                ## Can the points be made more legible for both poorly and well recorded species?
                layer(sp.polygons(aus)) +
                layer(sp.points(occ, pch = 19, cex = 0.15, 
                                col = c('red', 'transparent', 'transparent')[panel.number()]), data = list(occ = occ)))
        dev.off()
        
      } else {
        
        message(species, ' skipped - MESS map already run') 
        
      }
      
      
    })
    
  })
  
}