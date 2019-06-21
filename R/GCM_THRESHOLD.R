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
            
            
            
            #########################################################################################################################
            ## Then using this GCM consensus, calculate whether the species is likely to be present in each SUA. Decide on a threshold 
            ## of % area (10?) of the SUA that needs to be occupied, for each species to be considered present. 
            message('Running zonal stats for ', species, ' | 20', time_slice, ' combined suitability > ', thresh)
            names(SUA.RAS@data@attributes)
            
            ## Check the order of lists match, species, SUAs, areas need to match up ................................................
            unique(SUA.RAS)
            test = SUA.RAS * combo_suit_4GCM
 
            
            ## Should we also create a simple count of the no. of cells per SUA with where > 4 GCMs met the suitability threshold?
            ## Need to fix this so that it has the same order as the shapefile
            #sp.count <- spatialEco::zonal.stats(x = areal_unit, y = combo_suit_4GCM, stat = sum, trace = TRUE, plot = TRUE) 
            # z.med   <- spatialEco::zonal.stats(x = areal_unit, y = combo_suit_4GCM, stat = median, trace = TRUE, plot = TRUE)
            
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