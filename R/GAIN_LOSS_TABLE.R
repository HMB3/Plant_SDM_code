gain_loss_table = function(DIR_list, species_list, maxent_path, thresholds, percentiles, time_slice) {
  
  ## How can the shapefiles be read in once, not for each species?...................................................
  
  ###################################################################################################################
  ## Read in shapefiles :: this should be done outside the loop
  aus        = readRDS("F:/green_cities_sdm/data/base/CONTEXTUAL/aus_states.rds") %>%
    spTransform(ALB.CONICAL)
  
  LAND       = readRDS("F:/green_cities_sdm/data/base/CONTEXTUAL/LAND_world.rds")
  
  areal_unit = readRDS("F:/green_cities_sdm/data/base/CONTEXTUAL/SUA.rds") %>%
    spTransform(ALB.CONICAL)
  
  areal_unit = areal_unit[order(areal_unit$SUA_NAME11),]
  
  ## Loop over each directory
  lapply(DIR_list, function(DIR) { 
    
    ## And each species - although we don't want all possible combinations. use mapply in the function call
    lapply(species_list, function(species) {
      
      # for (slice in time_slice) {
      
      ###################################################################################################################
      ## Create a list of the rasters in each directory, then take the mean. Tidy this up with %, etc
      message('Calcualting mean of GCMs for ', species)
      
      ## Check if the mean GCM raster exists
      f_mean = sprintf('%s/%s/full/%s_20%s_suitability_mean.tif', maxent_path, species, species, time_slice)
      
      ## The mean of the GCMs doesn't exist, create it
      if(!file.exists(f_mean)) { 
        
        raster.list       = list.files(as.character(DIR), pattern = sprintf('bi%s.tif$', time_slice), full.names = TRUE)  
        suit              = stack(raster.list)
        suit.list         = unstack(suit)
        combo_suit_mean   = mean(suit)                            ## plot(mean.suit)
        
        writeRaster(combo_suit_mean , sprintf('%s/%s/full/%s_20%s_suitability_mean.tif', 
                                              maxent_path, species, species, time_slice), overwrite = TRUE)
        
      } else {
        
        ## Create another level without the mean calculation
        raster.list = list.files(as.character(DIR), pattern = sprintf('bi%s.tif$', time_slice), full.names = TRUE)  
        suit        = stack(raster.list)
        suit.list   = unstack(suit)
        
      }
      
      #########################################################################################################################
      ## Then, create rasters that meet habitat suitability criteria thresholds, determined by the rmaxent function
      for (thresh in thresholds) {
        
        for (percent in percentiles) { 
          
          ## Can we set 0 to NA in the rasters before running the calculations?
          ## Check if the combined suitability raster exists
          f_max_train_suit <- sprintf('%s/%s/full/%s_20%s%s%s.tif', maxent_path,
                                      species, species, time_slice, "_Max_train_sensit_above_", thresh)
          
          ## If it doesn't exist, create the suitability raster
          #if(!file.exists(f_max_train_suit)) {
          
          ## Print the species being analysed
          message('doing ', species, ' | Max train sensit > ', thresh, ' for 20', time_slice)
          
          ## Read in the current suitability raster
          f_current <- raster(sprintf('%s/%s/full/%s_current.tif', 
                                      maxent_path, species, species))
          
          ## First, create a simple function to threshold each of the rasters in raster.list
          thresh_greater  = function (x) {x > thresh}
          percent_greater = function (x) {x > percent}
          
          ## Then apply this to just the current suitability raster. These functions use the : 
          ## Maximum training sensitivity plus specificity Logistic threshold
          ## 10th percentile training presence training omission
          current_suit_thresh  = thresh_greater(f_current)
          current_suit_percent = percent_greater(f_current) 
          
          ## Now, apply these functions to the list of rasters (6 GCMs) for each species
          ## Also, the 'init' function initializes a raster object with values
          ## suit_ras_greater    = reduce(suit.list, thresh_above_fun, .init = suit.list[[1]] > thresh)
          
          ## Check the logic of removing the zeros.....................................................
          ## If we are just using the combined rasters without calculating the difference, don't worry about the zeros   
          
          #########################################################################################################################
          ## First, calculate the cells which are greater that the: 
          ## Maximum training sensitivity plus specificity Logistic threshold
          message('Running thresholds for ', species, ' | 20', time_slice, ' combined suitability > ', thresh)
          
          suit_ras1_thresh   = thresh_greater(suit.list[[1]])   ## Abbreviate these...
          suit_ras2_thresh   = thresh_greater(suit.list[[2]])
          suit_ras3_thresh   = thresh_greater(suit.list[[3]])
          suit_ras4_thresh   = thresh_greater(suit.list[[4]])   
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
          
          ## Another problem is that there is no current OR future 
          #if(summary(combo_suit_4GCM)[3, "layer"] > 0) {
          
          #########################################################################################################################
          ## Now create a raster of the gain, loss and stable
          ## Create a raster stack
          d <- stack(current_suit_thresh, combo_suit_4GCM)[]
          r <- raster(current_suit_thresh)
          z <- as.data.frame(d)
          
          ## Then classify the raster stack to make each value (i.e. outcome) unique
          r[d[, 1]==1 & d[, 2]==0] <- 1  ## 1 in current raster and 0 in future = LOSS
          r[d[, 1]==0 & d[, 2]==1] <- 2  ## 0 in current raster and 1 in future = GAIN
          r[d[, 1]==1 & d[, 2]==1] <- 3  ## 1 in current raster and 1 in future = STABLE
          
          ## Now convert the raster to a factor and assign lables to the levels
          gain_loss <- as.factor(r)
          levels(gain_loss)[[1]] <- data.frame(ID = 1:3, label = c('Lost', 'Gained', 'Stable'))
          
          ## Create a table of the gain/loss/stable :: write this to file as well
          gain_loss_table  = table(z[, 1], z[, 2])
          gain_loss_df     = as.data.frame(raster::freq(gain_loss))
          
          ## Save the gain/loss table
          saveRDS(gain_loss_df, sprintf('%s/%s/full/%s_20%s_%s%s.rds', maxent_path,
                                        species, species, time_slice, "gain_loss_table_", thresh))
          
        }
        
      }
      
    })
    
  })
  
}





