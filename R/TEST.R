#########################################################################################################################
## Loop over directories, species and one threshold for each, also taking a time_slice argument. Next, make the lists generic too
## pervious version in R/old/model_combine.R
SUA_cell_count = function(DIR_list, species_list, maxent_path, thresholds, 
                          percentiles, time_slice, write_rasters) {
  
  ## How can the shapefiles be read in once, not for each species?...................................................
  
  ###################################################################################################################
  ## Read in shapefiles :: this should be done outside the loop
  aus        = readRDS('data/base/CONTEXTUAL/aus_states.rds') %>%
    spTransform(ALB.CONICAL)
  
  LAND       = readRDS('data/base/CONTEXTUAL/LAND_world.rds')
  
  areal_unit = readRDS("./data/base/CONTEXTUAL/SUA/SUA_2016_AUST.rds") %>%
    spTransform(ALB.CONICAL)
  
  ##  Rasterize shapefile
  message('rasterizing SUA shapefile')
  areal_unit = areal_unit[order(areal_unit$SUA_NAME16),]
  f <- tempfile()
  
  ###################################################################################################################
  ## Rasterize the SUA shapefile
  writeOGR(areal_unit, tempdir(), basename(f), 'ESRI Shapefile')
  template <- raster('./output/maxent/SUA_TREES_ANALYSIS/Acacia_retinodes/full/Acacia_retinodes_ac85bi30.tif')
  
  areal_unit_rast <- gdal_rasterize(
    normalizePath(paste0(f, '.shp')), 
    
    'H:/green_cities_sdm/data/base/CONTEXTUAL/SUA/SUA_2016_AUST.tif', tr=res(template),
    te = c(bbox(template)), a = 'SUA_CODE16', a_nodata = 0, init = 0, ot = 'UInt16', output_Raster = TRUE)
  
  areal_unit_vec <- c(areal_unit_rast[])
  summary(areal_unit_vec)
  
  ## Loop over each directory
  lapply(DIR_list, function(DIR) { 
    
    ## And each species - although we don't want all possible combinations. use mapply in the function call
    lapply(species_list, function(species) {
      
      ###################################################################################################################
      ## Create a list of the rasters in each directory, then take the mean. Tidy this up with %, etc
      message('Running summary of SDM predictions within SUAs for ', species, ' using ', names(areal_unit)[1], " shapefile")
      message('Calcualting mean of 20', time_slice, ' GCMs for ', species)
      
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
        
        ## Check if the SUA summary table exists
        SUA_file =   sprintf('%s/%s/full/%s_20%s_%s%s.csv', maxent_path,
                             species, species, time_slice, "SUA_cell_count_", thresh)
        
        ## The mean of the GCMs doesn't exist, create it
        if(file.exists(SUA_file)) { 
          
          message(species, ' SUA_table already exists, skip')
          next
          
        }
        
        for (percent in percentiles) {
          
          ## Print the species being analysed
          message('doing ', species, ' | Logistic > ', thresh, ' for 20', time_slice)
          
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
          band_4           <- function(x) {ifelse(x >=  4, 1, 0) }
          combo_suit_4GCM  <- calc(combo_suit_thresh, fun = band_4)
          
          #########################################################################################################################
          ## Now create a raster of the gain, loss and stable
          ## Create a raster stack of the current and future rasters
          message ("Counting cells lost/gained/stable/never suitable per SUA")
          
          ## Create a table of cell counts using a raster stack of current and future data
          d <- as.data.frame(stack(current_suit_thresh, combo_suit_4GCM)[]) %>% 
            setNames(c('current', 'future')) %>% 
            mutate(SUA_CODE16 = areal_unit_vec,
                   cell_number = seq_len(ncell(current_suit_thresh))) %>% 
            as.tbl
          dim(d);summary(d)
          
          ## Then classify the cells of the raster stack into lost, gained, stable and never suitable
          d2 <- d %>% 
            na.omit %>% 
            
            mutate(lost   = current == 1 & future == 0,
                   gained = current == 0 & future == 1,
                   stable = current == 1 & future == 1,
                   never  = current == 0 & future == 0,
                   nodata = is.na(current) | is.na(future)) 
          d2$class <- apply(select(d2, lost:never), 1, which)
          dim(d2)
          
          ## Then group the cell counts by SUA
          d3 <- d2 %>% 
            group_by(SUA_CODE16) %>%
            
            summarize(CURRENT_SUITABLE = sum(current, na.rm = TRUE),
                      FUTURE_SUITABLE  = sum(future,  na.rm = TRUE),
                      LOST             = sum(lost,    na.rm = TRUE),
                      GAINED           = sum(gained,  na.rm = TRUE),
                      STABLE           = sum(stable,  na.rm = TRUE),
                      NEVER            = sum(never,   na.rm = TRUE),
                      NODAT            = sum(nodata,  na.rm = TRUE),
                      n_cells = n()) %>% 
            
            ## Then calculate change between current and future
            mutate(CHANGE    = FUTURE_SUITABLE - CURRENT_SUITABLE,
                   GAIN_LOSS = ifelse(CHANGE < 0, 'LOSS', ifelse(CHANGE > 0, 'GAIN', 'STABLE')),
                   GAIN_LOSS = ifelse(CURRENT_SUITABLE == 0 & FUTURE_SUITABLE == 0, 'NEVER', GAIN_LOSS))
          dim(d3)
          
          ## Add the species column
          d4 = d3 %>% 
            join(areal_unit@data, .) %>%
            add_column(., SPECIES = species,    .after = "AREASQKM16") %>%
            add_column(., PERIOD  = time_slice, .after = "SPECIES")    %>%
            add_column(., THRESH  = thresh,     .after = "PERIOD")
          View(d4)
          
          #r[d2$cell_number] <- d2$class
          
          #########################################################################################################################
          ## Now calculate the number of cells lost/gained/stable across Australia
          message ("Counting cells lost/gained/stable/never suitable across Australia")
          d5 <- stack(current_suit_thresh, combo_suit_4GCM)[]
          r <- raster(current_suit_thresh)
          z <- as.data.frame(d4)
          
          ## Then classify the raster stack to make each value (i.e. outcome) unique
          r[d5[, 1]==1 & d5[, 2]==0] <- 1  ## 1 in current raster and 0 in future = LOSS
          r[d5[, 1]==0 & d5[, 2]==1] <- 2  ## 0 in current raster and 1 in future = GAIN
          r[d5[, 1]==1 & d5[, 2]==1] <- 3  ## 1 in current raster and 1 in future = STABLE
          r[d5[, 1]==0 & d5[, 2]==0] <- 4  ## 0 in current raster and 0 in future = NEVER_SUIT
          
          
          ## Now convert the raster to a factor and assign lables to the levels
          gain_loss <- as.factor(r)
          levels(gain_loss)[[1]] <- data.frame(ID = 1:4, label = c('Lost', 'Gained', 'Stable', 'Never_Suitable'))
          z <- as.data.frame(d5)
          
          ## The gain/loss raster could be intersected with the SUAs, instead of the combo_suit_4GCM.................................
          ## Create a table of the gain/loss/stable :: write this to file as well
          gain_loss_table      = table(z[, 1], z[, 2])
          gain_loss_df         = as.data.frame(raster::freq(gain_loss))
          gain_loss_df$SPECIES = species
          gain_loss_df$PERIOD  = time_slice
          
          names(gain_loss_df)  = c("CHANGE", "COUNT", "SPECIES", "PERIOD")
          gain_loss_df         = gain_loss_df[, c("SPECIES", "PERIOD", "CHANGE", "COUNT")]
          
          ## Change values and remove the NA row
          gain_loss_df$CHANGE[gain_loss_df$CHANGE == 1] <- "LOST"
          gain_loss_df$CHANGE[gain_loss_df$CHANGE == 2] <- "GAINED"
          gain_loss_df$CHANGE[gain_loss_df$CHANGE == 3] <- "STABLE"
          gain_loss_df$CHANGE[gain_loss_df$CHANGE == 4] <- "NEVER_SUIT"
          gain_loss_df = head(gain_loss_df, 4)
          head(gain_loss_df)
          
          #########################################################################################################################
          ## Save the continental gain/loss table
          write.csv(gain_loss_df, sprintf('%s/%s/full/%s_20%s_%s%s.csv', maxent_path,
                                          species, species, time_slice, "gain_loss_table_", thresh), row.names = FALSE)
          
          #########################################################################################################################
          ## Save the SUA gain/loss table
          write.csv(d4, sprintf('%s/%s/full/%s_20%s_%s%s.csv', maxent_path,
                                species, species, time_slice, "SUA_cell_count_", thresh), row.names = FALSE)
          
          #########################################################################################################################
          #########################################################################################################################
          ## Now write the rasters
          ## If the rasters don't exist, write them for each species/threshold
          if(write_rasters == "TRUE") {
            
            ## Write the current suitability raster, thresholded using the Maximum training sensitivity plus specificity Logistic threshold
            message('Writing ', species, ' current', ' max train > ', thresh)
            writeRaster(current_suit_thresh, sprintf('%s/%s/full/%s_%s%s.tif', maxent_path,
                                                     species, species, "current_suit_above_", thresh), overwrite = TRUE)
            
            ## Write the combined suitability raster, thresholded using the maximum training value
            message('Writing ', species, ' | 20', time_slice, ' max train > ', thresh)
            writeRaster(combo_suit_thresh, sprintf('%s/%s/full/%s_20%s%s%s.tif', maxent_path,
                                                   species, species, time_slice, "_Max_train_sensit_above_", thresh), overwrite = TRUE)
            
            ## Write the combined suitability raster, thresholded using the percentile value
            message('Writing ', species, ' | 20', time_slice, ' 10th percentile > ', percent)
            writeRaster(combo_suit_percent, sprintf('%s/%s/full/%s_20%s%s%s.tif', maxent_path,
                                                    species, species, time_slice, "_10_percentile_omiss_above_", percent), overwrite = TRUE)
            
            ## Write the combined future raster with > 4 GCMs above the maximum training value
            message('Writing ', species, ' | 20', time_slice, ' 4 GCMs > ', percent)
            writeRaster(combo_suit_4GCM, sprintf('%s/%s/full/%s_20%s%s%s.tif', maxent_path,
                                                 species, species, time_slice, "_4GCMs_above_", thresh), overwrite = TRUE)
            
            ## Write out the gain/loss raster
            writeRaster(gain_loss, sprintf('%s/%s/full/%s_20%s%s%s.tif', maxent_path,
                                           species, species, time_slice, "_gain_loss_", thresh), datatype = 'INT2U', overwrite = TRUE)
            
            
            #########################################################################################################################
            ##  Write PNG files too
            
          } else {
            
            message(' skip raster writing')   ##
            
          }
          
        }
        
      }
      
    })
    
  })
  
}