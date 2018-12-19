#########################################################################################################################
############################################# FUNCTIONS FOR MAPPING SDMS ################################################ 
#########################################################################################################################


#########################################################################################################################
## PROJECT MAXENT MODELS
#########################################################################################################################


## flag isses with.......................................................................................................


#########################################################################################################################
## Create maxent maps for a given time period 
# x             = scen_2030[1]
# species       = map_spp_list[1]
# time_slice    = 30
# maxent_path   = "./output/maxent/WITHOUT_INV"
# climate_path  = "./data/base/worldclim/aus/1km/bio"
# grid_names    = grid.names
# current_grids = aus.grids.current
# scen_list     = scen_2030


## One function for all time periods
project_maxent_grids = function(shp, species_list, maxent_path, 
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
    lapply(species_list, function(species) {
      
      ## First, check if the maxent model exists
      ## Can we skip the species before dividing the rasters?................................................................
      save_name = gsub(' ', '_', species)
      if(file.exists(sprintf('%s/%s/full/maxent_fitted.rds', maxent_path, species))) {
        message('Doing ', species)
        
        ## Then, check if the species projection has already been run...
        if(!file.exists(sprintf('%s/%s/full/%s_%s.tif', maxent_path, species, species, x))) {
          
          ## Assign the scenario name (to use later in the plot)
          scen_name = eval(parse(text = sprintf('gcms.%s$GCM[gcms.%s$id == x]', time_slice, time_slice)))           
          
          ########################################################################################################################
          ## Now read in the SDM model calibrated on current conditions
          m <- readRDS(sprintf('%s/%s/full/maxent_fitted.rds', maxent_path, species)) 
          m <- m$me_full  ##
          
          ## Read in the occurrence points used to create the SDM :: need the transform to plot later
          occ <- readRDS(sprintf('%s/%s/%s_occ.rds', maxent_path, species, save_name)) %>%
            spTransform(ALB.CONICAL)  
          
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
                              maxent_path, species, species, x)
          
          if(!file.exists(f_future)) {
            
            ## Report which prediction is in progress
            message('Running future prediction for ', species, ' ', x) 
            
            ## Create the future raster
            pred.future <- rmaxent::project(
              m, s[[colnames(m@presence)]])$prediction_logistic
            writeRaster(pred.future, f_future, overwrite = TRUE)
            
            ## Re-project the current raster here
            #pred.current <- projectRaster(pred.current, crs = ALB.CONICAL)
            
            ## Now create the empty panel just before plotting
            empty_ras <- init(pred.current, function(x) NA) 
            
            ## Check exents, what is going wromg here?...........................................................................
            projection(aus);projection(occ);projection(empty_ras)
            projection(pred.current);projection(pred.future)
            
            identical(extent(pred.current), extent(pred.future))
            
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
      
    })
    
  })
  
}





#########################################################################################################################
## COMBINE GCM FUNCTIONS
#########################################################################################################################


#########################################################################################################################
## Loop over directories, species and one threshold for each, also taking a time_slice argument. Next, make the lists generic too
## pervious version in R/old/model_combine.R
SUA_cell_count = function(shp, vec, DIR_list, species_list, 
                          maxent_path, thresholds, percentiles, 
                          time_slice, write_rasters) {
  
  ## How can the shapefiles be read in once, not for each species?...................................................
  
  ###################################################################################################################
  ## Read in shapefiles :: this should be done outside the loop
  # aus        = readRDS('data/base/CONTEXTUAL/aus_states.rds') %>%
  #   spTransform(ALB.CONICAL)
  # 
  # LAND       = readRDS('data/base/CONTEXTUAL/LAND_world.rds')
  
  areal_unit      = shp 
  areal_unit_vec  = vec
  
  
  ###################################################################################################################
  ## Read in rasters of the SUA shapefile
  # areal_unit_rast = rast
  # areal_unit_vec  = vec

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
        # SUA_file =   sprintf('%s/%s/full/%s_20%s_%s%s.csv', maxent_path,
        #                      species, species, time_slice, "SUA_cell_count_", thresh)
        
        SUA_file =   sprintf('%s/%s/full/%s_20%s_%s%s.tif', maxent_path,
                             species, species, time_slice, "gain_loss_", thresh)

        
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
          #View(d4)
          
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
            
            ## Write the current suitability raster, thresholded using the Maximum training 
            ## sensitivity plus specificity Logistic threshold
            message('Writing ', species, ' current', ' max train > ', thresh)
            writeRaster(current_suit_thresh, sprintf('%s/%s/full/%s_%s%s.tif', maxent_path,
                                                     species, species, "current_suit_above_", thresh), 
                        overwrite = TRUE)
            
            ## Write the combined suitability raster, thresholded using the maximum training value
            message('Writing ', species, ' | 20', time_slice, ' max train > ', thresh)
            writeRaster(combo_suit_thresh, sprintf('%s/%s/full/%s_20%s%s%s.tif', maxent_path,
                                                   species, species, time_slice, "_Max_train_sensit_above_", thresh), 
                        overwrite = TRUE)
            
            ## Write the combined suitability raster, thresholded using the percentile value
            message('Writing ', species, ' | 20', time_slice, ' 10th percentile > ', percent)
            writeRaster(combo_suit_percent, sprintf('%s/%s/full/%s_20%s%s%s.tif', maxent_path,
                                                    species, species, time_slice, "_10_percentile_omiss_above_", percent), 
                        overwrite = TRUE)
            
            ## Write the combined future raster with > 4 GCMs above the maximum training value
            message('Writing ', species, ' | 20', time_slice, ' 4 GCMs > ', percent)
            writeRaster(combo_suit_4GCM, sprintf('%s/%s/full/%s_20%s%s%s.tif', maxent_path,
                                                 species, species, time_slice, "_4GCMs_above_", thresh), 
                        overwrite = TRUE)
            
            ## Write out the gain/loss raster
            writeRaster(gain_loss, sprintf('%s/%s/full/%s_20%s%s%s.tif', maxent_path,
                                           species, species, time_slice, "_gain_loss_", thresh), datatype = 'INT2U', 
                        overwrite = TRUE)
            
            #########################################################################################################################
            ## Write gain/loss PNG too
            # save_name = gsub(' ', '_', species)
            # empty_ras <- init(current_suit_thresh, function(x) NA)
            # 
            # png(sprintf('%s/%s/full/%s_%s.png', maxent_path, save_name, save_name, "current_suit"),
            #     11, 4, units = 'in', res = 300)
            # # 
            # # message('writing thresholded map for ', 'species')
            # # 
            # # occ <- readRDS(sprintf('%s/%s/%s_occ.rds', maxent_path, species, save_name)) %>%
            # #   spTransform(ALB.CONICAL)  
            # # 
            # ## Need an empty frame
            # print(levelplot(stack(f_current,
            #                       combo_suit_thresh,
            #                       #gain_loss,
            #                       quick = TRUE), margin = FALSE,
            # 
            #                 ## Create a colour scheme using colbrewer: 100 is to make it continuos
            #                 ## Also, make it a one-directional colour scheme
            #                 #scales      = list(draw = FALSE),
            #                 #at = seq(0, 1, length = 100),
            #                 #col.regions = colorRampPalette(rev(brewer.pal(11, 'Spectral'))),
            #                 # str(gain_loss@data@attributes)
            # 
            #                 ## Give each plot a name: the third panel is the GCM
            #                 names.attr = c('Current climatic suitability', paste0('Suitability in 20', 
            #                 time_slice)),#, 'Gain/loss raster'),
            #                 colorkey   = list(height = 0.5, width = 3), xlab = '', ylab = '',
            #                 main       = list(gsub('_current_suit_above_', '_', names(spp.thresh)), font = 4, cex = 2)) +
            # 
            #         ## Plot the Aus shapefile with the occurrence points for reference
            #         ## Can the points be made more legible for both poorly and well recorded species?
            #         latticeExtra::layer(sp.polygons(aus), data = list(aus = aus)))
            #         #layer(sp.polygons(aus), data = list(aus = aus))) #+
            #         # layer(sp.points(occ, pch = 19, cex = 0.15,
            #         #                 col = c('red', 'transparent', 'transparent')[panel.number()]), data = list(spp.points = occ)))
            # 
            # dev.off()
            
            
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











#########################################################################################################################
## MESS MAP FUNCTIONS
#########################################################################################################################


# species       = mess_spp[1]
# threshold     = MESS.THRESH[1]
# maxent_path   = maxent_path 
# climate_path  = "./data/base/worldclim/aus/1km/bio"
# grid_names    = grid.names
# current_grids = aus.grids.current


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
    swd  <- readRDS(sprintf('%s%s/swd.rds', maxent_path, species, species))
    occ  <- readRDS(sprintf('%s%s/swd.rds', maxent_path, species, species))
    occ  <- readRDS(sprintf('%s%s/%s_occ.rds', maxent_path, species, species)) %>%
      spTransform(ALB.CONICAL) 
    #ras_names <- names(m@presence)
    
    ## Select only the eight bioclim variables used in the maxent models
    #current_vars <- dropLayer(s_current, setdiff(names(s_current), vars))
    # names(current_grids) <- grid_names  ## Error in `names<-`(`*tmp*`, value = grid_names) : incorrect number of layer names
    # current_grids        <- subset(current_grids, intersect(names(current_grids), ras_names))
    
    ## First, check if the mess map has already been run
    if(!file.exists(sprintf('%s%s/full/%s_%s%s%s.tif', maxent_path, species, species, 
                            "current_suit_above_", thresh, "_notNovel"))) {
      
      ## Read in current continuous raster, and the binary raster thresholded raster (0-1)
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
      novel_current[novel_current==0] <- NA               ##  0 values are NA
      
      ##################################################################
      ## Write out the current mess maps 
      writeRaster(mess_current$similarity_min, sprintf('%s%s/full/%s_%s.tif', 
                                                       maxent_path, species, species, "current_mess_map"), 
                  overwrite = TRUE, datatype = 'INT2S')
      
      ##################################################################
      ## Create a PNG file of all the MESS output
      mapply(function(r, name) {
        
        p <- levelplot(r, margin = FALSE, scales = list(draw = FALSE), 
                       at = seq(minValue(r), maxValue(r), len = 100), 
                       colokey = list(height = 0.6), main = gsub('_', ' ', sprintf('%s (%s)', name, species))) + 
          
          layer(sp.polygons(aus_albers), data = list(aus_albers = aus_albers))  ## Use this in previous functions
        
        p <- diverge0(p, 'RdBu')
        f <- sprintf('%s%s/full/%s_messCurrent__%s.png', maxent_path, species, species, name)
        
        png(f, 8, 8, units = 'in', res = 300, type = 'cairo')
        print(p)
        dev.off()
        
      }, unstack(mess_current$similarity), names(mess_current$similarity))
      
      writeRaster(novel_current, sprintf('%s%s/full/%s_%s.tif', 
                                         maxent_path, species, species, "current_novel_map"), 
                  overwrite = TRUE, datatype = 'INT2S')
      
      hs_current_notNovel <- hs_current * is.na(novel_current) 
      
      ##################################################################
      # mask out novel environments 
      # is.na(novel_current) is a binary layer showing 
      # not novel [=1] vs novel [=0], 
      # so multiplying with hs_current will mask out novel
      writeRaster(hs_current_notNovel, sub('\\.tif', '_notNovel.tif', hs_current@file@name), 
                  overwrite = TRUE, datatype = 'INT2S')
      
    } else {
      
      message(species, ' skipped - MESS map already run') 
      
    }
    
  }, species_list, threshold_list, SIMPLIFY = FALSE)
  
}





#########################################################################################################################
## For each scenario, calcualte the mean annual temperature and annual rainfall anomaly
calculate.anomaly = function(scen_list, time_slice, climate_path) {
  
  #########################################################################################################################
  ## Create current rasters :: can the raster creation happen just once?
  ## And read in the current data. Can the specific "current" be removed withouth makin a seconf path argument?
  message('current rasters / 10')
  
  env.grids.current       = stack(file.path(sprintf('%s/current/bio_%02d.tif', climate_path, 1:19)))
  env.grids.current[[1]]  = env.grids.current[[1]]/10
  current.bio1            = env.grids.current[[1]]
  current.bio12           = env.grids.current[[12]]
  
  ## First, run a loop over each scenario :: 
  lapply(scen_list, function(x) {
    
    ## Assign the scenario name (to use later in the plot)
    scen_name = eval(parse(text = sprintf('gcms.%s$GCM[gcms.%s$id == x]', time_slice, time_slice)))
    
    #########################################################################################################################
    ## Create a raster stack for each 2050 GCM - also an empty raster for the final plot
    s <- stack(sprintf('%s/20%s/%s/%s%s.tif', climate_path, time_slice, x, x, 1:19))
    
    ########################################################################################################################
    ## Create future rasters
    message('20', time_slice, ' rasters / 10')
    
    s[[1]]       = s[[1]]/10
    future.bio1  = s[[1]]
    future.bio12 = s[[12]]
    
    ## Now subtract the current from the future...
    temp.anomaly = future.bio1  - current.bio1
    rain.anomaly = future.bio12 - current.bio12
    
    plot(temp.anomaly, main = paste0("BIO1  anomaly ", x))
    plot(rain.anomaly, main = paste0("BIO12 anomaly ", x))
    
    ########################################################################################################################
    ## Write the rasters for each species/threshold
    message('Writing BIO1/12 anomaly rasters for 20', time_slice, ' ', x)
    
    if(!file.exists(sprintf('%s/anomalies/%s_%s.tif', climate_path, "BIO1_anomaly",  x))) {
      
      ## Temp anomalies
      writeRaster(temp.anomaly, sprintf('%s/anomalies/%s_%s.tif', climate_path, "BIO1_anomaly",  x), overwrite = TRUE)
      
    } else {
      
      message('Skipped ', "BIO1 anomaly for ", x, ', already exists')
      
    }
    
    if(!file.exists(sprintf('%s/anomalies/%s_%s.tif', climate_path, "BIO12_anomaly",  x))) {
      
      ## Rain anomalies
      writeRaster(rain.anomaly, sprintf('%s/anomalies/%s_%s.tif', climate_path, "BIO12_anomaly", x), overwrite = TRUE)
      
    } else {
      
      message('Skipped ', "BIO12 anomaly for ", x, ', already exists')
      
    }
    
  })
  
} 





#########################################################################################################################
## For each scenario, calcualte the mean annual temperature and annual rainfall anomaly
gcm.value = function(scen_list, time_slice, climate_path) {
  
  #########################################################################################################################
  ## Create current rasters :: can the raster creation happen just once?
  ## And read in the current data. Can the specific "current" be removed withouth makin a seconf path argument?
  message('current rasters divided by 10')
  
  env.grids.current       = stack(file.path(sprintf('%s/current/bio_%02d.tif', climate_path, 1:19)))
  env.grids.current[[1]]  = env.grids.current[[1]]/10
  current.bio1            = env.grids.current[[1]]
  current.bio12           = env.grids.current[[12]]
  
  ## First, run a loop over each scenario :: 
  lapply(scen_list, function(x) {
    
    ## Assign the scenario name (to use later in the plot)
    scen_name = eval(parse(text = sprintf('gcms.%s$GCM[gcms.%s$id == x]', time_slice, time_slice)))
    
    #########################################################################################################################
    ## Create a raster stack for each 2050 GCM - also an empty raster for the final plot
    s <- stack(sprintf('%s/20%s/%s/%s%s.tif', climate_path, time_slice, x, x, 1:19))
    
    ########################################################################################################################
    ## Create future rasters
    message('20', time_slice, ' rasters divided by 10')
    
    s[[1]]       = s[[1]]/10
    future.bio1  = s[[1]]
    future.bio12 = s[[12]]
    
    ## Now subtract the current from the future...
    temp.anomaly = future.bio1  - current.bio1
    rain.anomaly = future.bio12 - current.bio12
    
    ## Also can we simply indicate the range or rainfall and temperautre values for each ?
    #temp.current = 
    
    plot(temp.anomaly, main = paste0("BIO1  anomaly ", x))
    plot(rain.anomaly, main = paste0("BIO12 anomaly ", x))
    
    ########################################################################################################################
    ## Write the rasters for each species/threshold
    message('Writing BIO1/12 anomaly rasters for 20', time_slice, ' ', x)
    
    if(!file.exists(sprintf('%s/anomalies/%s_%s.tif', climate_path, "BIO1_anomaly",  x))) {
      
      ## Temp anomalies
      writeRaster(temp.anomaly, sprintf('%s/anomalies/%s_%s.tif', climate_path, "BIO1_anomaly",  x), overwrite = TRUE)
      
    } else {
      
      message('Skipped ', "BIO1 anomaly for ", x, ', already exists')
      
    }
    
    if(!file.exists(sprintf('%s/anomalies/%s_%s.tif', climate_path, "BIO12_anomaly",  x))) {
      
      ## Rain anomalies
      writeRaster(rain.anomaly, sprintf('%s/anomalies/%s_%s.tif', climate_path, "BIO12_anomaly", x), overwrite = TRUE)
      
    } else {
      
      message('Skipped ', "BIO12 anomaly for ", x, ', already exists')
      
    }
    
  })
  
} 





#########################################################################################################################
## STACK FUNCTIONS
#########################################################################################################################


# dist_change_binary <- function(from, to) {
#   
#   ## Create a raster stack
#   d <- stack(from, to)[]
#   r <- raster(from)
#   
#   ## Classify the raster stack to make each value (i.e. outcome) unique
#   r[d[, 1]==1 & d[, 2]==0] <- 1  ## 1 in current raster and 0 in future = LOSS
#   r[d[, 1]==0 & d[, 2]==1] <- 2  ## 0 in current raster and 1 in future = GAIN
#   r[d[, 1]==1 & d[, 2]==1] <- 3  ## 1 in current raster and 1 in future = STABLE
#   
#   r_f <- as.factor(r)
#   
#   levels(r_f)[[1]] <- data.frame(ID=1:3, label=c('Lost', 'Gained', 'Stable'))
#   
#   ## Now create a level plot
#   png(sprintf('%s.png', time_slice), 6, 6, units='in', res=300)
#   
#   p <- levelplot(r_f, col.regions=c('purple', 'red', 'green'), scales = list(draw = FALSE)) +
#     layer(sp.polygons(aus))
#   print(p)
#   dev.off()
#   
#   #writeRaster(r_f, sprintf('%s.tif', ), datatype='INT2U', overwrite=T)
#   
#   r_f
#   
# }


# ## Need an empty frame
# print(levelplot(stack(empty,                ## needs to have a different colour scale,
#                       combo_suit_percent,
#                       combo_suit_thresh), margin = FALSE,
#                 
#                 ## Create a colour scheme using colbrewer: 100 is to make it continuos
#                 ## Also, make it a one-directional colour scheme
#                 scales      = list(draw = FALSE), 
#                 at = seq(0, 6, length = 100),
#                 col.regions = colorRampPalette(rev(brewer.pal(11, 'Spectral'))),
#                 
#                 ## Give each plot a name
#                 names.attr = c('Aus occurrences', 
#                                sprintf('20%s GCM 10thp train omission > %s', time_slice, percent), 
#                                sprintf('20%s GCM Max train logis > %s',      time_slice, thresh)),
#                 
#                 colorkey   = list(height = 0.5, width = 3), xlab = '', ylab = '',
#                 main       = list(gsub('_', ' ', species), font = 4, cex = 2)) +
#         
#         ## Plot the Aus shapefile with the occurrence points for reference
#         ## Can we assign different shapefiles to different panels, rather than to them all?
#         
#         layer(sp.polygons(aus)) + ## sp.polygons(aus))
#         layer(sp.points(occ, pch = 20, cex = 0.4, 
#                         col = c('red', 'transparent', 'transparent')[panel.number()]), data = list(occ = occ)))
# 
# ## Finish the device
# dev.off()


#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################
