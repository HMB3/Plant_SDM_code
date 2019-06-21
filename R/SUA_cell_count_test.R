#########################################################################################################################
## Loop over directories, species and one threshold for each, also taking a time_slice argument
SUA_cell_count = function(unit_path, unit_shp, unit_vec, aus_shp, world_shp,
                          DIR_list, species_list, 
                          maxent_path, thresholds, 
                          time_slice, write_rasters, nclust) {
  
  ###################################################################################################################
  ## Read in shapefiles: clunky, but how else will can you read in shapefiles as arguments?  
  areal_unit = readRDS(paste0(unit_path,  unit_shp)) %>%
    spTransform(ALB.CONICAL)
  
  areal_unit      = areal_unit[order(areal_unit$SUA_NAME16),] 
  areal_unit_vec  = readRDS(paste0(unit_path, unit_vec))
  
  aus_poly = readRDS(paste0(unit_path, aus_shp)) %>%
    spTransform(ALB.CONICAL)
  
  world_poly = readRDS(paste0(unit_path, world_shp)) %>%
    spTransform(CRS.WGS.84)
  
  ###################################################################################################################
  ## Loop over each directory
  lapply(DIR_list, function(DIR) { 
    
    ## And each species - although we don't want all possible combinations. use mapply in the function call
    lapply(species_list, function(species) {
      
      ###################################################################################################################
      ## Create a list of the rasters in each species directory for each time period, then take the mean
      message('Running summary of SDM predictions within SUAs for ', species, ' using ', names(areal_unit)[1], " shapefile")
      message('Calcualting mean of 20', time_slice, ' GCMs for ', species)
      
      ## Check if the mean GCM raster exists
      f_mean = sprintf('%s/%s/full/%s_20%s_suitability_mean.tif', maxent_path, species, species, time_slice)
      
      ## The mean of the GCMs doesn't exist, create it
      if(!file.exists(f_mean)) { 
        
        ## Read in all the habitat suitability rasters for each time period which are _not_ novel
        #raster.list       = list.files(as.character(DIR), pattern = sprintf('bi%s.tif$', time_slice), full.names = TRUE)  
        raster.list       = list.files(as.character(DIR), pattern = 'future_not_novel', full.names = TRUE) 
        suit              = stack(raster.list)
        suit.list         = unstack(suit)
        combo_suit_mean   = mean(suit)
        
        writeRaster(combo_suit_mean , sprintf('%s/%s/full/%s_20%s_suitability_mean.tif', 
                                              maxent_path, species, species, time_slice), overwrite = TRUE)
        
      } else {
        
        ## Create another level, without the mean calculation
        raster.list = list.files(as.character(DIR), pattern = sprintf('bi%s.tif$', time_slice), full.names = TRUE)  
        suit        = stack(raster.list)
        suit.list   = unstack(suit)
        
      }
      
      #########################################################################################################################
      ## Then, create rasters that meet habitat suitability criteria thresholds, determined by the rmaxent function
      #for (thresh in thresholds) {
      sua_combine_fun <- function(thresh) {
        
        ## Check if the SUA summary table exists
        SUA_file =   sprintf('%s/%s/full/%s_20%s_%s%s.tif', maxent_path,
                             species, species, time_slice, "gain_loss_", thresh)
        
        ## If the gain/loss raster already exists, skip the calculation
        if(file.exists(SUA_file)) { 
          
          message(species, ' ', 20, time_slice, ' SUA aggregation already exists, skip')
          next
          
        }
        
        ## Print the species being analysed
        message('doing ', species, ' | Logistic > ', thresh, ' for 20', time_slice)
        
        ## Read in the current suitability raster :: get the current_not_novel raster
        f_current <- raster(sprintf('%s/%s/full/%s_current_not_novel.tif', 
                                    maxent_path, species, species))
        
        ## First, create a simple function to threshold each of the rasters in raster.list,
        ## Then apply this to just the current suitability raster. 
        thresh_greater       = function (x) {x > thresh}
        current_suit_thresh  = thresh_greater(f_current)
        
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
        
        #########################################################################################################################
        ## Then sum them up: All the threshholds
        combo_suit_thresh   =  Reduce("+", list(suit_ras1_thresh, suit_ras2_thresh, suit_ras3_thresh,
                                                suit_ras4_thresh, suit_ras5_thresh, suit_ras6_thresh))
        
        #########################################################################################################################
        ## For each species, create a binary raster with cells > 4 GCMs above the maxent threshold = 1, and cells with < 4 GCMs = 0. 
        message('Calculating change for ', species, ' | 20', time_slice, ' combined suitability > ', thresh)
        
        ## Functions for thresholding rasters
        band_4           <- function(x) {ifelse(x >=  4, 1, 0) }
        combo_suit_4GCM  <- calc(combo_suit_thresh, fun = band_4)
        
        #########################################################################################################################
        ## Now create a raster of the gain, loss and stable
        ## Create a raster stack of the current and future rasters
        message ("Counting cells lost/gained/stable/never suitable, both across AUS and per SUA")
        
        #########################################################################################################################
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
        
        #########################################################################################################################
        ## Now calculate the number of cells lost/gained/stable across Australia
        message ("Counting cells lost/gained/stable/never suitable across Australia")
        d5 <- stack(current_suit_thresh, combo_suit_4GCM)[]
        r  <- raster(current_suit_thresh)                     ## rename, too cryptic.............................................
        z  <- as.data.frame(d4)                               ## rename, too cryptic.............................................
        
        ## Then classify the raster stack to make each value (i.e. outcome) unique
        r[d5[, 1]==1 & d5[, 2]==0] <- 1  ## 1 in current raster and 0 in future = LOSS
        r[d5[, 1]==0 & d5[, 2]==1] <- 2  ## 0 in current raster and 1 in future = GAIN
        r[d5[, 1]==1 & d5[, 2]==1] <- 3  ## 1 in current raster and 1 in future = STABLE
        r[d5[, 1]==0 & d5[, 2]==0] <- 4  ## 0 in current raster and 0 in future = NEVER_SUIT
        
        ## Create a table of these values, to merge with the levels later. This avoids the problem that not all the categories will
        ## be present for all species
        change_values = data.frame("ID" = 1:4, "CHANGE" = c("LOST", "GAINED", "STABLE", "NEVER"))
        
        
        ## Now convert the raster to a factor and assign lables to the levels
        gain_loss <- as.factor(r)
        levels(gain_loss)[[1]] <- data.frame(ID = 1:4, label = c('Lost', 'Gained', 'Stable', 'Never_Suitable'))
        z <- as.data.frame(d5)
        
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
        ## Save the gain/loss table for the whole of Australia
        message('Writing ', species, ' gain_loss tables for 20', time_slice)
        write.csv(gain_loss_df, sprintf('%s/%s/full/%s_20%s_%s%s.csv', maxent_path,
                                        species, species, time_slice, "gain_loss_table_", thresh), row.names = FALSE)
        
        ## Save the gain/loss table
        write.csv(d4, sprintf('%s/%s/full/%s_20%s_%s%s.csv', maxent_path,
                              species, species, time_slice, "SUA_cell_count_", thresh), row.names = FALSE)
        
        #########################################################################################################################
        ## Now write the rasters
        ## If the rasters don't exist, write them for each species/threshold
        if(write_rasters == "TRUE") {
          
          ## Write the current suitability raster, thresholded using the Maximum training 
          ## sensitivity plus specificity Logistic threshold
          message('Writing ', species, ' current', ' max train > ', thresh)
          writeRaster(current_suit_thresh, sprintf('%s/%s/full/%s_%s%s.tif', maxent_path,
                                                   species, species, "current_suit_not_novel_above_", thresh), 
                      overwrite = TRUE)
          
          ## Write the combined suitability raster, thresholded using the maximum training value
          message('Writing ', species, ' | 20', time_slice, ' max train > ', thresh)
          writeRaster(combo_suit_thresh, sprintf('%s/%s/full/%s_20%s%s%s.tif', maxent_path,
                                                 species, species, time_slice, "_log_thresh_above_", thresh), 
                      overwrite = TRUE)
          
          ## Write the combined future raster with > 4 GCMs above the maximum training value
          message('Writing ', species, ' | 20', time_slice, ' 4 GCMs > ', thresh)
          writeRaster(combo_suit_4GCM, sprintf('%s/%s/full/%s_20%s%s%s.tif', maxent_path,
                                               species, species, time_slice, "_4GCMs_above_", thresh), 
                      overwrite = TRUE)
          
          ## Write out the gain/loss raster
          writeRaster(gain_loss, sprintf('%s/%s/full/%s_20%s%s%s.tif', maxent_path,
                                         species, species, time_slice, "_gain_loss_", thresh), datatype = 'INT2U', 
                      overwrite = TRUE)
          
          #########################################################################################################################
          #########################################################################################################################
          ## Create a color scheme for the gain_loss plot
          SUA.plot.cols = brewer.pal(12, "Paired")
          ## Create dataframe of colors that match the categories :;
          # 1 <- "LOST"        SUA.plot.cols[8]
          # 2 <- "GAINED"      SUA.plot.cols[4]
          # 3 <- "STABLE"      SUA.plot.cols[1]
          # 4 <- "NEVER_SUIT"  SUA.plot.cols[9]
          ## 1 = STABLE, 4 = GAIN, 8 = LOSS, 9 = NEVER
          lc_id = c(1, 2, 3, 4)
          cover_palette <- c(SUA.plot.cols[8], SUA.plot.cols[4], SUA.plot.cols[1], SUA.plot.cols[11])
          
          colors        <- data.frame(id = lc_id, color = cover_palette, stringsAsFactors = FALSE)
          csort         <- colors[order(colors$id),] 
          
          ## Create Labels for the gain_loss plots 
          gain_plot <- ratify(gain_loss)
          rat <- levels(gain_plot)[[1]]
          
          ## Use an inner join, to accomodate the different categories. EG sometimes, there are 
          ## no cells which are gained for each species, etc.
          rat[["CHANGE"]]   <- join(rat, change_values, type = "inner")$CHANGE 
          levels(gain_plot) <- rat
          
          
          #########################################################################################################################
          ## Save the gain/loss raster to PNG
          save_name = gsub(' ', '_', species)
          identical(projection(aus_poly), projection(gain_plot))
          
          message('writing gain/loss png for ', 'species')
          png(sprintf('%s/%s/full/%s_%s_%s_20%s.png', maxent_path, save_name, save_name, "gain_loss", thresh, time_slice),
              16, 10, units = 'in', res = 500)
          
          ## Could add the SUA polygons as well
          print(levelplot(gain_plot, 
                          col.regions = csort$color, 
                          xlab = NULL, ylab = NULL,
                          main       = list(paste0(gsub('_', ' ', species), ' :: ',  20, 
                                                   time_slice, ' 4GCMs > ',  thresh), font = 4, cex = 2)) +
                  latticeExtra::layer(sp.polygons(aus_poly), data = list(aus_poly = aus_poly)))
          
          dev.off()
          
        } else {
          
          message(' skip raster writing')
          
        }
        
      }
      
      if (nclust==1) {
        
        lapply(thresholds, sua_combine_fun) 
        
      } else {
        
        cl <- makeCluster(nclust)
        clusterExport(cl, c(
          
          ## The number of objects that need exporting is full on - maybe there's a better way
          'unit_path',  'unit_shp',      'unit_vec',    'aus_shp', 'world_shp',
          'DIR_list',   'species_list',  'maxent_path', 'thresholds', 
          'time_slice', 'write_rasters', 'nclust'), envir = environment())

        clusterEvalQ(cl, {
          
          library(rmaxent)
          library(sp)
          library(raster)
          library(rasterVis)
          library(latticeExtra)
          library(magrittr)
          
        })
        
        message('Running project_maxent_grids_mess for ', length(species_list),
                ' species on ', nclust, ' cores for GCM ', x)
        
      }
      
    })
    
  })
  
}



