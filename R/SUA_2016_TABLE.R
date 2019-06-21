#########################################################################################################################
## Loop over directories, species and one threshold for each, also taking a time_slice argument. Next, make the lists generic too
## pervious version in R/old/model_combine.R
SUA_table = function(DIR_list, species_list, maxent_path, thresholds, percentiles, time_slice, area_occ) {
  
  ## How can the shapefiles be read in once, not for each species?...................................................
  
  ###################################################################################################################
  ## Read in shapefiles :: this should be done outside the loop
  aus        = readRDS("./data/base/CONTEXTUAL/aus_states.rds") %>%
    spTransform(ALB.CONICAL)
  
  LAND       = readRDS("./data/base/CONTEXTUAL/LAND_world.rds")
  
  areal_unit = readRDS("./data/base/CONTEXTUAL/SUA/SUA_2016_AUST.rds") %>%
    spTransform(ALB.CONICAL)
  
  areal_unit = areal_unit[order(areal_unit$SUA_NAME16),]
  
  ## Loop over each directory
  lapply(DIR_list, function(DIR) { 
    
    ## And each species - although we don't want all possible combinations. use mapply in the function call
    lapply(species_list, function(species) {
    
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
        
      } else {
        
        ## Create another level without the mean calculation
        raster.list = list.files(as.character(DIR), pattern = sprintf('bi%s.tif$', time_slice), full.names = TRUE)  
        suit        = stack(raster.list)
        suit.list   = unstack(suit)
        
      }
      
      #########################################################################################################################
      ## Then, create rasters that meet habitat suitability criteria thresholds, determined by the rmaxent function
      for (thresh in thresholds) {
        
        ## Can we set 0 to NA in the rasters before running the calculations?
        ## Check if the combined suitability raster exists
        f_max_train_suit <- sprintf('%s/%s/full/%s_20%s%s%s.tif', maxent_path,
                                    species, species, time_slice, "_Max_train_sensit_above_", thresh)
      
        
        ## Print the species being analysed
        message('doing ', species, ' | Max train sensit > ', thresh, ' for 20', time_slice)
        
        ## Read in the current suitability raster
        f_current <- raster(sprintf('%s/%s/full/%s_current.tif', 
                                    maxent_path, species, species))
        
        ## First, create a simple function to threshold each of the rasters in raster.list
        thresh_greater  = function (x) {x > thresh}
        
        ## Then apply this to just the current suitability raster. These functions use the : 
        ## Maximum training sensitivity plus specificity Logistic threshold
        ## 10th percentile training presence training omission
        current_suit_thresh  = thresh_greater(f_current)
        current_suit_percent = percent_greater(f_current) 
        
        ## Now, apply these functions to the list of rasters (6 GCMs) for each species
        ## Also, the 'init' function initializes a raster object with values
        ## suit_ras_greater    = reduce(suit.list, thresh_above_fun, .init = suit.list[[1]] > thresh)
        
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
        gain_loss_df = head(gain_loss_df, 3)
        head(gain_loss_df)
        
        ## Save the gain/loss table
        write.csv(gain_loss_df, sprintf('%s/%s/full/%s_20%s_%s%s.csv', maxent_path,
                                        species, species, time_slice, "gain_loss_table_", thresh), row.names = FALSE)
        
        ## Plot rasters to check
        plot(current_suit_thresh,   main = gsub('_', ' ', (sprintf('%s current suit > %s',             species, thresh))))
        plot(combo_suit_thresh,     main = gsub('_', ' ', (sprintf('%s future suit > %s',              species, thresh))))
        plot(combo_suit_4GCM,       main = gsub('_', ' ', (sprintf('%s 4+ GCMs > %s',                  species, thresh))))
        plot(gain_loss,             main = gsub('_', ' ', (sprintf('%s 4+ GCMs > %s plus current',     species, thresh))))
        
        #########################################################################################################################
        ## Then using this GCM consensus, calculate whether the species is likely to be present in each SUA.
        ## Decide on a threshold of % area (10?) of the SUA that needs to be occupied, for each species to be considered present. 
        message('Running SUA area calculations for ', species, ' | 20', time_slice, ' combined suitability > ', thresh)
        
        #########################################################################################################################
        ## Then, extract the values of the presence raster for each areal unit: generates a list
        ext.current  <- raster::extract(current_suit_thresh, areal_unit, method = 'simple')
        ext.future   <- raster::extract(combo_suit_4GCM,     areal_unit, method = 'simple')
        
        ## A function to tabulate the raster values by aerial unit, returning a data frame
        tabFunc <- function(indx, extracted, region, regname) {
          
          dat<-as.data.frame(table(extracted[[indx]]))
          dat$name<-region[[regname]][[indx]]
          return(dat)
          
        }
        
        ## Run through each areal unit and calculate a table of the count of raster cells
        tabs.current <- lapply(seq(ext.current), tabFunc, ext.current, areal_unit, "SUA_NAME16")
        tabs.current <- do.call("rbind", tabs.current )
        
        tabs.future  <- lapply(seq(ext.future), tabFunc, ext.future, areal_unit,   "SUA_NAME16")
        tabs.future  <- do.call("rbind", tabs.future )
        
        #########################################################################################################################
        ## Now mutate the current table
        PERECENT.AREA.CURRENT <- tabs.current %>%
          
          group_by(name) %>%                                          ## Group by region
          mutate(totcells = sum(Freq),                                ## How many cells overall?
                 percent.area = round(100 * Freq / totcells, 2)) %>%  ## Cells divided by total cells
          
          dplyr::select(-c(Freq, totcells)) %>%                       ## There is a select func in raster so need to specify
          spread(key = Var1, value = percent.area, fill = 0)     %>%  ## Make wide format
          as.data.frame()
        
        #########################################################################################################################
        ## Now mutate the future table
        PERECENT.AREA.FUTURE <- tabs.future %>%
          
          group_by(name) %>%                                          ## Group by region
          mutate(totcells = sum(Freq),                                ## How many cells overall?
                 percent.area = round(100 * Freq / totcells, 2)) %>%  ## Cells divided by total cells
          
          dplyr::select(-c(Freq, totcells)) %>%                       ## There is a select func in raster so need to specify
          spread(key = Var1, value = percent.area, fill = 0)     %>%  ## Make wide format
          as.data.frame()
        
        ## Rename and create a column for whether or not the species occupies each SUA according to an area threshold (area_occ)
        ## If both the current and future suitablility rasters are in the SUA
        ## Then name all three columns :: 0[col 1] = absent, 1[col 2] = present
        names(PERECENT.AREA.CURRENT) =  c('SUA_NAME16', 'Absent', 'Present') 
        names(PERECENT.AREA.FUTURE)  =  c('SUA_NAME16', 'Absent', 'Present')   ## has 0 length for some species
        
        PERECENT.AREA.CURRENT$species_present = ifelse(PERECENT.AREA.CURRENT$Present >= area_occ, 1, 0)
        PERECENT.AREA.FUTURE$species_present  = ifelse(PERECENT.AREA.FUTURE$Present  >= area_occ, 1, 0)
        
        head(PERECENT.AREA.CURRENT)
        head(PERECENT.AREA.FUTURE)
        
        ## Create a table with the columns: REGION, AREA SPECIES, STATS (for each time slice)
        ## ifelse(<condition>, <yes>, ifelse(<condition>, <yes>, <no>))
        GCM.AREA.SUMMARY <- data.frame(SUA          = areal_unit$SUA_NAME16, 
                                       AREA_SQKM    = areal_unit$AREA_SQKM,
                                       SPECIES      = species,
                                       PERIOD       = time_slice,
                                       AREA_THRESH  = area_occ,
                                       MAX_TRAIN    = thresh,
                                       CURRENT_AREA = PERECENT.AREA.CURRENT$Present,
                                       FUTURE_AREA  = PERECENT.AREA.FUTURE$Present,
                                       AREA_CHANGE  = PERECENT.AREA.FUTURE$Present - PERECENT.AREA.CURRENT$Present,
                                       PRESENT      = PERECENT.AREA.FUTURE$species_present)
        
        
        ## Now calculate :
        ## Gain (ie not suitable now but suitable in future)
        ## Loss (ie suitable now but not in future)
        ## Stable (suitable in both)
        GCM.AREA.SUMMARY$GAIN_LOSS  <- with(GCM.AREA.SUMMARY, ifelse(
          GCM.AREA.SUMMARY$AREA_CHANGE <= -1 & GCM.AREA.SUMMARY$AREA_CHANGE < 0, 'LOSS', ifelse(
            GCM.AREA.SUMMARY$AREA_CHANGE >= 1, 'GAIN', ifelse(
              GCM.AREA.SUMMARY$AREA_CHANGE == 0, 'STABLE', 'GAIN'))))
        
        ##
        GCM.AREA.SUMMARY$GAIN_LOSS 
        View(GCM.AREA.SUMMARY) 
        
        #########################################################################################################################
        ## Then save the table of SUA results for all species to a datafile...
        ## This would be the file to loop over to create a summary of species per SUA
        write.csv(GCM.AREA.SUMMARY, sprintf('%s/%s/full/%s_20%s_%s%s%s.csv', maxent_path,
                                            species, species, time_slice, area_occ, "pc_area_SUA_summary_", thresh), row.names = FALSE)
        
      }
    
    })
  
  })
  
}

