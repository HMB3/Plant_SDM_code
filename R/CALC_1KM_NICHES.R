#########################################################################################################################
## 3). INTERSECT SPECIES RECORDS WITH SIGNIFICANT URBAN AREAS
#########################################################################################################################


#########################################################################################################################
## Only run this code if we are saving niches 
if(calc_niche == "TRUE") {
  
  
  #########################################################################################################################
  ## Create a spatial points object
  COMBO.RASTER.1KM  <- CLEAN.INV
  
  
  #########################################################################################################################
  ## Create a spatial points object, and change to a projected system to calculate distance more accurately 
  coordinates(COMBO.RASTER.1KM)    <- ~lon+lat
  proj4string(COMBO.RASTER.1KM)    <- '+init=epsg:4326'
  COMBO.RASTER.1KM                 <- spTransform(COMBO.RASTER.1KM, CRS(sp_epsg54009))
  
  
  ## Now split using the data using the species column, and get the unique occurrence cells
  COMBO.RASTER.SPLIT.1KM <- split(COMBO.RASTER.1KM, COMBO.RASTER.1KM$searchTaxon)
  occurrence_cells_all   <- lapply(COMBO.RASTER.SPLIT.1KM, function(x) cellFromXY(template.raster.1km, x))
  length(occurrence_cells_all)   ## this is a list of dataframes, where the number of rows for each being the species table
  
  
  #########################################################################################################################
  ## Now get just one record within each 10*10km cell.
  ## Split the inventory data, and subsample inventory data to coarseer resoltion
  NICHE.1KM <- mapply(function(x, cells) {
    x[!duplicated(cells), ]
  }, COMBO.RASTER.SPLIT.1KM, occurrence_cells_all, SIMPLIFY = FALSE) %>% do.call(rbind, .)
  
  
  ## We want to know the count of species that occur in 'n' LGAs, across a range of climates. Read in LGA and SUA
  message(round(nrow(NICHE.1KM)/nrow(CLEAN.INV)*100, 2), " % records retained at 1km resolution")  
  length(unique(NICHE.1KM$searchTaxon))
  projection(LGA);projection(AUS);projection(SUA_2016)
  
  
  ## Use a projected, rather than geographic, coordinate system
  ## Not sure why, but this is needed for the spatial overlay step
  NICHE.1KM.84 = spTransform(NICHE.1KM, CRS.WGS.84)
  AUS.WGS      = spTransform(AUS,       CRS.WGS.84)
  POA.WGS      = spTransform(POA_2016,  CRS.WGS.84)
  SUA.WGS      = spTransform(SUA_2016 , CRS.WGS.84)
  LAND.WGS     = spTransform(LAND,      CRS.WGS.84)
  
  
  ## Create global niche and Australian niche.
  ## So we need a subset for Australia 
  NICHE.AUS <-  NICHE.1KM.84[AUS.WGS, ] 
  
  
  #########################################################################################################################
  ## Run join between species records and LGAs/SUAs :: Double check they are the same.
  ## See the ABS for details :: there are 563 LGAs
  ## http://www.abs.gov.au/ausstats/abs@.nsf/Lookup/by%20Subject/1270.0.55.003~July%202016~Main%20Features~Local%20Government%20Areas%20(LGA)~7
  message('Joining occurence data to SUAs for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")
  projection(NICHE.1KM.84);projection(SUA.WGS);projection(AUS.WGS);projection(POA.WGS)
  SUA.JOIN      = over(NICHE.1KM.84, SUA.WGS)
  POA.JOIN      = over(NICHE.1KM.84, POA.WGS)
  
  
  ## What are the lengths?
  length(unique(SUA.JOIN$SUA_NAME16))
  length(unique(POA.JOIN$POA_NAME16))
  
  
  ## Bind the SUA and POA data to the occurence data
  COMBO.SUA.POA = cbind.data.frame(NICHE.1KM.84, SUA.JOIN, POA.JOIN) 
  
  
  #########################################################################################################################
  ## Aggregate the number of SUAs and POAs each species is found in
  SUA.AGG   = tapply(COMBO.SUA.POA$SUA_NAME16, COMBO.SUA.POA$searchTaxon, function(x) length(unique(x))) ## group SUA by species name
  POA.AGG   = tapply(COMBO.SUA.POA$POA_NAME16, COMBO.SUA.POA$searchTaxon, function(x) length(unique(x))) ## group POA by species name
  SUA.AGG   = as.data.frame(SUA.AGG)
  POA.AGG   = as.data.frame(POA.AGG)
  
  SUA.AGG <- setDT(SUA.AGG, keep.rownames = TRUE)[]
  POA.AGG <- setDT(POA.AGG, keep.rownames = TRUE)[]
  names(SUA.AGG) = c("searchTaxon", "SUA_COUNT")
  names(POA.AGG) = c("searchTaxon", "POA_COUNT")
  
  
  ## Now create a table of all the SUA's that each species occurrs
  SUA.AGG = join(SUA.AGG, POA.AGG)
  saveRDS(SUA.AGG, paste0(DATA_path, 'SUA_SPP_COUNT_', save_run, '.rds'))
  
  
  
  
  
  #########################################################################################################################
  ## 4). CREATE NICHES FOR PROCESSED TAXA
  #########################################################################################################################
  
  
  #########################################################################################################################
  ## Now summarise the niches. But, figure out a cleaner way of doing this
  message('Estimating global niches for ', length(GBIF.spp), ' species across ', length(env.variables), ' climate variables')
  
  
  #########################################################################################################################
  ## Environmental variables
  
  
  # BIO1  = Annual Mean Temperature                                     ## 
  # BIO2  = Mean Diurnal Range (Mean of monthly (max temp - min temp))  ##
  # BIO3  = Isothermality (BIO2/BIO7) (* 100)
  # BIO4  = Temperature Seasonality (standard deviation *100)           ##
  # BIO5  = Max Temperature of Warmest Month                            ## 
  # BIO6  = Min Temperature of Coldest Month                            ## 
  # BIO7  = Temperature Annual Range (BIO5-BIO6)
  # BIO8  = Mean Temperature of Wettest Quarter
  # BIO9  = Mean Temperature of Driest Quarter
  # BIO10 = Mean Temperature of Warmest Quarter
  # BIO11 = Mean Temperature of Coldest Quarter
  # BIO12 = Annual Precipitation                                        ##
  # BIO13 = Precipitation of Wettest Month                              ##
  # BIO14 = Precipitation of Driest Month                               ##
  # BIO15 = Precipitation Seasonality (Coefficient of Variation)        ##
  # BIO16 = Precipitation of Wettest Quarter
  # BIO17 = Precipitation of Driest Quarter
  # BIO18 = Precipitation of Warmest Quarter
  # BIO19 = Precipitation of Coldest Quarter
  
  
  #########################################################################################################################
  ## Create niche summaries for each environmental condition like this...
  ## Here's what the function will produce :
  NICHE.AUS.DF = as.data.frame(NICHE.AUS)
  NICHE.GLO.DF = as.data.frame(NICHE.1KM.84) 
  
  head(niche_estimate (DF = NICHE.AUS.DF, colname = "Annual_mean_temp"))  ## Including the q05 and q95
  head(niche_estimate (DF = NICHE.GLO.DF, colname = "Annual_mean_temp"))  ## Including the q05 and q95
  
  
  #########################################################################################################################
  message('Estimating global niches for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")
  GLOB.NICHE <- env.variables %>% 
    
    ## Pipe the list into lapply
    lapply(function(x) {
      
      ## Now, use the niche width function on each colname
      ## Also, need to figure out how to make the aggregating column generic (species, genus, etc.)
      ## currently it only works hard-wired..............................................................................
      niche_estimate(DF = NICHE.GLO.DF, colname = x)
      
      ## Would be good to remove the duplicate columns here..............................................................
      
    }) %>% 
    
    ## finally, create one dataframe for all niches
    as.data.frame
  
  
  ## Remove duplicate Taxon columns and check the output :: would be great to skip these columns when running the function
  names(GLOB.NICHE)
  GLOB.NICHE <- subset(GLOB.NICHE, select = -c(searchTaxon.1,  searchTaxon.2,  searchTaxon.3,  searchTaxon.4,
                                               searchTaxon.5,  searchTaxon.6,  searchTaxon.7,  searchTaxon.7,
                                               searchTaxon.8,  searchTaxon.9,  searchTaxon.10, searchTaxon.11,
                                               searchTaxon.12, searchTaxon.13, searchTaxon.14, searchTaxon.15,
                                               searchTaxon.16, searchTaxon.17, searchTaxon.18, searchTaxon.19))
  
  
  #########################################################################################################################
  message('Estimating global niches for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")
  AUS.NICHE <- env.variables %>% 
    
    ## Pipe the list into lapply
    lapply(function(x) {
      
      ## Now, use the niche width function on each colname
      ## Also, need to figure out how to make the aggregating column generic (species, genus, etc.)
      ## currently it only works hard-wired..............................................................................
      niche_estimate(DF = NICHE.AUS.DF, colname = x)
      
      ## Would be good to remove the duplicate columns here..............................................................
      
    }) %>% 
    
    ## finally, create one dataframe for all niches
    as.data.frame
  
  
  ## Remove duplicate Taxon columns and check the output :: would be great to skip these columns when running the function
  names(AUS.NICHE)
  AUS.NICHE <- subset(AUS.NICHE, select = -c(searchTaxon.1,  searchTaxon.2,  searchTaxon.3,  searchTaxon.4,
                                             searchTaxon.5,  searchTaxon.6,  searchTaxon.7,  searchTaxon.7,
                                             searchTaxon.8,  searchTaxon.9,  searchTaxon.10, searchTaxon.11,
                                             searchTaxon.12, searchTaxon.13, searchTaxon.14, searchTaxon.15,
                                             searchTaxon.16, searchTaxon.17, searchTaxon.18, searchTaxon.19))
  
  
  #########################################################################################################################
  ## How are the AUS and GLOB niches related? Many species won't have both Australian and Global niches.
  ## So best to calculate the AUS niche as a separate table. Then, just use the global niche table for the rest of the code
  length(AUS.NICHE$searchTaxon); length(GLOB.NICHE$searchTaxon)
  setdiff(AUS.NICHE$searchTaxon, GLOB.NICHE$searchTaxon)
  setdiff(GLOB.NICHE$searchTaxon, AUS.NICHE$searchTaxon)
  
  
  #########################################################################################################################
  ## Add counts for each species, and record the total number of taxa processed
  GLOBAL_RECORDS = as.data.frame(table(COMBO.SUA.POA$searchTaxon))
  names(GLOBAL_RECORDS) = c("searchTaxon", "GLOBAL_RECORDS")
  identical(nrow(GLOBAL_RECORDS), nrow(GLOB.NICHE))
  
  
  #########################################################################################################################
  ## Add the counts of Australian records for each species to the niche database
  GLOB.NICHE = join(GLOBAL_RECORDS, GLOB.NICHE,  type = "right")
  GLOB.NICHE = join(SUA.AGG,        GLOB.NICHE,  type = "right")
  
  
  ## Check
  head(GLOB.NICHE$GLOBAL_RECORDS)
  head(GLOB.NICHE$SUA_COUNT)
  
  names(COMBO.NICHE)
  dim(COMBO.NICHE)
  
  
  
  
  
  #########################################################################################################################
  ## 5). CALCULATE AREA OF OCCUPANCY
  #########################################################################################################################
  
  
  ## Given a dataframe of georeferenced occurrences of one, or more, taxa, this function provide statistics values 
  ## (Extent of Occurrence, Area of Occupancy, number of locations, number of subpopulations) and provide a preliminary 
  ## conservation status following Criterion B of IUCN.
  
  
  #########################################################################################################################
  ## AREA OF OCCUPANCY (AOO).
  ## For every species in the list: calculate the AOO
  ## x = spp.geo[1]
  spp.geo = as.character(unique(COMBO.SUA.POA$searchTaxon))
  
  GBIF.AOO <- spp.geo %>%
    
    ## Pipe the list into lapply
    lapply(function(x) {
      
      ## Subset the the data frame to calculate area of occupancy according the IUCN.eval
      DF   = subset(COMBO.SUA.POA, searchTaxon == x)[, c("lat", "lon", "searchTaxon")]
      
      message('Calcualting geographic ranges for ', x, ', ', nrow(DF), ' records')
      AOO  = AOO.computing(XY = DF, Cell_size_AOO = 2)  ## Grid size in decimal degrees. Changes the results
      AOO  = as.data.frame(AOO)
      AOO$searchTaxon  <- rownames(AOO)
      rownames(AOO)    <- NULL
      
      ## Remeber to explicitly return the df at the end of loop, so we can bind
      return(AOO)
      
    }) %>%
    
    ## Finally, create one dataframe for all niches
    bind_rows
  
  head(GBIF.AOO)
  
  
  #########################################################################################################################
  ## Now join on the geographic range and glasshouse data
  identical(nrow(GBIF.AOO), length(GBIF.spp))
  GLOB.NICHE = join(GBIF.AOO, GLOB.NICHE, type = "right")
  
  
  #########################################################################################################################
  ## 6). CALCULATE THE ZONAL STATS FOR POSTAL AREAS
  #########################################################################################################################
  
  
  ## This code can be used to re-do the zonal stats for the SUAs in step 9................................................. 
  
  
  #########################################################################################################################
  ## Just use the two environmental conditions likely  to be used for ranges
  POA_temp = aus.grids.current[[1]]
  POA_rain = aus.grids.current[[12]]
  
  ## Make sure the raster extents match with the POA 
  POA_SF <- POA.WGS %>%
    spTransform(., ALB.CON) %>%
    crop(., extent(POA_temp)) %>%
    st_as_sf()
  
  
  #########################################################################################################################
  ## Calculate the mean temp and rain in each POA - doesn't need a loop
  POA_SF$Annual_mean_temp    <- exact_extract(POA_temp, POA_SF, weighted.mean, na.rm = TRUE)
  POA_SF$Annual_precip       <- exact_extract(POA_rain, POA_SF, weighted.mean, na.rm = TRUE)
  
  POA_SF$Annual_mean_temp_30 <- exact_extract(aus.temp.2030, POA_SF, weighted.mean, na.rm = TRUE)
  POA_SF$Annual_precip_30    <- exact_extract(aus.rain.2030, POA_SF, weighted.mean, na.rm = TRUE)
  
  POA_SF$Annual_mean_temp_50 <- exact_extract(aus.temp.2050, POA_SF, weighted.mean, na.rm = TRUE)
  POA_SF$Annual_precip_50    <- exact_extract(aus.rain.2050, POA_SF, weighted.mean, na.rm = TRUE)
  
  POA_SF$Annual_mean_temp_70 <- exact_extract(aus.temp.2070, POA_SF, weighted.mean, na.rm = TRUE)
  POA_SF$Annual_precip_70    <- exact_extract(aus.rain.2070, POA_SF, weighted.mean, na.rm = TRUE)
  
  length(POA_SF$Annual_mean_temp_30);length(POA_SF$Annual_precip_30);
  length(POA_SF$Annual_mean_temp_50);length(POA_SF$Annual_precip_50);
  length(POA_SF$Annual_mean_temp_70);length(POA_SF$Annual_precip_70);
  
  
  
  ## Create a dataframe of the temperature and rainfal
  POA_climate          = as.data.frame(POA_SF)
  POA_climate$geometry = NULL
  head(POA_climate$POA_CODE16)
  
  #########################################################################################################################
  ## How to include the POA? Could add them all 
  POA.SYD = POA_climate[POA_climate$POA_CODE16 %in% 2000 , ] 
  POA.BRI = POA_climate[POA_climate$POA_CODE16 %in% 4000 , ] 
  POA.MEL = POA_climate[POA_climate$POA_CODE16 %in% 3000 , ] 
  POA.PER = POA_climate[POA_climate$POA_CODE16 %in% 6000 , ]
  POA.ADE = POA_climate[POA_climate$POA_CODE16 %in% 5000 , ] 
  POA.DAR = POA_climate[POA_climate$POA_CODE16 %in% 0800 , ] 
  POA.HOB = POA_climate[POA_climate$POA_CODE16 %in% 7000 , ] 
  POA.CAN = POA_climate[POA_climate$POA_CODE16 %in% 2601 , ] 
  
  
  
  
  
  #########################################################################################################################
  ## 7). PLOT HISTOGRAMS AND BAR CHARTS FOR EACH SPECIES AT 1KM
  #########################################################################################################################
  
  
  ##############################################################################################
  ## Plot histograms of temperature and rainfall
  ## species = spp.geo[2]
  for (species in spp.geo) {
    
    ## Subset the spatial dataframe into records for each species
    SP.DF     <- NICHE.1KM.84[NICHE.1KM.84$searchTaxon %in% species , ] 
    
    TMP.GLO   <- subset(GLOB.NICHE,   searchTaxon == species)[c("searchTaxon", "Annual_mean_temp_q95_q05", 
                                                                "Annual_mean_temp_q05", "Annual_mean_temp_q95")]
    
    TMP.AUS   <- subset(AUS.NICHE,    searchTaxon == species)[c("searchTaxon", "Annual_mean_temp_q95_q05", 
                                                                "Annual_mean_temp_q05", "Annual_mean_temp_q95")]
    
    #############################################################
    ## Now, build a df of the temperature vectors 
    if(nrow(TMP.GLO) > 0){
      TMP.GLO$RANGE = "GLOBAL"
    } else {
      message("No global data for ", species)
    }
    
    if(nrow(TMP.AUS) > 0){
      TMP.AUS$RANGE = "AUS"
    } else {
      message("No Australian data for ", species)
    }
    
    TMP.RANGE <- rbind(TMP.GLO, TMP.AUS)
    names(TMP.RANGE)[2] = c("Temperature_range")
    
    ## Subset DF into records for each species
    DF     <- subset(COMBO.SUA.POA, searchTaxon == species)
    DF.OCC <- subset(COMBO.SUA.POA, searchTaxon == species & SOURCE != "INVENTORY")
    DF.INV <- subset(COMBO.SUA.POA, searchTaxon == species & SOURCE == "INVENTORY")
    
    #############################################################
    ## Plot occurrence points by source for the world
    message('Writing global occ sources for ', species)
    png(sprintf("./data/ANALYSIS/SPECIES_RANGES/%s_%s", species, "1km_occ_points_source.png"),
        16, 10, units = 'in', res = 500)
    
    plot(LAND.84, main = paste0("Global points for ", species),
         lwd = 0.01, asp = 1, col = 'grey', bg = 'sky blue')
    
    points(SP.DF,
           pch = ".", cex = 3.3, cex.lab = 3, cex.main = 4, cex.axis = 2,
           xlab = "", ylab = "", asp = 1,
           col = factor(SP.DF$SOURCE))
    
    dev.off()
    
    #############################################################
    ## Plot temperature barchart
    message('Writing global temp histograms for ', species)
    png(sprintf("./data/ANALYSIS/SPECIES_RANGES/%s_%s", species, "temp_barchart_1km_records.png"),
        16, 10, units = 'in', res = 500)
    
    ## Use the 'SOURCE' column to create a histogram for each source.
    max.temp = max(TMP.RANGE$Annual_mean_temp_q95)+1
    min.temp = min(TMP.RANGE$Annual_mean_temp_q05)-1
    
    temp.bar = 
      ggplot(TMP.RANGE, aes(y = Temperature_range, x = RANGE, fill = RANGE)) +
      
      # scale_x_discrete(limits = c(min(TMP.RANGE$Annual_mean_temp_q05), 
      #                             max(TMP.RANGE$Annual_mean_temp_q95))) +
      
      geom_bar(stat = "identity", position = "identity") +
      coord_flip() +
      
      ## Add some median lines : overall, ALA and GBIF
      geom_vline(aes(xintercept = POA.SYD$Annual_mean_temp),
                 col = 'blue', size = 1) +
      geom_vline(aes(xintercept = POA.SYD$Annual_mean_temp_50),
                 col = 'red', size = 1) +
      geom_vline(aes(xintercept = POA.SYD$Annual_mean_temp_70),
                 col = 'green', size = 1) +
      
      
      ggtitle(paste0("Worldclim temperature ranges for ", species)) +
      
      ## Add themes
      theme(axis.title.x     = element_text(colour = "black", size = 35),
            axis.text.x      = element_text(size = 25),
            
            axis.title.y     = element_text(colour = "black", size = 35),
            axis.text.y      = element_text(size = 25),
            
            panel.background = element_blank(),
            panel.border     = element_rect(colour = "black", fill = NA, size = 3),
            plot.title       = element_text(size   = 40, face = "bold"),
            legend.text      = element_text(size   = 20),
            legend.title     = element_text(size   = 20),
            legend.key.size  = unit(1.5, "cm"))
    
    ## Print the plot and close the device 
    print(temp.hist + ggtitle(paste0("Worldclim temp niches for ", species)))
    dev.off()
    
    
    
    
    #############################################################
    ## Plot temperature histograms
    message('Writing global temp histograms for ', species)
    png(sprintf("./data/ANALYSIS/SPECIES_RANGES/%s_%s", species, "temp_niche_histograms_1km_records.png"),
        16, 10, units = 'in', res = 500)
    
    ## Use the 'SOURCE' column to create a histogram for each source.
    temp.hist = ggplot(DF, aes(x = Annual_mean_temp, group = SOURCE, fill = SOURCE)) +
      
      geom_histogram(position = "identity", alpha = 0.5, binwidth = 0.25,
                     aes(y =..density..))  + 
      geom_density(col = 4, alpha = 0.5) +
      
      ## Add some median lines : overall, ALA and GBIF
      geom_vline(aes(xintercept = median(DF.INV$Annual_mean_temp)),
                 col = 'blue', size = 1) +
      geom_vline(aes(xintercept = median(DF.OCC$Annual_mean_temp)),
                 col = 'red', size = 1) +
      
      ggtitle(paste0("Worldclim temp niches for ", species)) +
      
      ## Add themes
      theme(axis.title.x     = element_text(colour = "black", size = 35),
            axis.text.x      = element_text(size = 25),
            
            axis.title.y     = element_text(colour = "black", size = 35),
            axis.text.y      = element_text(size = 25),
            
            panel.background = element_blank(),
            panel.border     = element_rect(colour = "black", fill = NA, size = 3),
            plot.title       = element_text(size   = 40, face = "bold"),
            legend.text      = element_text(size   = 20),
            legend.title     = element_text(size   = 20),
            legend.key.size  = unit(1.5, "cm"))
    
    ## Print the plot and close the device 
    print(temp.hist + ggtitle(paste0("Worldclim temp niches for ", species)))
    dev.off()
    
    #############################################################
    ## Plot rainfall histograms
    message('Writing global rain histograms for ', species)
    png(sprintf("./data/ANALYSIS/SPECIES_RANGES/%s_%s", species, "rain_niche_histograms_1km_records.png"),
        16, 10, units = 'in', res = 500)
    
    ## Use the 'SOURCE' column to create a histogram for each source.
    rain.hist = ggplot(DF, aes(x = Annual_precip, group = SOURCE, fill = SOURCE)) +
      
      geom_histogram(position = "identity", alpha = 0.5, binwidth = 15,
                     aes(y =..density..))  + 
      geom_density(col = 4, alpha = 0.5) +
      
      ## Add some median lines : overall, ALA and GBIF
      geom_vline(aes(xintercept = median(DF.INV$Annual_precip)),
                 col = 'blue', size = 1) +
      geom_vline(aes(xintercept = median(DF.OCC$Annual_precip)),
                 col = 'red', size = 1) +
      
      ggtitle(paste0("Worldclim rain niches for ", species)) +
      
      ## Add themes
      theme(axis.title.x     = element_text(colour = "black", size = 35),
            axis.text.x      = element_text(size = 25),
            
            axis.title.y     = element_text(colour = "black", size = 35),
            axis.text.y      = element_text(size = 25),
            
            panel.background = element_blank(),
            panel.border     = element_rect(colour = "black", fill = NA, size = 3),
            plot.title       = element_text(size   = 40, face = "bold"),
            legend.text      = element_text(size   = 20),
            legend.title     = element_text(size   = 20),
            legend.key.size  = unit(1.5, "cm"))
    
    ## Print the plot and close the device 
    print(rain.hist)
    dev.off()
    
  }
  
  
  #########################################################################################################################
  ## 7). JOIN ON CONTEXTUAL DATA
  #########################################################################################################################
  
  
  #########################################################################################################################
  ## Now join the horticultural contextual data onto one or both tables ()
  message('Joining contextual data for raster and niche files', length(GBIF.spp), ' species in the set ', "'", save_run, "'")
  COMBO.RASTER.CONTEXT = CLEAN.INV
  names(COMBO.RASTER.CONTEXT)
  
  
  #########################################################################################################################
  ## Now join the hort context data and the tree inventory plantings to the niche
  COMBO.NICHE.CONTEXT = join(CLEAN.GROW, COMBO.NICHE,         type = "right") 
  COMBO.NICHE.CONTEXT = join(TI.LUT,     COMBO.NICHE.CONTEXT, type = "right")
  
  
  ## Check the columns are still there
  head(COMBO.NICHE.CONTEXT$AUS_RECORDS)
  head(COMBO.NICHE.CONTEXT$GLOBAL_RECORDS)
  head(COMBO.NICHE.CONTEXT$Plantings)
  head(COMBO.NICHE.CONTEXT$SUA_COUNT)
  
  
  #########################################################################################################################
  ## Is there a relationship between the number of records in each category of records?
  ## Aus records now gets removed, because we are using maxent records
  #windows();
  plot(COMBO.NICHE.CONTEXT$AUS_RECORDS, COMBO.NICHE.CONTEXT$GLOBAL_RECORDS)
  plot(COMBO.NICHE.CONTEXT$AUS_RECORDS, COMBO.NICHE.CONTEXT$Plantings)
  
  
  lm.glob  = lm(COMBO.NICHE.CONTEXT$AUS_RECORDS ~ COMBO.NICHE.CONTEXT$GLOBAL_RECORDS)
  lm.plant = lm(COMBO.NICHE.CONTEXT$AUS_RECORDS ~ COMBO.NICHE.CONTEXT$Plantings)
  
  layout(matrix(c(1,1,2,3), 2, 1, byrow = TRUE))
  
  plot(COMBO.NICHE.CONTEXT$GLOBAL_RECORDS, COMBO.NICHE.CONTEXT$AUS_RECORDS, 
       pch = 19, col  = "blue",
       xlab = "Global records", ylab = "Australian records", 
       abline(lm.glob), 
       main = save_run, cex = 2)
  
  legend("topleft", bty = "n", 
         legend = paste("R2 is", format(summary(lm.glob)$adj.r.squared, digits = 4)))
  
  plot(COMBO.NICHE.CONTEXT$AUS_RECORDS, COMBO.NICHE.CONTEXT$Plantings, 
       pch = 19, col  = "blue",
       xlab = "Aus records", ylab = "Plantings", 
       abline(lm.plant), 
       main = save_run, cex = 2)
  
  legend("topleft", bty = "n", 
         legend = paste("R2 is", format(summary(lm.plant)$adj.r.squared, digits = 4)))
  
  #########################################################################################################################
  ## Now combine the SDM output with the niche context data 
  ## CLEAN.GROW needs to be put through GBIF and TPL........................................................................
  COMBO.NICHE.CONTEXT = COMBO.NICHE.CONTEXT[order(COMBO.NICHE.CONTEXT$searchTaxon),]
  
  
  ## Set NA to blank, then sort by no. of growers
  COMBO.NICHE.CONTEXT$Total_growers[is.na(COMBO.NICHE.CONTEXT$Total.growers)] <- 0
  COMBO.NICHE.CONTEXT = COMBO.NICHE.CONTEXT[with(COMBO.NICHE.CONTEXT, rev(order(Total_growers))), ]
  
  
  ## View the data
  dim(COMBO.RASTER.CONTEXT)
  dim(COMBO.NICHE.CONTEXT)
  
  
  ## Print the dataframe dimensions to screen
  dim(CLEAN.INV)
  dim(COMBO.NICHE.CONTEXT)
  dim(COMBO.RASTER.CONTEXT)
  
  length(unique(CLEAN.INV$searchTaxon))
  length(COMBO.NICHE.CONTEXT$searchTaxon)
  length(unique(COMBO.RASTER.CONTEXT$searchTaxon))
  unique(COMBO.RASTER.CONTEXT$SOURCE)
  
  
  ## Check the column order for the niche and the record tables
  ## Note that for records with catalogue numbers, you can actually search GBIF and ALA for those occurrences
  head(COMBO.NICHE.CONTEXT)[1:12]
  head(COMBO.RASTER.CONTEXT)[1:16]
  
  
  #########################################################################################################################
  ## save .rds file for the next session
  message('Writing 1km resolution niche and raster data for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")
  saveRDS(COMBO.NICHE.CONTEXT,    paste0(DATA_path, 'COMBO_NICHE_CONTEXT_',  OCC_SOURCE, '_1KM_', save_run, '.rds'))
  saveRDS(COMBO.RASTER.CONTEXT,   paste0(DATA_path, 'COMBO_RASTER_CONTEXT_', OCC_SOURCE, '_1KM_', save_run, '.rds'))
  write.csv(COMBO.NICHE.CONTEXT,  paste0(DATA_path, 'COMBO_NICHE_CONTEXT_',  OCC_SOURCE, '_1KM_', save_run, '.csv'), row.names = FALSE)
  
  
} else {
  
  message(' skip file saving, ', length(GBIF.spp), ' species analysed')   ##
  
}