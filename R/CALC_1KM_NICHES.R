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

  
  ##
  unique(SUA.JOIN$SUA_NAME16)
  unique(POA.JOIN$POA_NAME16)
  
  
  ##
  COMBO.SUA.POA = cbind.data.frame(NICHE.1KM.84, SUA.JOIN, POA.JOIN) 
  
  
  #########################################################################################################################
  ## AGGREGATE THE NUMBER OF SUAs EACH SPECIES IS FOUND IN
  SUA.AGG   = tapply(COMBO.SUA.POA$SUA_NAME16, COMBO.SUA.POA$searchTaxon, function(x) length(unique(x))) ## group SUA by species name
  SUA.AGG   = as.data.frame(SUA.AGG)
  
  SUA.AGG <- setDT(SUA.AGG, keep.rownames = TRUE)[]
  names(SUA.AGG) = c("searchTaxon", "SUA_COUNT")

  
  ## Now create a table of all the SUA's that each species occurrs
  SUA.SPP.COUNT = as.data.frame(table(COMBO.SUA.POA[["SUA_NAME16"]], COMBO.SUA.POA[["searchTaxon"]]))
  names(SUA.SPP.COUNT) = c("SUA", "SPECIES", "SUA_RECORDS")
  
  
  ## Save .rds file for the next session
  saveRDS(SUA.SPP.COUNT, paste0(DATA_path, 'SUA_SPP_COUNT_', save_run, '.rds'))
  
  
  ## Check : That's ok, but we want a table of which SUA each species is actually in.
  dim(SUA.AGG)
  head(SUA.AGG)
  
  
  ## Remove duplicate coordinates
  drops <- c("lon.1", "lat.1")
  COMBO.SUA.POA <- COMBO.SUA.POA[ , !(names(COMBO.SUA.POA) %in% drops)]
  names(COMBO.SUA.POA)
  dim(COMBO.SUA.POA)
  str(unique(COMBO.SUA.POA$searchTaxon))
  unique(COMBO.SUA.POA$SOURCE)
  
  
  
  
  
  #########################################################################################################################
  ## 4). CREATE NICHES FOR SELECTED TAXA
  #########################################################################################################################
  
  
  #########################################################################################################################
  ## Now summarise the niches. But, figure out a cleaner way of doing this
  message('Estimating global niches for ', length(GBIF.spp), ' species across ', length(env.variables), ' climate variables')
  
  
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
  names(COMBO.NICHE)
  COMBO.NICHE <- subset(COMBO.NICHE, select = -c(searchTaxon.1,  searchTaxon.2,  searchTaxon.3,  searchTaxon.4,
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
  names(COMBO.NICHE)
  COMBO.NICHE <- subset(COMBO.NICHE, select = -c(searchTaxon.1,  searchTaxon.2,  searchTaxon.3,  searchTaxon.4,
                                                 searchTaxon.5,  searchTaxon.6,  searchTaxon.7,  searchTaxon.7,
                                                 searchTaxon.8,  searchTaxon.9,  searchTaxon.10, searchTaxon.11,
                                                 searchTaxon.12, searchTaxon.13, searchTaxon.14, searchTaxon.15,
                                                 searchTaxon.16, searchTaxon.17, searchTaxon.18, searchTaxon.19))
  
  
  
  #########################################################################################################################
  ## Add counts for each species, and record the total number of taxa processed
  GLOBAL_RECORDS = as.data.frame(table(COMBO.SUA.POA$searchTaxon))
  names(GLOBAL_RECORDS) = c("searchTaxon", "GLOBAL_RECORDS")
  identical(nrow(GLOBAL_RECORDS), nrow(COMBO.NICHE))
  
  
  #########################################################################################################################
  ## Add the counts of Australian records for each species to the niche database
  Total.taxa.processed = nrow(COMBO.NICHE)
  COMBO.NICHE  = join(GLOBAL_RECORDS, COMBO.NICHE, type = "right")
  COMBO.NICHE  = join(SUA.AGG, COMBO.NICHE,        type = "right")
  
  
  head(COMBO.NICHE$AUS_RECORDS)
  head(COMBO.NICHE$SUA_COUNT)
  
  
  names(COMBO.NICHE)
  dim(COMBO.NICHE)
  
  
  
  
  
  #########################################################################################################################
  ## 5). CALCULATE AREA OF OCCUPANCY
  #########################################################################################################################
  
  
  ## Given a dataframe of georeferenced occurrences of one, or more, taxa, this function provide statistics values 
  ## (Extent of Occurrence, Area of Occupancy, number of locations, number of subpopulations) and provide a preliminary 
  ## conservation status following Criterion B of IUCN.
  
  ## The IUCN.eval function needs a data frame, not a spdf
  AUS.AOO.DAT = COMBO.SUA.POA 
  
  
  #########################################################################################################################
  ## AREA OF OCCUPANCY (AOO).
  ## For every species in the list: calculate the AOO
  ## x = spp.geo[1]
  spp.geo = as.character(unique(AUS.AOO.DAT$searchTaxon))
  
  GBIF.AOO <- spp.geo %>%
    
    ## Pipe the list into lapply
    lapply(function(x) {
      
      ## Subset the the data frame to calculate area of occupancy according the IUCN.eval
      DF   = subset(AUS.AOO.DAT, searchTaxon == x)[, c("lat", "lon", "searchTaxon")]
      
      message('Calcualting geographic ranges for ', x, ', ', nrow(DF), ' records')
      #IUCN.eval(DF, Cell_size_AOO = 10, DrawMap = FALSE, country_map = land, SubPop = FALSE) 
      AOO  = AOO.computing(XY = DF)
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
  COMBO.NICHE = join(GBIF.AOO, COMBO.NICHE, type = "right")
  
  
  
  
  
  #########################################################################################################################
  ## 6). PLOT HISTOGRAMS FOR EACH SPECIES AT 1KM
  #########################################################################################################################

    
  ##############################################################################################
  ## Plot histograms of temperature and rainfall
  ## species = plot.taxa[9]
  for (species in plot.taxa) {
    
    ## Subset the spatial dataframe into records for each species
    SP.DF  <- subset(NICHE.1KM.84, searchTaxon == species)
    
    ## Subset DF into records for each species
    DF     <- subset(COMBO.SUA.POA, searchTaxon == species)
    DF.OCC <- subset(COMBO.SUA.POA, searchTaxon == species & SOURCE != "INVENTORY")
    DF.INV <- subset(COMBO.SUA.POA, searchTaxon == species & SOURCE == "INVENTORY")
    
    #############################################################
    ## Plot occurrence points by source for the world
    message('Writing global occ sources for ', species)
    png(sprintf("./data/ANALYSIS/CLEAN_GBIF/%s_%s", species, "occ_points_SOURCE.png"),
        16, 10, units = 'in', res = 500)
    
    plot(AUS.84, main = paste0("Global points for ", species),
         lwd = 0.01, asp = 1, col = 'grey', bg = 'sky blue')
    
    points(SP.DF,
           pch = ".", cex = 3.3, cex.lab = 3, cex.main = 4, cex.axis = 2,
           xlab = "", ylab = "", asp = 1,
           col = factor(SP.DF$SOURCE))
    
    dev.off()
    
    #############################################################
    ## Plot temperature histograms
    message('Writing global temp histograms for ', species)
    png(sprintf("./data/ANALYSIS/CLEAN_GBIF/%s_%s", species, "temp_niche_histograms_all_records.png"),
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
    png(sprintf("./data/ANALYSIS/CLEAN_GBIF/%s_%s", species, "rain_niche_histograms_all_records.png"),
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