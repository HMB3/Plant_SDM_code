#########################################################################################################################
################################################# CALCULATE 1KM NICHES ##################################################
#########################################################################################################################


## This code estimates realised environmental and geographic ranges (i.e. niches) for a list of species


#########################################################################################################################
## Only run this code if we are saving niches 
if(calc_niche == "TRUE") {
  
  
  ## If we need to read the data in? 
  if(read_data == "TRUE") {
    
    ## read in RDS files from previous step
    CLEAN.INV = readRDS(paste0(DATA_path, 'CLEAN_INV_', save_run, '.rds'))
    message('Species overlap ', length(intersect(GBIF.spp, unique(CLEAN.INV$searchTaxon))))
    rasterTmpFile()
    
  } else {
    
    message(' skip file reading, not many species analysed')   ##
    
  }
  
  
  
  #########################################################################################################################
  ## Create a spatial points object
  ## We want to summarise the niches at 1km, but including all the environmental variables (EG PET, etc,),
  ## Not just those used in the SDM table (i.e. worldclim so far)
  ## CLEAN.INV = CLEAN.INV[CLEAN.INV$searchTaxon %in% GBIF.spp, ]
  NICHE.1KM    <- CLEAN.INV
  NICHE.1KM.84 <- SpatialPointsDataFrame(coords      = NICHE.1KM[c("lon", "lat")],
                                         data        = NICHE.1KM,
                                         proj4string = CRS.WGS.84)
  
  
  ## Use a projected, rather than geographic, coordinate system
  ## Not sure why, but this is needed for the spatial overlay step
  #NICHE.1KM.84 = spTransform(NICHE.1KM, CRS.WGS.84)
  AUS.WGS      = spTransform(AUS,       CRS.WGS.84)
  POA.WGS      = spTransform(POA_2016,  CRS.WGS.84)
  SUA.WGS      = spTransform(SUA_2016 , CRS.WGS.84)
  LAND.WGS     = spTransform(LAND,      CRS.WGS.84)
  
  
  ## Create global niche and Australian niche.
  ## So we need a subset for Australia
  ## The ,] acts like a clip in ArcMap 
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
  ## Can this be done using 
  COMBO.SUA.POA = cbind.data.frame(NICHE.1KM.84, SUA.JOIN, POA.JOIN) 
  
  
  #########################################################################################################################
  ## Aggregate the number of SUAs and POAs each species is found in
  SUA.AGG   = tapply(COMBO.SUA.POA$SUA_NAME16, COMBO.SUA.POA$searchTaxon, function(x) length(unique(x))) ## group SUA by species name
  POA.AGG   = tapply(COMBO.SUA.POA$POA_NAME16, COMBO.SUA.POA$searchTaxon, function(x) length(unique(x))) ## group POA by species name
  SUA.AGG   = as.data.frame(SUA.AGG)
  POA.AGG   = as.data.frame(POA.AGG)
  
  
  ## Create tables
  SUA.AGG <- setDT(SUA.AGG, keep.rownames = TRUE)[]
  POA.AGG <- setDT(POA.AGG, keep.rownames = TRUE)[]
  names(SUA.AGG) = c("searchTaxon", "SUA_count")
  names(POA.AGG) = c("searchTaxon", "POA_count")
  
  
  ## Now create a table of all the SUA's that each species occurrs
  SUA.AGG = join(SUA.AGG, POA.AGG)
  saveRDS(SUA.AGG, paste0(DATA_path, 'SUA_SPP_COUNT_', save_run, '.rds'))
  
  
  
  
  
  #########################################################################################################################
  ## 4). CREATE NICHES FOR PROCESSED TAXA
  #########################################################################################################################
  
  
  #########################################################################################################################
  ## Now summarise the niches. But, figure out a cleaner way of doing this
  message('Estimating global niches for ', length(unique(CLEAN.INV$searchTaxon)), ' species across ', 
          length(env.variables), ' climate variables')
  
  
  #########################################################################################################################
  ## Create niche summaries for each environmental condition like this...commit
  ## Here's what the function will produce :
  NICHE.AUS.DF = NICHE.AUS %>% 
    as.data.frame() %>% 
    dplyr::select(., searchTaxon, one_of(env.variables))
  
  
  NICHE.GLO.DF = NICHE.1KM.84 %>% 
    as.data.frame() %>% 
    dplyr::select(., searchTaxon, one_of(env.variables))
  
  
  ## new_DF <- NICHE.AUS.DF[rowSums(is.na(NICHE.AUS.DF)) > 0,]
  NICHE.AUS.DF = completeFun(NICHE.AUS.DF, "PET")
  NICHE.GLO.DF = completeFun(NICHE.GLO.DF, "PET")
  
  
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
  ## Add global counts for each species, and record the total number of taxa processed
  ## This is actually redundant - the records will the be the same as the maxent_records
  Global_records = as.data.frame(table(COMBO.SUA.POA$searchTaxon))
  names(Global_records) = c("searchTaxon", "Global_records")
  identical(nrow(Global_records), nrow(GLOB.NICHE))
  
  
  ## Add the count of Australian records - this is not necessarily the same as maxent_records 
  Aus_records = as.data.frame(table(NICHE.AUS.DF$searchTaxon))
  names(Aus_records) = c("searchTaxon", "Aus_records")
  identical(nrow(Aus_records), nrow(AUS.NICHE))
  
  
  #########################################################################################################################
  ## Add the counts of Australian records for each species to the niche database
  GLOB.NICHE = join(Aus_records,    GLOB.NICHE,  type = "right")
  GLOB.NICHE = join(Global_records, GLOB.NICHE,  type = "right")
  GLOB.NICHE = join(SUA.AGG,        GLOB.NICHE,  type = "right")
  
  
  ## Check the record and POA counts
  head(GLOB.NICHE$Aus_records,    20)
  head(GLOB.NICHE$Global_records, 20)
  head(GLOB.NICHE$SUA_count, 20)
  head(GLOB.NICHE$POA_count, 20)
  
  names(GLOB.NICHE)
  dim(GLOB.NICHE)
  dim(AUS.NICHE)
  
  
  
  
  
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
  ## 6). PLOT CONVEX HULL OF DATA
  #########################################################################################################################
  
  
  #########################################################################################################################
  ## Create a convex hull plot of points for all specues, according to source :: OCC and GBIF
  # CLEAN.INV <- mutate(CLEAN.INV, OCC_TYPE = ifelse(grepl("INVENTORY", SOURCE), "INV", "OCC"))
  # unique(CLEAN.INV$SOURCE)
  # unique(CLEAN.INV$OCC_TYPE)
  
  
  #########################################################################################################################
  ## Calculate the hulls for each group
  #test = subset(CLEAN.INV, searchTaxon == CLEAN.INV$searchTaxon[1])
  
  
  # ## Create PNG file
  # png(sprintf("./data/ANALYSIS/SPECIES_RANGES/%s", "all_species_1km_convex_hull.png"),
  #     16, 10, units = 'in', res = 500)
  # 
  # p <- ggplot(CLEAN.INV, aes(Annual_mean_temp, Annual_precip, fill = OCC_TYPE, color = OCC_TYPE)) 
  # 
  # hull_occ_source <- CLEAN.INV %>%
  #   group_by(OCC_TYPE) %>%
  #   slice(chull(Annual_mean_temp, Annual_precip))
  # 
  # ## Update the plot with a fill group, and overlay the new hulls
  # p + geom_polygon(data = hull_occ_source, alpha = 0.3) +
  #   geom_point(shape = 21, size = 2.5) + ## geom_density_2d
  #   
  #   ## Add x,y, and title
  #   labs(x = "Mean annual temp", 
  #        y = "Annual precipitation",
  #        title = "Convex Hull for") +
  #   
  #   ## Add themes
  #   theme(axis.title.x     = element_text(colour = "black", size = 35),
  #         axis.text.x      = element_text(size = 20),
  #         
  #         axis.title.y     = element_text(colour = "black", size = 35),
  #         axis.text.y      = element_text(size = 20),
  # 
  #         panel.background = element_blank(),
  #         panel.border     = element_rect(colour = "black", fill = NA, size = 1.5),
  #         plot.title       = element_text(size   = 40, face = "bold"),
  #         legend.text      = element_text(size   = 20),
  #         legend.title     = element_text(size   = 20),
  #         legend.key.size  = unit(1.5, "cm"))
  # 
  # ## close device
  # dev.off()
  
  
  
  
  
  #########################################################################################################################
  ## 8). JOIN NICHE DATE TO HORTICULTURAL CONTEXTUAL DATA
  #########################################################################################################################
  
  
  ## Need to figure out what make sense to include in table - Global niche, Australian niche? Etc.
  
  
  #########################################################################################################################
  ## Now join the horticultural contextual data onto one or both tables ()
  message('Joining contextual data for raster and niche files ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")
  COMBO.RASTER.CONTEXT = CLEAN.INV
  names(COMBO.RASTER.CONTEXT)
  
  
  #########################################################################################################################
  ## Now join the hort context data and the tree inventory plantings to the niche
  COMBO.NICHE.CONTEXT = join(CLEAN.GROW, GLOB.NICHE,           type = "right") 
  COMBO.NICHE.CONTEXT = join(TI.LUT,     COMBO.NICHE.CONTEXT,  type = "right")
  COMBO.NICHE.CONTEXT = COMBO.NICHE.CONTEXT[c(1, 3:14, 2, 15:ncol(COMBO.NICHE.CONTEXT))]
  
  
  ## Check the columns are still there
  head(COMBO.NICHE.CONTEXT$Aus_records,    20)
  head(COMBO.NICHE.CONTEXT$Global_records, 20)
  head(COMBO.NICHE.CONTEXT$Plantings,      20)
  head(COMBO.NICHE.CONTEXT$SUA_count,      20)
  
  
  #########################################################################################################################
  ## Is there a relationship between the number of records in each category of records?
  ## Aus records now gets removed, because we are using maxent records
  #windows();
  # plot(COMBO.NICHE.CONTEXT$Aus_records, COMBO.NICHE.CONTEXT$Global_records)
  # plot(COMBO.NICHE.CONTEXT$Aus_records, COMBO.NICHE.CONTEXT$Plantings)
  
  
  lm.glob  = lm(COMBO.NICHE.CONTEXT$Aus_records ~ COMBO.NICHE.CONTEXT$Global_records)
  lm.plant = lm(COMBO.NICHE.CONTEXT$Aus_records ~ COMBO.NICHE.CONTEXT$Plantings)
  
  layout(matrix(c(1,1,2,3), 2, 1, byrow = TRUE))
  
  plot(COMBO.NICHE.CONTEXT$Global_records, COMBO.NICHE.CONTEXT$Aus_records,
       pch = 19, col  = "blue",
       xlab = "Global records", ylab = "Australian records",
       abline(lm.glob),
       main = save_run, cex = 2)
  
  legend("topleft", bty = "n",
         legend = paste("R2 is", format(summary(lm.glob)$adj.r.squared, digits = 4)))
  
  plot(COMBO.NICHE.CONTEXT$Plantings, COMBO.NICHE.CONTEXT$Aus_records,
       pch = 19, col  = "blue",
       xlab = "Plantings", 
       ylab = "Australian records", 
       abline(lm.plant),
       main = save_run, cex = 2)
  
  legend("topright", bty = "n",
         legend = paste("R2 is", format(summary(lm.plant)$adj.r.squared, digits = 4)))
  
  
  #########################################################################################################################
  ## Now combine the SDM output with the niche context data 
  COMBO.NICHE.CONTEXT = COMBO.NICHE.CONTEXT[order(COMBO.NICHE.CONTEXT$searchTaxon),]
  
  
  ## Set Total.growers NA to blank, then sort by no. of growers
  COMBO.NICHE.CONTEXT$Total_growers[is.na(COMBO.NICHE.CONTEXT$Total_growers)] <- 0
  COMBO.NICHE.CONTEXT = COMBO.NICHE.CONTEXT[with(COMBO.NICHE.CONTEXT, rev(order(Total_growers))), ]
  
  
  ## Print the dataframe dimensions to screen
  dim(COMBO.RASTER.CONTEXT)
  dim(COMBO.NICHE.CONTEXT)
  
  length(unique(CLEAN.INV$searchTaxon))
  length(unique(COMBO.RASTER.CONTEXT$searchTaxon))
  length(COMBO.NICHE.CONTEXT$searchTaxon)
  unique(COMBO.RASTER.CONTEXT$SOURCE)
  
  
  ## Check the column order for the niche and the record tables
  ## Note that for records with catalogue numbers, you can actually search GBIF and ALA for those occurrences
  head(COMBO.NICHE.CONTEXT)[1:15]
  head(COMBO.RASTER.CONTEXT)[1:16]
  
  
  #########################################################################################################################
  ## save .rds file for the next session
  message('Writing 1km resolution niche and raster data for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")
  saveRDS(COMBO.NICHE.CONTEXT,    paste0(DATA_path, 'COMBO_NICHE_CONTEXT_',  save_run, '.rds'))
  saveRDS(COMBO.RASTER.CONTEXT,   paste0(DATA_path, 'COMBO_RASTER_CONTEXT_', save_run, '.rds'))
  write.csv(COMBO.NICHE.CONTEXT,  paste0(DATA_path, 'COMBO_NICHE_CONTEXT_',  save_run, '.csv'), row.names = FALSE)
  
  
} else {
  
  message(' skip file saving, ', length(GBIF.spp), ' species analysed')   ##
  
}





#########################################################################################################################
#####################################################  TBC ##############################################################
#########################################################################################################################