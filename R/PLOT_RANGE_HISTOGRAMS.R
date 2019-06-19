#########################################################################################################################
################################################# PLOT SPECIES NICHES ###################################################
#########################################################################################################################


## This code creates barplots, histograms, and convex-hull plots for a set of species
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
## Calculate the mean temp and rain in each POA - don't do inside a loop
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
## species = spp.geo[1]
## geom_rect needs xmin ymax, ymin ymax 
for (species in spp.geo) {
  
  ## Subset the spatial dataframe into records for each species
  SP.DF     <- NICHE.1KM.84[NICHE.1KM.84$searchTaxon %in% species , ]
  DF        <- CLEAN.INV[CLEAN.INV$searchTaxon %in% species , ]
  
  TMP.GLO   <- subset(GLOB.NICHE,   searchTaxon == species)
  TMP.AUS   <- subset(AUS.NICHE,    searchTaxon == species)

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
    message("No Australian data for ", species, " don't plot the range")
  }

  TMP.RANGE <- rbind(TMP.GLO, TMP.AUS)
  names(TMP.RANGE)[2] = c("Temperature_range")

  ## Subset DF into records for each species
  DF     <- subset(COMBO.SUA.POA, searchTaxon == species)
  DF.OCC <- subset(COMBO.SUA.POA, searchTaxon == species & SOURCE != "INVENTORY")
  DF.INV <- subset(COMBO.SUA.POA, searchTaxon == species & SOURCE == "INVENTORY")
   
  # #############################################################
  # ## Plot occurrence points by source for the world
  # message('Writing global occ sources for ', species)
  # png(sprintf("./data/ANALYSIS/SPECIES_RANGES/%s_%s", species, "1km_occ_points_source.png"),
  #     16, 10, units = 'in', res = 500)
  # 
  plot(LAND.WGS, main = paste0("Global points for ", species),
       lwd = 0.01, asp = 1, col = 'grey', bg = 'sky blue')

  points(SP.DF,
         pch = ".", cex = 3.3, cex.lab = 3, cex.main = 4, cex.axis = 2,
         xlab = "", ylab = "", asp = 1,
         col = factor(SP.DF$SOURCE))
  # 
  # dev.off()
  # 
  # #############################################################
  # ## Plot temperature barchart
  # message('Writing global temp histograms for ', species)
  # png(sprintf("./data/ANALYSIS/SPECIES_RANGES/%s_%s", species, "temp_barchart_1km_records.png"),
  #     16, 10, units = 'in', res = 500)
  # 
  # ## Use the 'SOURCE' column to create a histogram for each source.
  # ## Back to here..............................................................................
  max.temp = max(TMP.RANGE$Annual_mean_temp_q95)+5
  min.temp = min(TMP.RANGE$Annual_mean_temp_q05)-5

  temp.bar =
    ggplot(TMP.RANGE, aes(y = Temperature_range, x = RANGE, fill = RANGE)) +   ## supply xmin, etc in aes

    # scale_y_discrete(limits = c(min.temp,
    #                             max.temp)) +

    geom_bar(stat = "identity", position = "identity", width = 0.1) +          ## use geom_rect here   
    coord_flip() +                                                             ## geom_segment

    ## Add some median lines : overall, ALA and GBIF
    ## This will only work if we plot the full range of temperatures on the x-axis
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

          panel.background = element_blank(),
          panel.border     = element_rect(colour = "black", fill = NA, size = 3),
          plot.title       = element_text(size   = 40, face = "bold"),
          legend.text      = element_text(size   = 20),
          legend.title     = element_text(size   = 20),
          legend.key.size  = unit(1.5, "cm"))
  # 
  # ## Print the plot and close the device
  # print(temp.hist + ggtitle(paste0("Worldclim temp niches for ", species)))
  # dev.off()
  
  
  #############################################################
  ## PLot convex Hull
  message('Writing global convex hulls for ', species)
  DF.CONV <- mutate(DF, OCC_TYPE = ifelse(grepl("INVENTORY", SOURCE), "INV", "OCC"))
  unique(DF.CONV)
  unique(DF.CONV$OCC_TYPE)
  
  
  png(sprintf("./data/ANALYSIS/SPECIES_RANGES/%s_%s", species, "_1km_convex_hull.png"),
      16, 10, units = 'in', res = 500)
  
  convex <- ggplot(DF.CONV, aes(Annual_mean_temp, Annual_precip, fill = OCC_TYPE, color = OCC_TYPE)) 
  
  hull_occ_source <- DF.CONV %>%
    group_by(OCC_TYPE) %>%
    slice(chull(Annual_mean_temp, Annual_precip))
  
  ## Update the plot with a fill group, and overlay the new hulls
  convex + geom_polygon(data = hull_occ_source, alpha = 0.3) +
    geom_point(shape = 21, size = 2.5) + ## geom_density_2d
    
    ## Add x,y, and title
    labs(x = "Mean annual temp", 
         y = "Annual precipitation",
         title = "Convex Hull for ", species) +
    
    ## Add themes
    theme(axis.title.x     = element_text(colour = "black", size = 35),
          axis.text.x      = element_text(size = 20),
          
          axis.title.y     = element_text(colour = "black", size = 35),
          axis.text.y      = element_text(size = 20),
          
          panel.background = element_blank(),
          panel.border     = element_rect(colour = "black", fill = NA, size = 1.5),
          plot.title       = element_text(size   = 40, face = "bold"),
          legend.text      = element_text(size   = 20),
          legend.title     = element_text(size   = 20),
          legend.key.size  = unit(1.5, "cm"))
  
  
  ##
  print(convex + ggtitle(paste0("Worldclim temp niches for ", species)))
  
  ## close device
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
    geom_vline(aes(xintercept = median(DF$Annual_mean_temp)),
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
    geom_vline(aes(xintercept = median(DF$Annual_precip)),
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
  print(temp.hist + ggtitle(paste0("Worldclim rain niches for ", species)))
  dev.off()
  
}





#########################################################################################################################
#####################################################  TBC ##############################################################
#########################################################################################################################