#########################################################################################################################
################################################# PLOT SPECIES NICHES ###################################################
#########################################################################################################################


## This code creates barplots, histograms, and convex-hull plots for a set of species
## This code can be used to re-do the zonal stats for the SUAs in step 9................................................. 
## If we need to read the data in? 
if(read_data == "TRUE") {
  
  ## read in RDS files from previous step
  ## CLEAN.INV = CLEAN.INV[CLEAN.INV$searchTaxon %in% na.error$searchTaxon[1:50], ]
  CLEAN.INV = readRDS(paste0(DATA_path, 'CLEAN_INV_', save_run, '.rds'))
  message('Species overlap ', length(intersect(GBIF.spp, unique(CLEAN.INV$searchTaxon))))
  rasterTmpFile()
  
} else {
  
  message(' skip file reading, not many species analysed')   ##
  
}


## Create a list of species to loop over
spp.geo = as.character(unique(CLEAN.INV$searchTaxon))


#########################################################################################################################
## Just use the two environmental conditions likely  to be used for ranges
# POA_temp = aus.grids.current[[1]]
# POA_rain = aus.grids.current[[12]]
# 
# ## Make sure the raster extents match with the POA
# POA_SF <- POA.WGS %>%
#   spTransform(., ALB.CON) %>%
#   crop(., extent(POA_temp)) %>%
#   st_as_sf()


#########################################################################################################################
## Calculate the mean temp and rain in each POA - don't do inside a loop
# POA_SF$Annual_mean_temp    <- exact_extract(POA_temp, POA_SF, weighted.mean, na.rm = TRUE)
# POA_SF$Annual_precip       <- exact_extract(POA_rain, POA_SF, weighted.mean, na.rm = TRUE)
# 
# POA_SF$Annual_mean_temp_30 <- exact_extract(aus.temp.2030, POA_SF, weighted.mean, na.rm = TRUE)
# POA_SF$Annual_precip_30    <- exact_extract(aus.rain.2030, POA_SF, weighted.mean, na.rm = TRUE)
# 
# POA_SF$Annual_mean_temp_50 <- exact_extract(aus.temp.2050, POA_SF, weighted.mean, na.rm = TRUE)
# POA_SF$Annual_precip_50    <- exact_extract(aus.rain.2050, POA_SF, weighted.mean, na.rm = TRUE)
# 
# POA_SF$Annual_mean_temp_70 <- exact_extract(aus.temp.2070, POA_SF, weighted.mean, na.rm = TRUE)
# POA_SF$Annual_precip_70    <- exact_extract(aus.rain.2070, POA_SF, weighted.mean, na.rm = TRUE)
# 
# length(POA_SF$Annual_mean_temp_30);length(POA_SF$Annual_precip_30);
# length(POA_SF$Annual_mean_temp_50);length(POA_SF$Annual_precip_50);
# length(POA_SF$Annual_mean_temp_70);length(POA_SF$Annual_precip_70);


## Create a dataframe of the temperature and rainfal
# POA_climate          = as.data.frame(POA_SF)
# POA_climate$geometry = NULL
# head(POA_climate$POA_CODE16)


#########################################################################################################################
## How to include the POA? Could add them all
# POA.SYD = POA_climate[POA_climate$POA_CODE16 %in% 2000 , ]
# POA.BRI = POA_climate[POA_climate$POA_CODE16 %in% 4000 , ]
# POA.MEL = POA_climate[POA_climate$POA_CODE16 %in% 3000 , ]
# POA.PER = POA_climate[POA_climate$POA_CODE16 %in% 6000 , ]
# POA.ADE = POA_climate[POA_climate$POA_CODE16 %in% 5000 , ]
# POA.DAR = POA_climate[POA_climate$POA_CODE16 %in% 0800 , ]
# POA.HOB = POA_climate[POA_climate$POA_CODE16 %in% 7000 , ]
# POA.CAN = POA_climate[POA_climate$POA_CODE16 %in% 2601 , ]





#########################################################################################################################
## 7). PLOT HISTOGRAMS AND BAR CHARTS FOR EACH SPECIES AT 1KM
#########################################################################################################################


##############################################################################################
## Plot histograms of temperature and rainfall
##  spp = spp.geo[1]
## geom_rect needs xmin ymax, ymin ymax 
for (spp in spp.geo) {
  
  ## Subset the spatial dataframe into records for each spp
  SP.DF     <- NICHE.1KM.84[NICHE.1KM.84$searchTaxon %in% spp , ]
  DF        <- CLEAN.INV[CLEAN.INV$searchTaxon %in% spp , ]
  
  # TMP.GLO   <- subset(GLOB.NICHE,   searchTaxon == spp)
  # TMP.AUS   <- subset(AUS.NICHE,    searchTaxon == spp)

  #############################################################
  ## Now, build a df of the temperature vectors
  # if(nrow(TMP.GLO) > 0){
  #   TMP.GLO$RANGE = "GLOBAL"
  # } else {
  #   message("No global data for ", spp)
  # }
  # 
  # if(nrow(TMP.AUS) > 0){
  #   TMP.AUS$RANGE = "AUS"
  # } else {
  #   message("No Australian data for ", spp, " don't plot the range")
  # }
  # 
  # TMP.RANGE <- rbind(TMP.GLO, TMP.AUS)
  # names(TMP.RANGE)[2] = c("Temperature_range")

  ## Subset DF into records for each spp
  # DF     <- subset(COMBO.SUA.POA, searchTaxon == spp)
  # DF.OCC <- subset(COMBO.SUA.POA, searchTaxon == spp & SOURCE != "INVENTORY")
  # DF.INV <- subset(COMBO.SUA.POA, searchTaxon == spp & SOURCE == "INVENTORY")
   
  # #############################################################
  # ## Plot occurrence points by source for the world
  # message('Writing global occ sources for ', spp)
  # png(sprintf("./data/ANALYSIS/SPECIES_RANGES/%s_%s", spp, "1km_occ_points_source.png"),
  #     16, 10, units = 'in', res = 500)
  # 
  # plot(LAND.WGS, main = paste0("Global points for ", spp),
  #      lwd = 0.01, asp = 1, col = 'grey', bg = 'sky blue')
  # 
  # points(SP.DF,
  #        pch = ".", cex = 3.3, cex.lab = 3, cex.main = 4, cex.axis = 2,
  #        xlab = "", ylab = "", asp = 1,
  #        col = factor(SP.DF$SOURCE))
  # 
  # dev.off()
  # 
  # #############################################################
  # ## Plot temperature barchart
  # message('Writing global temp histograms for ', spp)
  # png(sprintf("./data/ANALYSIS/SPECIES_RANGES/%s_%s", spp, "temp_barchart_1km_records.png"),
  #     16, 10, units = 'in', res = 500)
  # 
  # ## Use the 'SOURCE' column to create a histogram for each source.
  # ## Back to here..............................................................................
  # max.temp = max(TMP.RANGE$Annual_mean_temp_q95)+5
  # min.temp = min(TMP.RANGE$Annual_mean_temp_q05)-5
  # 
  # temp.bar =
  #   ggplot(TMP.RANGE, aes(y = Temperature_range, x = RANGE, fill = RANGE)) +   ## supply xmin, etc in aes

    # scale_y_discrete(limits = c(min.temp,
    #                             max.temp)) +

    # geom_bar(stat = "identity", position = "identity", width = 0.1) +          ## use geom_rect here   
    # coord_flip() +                                                             ## geom_segment
    # 
    # ## Add some median lines : overall, ALA and GBIF
    # ## This will only work if we plot the full range of temperatures on the x-axis
    # geom_vline(aes(xintercept = POA.SYD$Annual_mean_temp),
    #            col = 'blue', size = 1) +
    # geom_vline(aes(xintercept = POA.SYD$Annual_mean_temp_50),
    #            col = 'red', size = 1) +
    # geom_vline(aes(xintercept = POA.SYD$Annual_mean_temp_70),
    #            col = 'green', size = 1) +
    # 
    # 
    # ggtitle(paste0("Worldclim temperature ranges for ", spp)) +
    # 
    # ## Add themes
    # theme(axis.title.x     = element_text(colour = "black", size = 35),
    #       axis.text.x      = element_text(size = 25),
    # 
    #       panel.background = element_blank(),
    #       panel.border     = element_rect(colour = "black", fill = NA, size = 3),
    #       plot.title       = element_text(size   = 40, face = "bold"),
    #       legend.text      = element_text(size   = 20),
    #       legend.title     = element_text(size   = 20),
    #       legend.key.size  = unit(1.5, "cm"))
  # 
  # ## Print the plot and close the device
  # print(temp.hist + ggtitle(paste0("Worldclim temp niches for ", spp)))
  # dev.off()
  
  
  #############################################################
  ## Plot convex Hull
  message('Writing global convex hulls for ', spp)
  #DF.CONV <- plyr::mutate(DF, OCC_TYPE = ifelse(grepl("INVENTORY", SOURCE), "INV", "OCC"))
  
  ## Start PNG device
  png(sprintf("./data/ANALYSIS/SPECIES_RANGES/%s_%s", spp, "_1km_convex_hull.png"),
      16, 10, units = 'in', res = 500)
  
  find_hull <- function(DF) DF[chull(DF$Annual_mean_temp, DF$Annual_precip), ]
  hulls     <- ddply(DF, "SOURCE", find_hull)
  
  plot      <- ggplot(data = DF, aes(x = Annual_mean_temp, 
                                     y = Annual_precip, colour = SOURCE, fill = SOURCE)) +
    geom_point() + 
    geom_polygon(data = hulls, alpha = 0.5) +
    labs(x = "Annual_mean_temp", y = "Annual_precip") +
    
    theme(axis.title.x     = element_text(colour = "black", size = 35),
          axis.text.x      = element_text(size = 20),
          
          axis.title.y     = element_text(colour = "black", size = 35),
          axis.text.y      = element_text(size = 20),
          
          panel.background = element_blank(),
          panel.border     = element_rect(colour = "black", fill = NA, size = 1.5),
          plot.title       = element_text(size   = 40, face = "bold"),
          legend.text      = element_text(size   = 20),
          legend.title     = element_text(size   = 20),
          legend.key.size  = unit(1.5, "cm")) + 
    
  ggtitle(paste0("Convex Hull for ", spp))
  print(plot)
  
  ## close device
  dev.off()
  
  #############################################################
  ## Plot temperature histograms
  message('Writing global temp histograms for ', spp)
  png(sprintf("./data/ANALYSIS/SPECIES_RANGES/%s_%s", spp, "temp_niche_histograms_1km_records.png"),
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
    
    ggtitle(paste0("Worldclim temp niches for ", spp)) +
    
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
  print(temp.hist + ggtitle(paste0("Worldclim temp niches for ", spp)))
  dev.off()
  
  #############################################################
  ## Plot rainfall histograms
  message('Writing global rain histograms for ', spp)
  png(sprintf("./data/ANALYSIS/SPECIES_RANGES/%s_%s", spp, "rain_niche_histograms_1km_records.png"),
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
    
    ggtitle(paste0("Worldclim rain niches for ", spp)) +
    
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
  print(rain.hist + ggtitle(paste0("Worldclim rain niches for ", spp)))
  dev.off()
  
}





#########################################################################################################################
#####################################################  TBC ##############################################################
#########################################################################################################################