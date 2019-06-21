#########################################################################################################################
############################################# CHECK COORDCLEAN SETTINGS #################################################
#########################################################################################################################


#########################################################################################################################
##  PLOT GBIF ERRORS TO CHECK
#########################################################################################################################


#########################################################################################################################
## Try plotting the points which are outliers for a subset of species and label them
CLEAN.PLOT = SpatialPointsDataFrame(coords      = TEST.GEO[c("lon", "lat")],
                                    data        = TEST.GEO,
                                    proj4string = CRS.WGS.84)

ALL.PLOT = SpatialPointsDataFrame(coords      = CLEAN.TRUE[c("lon", "lat")],
                                  data        = CLEAN.TRUE,
                                  proj4string = CRS.WGS.84)

## Create global and australian shapefile in the local coordinate system
LAND.84 = LAND %>%
  spTransform(CRS.WGS.84)
AUS.84 = AUS %>%
  spTransform(CRS.WGS.84)


#########################################################################################################################
## Get the first 10 species
## species = plot.taxa[1]
## species = plot.taxa[1]
plot.taxa <- as.character(unique(CLEAN.PLOT$searchTaxon))
for (species in plot.taxa) {
  
  ## Plot a subset of taxa
  CLEAN.PLOT.PI = CLEAN.PLOT[ which(CLEAN.PLOT$searchTaxon == species), ]
  
  message("plotting occ data for ", species, ", ", 
          nrow(CLEAN.PLOT.PI), " records flagged as either ",
          unique(CLEAN.PLOT.PI$coord_summary))
  
  #############################################################
  ## Plot true and false points for the world
  ## Black == FALSE
  ## Red   == TRUE
  message('Writing map of global coord clean records for ', species)
  png(sprintf("./data/ANALYSIS/CLEAN_GBIF/%s_%s", species, "global_check_cc.png"),
      16, 10, units = 'in', res = 500)
  
  par(mfrow = c(1,2))
  plot(LAND.84, main = paste0(nrow(subset(CLEAN.PLOT.PI, coord_summary == "FALSE")),
                              " Global clean_coord 'FALSE' points for ", species),
       lwd = 0.01, asp = 1, col = 'grey', bg = 'sky blue')
  
  points(CLEAN.PLOT.PI,
         pch = ".", cex = 3.3, cex.lab = 3, cex.main = 4, cex.axis = 2, 
         xlab = "", ylab = "", asp = 1,
         col = factor(CLEAN.PLOT.PI$coord_summary))
  
  # legend("topleft", legend = levels(as.factor(CLEAN.PLOT.PI$coord_summary)), 
  #        pch = 16,  col = unique(CLEAN.PLOT.PI$coord_summary))
  
  #############################################################
  ## Plot true and false points for the world
  ## Black == FALSE
  ## Red   == TRUE
  plot(AUS.84, main = paste0(nrow(subset(CLEAN.PLOT.PI, coord_summary == "FALSE")),
                             " Global clean_coord 'FALSE' points for ", species),
       lwd = 0.01, asp = 1, bg = 'sky blue', col = 'grey')
  
  points(CLEAN.PLOT.PI,
         pch = ".", cex = 3.3, cex.lab = 3, cex.main = 4, cex.axis = 2, 
         xlab = "", ylab = "", asp = 1,
         col = factor(CLEAN.PLOT.PI$coord_summary))
  
  # legend("topleft", legend = levels(as.factor(CLEAN.PLOT.PI$coord_summary)), 
  #        pch = 16,  col = unique(CLEAN.PLOT.PI$coord_summary))
  
  dev.off()
  
}


#########################################################################################################################
## This code could be used to manually clean outlying records, using the OBS column......................................
## Note that for records with catalogue numbers, you can actually search GBIF and ALA for those occurrences. This seems
## like overkill for so many species. But you could make a list of the OBS column, then exclude these from the larger
## dataframe. Seems like you would only do this on the final dataset, filtered to 1km resolution.


## Select columns to send to the shapefile
GBIF.ALA.CHECK  = dplyr::select(TEST.GEO,     CC.OBS, searchTaxon, scientificName, lat, lon, SOURCE, year, 
                                country, locality, basisOfRecord, institutionCode,
                                coord_spp, coord_val,  coord_equ,  coord_zer,  coord_cap, 
                                coord_cen, coord_gbf, coord_inst, coord_summary)


## Rename the fields so that ArcMap can handle them
GBIF.ALA.CHECK  = dplyr::rename(GBIF.ALA.CHECK, 
                                OBS       = CC.OBS,
                                TAXON     = searchTaxon,
                                LAT       = lat,
                                LON       = lon,
                                SC_NAME   = scientificName,
                                BASIS     = basisOfRecord,                
                                LOCAL     = locality,                      
                                INSTIT    = institutionCode,                
                                COUNTRY   = country,                
                                YEAR      = year,
                                COORD_CL  = coord_summary)
names(GBIF.ALA.CHECK)


## Then create a SPDF
GBIF.ALA.SPDF    = SpatialPointsDataFrame(coords      = GBIF.ALA.CHECK[c("LON", "LAT")],
                                          data        = GBIF.ALA.CHECK,
                                          proj4string = CRS.WGS.84)






#########################################################################################################################
## CREATE HISTOGRAMS FOR EACH SPECIES, USING INVENTORY DATA
#########################################################################################################################


## Create a DF for the histograms
HIST.PLOT = SpatialPointsDataFrame(coords      = CLEAN.INV[c("lon", "lat")],
                                   data        = CLEAN.INV,
                                   proj4string = CRS.WGS.84)


##############################################################################################
## Plot histograms of temperature and rainfall
## species = plot.taxa[5]
for (species in plot.taxa) {
  
  ## Subset the spatial dataframe into records for each species
  #SP.DF  <- subset(HIST.PLOT, searchTaxon == i)
  SP.DF  <- HIST.PLOT[ which(HIST.PLOT$searchTaxon == species), ]
  
  ## Subset DF into records for each species
  DF     <- CLEAN.INV[ which(CLEAN.INV$searchTaxon == species), ] #subset(CLEAN.TRUE, searchTaxon == i)
  DF.OCC <- CLEAN.INV[ which(CLEAN.INV$searchTaxon == species & CLEAN.INV$SOURCE != "INVENTORY") , ]
  DF.INV <- CLEAN.INV[ which(CLEAN.INV$searchTaxon == species & CLEAN.INV$SOURCE != "INVENTORY") , ]
  
  #############################################################
  ## Plot occurrence points by source for the world
  message('Writing global occ sources for ', species)
  png(sprintf("./data/ANALYSIS/CLEAN_GBIF/%s_%s", species, "occ_points_SOURCE.png"),
      16, 10, units = 'in', res = 500)
  
  plot(LAND.84, main = paste0("Global points for ", species),
       lwd = 0.01, asp = 1, col = 'grey', bg = 'sky blue')
  
  points(SP.DF,
         pch = ".", cex = 3.3, cex.lab = 3, cex.main = 4, cex.axis = 2,
         xlab = "", ylab = "", asp = 1,
         col = factor(SP.DF$SOURCE))
  
  dev.off()
  
  ## Write out shapefile, to check if the inventory data is being used
  # message("Writing coord_clean shapefile for ", species)
  # tmp <- HIST.PLOT[HIST.PLOT$searchTaxon == species, ]
  # writeOGR(tmp, dsn = "./data/ANALYSIS/CLEAN_GBIF", paste0(species, '_occ_points_source'), 
  #          driver = "ESRI Shapefile", overwrite_layer = TRUE)
  
  #############################################################
  ## Plot temperature histograms
  # message('Writing global temp histograms for ', species)
  # png(sprintf("./data/ANALYSIS/CLEAN_GBIF/%s_%s", species, "temp_niche_histograms_all_records.png"),
  #     16, 10, units = 'in', res = 500)
  # 
  # ## Use the 'SOURCE' column to create a histogram for each source.
  # temp.hist = ggplot(DF, aes(x = Annual_mean_temp, group = SOURCE, fill = SOURCE)) +
  #   
  #   geom_histogram(position = "identity", alpha = 0.5, binwidth = 0.25,
  #                  aes(y =..density..))  + 
  #   geom_density(col = 4, alpha = 0.5) +
  #   
  #   ## Add some median lines : overall, ALA and GBIF
  #   geom_vline(aes(xintercept = median(DF.INV$Annual_mean_temp)),
  #              col = 'blue', size = 1) +
  #   geom_vline(aes(xintercept = median(DF.OCC$Annual_mean_temp)),
  #              col = 'red', size = 1) +
  #   
  #   ggtitle(paste0("Worldclim temp niches for ", species)) +
  #   
  #   ## Add themes
  #   theme(axis.title.x     = element_text(colour = "black", size = 35),
  #         axis.text.x      = element_text(size = 25),
  #         
  #         axis.title.y     = element_text(colour = "black", size = 35),
  #         axis.text.y      = element_text(size = 25),
  #         
  #         panel.background = element_blank(),
  #         panel.border     = element_rect(colour = "black", fill = NA, size = 3),
  #         plot.title       = element_text(size   = 40, face = "bold"),
  #         legend.text      = element_text(size   = 20),
  #         legend.title     = element_text(size   = 20),
  #         legend.key.size  = unit(1.5, "cm"))
  # 
  # ## Print the plot and close the device 
  # print(temp.hist + ggtitle(paste0("Worldclim temp niches for ", species)))
  # dev.off()
  # 
  # #############################################################
  # ## Plot rainfall histograms
  # message('Writing global rain histograms for ', species)
  # png(sprintf("./data/ANALYSIS/CLEAN_GBIF/%s_%s", species, "rain_niche_histograms_all_records.png"),
  #     16, 10, units = 'in', res = 500)
  # 
  # ## Use the 'SOURCE' column to create a histogram for each source.
  # rain.hist = ggplot(DF, aes(x = Annual_precip, group = SOURCE, fill = SOURCE)) +
  #   
  #   geom_histogram(position = "identity", alpha = 0.5, binwidth = 15,
  #                  aes(y =..density..))  + 
  #   geom_density(col = 4, alpha = 0.5) +
  #   
  #   ## Add some median lines : overall, ALA and GBIF
  #   geom_vline(aes(xintercept = median(DF.INV$Annual_precip)),
  #              col = 'blue', size = 1) +
  #   geom_vline(aes(xintercept = median(DF.OCC$Annual_precip)),
  #              col = 'red', size = 1) +
  #   
  #   ggtitle(paste0("Worldclim rain niches for ", species)) +
  #   
  #   ## Add themes
  #   theme(axis.title.x     = element_text(colour = "black", size = 35),
  #         axis.text.x      = element_text(size = 25),
  #         
  #         axis.title.y     = element_text(colour = "black", size = 35),
  #         axis.text.y      = element_text(size = 25),
  #         
  #         panel.background = element_blank(),
  #         panel.border     = element_rect(colour = "black", fill = NA, size = 3),
  #         plot.title       = element_text(size   = 40, face = "bold"),
  #         legend.text      = element_text(size   = 20),
  #         legend.title     = element_text(size   = 20),
  #         legend.key.size  = unit(1.5, "cm"))
  # 
  # ## Print the plot and close the device 
  # print(rain.hist)
  # dev.off()
  
}





#########################################################################################################################
#####################################################  TBC ##############################################################
#########################################################################################################################