#########################################################################################################################
############################################# CHECK COORDCLEAN SETTINGS #################################################
#########################################################################################################################


## Rename the columns to fit the CleanCoordinates format and create a tibble. 
# TIB.GBIF <- COMBO.RASTER.CONVERT %>% dplyr::rename(species          = searchTaxon,
#                                                    decimallongitude = lon, 
#                                                    decimallatitude  = lat) %>%
#   
#   ## The create a tibble for running the spatial outlier cleaning
#   timetk::tk_tbl() %>% 
#   
#   ## Consider the arguments. We've already stripped out the records that fall outside
#   ## the worldclim raster boundaries, so the sea test is probably not the most important
#   ## Study area is the globe, but we are only projecting models onto Australia
#   
#   ## The geographic outlier detection is not working here. Try it in step 6
#   ## Also, the duplicates step is not working, flagging too many records
#   clean_coordinates(.,
#                     verbose         = TRUE,
#                     tests = c("capitals",     "centroids", "equal", "gbif", 
#                               "institutions", "zeros"), ## duplicates flagged too many
#                     
#                     capitals_rad    = 10000,  ## remove records within 10km  of capitals
#                     centroids_rad   = 5000    ## remove records within 5km of country centroids
#                     # outliers_method = "distance", ## The other checks are not producing reliable results
#                     # outliers_td     = 800,
#                     # outliers_mtp    = 5
#   ) %>%
#   
#   ## The select the relevant columns and rename
#   select(., species, CC.OBS, .val,  .equ, .zer,  .cap,   
#          .cen,    .gbf,   .inst, .summary) 
# 
# ## Then rename
# summary(TIB.GBIF)
# names(TIB.GBIF) = c("coord_spp", "CC.OBS",    "coord_val",  "coord_equ",  "coord_zer",  "coord_cap", 
#                     "coord_cen", "coord_gbf", "coord_inst", "coord_summary")
# 
# 
# ## Flagging ~ x%, excluding the spatial outliers. Seems reasonable?
# message(round(with(TIB.GBIF, table(coord_summary)/sum(table(coord_summary))*100), 2), " % records removed")





#########################################################################################################################
##  PLOT GBIF ERRORS TO CHECK
#########################################################################################################################


#########################################################################################################################
## Try plotting the points which are outliers for a subset of species and label them
CLEAN.PLOT = SpatialPointsDataFrame(coords      = TEST.GEO[c("lon", "lat")],
                                    data        = TEST.GEO,
                                    proj4string = CRS.WGS.84)

## Create global and australian shapefile in the local coordinate system
LAND.84 = LAND %>%
  spTransform(CRS.WGS.84)
AUS.84 = AUS %>%
  spTransform(CRS.WGS.84)


#########################################################################################################################
## Get the first 10 species
## species = plot.taxa[9]
plot.taxa <- as.character(unique(CLEAN.PLOT$searchTaxon))
for (species in plot.taxa) {
  
  ## Plot a subset of taxa
  CLEAN.PLOT.PI = subset(CLEAN.PLOT, searchTaxon == species)
  
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
                                Taxonomic.status, New.Taxonomic.status,
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
                                TAX_STAT  = Taxonomic.status,
                                NEW_STA   = New.Taxonomic.status,
                                COORD_CL  = coord_summary)
names(GBIF.ALA.CHECK)


## Then create a SPDF
GBIF.ALA.SPDF    = SpatialPointsDataFrame(coords      = GBIF.ALA.CHECK[c("LON", "LAT")],
                                          data        = GBIF.ALA.CHECK,
                                          proj4string = CRS.WGS.84)


## Create list of species 
plot.taxa <- unique(GBIF.ALA.SPDF$TAXON)


#########################################################################################################################
## Then, loop over the species list and create a shapefile for each 
## For lots of taxa, this would be needed to 
for (species in plot.taxa) {
  
  ## Need to check the OBS column matches up - or do we not need this again?
  message("Writing coord_clean shapefile for ", species)
  tmp <- GBIF.ALA.SPDF[GBIF.ALA.SPDF$TAXON == species, ]
  writeOGR(tmp, dsn = "./data/ANALYSIS/CLEAN_GBIF", paste0(species, '_coord_clean'), driver = "ESRI Shapefile", overwrite_layer = TRUE)
  
}


## Also write out the whole shapefile, just in case
# writeOGR(obj = GBIF.ALA.SPDF, dsn = "./data/ANALYSIS/CLEAN_GBIF", 
#          layer = paste0("CLEAN_SPDF_", save_run), driver = "ESRI Shapefile")





#########################################################################################################################
##  PLOT GBIF SPATIAL OUTLIERS TO CHECK
#########################################################################################################################


#########################################################################################################################
## Create a tibble to supply to coordinate cleaner
# test.geo = SpatialPointsDataFrame(coords      = TEST.GEO[c("lon", "lat")],
#                                   data        = TEST.GEO,
#                                   proj4string = CRS.WGS.84)
# 
# SDM.COORDS  <- test.geo %>% 
#   
#   spTransform(., CRS.WGS.84) %>% 
#   
#   as.data.frame() %>% 
#   
#   select(searchTaxon, lon, lat, CC.OBS, SOURCE) %>%
#   
#   dplyr::rename(species          = searchTaxon,
#                 decimallongitude = lon, 
#                 decimallatitude  = lat) %>%
#   
#   timetk::tk_tbl()
# 
# 
# ## Check
# dim(SDM.COORDS)
# head(SDM.COORDS)
# class(SDM.COORDS)
# summary(SDM.COORDS$decimallongitude)
# identical(SDM.COORDS$index, SDM.COORDS$CC.OBS)
# length(unique(SDM.COORDS$species))
# 
# 
# #########################################################################################################################
# ## Check how many records each species has
# COMBO.LUT <- SDM.COORDS %>% 
#   as.data.frame() %>%
#   select(species) %>%
#   table() %>%
#   as.data.frame() 
# COMBO.LUT <- setDT(COMBO.LUT, keep.rownames = FALSE)[]
# names(COMBO.LUT) = c("species", "FREQUENCY")
# COMBO.LUT = COMBO.LUT[with(COMBO.LUT, rev(order(FREQUENCY))), ]
# 
# 
# ## Watch out here - this sorting could cause problems for the order of the data frame once it's stitched back together
# ## If we we use species to join the data back together, will it preserve the order? 
# LUT.100K = as.character(subset(COMBO.LUT, FREQUENCY < 65000)$species)
# LUT.100K = trimws(LUT.100K [order(LUT.100K)])
# length(LUT.100K)
# 
# 
# ## Unfortunately, the cc_outl function can't handle vectors of a certain size - over 40 GB at least.
# ## So we have to run this after the SDM step
# ## The current settings are not getting enough spatial outliers..........................................................
# 
# ## Create a data frame of species name and spatial outlier
# SPAT.OUT <- LUT.100K  %>%
#   
#   ## pipe the list of species into lapply
#   lapply(function(x) {
#     
#     ## Create the species df by subsetting by species
#     f <- subset(SDM.COORDS, species == x)
#     
#     ## Run the spatial outlier detection
#     message("Running spatial outlier detection for ", nrow(f), " records for ", x)
#     sp.flag <- cc_outl(f,
#                        lon     = "decimallongitude",
#                        lat     = "decimallatitude",
#                        species = "species",
#                        method  = "quantile", #"distance",
#                        #mltpl   = 4,
#                        #tdi     = 300,
#                        value   = "flagged",
#                        verbose = "TRUE")
#     
#     ## Now add attache column for species, and the flag for each record
#     #d = data.frame(searchTaxon = x, SPAT_OUT = sp.flag)
#     d = cbind(searchTaxon = x,
#               SPAT_OUT = sp.flag, f)[c("searchTaxon", "SPAT_OUT", "CC.OBS")]
#     
#     ## Remeber to explicitly return the df at the end of loop, so we can bind
#     return(d)
#     
#   }) %>%
#   
#   ## Finally, bind all the rows together
#   bind_rows
# 
# gc()


# #########################################################################################################################
# ## Try plotting the points which are outliers for a subset of species and label them
# SPAT.FLAG = join(as.data.frame(test.geo), SPAT.OUT)    ## Join means the skipped species are left out
# identical(SPAT.FLAG$searchTaxon, test.geo$coord_spp)                                       ## order matches 
# identical(SPAT.FLAG$CC.OBS,      test.geo$CC.OBS)                                          ## order matches
# dim(SPAT.FLAG)
# 
# 
# CLEAN.OUT = SpatialPointsDataFrame(coords      = SPAT.FLAG[c("lon", "lat")],
#                                    data        = SPAT.FLAG,
#                                    proj4string = CRS.WGS.84)
# 
# 
# #########################################################################################################################
# ## Get the first 10 species
# ## species = plot.taxa[4]
# for (species in plot.taxa[2]) {
#   
#   ## Plot a subset of taxa
#   CLEAN.OUT.PI      = subset(CLEAN.OUT, searchTaxon == species)
#   
#   ## Create a color scheme that will have the same colour for true/false across all maps
#   # cc.out.col           <- brewer.pal(3,"Set1")
#   # names(cc.out.col)    <- levels(as.factor(CLEAN.OUT$SPAT_OUT))
#   # cc.out.col           <- as.data.frame(cc.out.col)
#   # cc.out.col$SPAT_OUT  <- rownames(cc.out.col)
#   # rownames(cc.out.col) <- NULL
#   # names(cc.out.col)    <- c("map_col", "SPAT_OUT")
#   
#   # CLEAN.OUT@data = data.frame(CLEAN.OUT@data, 
#   #                             cc.out.col[match(CLEAN.OUT@data[,"SPAT_OUT"], 
#   #                                              cc.out.col[,"SPAT_OUT"]),])
#   # colcheck          = CLEAN.OUT.PI[c("searchTaxon", "lat", "lon", "SPAT_OUT", "map_col")]
#   
#   message("plotting occ data for ", species, ", ", 
#           nrow(CLEAN.OUT.PI), " records flagged as either ",
#           unique(CLEAN.OUT.PI$SPAT_OUT))
#   
#   #############################################################
#   ## Plot true and false points for the world
#   ## Black == FALSE
#   ## Red   == TRUE
#   message('Writing map of global coord clean records for ', species)
#   png(sprintf("./data/ANALYSIS/CLEAN_GBIF/%s_%s", species, "global_check_cc_out.png"),
#       16, 10, units = 'in', res = 500)
#   
#   par(mfrow = c(1,2))
#   plot(LAND.84, main = paste0(nrow(subset(CLEAN.OUT.PI, SPAT_OUT == "FALSE")),
#                               " Global cc_outl 'FALSE' points for ", species),
#        lwd = 0.01, asp = 1, bg = 'sky blue', col = 'grey')
#   
#   points(CLEAN.OUT.PI,
#          pch  = ".", cex = 3.3, cex.lab = 3, cex.main = 4, cex.axis = 2, 
#          xlab = "", ylab = "", asp = 1,
#          col  = factor(CLEAN.OUT.PI$SPAT_OUT))
#   
#   #############################################################
#   ## Plot true and false points for the world
#   ## Black == FALSE
#   ## Red   == TRUE
#   plot(AUS.84, main = paste0(nrow(subset(CLEAN.OUT.PI, SPAT_OUT == "FALSE")),
#                              " Global cc_outl 'FALSE' points for ", species),
#        lwd = 0.01, asp = 1, bg = 'sky blue', col = 'grey')
#   
#   points(CLEAN.OUT.PI,
#          pch = ".", cex = 3.3, cex.lab = 3, cex.main = 4, cex.axis = 2, 
#          xlab = "", ylab = "", asp = 1,
#          col = factor(CLEAN.OUT.PI$SPAT_OUT)) # CLEAN.OUT.PI$map_col
#   
#   # legend(x      = 'bottomleft', 
#   #        legend = as.character(cc.out.col$SPAT_OUT),
#   #        col    = cc.out.col$map_col, pch = par("pch"), bty = 'n', xjust = 1)
#   
#   dev.off()
#   
# }
# 
# 
# #########################################################################################################################
# ## Save shapefiles data
# ## Select columns to send to the shapefile
# CC.OUT.CHECK  = dplyr::select(SPAT.FLAG,     CC.OBS, searchTaxon, scientificName, lat, lon, SOURCE, year, 
#                               country, locality, basisOfRecord, institutionCode, 
#                               SPAT_OUT)
# 
# 
# ## Rename the fields so that ArcMap can handle them
# CC.OUT.CHECK  = dplyr::rename(CC.OUT.CHECK , 
#                               OBS       = CC.OBS,
#                               TAXON     = searchTaxon,
#                               LAT       = lat,
#                               LON       = lon,
#                               SC_NAME   = scientificName,
#                               BASIS     = basisOfRecord,                
#                               LOCAL     = locality,                      
#                               INSTIT    = institutionCode,                
#                               COUNTRY   = country,                
#                               YEAR      = year,
#                               CC_OUTL   = SPAT_OUT)
# names(GBIF.ALA.CHECK)
# 
# 
# 
# for (species in plot.taxa[2]) {
#   
#   ## Need to check the OBS column matches up - or do we not need this again?
#   message("Writing spatial outlier shapefile for ", species )
#   tmp <- CLEAN.OUT[CLEAN.OUT$searchTaxon == species , ]
#   writeOGR(tmp, dsn = "./data/ANALYSIS/CLEAN_GBIF", paste0(species, '_cc_outl'), 
#            driver = "ESRI Shapefile", overwrite_layer = TRUE)
#   
# }





#########################################################################################################################
## CREATE HISTOGRAMS FOR EACH SPECIES
#########################################################################################################################


## Create a DF for the histograms
HIST.PLOT = SpatialPointsDataFrame(coords      = CLEAN.TRUE[c("lon", "lat")],
                                   data        = CLEAN.TRUE,
                                   proj4string = CRS.WGS.84)


##############################################################################################
## Plot histograms of temperature and rainfall
## species = plot.taxa[9]
for (species in plot.taxa) {
  
  ## Subset the spatial dataframe into records for each species
  SP.DF  <- subset(HIST.PLOT, searchTaxon == species)
  
  ## Subset DF into records for each species
  DF     <- subset(CLEAN.TRUE, searchTaxon == species)
  DF.OCC <- subset(CLEAN.TRUE, searchTaxon == species & SOURCE != "INVENTORY")
  DF.INV <- subset(CLEAN.TRUE, searchTaxon == species & SOURCE == "INVENTORY")
  
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
#####################################################  TBC ##############################################################
#########################################################################################################################