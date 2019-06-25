#########################################################################################################################
############################################# CHECK COORDCLEAN SETTINGS #################################################
#########################################################################################################################


#########################################################################################################################
##  PLOT GBIF ERRORS TO CHECK
#########################################################################################################################


#########################################################################################################################
## Try plotting the points which are outliers for a subset of spp and label them
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
## Get the first 10 spp
## spp = plot.taxa[1]
## spp = plot.taxa[1]
plot.taxa <- as.character(unique(CLEAN.PLOT$searchTaxon))
for (spp in plot.taxa) {
  
  ## Plot a subset of taxa
  CLEAN.PLOT.PI = CLEAN.PLOT[ which(CLEAN.PLOT$searchTaxon == spp), ]
  
  message("plotting occ data for ", spp, ", ", 
          nrow(CLEAN.PLOT.PI), " records flagged as either ",
          unique(CLEAN.PLOT.PI$coord_summary))
  
  #############################################################
  ## Plot true and false points for the world
  ## Black == FALSE
  ## Red   == TRUE
  message('Writing map of global coord clean records for ', spp)
  png(sprintf("./data/ANALYSIS/CLEAN_GBIF/%s_%s", spp, "global_check_cc.png"),
      16, 10, units = 'in', res = 500)
  
  par(mfrow = c(1,2))
  plot(LAND.84, main = paste0(nrow(subset(CLEAN.PLOT.PI, coord_summary == "FALSE")),
                              " Global clean_coord 'FALSE' points for ", spp),
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
                             " Global clean_coord 'FALSE' points for ", spp),
       lwd = 0.01, asp = 1, bg = 'sky blue', col = 'grey')
  
  points(CLEAN.PLOT.PI,
         pch = ".", cex = 3.3, cex.lab = 3, cex.main = 4, cex.axis = 2, 
         xlab = "", ylab = "", asp = 1,
         col = factor(CLEAN.PLOT.PI$coord_summary))
  
  dev.off()
  
}





#########################################################################################################################
##  PLOT SPATIAL OUTLIERS TO CHECK
#########################################################################################################################


#########################################################################################################################
## Create a tibble to supply to coordinate cleaner
test.geo = SpatialPointsDataFrame(coords      = TEST.GEO[c("lon", "lat")],
                                  data        = TEST.GEO,
                                  proj4string = CRS.WGS.84)

SDM.COORDS  <- test.geo %>% 
  
  spTransform(., CRS.WGS.84) %>% 
  
  as.data.frame() %>% 
  
  select(searchTaxon, lon, lat, CC.OBS, SOURCE) %>%
  
  dplyr::rename(species          = searchTaxon,
                decimallongitude = lon, 
                decimallatitude  = lat) %>%
  
  timetk::tk_tbl()


## Check
dim(SDM.COORDS)
head(SDM.COORDS)
class(SDM.COORDS)
summary(SDM.COORDS$decimallongitude)
identical(SDM.COORDS$index, SDM.COORDS$CC.OBS)
length(unique(SDM.COORDS$species))


#########################################################################################################################
## Check how many records each spp has
COMBO.LUT <- SDM.COORDS %>% 
  as.data.frame() %>%
  select(species) %>%
  table() %>%
  as.data.frame() 
COMBO.LUT <- setDT(COMBO.LUT, keep.rownames = FALSE)[]
names(COMBO.LUT) = c("species", "FREQUENCY")
COMBO.LUT = COMBO.LUT[with(COMBO.LUT, rev(order(FREQUENCY))), ]
View(COMBO.LUT)


## Watch out here - this sorting could cause problems for the order of the data frame once it's stitched back together
## If we we use spp to join the data back together, will it preserve the order? 
LUT.100K = as.character(subset(COMBO.LUT, FREQUENCY < 100000)$species)
LUT.100K = trimws(LUT.100K [order(LUT.100K)])
length(LUT.100K)


## Unfortunately, the cc_outl function can't handle vectors of a certain size - over 40 GB at least.
## So we have to run this afterthe SDM step
## The current settings are not getting enough spatial outliers..........................................................

## Create a data frame of spp name and spatial outlier
SPAT.OUT <- LUT.100K  %>%
  
  ## pipe the list of spp into lapply
  lapply(function(x) {
    
    ## Create the spp df by subsetting by spp
    f <- subset(SDM.COORDS, species == x)
    
    ## Run the spatial outlier detection
    message("Running spatial outlier detection for ", x)
    message(dim(f)[1], " records for ", x)
    sp.flag <- cc_outl(f,
                       lon     = "decimallongitude",
                       lat     = "decimallatitude",
                       species = "species",
                       method  = "quantile", 
                       mltpl   = 10,
                       value   = "flagged",
                       verbose = "TRUE")
    
    ## Now add attache column for spp, and the flag for each record
    d = cbind(searchTaxon = x,
              SPAT_OUT = sp.flag, f)[c("searchTaxon", "SPAT_OUT", "CC.OBS")]
    
    ## Remeber to explicitly return the df at the end of loop, so we can bind
    return(d)
    
  }) %>%
  
  ## Finally, bind all the rows together
  bind_rows

gc()


#########################################################################################################################
## Join the data back on
SPAT.FLAG = join(as.data.frame(test.geo), SPAT.OUT)    ## Join means the skipped spp are left out
dim(SPAT.FLAG)


#########################################################################################################################
## Try plotting the points which are outliers for a subset of spp and label them
SPAT.FLAG = SpatialPointsDataFrame(coords      = SPAT.FLAG[c("lon", "lat")],
                                    data        = SPAT.FLAG,
                                    proj4string = CRS.WGS.84)


#########################################################################################################################
## Get the first 10 spp
## spp = plot.taxa[1]
plot.taxa <- as.character(unique(SPAT.FLAG$searchTaxon))
for (spp in plot.taxa) {
  
  ## Plot a subset of taxa
  CLEAN.PLOT.PI   = subset(SPAT.FLAG, searchTaxon == spp)
  
  message("plotting occ data for ", spp, ", ", 
          nrow(CLEAN.PLOT.PI ), " clean records")
  
  #############################################################
  ## Plot true and false points for the world
  ## Black == FALSE
  ## Red   == TRUE
  message('Writing map of global coord clean records for ', spp)
  png(sprintf("./data/ANALYSIS/CLEAN_GBIF/%s_%s", spp, "global_spatial_outlier_check.png"),
      16, 10, units = 'in', res = 500)
  
  par(mfrow = c(1,2))
  plot(LAND.84, main = paste0(nrow(subset(CLEAN.PLOT.PI, SPAT_OUT == "FALSE")),
                              "Spatial outlier 'FALSE' points for ", spp),
       lwd = 0.01, asp = 1, col = 'grey', bg = 'sky blue')
  
  points(CLEAN.PLOT.PI,
         pch = ".", cex = 3.3, cex.lab = 3, cex.main = 4, cex.axis = 2, 
         xlab = "", ylab = "", asp = 1,
         col = factor(CLEAN.PLOT.PI$SPAT_OUT))
  
  #############################################################
  ## Plot true and false points for the world
  ## Black == FALSE
  ## Red   == TRUE
  plot(AUS.84, main = paste0("Australian points for ", spp),
       lwd = 0.01, asp = 1, bg = 'sky blue', col = 'grey')
  
  points(CLEAN.PLOT.PI,
         pch = ".", cex = 3.3, cex.lab = 3, cex.main = 4, cex.axis = 2, 
         xlab = "", ylab = "", asp = 1,
         col = factor(CLEAN.PLOT.PI$SPAT_OUT))
  
  dev.off()
  
}





#########################################################################################################################
#####################################################  TBC ##############################################################
#########################################################################################################################