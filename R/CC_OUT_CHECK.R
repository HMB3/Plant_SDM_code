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
## Check how many records each species has
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
## If we we use species to join the data back together, will it preserve the order? 
LUT.100K = as.character(subset(COMBO.LUT, FREQUENCY < 65000)$species)
LUT.100K = trimws(LUT.100K [order(LUT.100K)])
length(LUT.100K)


## Unfortunately, the cc_outl function can't handle vectors of a certain size - over 40 GB at least.
## So we have to run this afterthe SDM step
## The current settings are not getting enough spatial outliers..........................................................

## Create a data frame of species name and spatial outlier
SPAT.OUT <- LUT.100K  %>%
  
  ## pipe the list of species into lapply
  lapply(function(x) {
    
    ## Create the species df by subsetting by species
    f <- subset(SDM.COORDS, species == x)
    
    ## Run the spatial outlier detection
    message("Running spatial outlier detection for ", x)
    message(dim(f)[1], " records for ", x)
    sp.flag <- cc_outl(f,
                       lon     = "decimallongitude",
                       lat     = "decimallatitude",
                       species = "species",
                       method  = "quantile", #"distance",
                       #mltpl   = 4,
                       #tdi     = 300,
                       value   = "flagged",
                       verbose = "FALSE")
    
    ## Now add attache column for species, and the flag for each record
    #d = data.frame(searchTaxon = x, SPAT_OUT = sp.flag)
    d = cbind(searchTaxon = x,
              SPAT_OUT = sp.flag, f)[c("searchTaxon", "SPAT_OUT", "CC.OBS")]
    
    ## Remeber to explicitly return the df at the end of loop, so we can bind
    return(d)
    
  }) %>%
  
  ## Finally, bind all the rows together
  bind_rows

gc()


#########################################################################################################################
##
SPAT.FLAG = join(as.data.frame(test.geo), SPAT.OUT)    ## Join means the skipped species are left out
dim(SPAT.FLAG)


#########################################################################################################################
## Try plotting the points which are outliers for a subset of species and label them
CLEAN.PLOT = SpatialPointsDataFrame(coords      = SPAT.FLAG[c("lon", "lat")],
                                    data        = SPAT.FLAG,
                                    proj4string = CRS.WGS.84)

## Create global and australian shapefile in the local coordinate system
LAND.84 = LAND %>%
  spTransform(CRS.WGS.84)
AUS.84 = AUS %>%
  spTransform(CRS.WGS.84)

#########################################################################################################################
## Get the first 10 species
## species = plot.taxa[9]
plot.taxa <- as.character(unique(CLEAN.TRUE$searchTaxon))
for (species in plot.taxa) {
  
  ## Plot a subset of taxa
  CLEAN.PLOT.PI   = subset(CLEAN.PLOT,       searchTaxon == species)
  
  message("plotting occ data for ", species, ", ", 
          nrow(CLEAN.TRUE.P.I), " clean records")
  
  #############################################################
  ## Plot true and false points for the world
  ## Black == FALSE
  ## Red   == TRUE
  message('Writing map of global coord clean records for ', species)
  png(sprintf("./data/ANALYSIS/CLEAN_GBIF/CHECK/%s_%s", species, "global_coord_check_out.png"),
      16, 10, units = 'in', res = 500)
  
  par(mfrow = c(1,2))
  plot(LAND.84, main = paste0("Global points for ", species),
       lwd = 0.01, asp = 1, col = 'grey', bg = 'sky blue')
  
  points(CLEAN.PLOT.PI,
         pch = ".", cex = 3.3, cex.lab = 3, cex.main = 4, cex.axis = 2, 
         xlab = "", ylab = "", asp = 1,
         col = factor(CLEAN.PLOT.PI$SPAT_OUT))
  
  #############################################################
  ## Plot true and false points for the world
  ## Black == FALSE
  ## Red   == TRUE
  plot(AUS.84, main = paste0("Australian points for ", species),
       lwd = 0.01, asp = 1, bg = 'sky blue', col = 'grey')
  
  points(CLEAN.PLOT.PI,
         pch = ".", cex = 3.3, cex.lab = 3, cex.main = 4, cex.axis = 2, 
         xlab = "", ylab = "", asp = 1,
         col = factor(CLEAN.PLOT.PI$SPAT_OUT))
  
  dev.off()
  
}



