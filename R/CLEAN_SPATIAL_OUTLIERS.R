#########################################################################################################################
################################################# FLAG SPATIAL OUTLIERS #################################################
#########################################################################################################################


## This code uses the CoordinateCleaner package to Automatcially screen out the dodgy spatial records. See ::
## https://github.com/azizka/CoordinateCleaner


#########################################################################################################################
## 1). TRY CLEANING FILTERED DATA FOR SPATIAL OUTLIERS
#########################################################################################################################


#########################################################################################################################
## Create a tibble to supply to coordinate cleaner
SDM.COORDS  <- SDM.DATA.ALL %>% 
  
  spTransform(., CRS.WGS.84) %>% 
  
  as.data.frame() %>% 
  
  select(searchTaxon, lon, lat, OBS) %>%
  
  dplyr::rename(species          = searchTaxon,
                decimallongitude = lon, 
                decimallatitude  = lat) %>%
  
  timetk::tk_tbl()


## Check
head(SDM.COORDS)
class(SDM.COORDS)
summary(SDM.COORDS$decimallongitude)
identical(SDM.COORDS$index, SDM.COORDS$OBS)


#########################################################################################################################
## Check how big the table is
COMBO.LUT <- SDM.COORDS %>% 
  as.data.frame() %>%
  select(species) %>%
  table() %>%
  as.data.frame(row.names = TRUE) 
COMBO.LUT <- setDT(COMBO.LUT, keep.rownames = TRUE)[]
names(COMBO.LUT) = c("species", "FREQUENCY")
COMBO.LUT = COMBO.LUT[with(COMBO.LUT, rev(order(FREQUENCY))), ] 
View(COMBO.LUT)


LUT.100K = as.character(subset(COMBO.LUT, FREQUENCY < 65000)$species)
LUT.100K = LUT.100K [order(LUT.100K)]
length(LUT.100K)


## Unfortunately, the cc_outl function can't handle vectors of a certain size - over 40 GB at least.
## So we are problably better off moving this code to after the data has been thinned by the SDM step


## Create a data frame of species name and spatial outlier
SPAT.OUT <- LUT.100K %>%  ## unique(TIB.GBIF$species)
  
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
                       method  = "distance",
                       #mltpl   = 5,
                       tdi     = 300,
                       value   = "flags")
    
    ## Now add attache column for species, and the flag for each record
    #d = data.frame(searchTaxon = x, SPAT_OUT = sp.flag)
    d = cbind(searchTaxon = x,
              SPAT_OUT = sp.flag, f)[c("searchTaxon", "SPAT_OUT")]
    
    ## Remeber to explicitly return the df at the end of loop, so we can bind
    return(d)
    
  }) %>%
  
  ## Finally, bind all the rows together
  bind_rows


## These settings produce too many outliers. Try changing the settings..................................................
identical(dim(SPAT.OUT)[1], dim(SDM.COORDS)[1])
unique(SPAT.OUT$searchTaxon)
head(SPAT.OUT)
saveRDS(SPAT.OUT, paste0('data/base/HIA_LIST/COMBO/ALA_GBIF_SPAT_OUT_', save_run, '.rds'))





#########################################################################################################################
## 2). FILTER DATA TO THE CLEAN RECORDS
#########################################################################################################################


#########################################################################################################################
## Join data :: exclude the decimal lat/long, check the length 
identical(dim(SDM.COORDS)[1],dim(SPAT.OUT)[1])#;length(GBIF.SPAT.OUT)
names(SPAT.OUT)[1] = c("coord_spp")
identical(SDM.COORDS$searchTaxon, SPAT.OUT$coord_spp)                                                  ## order matches


#########################################################################################################################
## Is the species column the same as the searchTaxon column?
SPAT.FLAG = merge(SDM.COORDS, SPAT.OUT, by = "searchTaxon")
dim(SPAT.FLAG)
summary(SPAT.FLAG)
identical(SPAT.FLAG$searchTaxon, SPAT.FLAG$coord_spp)                                                    ## order matches


## Not sure why the inverse did not work :: get only the records which were not flagged as being dodgy.
dim(subset(SPAT.FLAG, summary == "TRUE"))         #  | GBIF.SPAT.OUT == "TRUE"))
SPAT.TRUE = subset(SPAT.FLAG, summary == "TRUE")  # & GBIF.SPAT.OUT == "TRUE")
unique(SPAT.FLAG$SPAT_OUT)   
  

## What percentage of records are retained?
length(unique(SPAT.TRUE$searchTaxon))
message(round(dim(SPAT.TRUE)[1]/dim(SDM.COORDS)[1]*100, 2), " % records retained")                                               


## Save the spatial flags
saveRDS(SPAT.FLAG, paste0('data/base/HIA_LIST/COMBO/SPAT_FLAG_', save_run, '.rds'))





#########################################################################################################################
## 4). CREATE SHAPEFILES TO CHECK OUTLIERS ARCMAP
#########################################################################################################################


## Rename the fields so that ArcMap can handle them
SPAT.OUT.CHECK     = dplyr::rename(SPAT.FLAG, 
                                   TAXON     = searchTaxon,
                                   LAT       = lat,
                                   LON       = lon)
names(SPAT.OUT.CHECK)


#########################################################################################################################
## Then create a SPDF
SPAT.OUT.SPDF    = SpatialPointsDataFrame(coords      = SPAT.OUT.CHECK[c("LON", "LAT")],
                                          data        = SPAT.OUT.CHECK,
                                          proj4string = CRS.WGS.84)


## Write the shapefile out
writeOGR(obj    = SPAT.OUT.SPDF, 
         dsn    = "./data/base/HIA_LIST/COMBO", 
         layer  = paste0('SPAT_OUT_CHECK_', save_run),
         driver = "ESRI Shapefile")





#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################