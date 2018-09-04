#########################################################################################################################
## 1). TRY CLEANING FILTERED DATA FOR SPATIAL OUTLIERS
#########################################################################################################################


#########################################################################################################################
## Create a tibble to supply to coordinate cleaner
SDM.COORDS  <- SDM.DATA.ALL %>% 
  
  spTransform(., CRS.WGS.84) %>% 
  
  as.data.frame() %>% 
  
  select(searchTaxon, lon, lat) %>%
  
  dplyr::rename(species          = searchTaxon,
                decimallongitude = lon, 
                decimallatitude  = lat) %>%
  
  timetk::tk_tbl()


## Get the 
head(SDM.COORDS)
class(SDM.COORDS)
summary(SDM.COORDS$decimallongitude)


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


LUT.100K = as.character(subset(COMBO.LUT, FREQUENCY < 100000)$species)
LUT.100K = LUT.100K [order(LUT.100K)]


## Unfortunately, the cc_outl function can't handle vectors of a certain size - over 40 GB at least.
## So we are problably better off moving this code to after the data has been thinned by the SDM step


## Create a data frame of species name and spatial outlier
SPAT.OUT <- LUT.100K %>%  ## unique(TIB.GBIF$species)
  
  ## pipe the list of species into lapply
  lapply(function(x) {
    
    ## Create the species df by subsetting by species
    f <- subset(SDM.COORDS, species == x)
    message("processing ", dim(d)[1], " records for ", x)
    
    ## Run the spatial outlier detection
    message("Running spatial outlier detection for ", x)
    sp.flag <- cc_outl(f,
                       lon     = "decimallongitude",
                       lat     = "decimallatitude",
                       species = "species",
                       method  = "quantile",
                       mltpl   = 5,
                       tdi     = 1000,
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


## Save data
# identical(dim(SPAT.OUT)[1], dim(SDM.COORDS)[1])
# unique(SPAT.OUT$searchTaxon)
# head(SPAT.OUT)
# saveRDS(SPAT.OUT, paste0('data/base/HIA_LIST/COMBO/ALA_GBIF_SPAT_OUT_', save_run, '.rds'))





#########################################################################################################################
## 2). FILTER DATA TO THE CLEAN RECORDS
#########################################################################################################################


#########################################################################################################################
## Join data :: exclude the decimal lat/long, check the length 
identical(dim(CLEAN.TRUE)[1],dim(SPAT.OUT)[1])#;length(GBIF.SPAT.OUT)
names(FLAGS)[1] = c("coord_spp")
identical(COMBO.RASTER.CONVERT$searchTaxon, FLAGS$coord_spp)                                            ## order matches
identical(dim(FLAGS)[1], dim(GBIF.TRIM.GEO)[1])


#########################################################################################################################
## Is the species column the same as the searchTaxon column?
TEST.GEO = cbind(COMBO.RASTER.CONVERT, SPAT.OUT)
summary(TEST.GEO)
identical(TEST.GEO$searchTaxon, TEST.GEO$coord_spp)                                                     ## order matches


## Not sure why the inverse did not work :: get only the records which were not flagged as being dodgy.
dim(subset(TEST.GEO, summary == "TRUE")) #  | GBIF.SPAT.OUT == "TRUE"))
CLEAN.TRUE = subset(TEST.GEO, summary == "TRUE") # & GBIF.SPAT.OUT == "TRUE")
unique(CLEAN.TRUE$summary)   
#unique(CLEAN.TRUE$GBIF.SPAT.OUT)   


## What percentage of records are retained?
identical(dim(CLEAN.TRUE)[1], (dim(COMBO.RASTER.CONVERT)[1] - dim(subset(TEST.GEO, summary == "FALSE"))[1]))
length(unique(CLEAN.TRUE$searchTaxon))
message(round(dim(CLEAN.TRUE)[1]/dim(TEST.GEO)[1]*100, 2), " % records retained")                                               


#########################################################################################################################
## Now bind on the urban tree inventory data. We are assuming this data is clean, after we manually fix the taxonomy
## Check the NAs
intersect(names(TI.RASTER.CONVERT), names(CLEAN.TRUE))
CLEAN.TRUE = bind_rows(CLEAN.TRUE, TI.RASTER.CONVERT)


names(CLEAN.TRUE)
unique(CLEAN.TRUE$SOURCE) 
unique(CLEAN.TRUE$INVENTORY) 
length(unique(CLEAN.TRUE$searchTaxon))
summary(CLEAN.TRUE$Annual_mean_temp)


## Then create a unique ID column which can be used to identify outlier records 
CLEAN.TRUE$OBS <- 1:nrow(CLEAN.TRUE)
dim(CLEAN.TRUE)[1];length(CLEAN.TRUE$OBS)  
identical(dim(CLEAN.TRUE)[1], length(CLEAN.TRUE$OBS))


## How many records are added by including the tree inventories?
message("Tree inventory data increases records by ", round(dim(CLEAN.TRUE)[1]/dim(TEST.GEO)[1]*100, 2), " % ")   


## save niches
saveRDS(CLEAN.TRUE, paste0('data/base/HIA_LIST/COMBO/CLEAN_TRUE_', save_run, '.rds'))





#########################################################################################################################
## 4). CREATE SHAPEFILES TO CHECK OUTLIERS ARCMAP
#########################################################################################################################


## Select columns
GBIF.ALA.CHECK  = select(CLEAN.TRUE,     OBS, searchTaxon, scientificName, lat, lon, SOURCE, INVENTORY, year, coordinateUncertaintyInMetres,
                         geodeticDatum,  country, locality, basisOfRecord, institutionCode, rank, Taxonomic.status, New.Taxonomic.status)


## Rename the fields so that ArcMap can handle them
GBIF.ALA.CHECK     = dplyr::rename(GBIF.ALA.CHECK, 
                                   TAXON     = searchTaxon,
                                   LAT       = lat,
                                   LON       = lon,
                                   SC_NAME   = scientificName,
                                   RANK      = rank,
                                   BASIS     = basisOfRecord,                
                                   LOCAL     = locality,                      
                                   INSTIT    = institutionCode,                
                                   COUNTRY   = country,                
                                   COORD_UN  = coordinateUncertaintyInMetres,
                                   DATUM     = geodeticDatum,                 
                                   YEAR      = year)
names(GBIF.ALA.CHECK)


#########################################################################################################################
## Then create SPDF
GBIF.ALA.SPDF    = SpatialPointsDataFrame(coords      = GBIF.ALA.CHECK[c("LON", "LAT")],
                                          data        = GBIF.ALA.CHECK,
                                          proj4string = CRS.WGS.84)




#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################






