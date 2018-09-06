#########################################################################################################################
#############################################  PREPARE DATA FOR MAXENT ################################################## 
#########################################################################################################################


#########################################################################################################################
## This code prepares the occurrence data for the SDM analyses


#########################################################################################################################
## Read in all data to run the SDM code :: species lists, shapefile, rasters & tables
#CLEAN.TRUE = readRDS('data/base/HIA_LIST/COMBO/CLEAN_TRUE_OLD_ALA.rds')
rasterTmpFile()


#########################################################################################################################
## 1). PREPARE DATA TABLE FOR SDM ANALYSIS
#########################################################################################################################


#########################################################################################################################
## Create a table with all the variables needed for SDM analysis
dim(CLEAN.TRUE)
length(unique(CLEAN.TRUE$searchTaxon))
length(unique(CLEAN.TRUE$OBS))
unique(CLEAN.TRUE$SOURCE)


## Select only the columns needed
COMBO.RASTER.ALL  <- dplyr::select(CLEAN.TRUE, searchTaxon, lon, lat, SOURCE, OBS,
                                   
                                   Annual_mean_temp,     Mean_diurnal_range,  Isothermality,     Temp_seasonality, 
                                   Max_temp_warm_month,  Min_temp_cold_month, Temp_annual_range, Mean_temp_wet_qu,
                                   Mean_temp_dry_qu,     Mean_temp_warm_qu,   Mean_temp_cold_qu, 
                                   
                                   Annual_precip,        Precip_wet_month,    Precip_dry_month,  Precip_seasonality,   
                                   Precip_wet_qu,        Precip_dry_qu,       Precip_warm_qu,    Precip_col_qu)


#########################################################################################################################
## Create a spatial points object, and change to a projected system to calculate distance more accurately 
coordinates(COMBO.RASTER.ALL)    <- ~lon+lat
proj4string(COMBO.RASTER.ALL)    <- '+init=epsg:4326'
COMBO.RASTER.ALL                 <- spTransform(COMBO.RASTER.ALL, CRS('+init=ESRI:54009'))


## Now split using the data using the species column, and get the unique occurrence cells
COMBO.RASTER.SPLIT.ALL <- split(COMBO.RASTER.ALL, COMBO.RASTER.ALL$searchTaxon)
occurrence_cells_all   <- lapply(COMBO.RASTER.SPLIT.ALL, function(x) cellFromXY(template.raster, x))
length(occurrence_cells_all)   ## this is a list of dataframes, where the number of rows for each being the species table


#########################################################################################################################
## Now get just one record within each 10*10km cell.
SDM.DATA.ALL <- mapply(function(x, cells) {
  x[!duplicated(cells), ]
}, COMBO.RASTER.SPLIT.ALL, occurrence_cells_all, SIMPLIFY = FALSE) %>% do.call(rbind, .)


## Check to see we have 19 variables + the species for the standard predictors, and 19 for all predictors
## Check data :: template, data table and species 
dim(template.raster)
dim(SDM.DATA.ALL)
names(SDM.DATA.ALL)
length(unique(SDM.DATA.ALL$searchTaxon))
class(Koppen_1975)


## What resolution is the template raster at?
xres(template.raster);yres(template.raster)





#########################################################################################################################
## 2). TRY CLEANING FILTERED DATA FOR SPATIAL OUTLIERS
#########################################################################################################################


## Check dimensions
dim(SDM.DATA.ALL)
length(unique(SDM.DATA.ALL$searchTaxon))
length(unique(SDM.DATA.ALL$OBS))

identical(head(SDM.DATA.ALL$OBS, 100), head(CLEAN.TRUE$OBS))
identical(tail(SDM.DATA.ALL$OBS, 100), tail(CLEAN.TRUE$OBS))
unique(SDM.DATA.ALL$SOURCE)


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
dim(SDM.COORDS)
head(SDM.COORDS)
class(SDM.COORDS)
summary(SDM.COORDS$decimallongitude)
identical(SDM.COORDS$index, SDM.COORDS$OBS)


#########################################################################################################################
## Check how many records each species has
COMBO.LUT <- SDM.COORDS %>% 
  as.data.frame() %>%
  select(species) %>%
  table() %>%
  as.data.frame(row.names = TRUE) 
COMBO.LUT <- setDT(COMBO.LUT, keep.rownames = TRUE)[]
names(COMBO.LUT) = c("species", "FREQUENCY")
COMBO.LUT = COMBO.LUT[with(COMBO.LUT, rev(order(FREQUENCY))), ] 


## Watch out here - this sorting could cause problems for the order of the data frame once it's stitched back together
## If we we use species to join the data back together, will it preserve the order? 
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
                       tdi     = 800,
                       value   = "flags")
    
    ## Now add attache column for species, and the flag for each record
    #d = data.frame(searchTaxon = x, SPAT_OUT = sp.flag)
    d = cbind(searchTaxon = x,
              SPAT_OUT = sp.flag, f)[c("searchTaxon", "SPAT_OUT", "index", "OBS")]
    
    ## Remeber to explicitly return the df at the end of loop, so we can bind
    return(d)
    
  }) %>%
  
  ## Finally, bind all the rows together
  bind_rows


## These settings produce too many outliers. Try changing the settings..................................................
identical(SPAT.OUT$index, SPAT.OUT$OBS)
length(unique(SPAT.OUT$searchTaxon))
head(SPAT.OUT)
table(SPAT.OUT$SPAT_OUT)
saveRDS(SPAT.OUT, paste0('data/base/HIA_LIST/COMBO/ALA_GBIF_SPAT_OUT_', save_run, '.rds'))
#SPAT.OUT = readRDS('data/base/HIA_LIST/COMBO/ALA_GBIF_SPAT_OUT_OLD_ALA.rds')





#########################################################################################################################
## 3). FILTER DATA TO THE CLEAN RECORDS
#########################################################################################################################


## Check the species which were not flagged as having spatial outliers are not missing by OUT = TRUE


#########################################################################################################################
## Join data :: Best to use the 'OBS' column here
## Not all species can be processed
intersect(names(SDM.DATA.ALL), names(SPAT.OUT))
SPAT.FLAG = merge(as.data.frame(SDM.DATA.ALL), SPAT.OUT, all = TRUE)
SPAT.FLAG$SPAT_OUT <- ifelse(is.na(SPAT.FLAG$SPAT_OUT), 'NOT_RUN', SPAT.FLAG$SPAT_OUT)


table(SPAT.FLAG$SPAT_OUT)
round(with(SPAT.FLAG, table(SPAT_OUT)/sum(table(SPAT_OUT))*100), 2)


identical(dim(SDM.COORDS)[1],dim(SPAT.FLAG)[1])
dim(SPAT.FLAG)
length(unique(SPAT.OUT$searchTaxon))
unique(SPAT.OUT$SPAT_OUT)
unique(SPAT.FLAG$SOURCE)


## Save the flags before subsetting
saveRDS(SPAT.FLAG, paste0('data/base/HIA_LIST/COMBO/SPAT_FLAG_', save_run, '.rds'))


## Just get the records that were not spatial outliers
SDM.SPAT.ALL = subset(SPAT.FLAG, SPAT_OUT == "TRUE" | SPAT_OUT == "NOT_RUN")  # & GBIF.SPAT.OUT == "TRUE")
unique(SDM.SPAT.ALL$SPAT_OUT)   


## What percentage of records are retained?
length(unique(SDM.SPAT.ALL$searchTaxon))
message(round(dim(SDM.SPAT.ALL)[1]/dim(SPAT.FLAG)[1]*100, 2), " % records retained")                                               





#########################################################################################################################
## 4). CREATE SHAPEFILES TO CHECK OUTLIERS ARCMAP
#########################################################################################################################


## Rename the fields so that ArcMap can handle them
SPAT.OUT.CHECK     = SPAT.FLAG %>% 
  select(OBS, searchTaxon, lat, lon, SOURCE, SPAT_OUT, index) %>%
  dplyr::rename(TAXON     = searchTaxon,
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
         driver = "ESRI Shapefile", overwrite_layer = TRUE)


## Convert back to format for SDMs
SDM.SPAT.ALL    = SpatialPointsDataFrame(coords       = SDM.SPAT.ALL[c("lon", "lat")],
                                         data         = SDM.SPAT.ALL,
                                         proj4string  = CRS('+init=ESRI:54009'))
projection(SDM.SPAT.ALL)





#########################################################################################################################
## 3). CREATE BACKGROUND POINTS AND VARIBALE NAMES
#########################################################################################################################


#########################################################################################################################
## Add in random records from previously saved runs
background = readRDS("./data/base/HIA_LIST/COMBO/SDM_DATA_CLEAN_052018.rds")
background = background [!background$searchTaxon %in% GBIF.spp, ]               ## Don't add records for other species
length(unique(background$searchTaxon));dim(background)
intersect(unique(background$searchTaxon), GBIF.spp)


## Make the columns match......................................................
names(background)
names(SDM.SPAT.ALL)
setdiff(names(SDM.SPAT.ALL), names(background))

background$OBS      = "BG"
background$SOURCE   = "BG"
background$SPAT_OUT = "BG"
background$index    = "BG"


drops <- c("lon", "lat") # list of col names
SDM.SPAT.ALL <- SDM.SPAT.ALL[,!(names(SDM.SPAT.ALL) %in% drops)]
setdiff(names(SDM.SPAT.ALL), names(background))


## Now bind on the background points............................................................................................
projection(SDM.SPAT.ALL);projection(background)
SDM.SPAT.ALL = rbind(SDM.SPAT.ALL, background)


#########################################################################################################################
## Save the final SDM dataset
saveRDS(SDM.SPAT.ALL, paste0('data/base/HIA_LIST/COMBO/SDM_SPAT_ALL_', save_run, '.rds'))








#########################################################################################################################
##################################################### TBC ############################################################### 
#########################################################################################################################