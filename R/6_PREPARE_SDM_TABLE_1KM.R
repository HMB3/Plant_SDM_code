#########################################################################################################################
#############################################  PREPARE DATA FOR MAXENT ################################################## 
#########################################################################################################################


#########################################################################################################################
## This code prepares the occurrence data for the SDM analyses, by selecting only one record per 1km grid cell, and 
## saving the output as spatial points data frame - the format required by the dismo fucntion


## The current cc_outl settings are not getting enough spatial outliers...................................................

## Print the species run to the screen
message('Preparing SDM table for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")


#########################################################################################################################
## Read in all data to run the SDM code :: species lists, shapefile, rasters & tables
if(read_data == "TRUE") {
  
  ## read in RDS files from previous step
  CLEAN.TRUE = readRDS(paste0(DATA_path, 'CLEAN_TRUE_', save_run, '.rds'))
  length(intersect(GBIF.spp, unique(CLEAN.TRUE$searchTaxon)))
  rasterTmpFile()
  
} else {
  
  message(' skip file reading, not many species analysed')   ##
  
}


#########################################################################################################################
##  some proj libs do not like the 'ESRI' string code - set this here
sp_epsg54009 = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs84=0,0,0"





#########################################################################################################################
## 1). RESTRICT DATA TABLE TO ONE RECORD PER 1KM GRID CELL
#########################################################################################################################


#########################################################################################################################
## Create a table with all the variables needed for SDM analysis
## CLEAN.TRUE = readRDS(paste0(DATA_path, 'COMBO_RASTER_ALL_WPW_TEST', '.rds'))
## CLEAN.TRUE = CLEAN.TRUE[CLEAN.TRUE$searchTaxon %in% GBIF.spp, ]
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
COMBO.RASTER.ALL                 <- spTransform(COMBO.RASTER.ALL, CRS(sp_epsg54009))


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
## 1). TRY CLEANING FILTERED DATA FOR SPATIAL OUTLIERS
#########################################################################################################################


## The cc_outl function is a bottleneck. I hoped to make the code run for one species at a time - this works locally,
## but can't get it working easily on katana for all the species....................................................


## Check dimensions
dim(SDM.DATA.ALL)
length(unique(SDM.DATA.ALL$searchTaxon))
length(unique(SDM.DATA.ALL$OBS))

identical(head(SDM.DATA.ALL$OBS, 100), head(CLEAN.TRUE$OBS, 100))   ## should be false - SDM.DATA is a subset of CLEAN.TRUE
identical(tail(SDM.DATA.ALL$OBS, 100), tail(CLEAN.TRUE$OBS, 100))
unique(SDM.DATA.ALL$SOURCE)


#########################################################################################################################
## Create a tibble to supply to coordinate cleaner
SDM.COORDS  <- SDM.DATA.ALL %>% 
  
  spTransform(., CRS.WGS.84) %>% 
  
  as.data.frame() %>% 
  
  select(searchTaxon, lon, lat, OBS, SOURCE) %>%
  
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
                       method  = "distance",
                       #mltpl   = 5,
                       tdi     = 300,
                       value   = "flags",
                       verbose = "FALSE")
    
    ## Now add attache column for species, and the flag for each record
    #d = data.frame(searchTaxon = x, SPAT_OUT = sp.flag)
    d = cbind(searchTaxon = x,
              SPAT_OUT = sp.flag, f)[c("searchTaxon", "SPAT_OUT", "index", "OBS")]
    
    ## Remeber to explicitly return the df at the end of loop, so we can bind
    return(d)
    
  }) %>%
  
  ## Finally, bind all the rows together
  bind_rows

gc()


## These settings produce too few outliers. Try changing the settings..................................................
print(table(SPAT.OUT$SPAT_OUT, exclude = NULL))
identical(SPAT.OUT$index, SPAT.OUT$OBS)
length(unique(SPAT.OUT$searchTaxon))
head(SPAT.OUT)


## Rename the spatial outlier species column, as another identifier
names(SPAT.OUT)[names(SPAT.OUT) == 'searchTaxon'] <- 'SPAT_SPP'


#########################################################################################################################
## save data
if(save_data == "TRUE") {
  
  ## save .rds file for the next session
  saveRDS(SPAT.OUT, paste0(DATA_path, 'ALA_GBIF_SPAT_OUT_', save_run, '.rds'))
  
} else {
  
  message(' skip file saving, not many species analysed')   ##
  
}





#########################################################################################################################
## 2). FILTER DATA TO THE CLEAN RECORDS
#########################################################################################################################


#########################################################################################################################
## Join data :: Best to use the 'OBS' column here
identical(dim(SDM.COORDS)[1],dim(SPAT.OUT)[1])
identical(dim(SDM.DATA.ALL)[1],dim(SPAT.OUT)[1])
identical(SDM.DATA.ALL$searchTaxon, SPAT.OUT$SPAT_SPP)
length(unique(SPAT.OUT$SPAT_SPP))


SPAT.FLAG = join(as.data.frame(SDM.DATA.ALL), SPAT.OUT)    ## Join means the skipped species are left out
dim(SPAT.FLAG)


## Check the join is working 
message('Checking spatial flags for ', length(unique(SPAT.FLAG$searchTaxon)), ' species in the set ', "'", save_run, "'")
print(table(SPAT.FLAG$SPAT_OUT, exclude = NULL))
length(unique(SPAT.FLAG$searchTaxon))
unique(SPAT.FLAG$SOURCE)
unique(SPAT.FLAG$SPAT_OUT)


## Just get the records that were not spatial outliers
SDM.SPAT.ALL = subset(SPAT.FLAG, SPAT_OUT == "TRUE")
unique(SDM.SPAT.ALL$SPAT_OUT)   
unique(SDM.SPAT.ALL$SOURCE) 
length(unique(SDM.SPAT.ALL$searchTaxon))


## What percentage of records are retained?
length(unique(SDM.SPAT.ALL$searchTaxon))
message(round(dim(SDM.SPAT.ALL)[1]/dim(SPAT.FLAG)[1]*100, 2), " % records retained")                                               


#########################################################################################################################
## Convert back to format for SDMs
SDM.SPAT.ALL    = SpatialPointsDataFrame(coords      = SDM.SPAT.ALL[c("lon", "lat")],
                                         data        = SDM.SPAT.ALL,
                                         proj4string = CRS(sp_epsg54009))
projection(SDM.SPAT.ALL)


#########################################################################################################################
## save data
if(save_data == "TRUE") {
  
  ## Save .rds file for the next session
  saveRDS(SPAT.FLAG, paste0(DATA_path, 'SPAT_FLAG_', save_run, '.rds'))
  
  writeOGR(obj    = SDM.SPAT.ALL, 
           dsn    = SHP_path, 
           layer  = paste0('SPAT_OUT_CHECK_', save_run),
           driver = "ESRI Shapefile", overwrite_layer = TRUE)
  
} else {
  
  message(' skip file saving, not many species analysed')   ##
  
}

## get rid of some memory
gc()





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


#########################################################################################################################
## Write the shapefile out
if(save_data == "TRUE") {
  
  ## save .shp for future refrence 
  writeOGR(obj    = SPAT.OUT.SPDF, 
           dsn    = DATA_path, 
           layer  = paste0('SPAT_OUT_CHECK_', save_run),
           driver = "ESRI Shapefile", overwrite_layer = TRUE)
  
} else {
  
  message(' skip file saving, not many species analysed')   ##
  
}







#########################################################################################################################
## 3). CREATE BACKGROUND POINTS AND VARIBALE NAMES
#########################################################################################################################


#########################################################################################################################
##  Dodgy, but need to make sure there are no ESRI: codes in the coord systems
projection(background) = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs84=0,0,0"


#########################################################################################################################
## Add in random records from previously saved runs
background = background [!background$searchTaxon %in% GBIF.spp, ]               ## Don't add records for other species
length(unique(background$searchTaxon));dim(background)
intersect(unique(background$searchTaxon), GBIF.spp)


SDM.SPAT.ALL.DF = as.data.frame(SDM.SPAT.ALL)
round(with(SDM.SPAT.ALL.DF, table(SOURCE)/sum(table(SOURCE))*100), 1)
round(with(background, table(SOURCE)/sum(table(SOURCE))*100), 1)


#########################################################################################################################
## Make the columns match......................................................
names(background)
names(SDM.SPAT.ALL)
setdiff(names(SDM.SPAT.ALL), names(background))

background$OBS      = "BG"
background$SOURCE   = "BG"
background$SPAT_OUT = "BG"
background$index    = "BG"
background$SPAT_SPP = "BG" 


drops <- c("lon", "lat") # list of col names
SDM.SPAT.ALL <- SDM.SPAT.ALL[,!(names(SDM.SPAT.ALL) %in% drops)]
setdiff(names(SDM.SPAT.ALL), names(background))


#########################################################################################################################
## Now bind on the background points
#browser()

SDM.SPAT.ALL = rbind(SDM.SPAT.ALL, background)
projection(SDM.SPAT.ALL);
projection(background)


## And check what the
unique(SDM.SPAT.ALL$SOURCE) 
length(unique(SDM.SPAT.ALL$searchTaxon))


#########################################################################################################################
## save data
if(save_data == "TRUE") {
  
  ## save .rds file for the next session
  saveRDS(SDM.SPAT.ALL, paste0(DATA_path, 'SDM_SPAT_ALL_',  save_run, '.rds'))
  
} else {
  
  message(' skip file saving, not many species analysed')   ##
  
}

## get rid of some memory
gc()





#########################################################################################################################
##################################################### TBC ############################################################### 
#########################################################################################################################