#########################################################################################################################
#############################################  PREPARE DATA FOR MAXENT ################################################## 
#########################################################################################################################


#########################################################################################################################
## This code prepares the occurrence data for the SDM analyses, saving the output as spatial points data frame, the format 
## required by the dismo fucntion. Background points are added to from both native australian trees and urban tree inventories.


#########################################################################################################################
## Read in all data to run the SDM code :: species lists, shapefile, rasters & tables
if(read_data == "TRUE") {
  
  ## read in RDS files from previous step
  CLEAN.INV = readRDS(paste0(DATA_path, 'CLEAN_INV_', save_run, '.rds'))
  message('Species overlap ', length(intersect(GBIF.spp, unique(CLEAN.INV$searchTaxon))))
  rasterTmpFile()
  
} else {
  
  message(' skip file reading, not many species analysed')   ##
  
}





#########################################################################################################################
## 1). RESTRICT DATA TABLE TO ONE RECORD PER 1KM GRID CELL
#########################################################################################################################


#########################################################################################################################
## Create a table with all the variables needed for SDM analysis
## This is the step where an ad/hoc version comes in 
## CLEAN.INV <- CLEAN.INV[CLEAN.INV$searchTaxon %in% GBIF.spp, ]
## Print the species run to the screen
message('Preparing SDM table for ', length(unique(CLEAN.INV$searchTaxon)), ' species in the set ', "'", save_run, "'",
        'using ', unique(CLEAN.INV$SOURCE), ' data')


## Select only the columns needed. This also needs to use the variable names
CLEAN.INV         <- CLEAN.INV[CLEAN.INV$searchTaxon %in% GBIF.spp, ]
length(unique(CLEAN.INV$searchTaxon))


COMBO.RASTER.ALL  <- dplyr::select(CLEAN.INV, searchTaxon, lon, lat, SOURCE, CC.OBS,
                                   
                                   Annual_mean_temp,     Mean_diurnal_range,  Isothermality,     Temp_seasonality, 
                                   Max_temp_warm_month,  Min_temp_cold_month, Temp_annual_range, Mean_temp_wet_qu,
                                   Mean_temp_dry_qu,     Mean_temp_warm_qu,   Mean_temp_cold_qu, 
                                   
                                   Annual_precip,        Precip_wet_month,    Precip_dry_month,  Precip_seasonality,   
                                   Precip_wet_qu,        Precip_dry_qu,       Precip_warm_qu,    Precip_col_qu)#,
                                   #PC1_WGS84, PC2_WGS84, PC3_WGS84, TWI, TPI)


#########################################################################################################################
## Create a spatial points object, and change to a projected system to calculate distance more accurately 
## This is the mollweide projection used for the SDMs
coordinates(COMBO.RASTER.ALL)    <- ~lon+lat
proj4string(COMBO.RASTER.ALL)    <- '+init=epsg:4326'
COMBO.RASTER.ALL                 <- spTransform(COMBO.RASTER.ALL, CRS(sp_epsg54009))


## Now split using the data using the species column, and get the unique occurrence cells
COMBO.RASTER.SPLIT.ALL <- split(COMBO.RASTER.ALL, COMBO.RASTER.ALL$searchTaxon)
occurrence_cells_all   <- lapply(COMBO.RASTER.SPLIT.ALL, function(x) cellFromXY(template.raster.1km, x))


## Check with a message, but could check with a fail 
message('Split prodcues ', length(occurrence_cells_all), ' data frames for ', length(GBIF.spp), ' species')   ## This is a list of dataframes 


#########################################################################################################################
## Now get just one record within each 1*1km cell.
SDM.DATA.ALL <- mapply(function(x, cells) {
  x[!duplicated(cells), ]
}, COMBO.RASTER.SPLIT.ALL, occurrence_cells_all, SIMPLIFY = FALSE) %>% do.call(rbind, .)


## Check to see we have 19 variables + the species for the standard predictors, and 19 for all predictors
message(round(nrow(SDM.DATA.ALL)/nrow(CLEAN.INV)*100, 2), " % records retained at 1km resolution")  


dim(template.raster.1km)
dim(SDM.DATA.ALL)
names(SDM.DATA.ALL)
unique(SDM.DATA.ALL$SOURCE)
length(unique(SDM.DATA.ALL$searchTaxon))
class(Koppen_1975_1km)


## What resolution is the template raster at?
xres(template.raster.1km);yres(template.raster.1km)





#########################################################################################################################
## 1). TRY CLEANING FILTERED DATA FOR SPATIAL OUTLIERS
#########################################################################################################################


## The cc_outl function has been tweaked and sped up.


#########################################################################################################################
## Create a unique identifier for spatial cleaning. This is used for automated cleaing of the records, and also saving shapefiles
## But this will not be run for all species linearly. So, it probably needs to be a combination of species and number
SDM.DATA.ALL$SPOUT.OBS <- 1:nrow(SDM.DATA.ALL)
SDM.DATA.ALL$SPOUT.OBS <- paste0(SDM.DATA.ALL$SPOUT.OBS, "_SPOUT_", SDM.DATA.ALL$searchTaxon)
SDM.DATA.ALL$SPOUT.OBS <- gsub(" ",     "_",  SDM.DATA.ALL$SPOUT.OBS, perl = TRUE)
length(SDM.DATA.ALL$SPOUT.OBS);length(unique(SDM.DATA.ALL$SPOUT.OBS))


## Check dimensions
dim(SDM.DATA.ALL)
length(unique(SDM.DATA.ALL$searchTaxon))
length(unique(SDM.DATA.ALL$SPOUT.OBS))
unique(SDM.DATA.ALL$SOURCE)


#########################################################################################################################
## Create a tibble to supply to coordinate cleaner
SDM.COORDS  <- SDM.DATA.ALL %>% 
  
  spTransform(., CRS.WGS.84) %>% 
  
  as.data.frame() %>% 
  
  select(searchTaxon, lon, lat, SPOUT.OBS, SOURCE) %>%
  
  dplyr::rename(species          = searchTaxon,
                decimallongitude = lon, 
                decimallatitude  = lat) %>%
  
  timetk::tk_tbl()


## Check
dim(SDM.COORDS)
head(SDM.COORDS)
class(SDM.COORDS)
summary(SDM.COORDS$decimallongitude)
identical(SDM.COORDS$index, SDM.COORDS$SPOUT.OBS)
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
LUT.100K = as.character(subset(COMBO.LUT, FREQUENCY < 100000)$species)
LUT.100K = trimws(LUT.100K [order(LUT.100K)])
length(LUT.100K)


## See communications with Alex Zizka
## Check the output for patterns - make the settings strict here, as outliers could be bogus after 1km thinning
## Generally, species with heaps of records, especially those with clumped/biased records, get more outliers

## Create a data frame of species name and spatial outlier
SPAT.OUT <- LUT.100K  %>%
  
  ## pipe the list of species into lapply
  lapply(function(x) {
    
    ## Create the species df by subsetting by species
    f <- subset(SDM.COORDS, species == x)
    
    ## Run the spatial outlier detection
    message("Running spatial outlier detection for ", nrow(f), " records for ", x)
    sp.flag <- cc_outl(f,
                       lon     = "decimallongitude",
                       lat     = "decimallatitude",
                       species = "species",
                       method  = "quantile", #"distance",
                       mltpl   = 10,
                       #tdi     = 300,
                       value   = "flagged",
                       verbose = "TRUE")
    
    ## Now add attache column for species, and the flag for each record
    d = cbind(searchTaxon = x,
              SPAT_OUT = sp.flag, f)[c("searchTaxon", "SPAT_OUT", "SPOUT.OBS")]
    
    ## Remeber to explicitly return the df at the end of loop, so we can bind
    return(d)
    
  }) %>%
  
  ## Finally, bind all the rows together
  bind_rows

gc()


## How many species are flagged as spatial outliers?
print(table(SPAT.OUT$SPAT_OUT, exclude = NULL))
length(unique(SPAT.OUT$searchTaxon))
head(SPAT.OUT)





#########################################################################################################################
## 2). FILTER DATA TO REMOVE SPATIAL OUTLIERS
#########################################################################################################################


#########################################################################################################################
## Join data :: Best to use the 'OBS' column here
identical(nrow(SDM.COORDS), nrow(SPAT.OUT))
identical(SDM.DATA.ALL$searchTaxon, SPAT.OUT$searchTaxon)
length(unique(SPAT.OUT$searchTaxon))


## This explicit join is required. Check the species have been analysed in exactly the same order
SPAT.FLAG = join(as.data.frame(SDM.DATA.ALL), SPAT.OUT, by = c("SPOUT.OBS", "searchTaxon") , type = "left", match = "first")    
identical(SDM.DATA.ALL$searchTaxon, SPAT.FLAG$searchTaxon)


## Check the join is working 
message('Checking spatial flags for ', length(unique(SPAT.FLAG$searchTaxon)), ' species in the set ', "'", save_run, "'")
print(table(SPAT.FLAG$SPAT_OUT, exclude = NULL))
length(unique(SPAT.FLAG$searchTaxon))
unique(SPAT.FLAG$SOURCE)
unique(SPAT.FLAG$SPAT_OUT)


## Just get the records that were not spatial outliers.
SDM.SPAT.ALL = subset(SPAT.FLAG, SPAT_OUT == "TRUE")
unique(SDM.SPAT.ALL$SPAT_OUT)   
unique(SDM.SPAT.ALL$SOURCE) 
length(unique(SDM.SPAT.ALL$searchTaxon))


## What percentage of records are retained?
message(round(nrow(SDM.SPAT.ALL)/nrow(SPAT.FLAG)*100, 2), " % records retained after spatial outlier detection")                                               


#########################################################################################################################
## Now calculate the niches here
## source('./R/CALC_1KM_NICHES.R')


#########################################################################################################################
## Convert back to format for SDMs :: use Mollweide projection
SDM.SPAT.ALL    = SpatialPointsDataFrame(coords      = SDM.SPAT.ALL[c("lon", "lat")],
                                         data        = SDM.SPAT.ALL,
                                         proj4string = CRS(sp_epsg54009))
projection(SDM.SPAT.ALL)
message(length(unique(SDM.SPAT.ALL$searchTaxon)), ' species processed through from download to SDM table')


#########################################################################################################################
## save data
if(save_data == "TRUE") {
  
  ## save .rds file for the next session
  saveRDS(SDM.SPAT.ALL, paste0(DATA_path, 'SDM_SPAT_ALL_',  save_run, '.rds'))
  
} else {
  
  message(' skip file saving, not many species analysed')   ##
  
}


## The here, read in the combined data..................................................................................


#########################################################################################################################
## save data
# if(save_data == "TRUE") {
#   
#   ## Save .rds file for the next session
#   saveRDS(SDM.SPAT.ALL, paste0(DATA_path, 'SDM_SPAT_ALL', save_run, '.rds'))
#   
#   writeOGR(obj    = SDM.SPAT.ALL, 
#            dsn    = SHP_path, 
#            layer  = paste0('SPAT_OUT_CHECK_', save_run),
#            driver = "ESRI Shapefile", overwrite_layer = TRUE)
#   
# } else {
#   
#   message(' skip file saving, not many species analysed')   ##
#   
# }
# 
# ## get rid of some memory
# gc()





#########################################################################################################################
## 4). CREATE SHAPEFILES TO CHECK OUTLIERS ARCMAP
#########################################################################################################################


## Rename the fields so that ArcMap can handle them
# SPAT.OUT.CHECK     = SPAT.FLAG %>%
#   select(SPOUT.OBS, searchTaxon, lat, lon, SOURCE, SPAT_OUT) %>%
#   dplyr::rename(TAXON     = searchTaxon,
#                 LAT       = lat,
#                 LON       = lon)
# names(SPAT.OUT.CHECK)
# 
# 
# #########################################################################################################################
# ## Then create a SPDF
# SPAT.OUT.SPDF    = SpatialPointsDataFrame(coords      = SPAT.OUT.CHECK[c("LON", "LAT")],
#                                           data        = SPAT.OUT.CHECK,
#                                           proj4string = CRS.WGS.84)


#########################################################################################################################
## Write the shapefile out
# if(save_data == "TRUE") {
#   
#   ## save .shp for future refrence 
#   writeOGR(obj    = SPAT.OUT.SPDF, 
#            dsn    = "./data/ANALYSIS/CLEAN_GBIF", 
#            layer  = paste0('SPAT_OUT_CHECK_', save_run),
#            driver = "ESRI Shapefile", overwrite_layer = TRUE)
#   
# } else {
#   
#   message(' skip file saving, not many species analysed')   ##
#   
# }
# 




#########################################################################################################################
## 5). FOR SMALL RUNS OF SPECIES(EG HOLLOWS),  CREATE BACKGROUND POINTS AND VARIBALE NAMES
#########################################################################################################################


## Use one data frame for all species analysis, to save mucking around with background points............................
# background.points = readRDS(paste0(DATA_path, 'SDM_SPAT_OCC_BG_ALL_EVREGREEN_JULY_2018.rds'))
# SDM.SPAT.OCC.BG   = rbind(SDM.SPAT.ALL, background.points)


#########################################################################################################################
## Re-run this after running steps 1-5 for all 4k species, and use that for background.points instead
## Add in random records from previously saved runs :: get all the species which have not
# background.points = background.points[!background.points$searchTaxon %in% GBIF.spp, ]   ## Don't add records for other species
# length(unique(background.points$searchTaxon));dim(background.points)
# intersect(unique(background.points$searchTaxon), GBIF.spp)


# setdiff(names(SDM.SPAT.ALL), names(background.points))
# setdiff(names(background.points), names(SDM.SPAT.ALL))
# setdiff(names(SDM.SPAT.ALL), names(background.points))


#########################################################################################################################
## Now bind on the background points
# projection(SDM.SPAT.ALL);
# projection(background.points)
# setdiff(names(SDM.SPAT.ALL), names(SDM.SPAT.ALL))
# identical(length(names(SDM.SPAT.ALL)), length(names(background.points)))


# SDM.SPAT.OCC.BG = rbind(SDM.SPAT.ALL, background.points)
# unique(SDM.SPAT.OCC.BG$SOURCE)
# table(SDM.SPAT.OCC.BG$SOURCE)
# length(unique(SDM.SPAT.OCC.BG$searchTaxon))


#########################################################################################################################
## save data
# if(save_data == "TRUE") {
#   
#   ## save .rds file for the next session
#   SDM.SPAT.OCC.BG = 
#   saveRDS(SDM.SPAT.OCC.BG, paste0(DATA_path, 'SDM_SPAT_OCC_BG_',  save_run, '.rds'))
# 
# } else {
#   
#   message(' skip file saving, not many species analysed')   ##
#   
# }

## get rid of some memory
gc()





#########################################################################################################################
##################################################### TBC ############################################################### 
#########################################################################################################################