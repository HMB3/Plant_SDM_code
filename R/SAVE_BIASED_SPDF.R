#########################################################################################################################
################################################# THIN RECORDS ##########################################################
#########################################################################################################################


## This code spatially thins species occurrence records to help address problems associated with spatial sampling biases. 
## Ideally, thinning removes the fewest records necessary to substantially reduce the effects of sampling bias, while 
## simultaneously retaining the greatest amount of useful information. 

## Use the SDM dataset - one record per 1km cell - rather than the original dataset......................................

## Use the stratified function: https://www.rdocumentation.org/packages/fifer/versions/1.0/topics/stratified
## spThin didn't really work https://cran.r-project.org/web/packages/spThin/vignettes/spThin_vignette.html


## Create lists
source('./R/HIA_LIST_UPDATE.R')


#########################################################################################################################
## 1). CREATE SPATIAL DATAFRAME FOR BIASED SPECIES
#########################################################################################################################


#########################################################################################################################
## Read in background data
background = readRDS("./data/base/HIA_LIST/COMBO/SDM_DATA_CLEAN_052018.rds")
background = background [!background$searchTaxon %in% GBIF.spp, ]               ## Don't add records for other species
length(unique(background$searchTaxon));dim(background)


#########################################################################################################################
## Load GBIF data and rain shapefile
BIAS.DATA.ALL        = readRDS("./data/base/HIA_LIST/COMBO/SDM_DATA_INV_HIA.rds") 
SPP.BIAS             = intersect(SPP.BIAS, GBIF.spp)    ## just re-run the models for species on the list
SPP_BIAS             = gsub(" ", "_", SPP.BIAS)
length(unique(BIAS.DATA.ALL$searchTaxon))
 
 
## Project the SDM data into WGS
BIAS.DATA.ALL <- spTransform(BIAS.DATA.ALL, CRS("+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
BG.DATA.ALL   <- spTransform(background, CRS("+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
projection(BIAS.DATA.ALL);projection(background)


## Get the coordinates
BIAS.COORDS = coordinates(BIAS.DATA.ALL)
BG.COORDS   = coordinates(background)
class(BIAS.COORDS)
tail(BIAS.COORDS)
tail(BIAS.DATA.ALL)[, c(1:3)]   ## indices match


## Bind the coordinates to the SDM table
BIAS.COORDS  = as.data.frame(BIAS.COORDS)
BIAS.DATA.DF = cbind(as.data.frame(BIAS.DATA.ALL), BIAS.COORDS)
BIAS.DATA.DF = BIAS.DATA.DF[, c(1:22)]

BG.COORDS    = as.data.frame(BG.COORDS)
BG.DATA.DF   = cbind(as.data.frame(BG.DATA.ALL), BG.COORDS)
BG.DATA.DF   = BG.DATA.DF[, c(1:22)]


names(BG.DATA.DF)
names(BIAS.DATA.DF)
dim(BIAS.DATA.DF)


#########################################################################################################################
## Intersect BOM with rainfall data
BIAS.DATA.SP   = SpatialPointsDataFrame(coords      = BIAS.DATA.DF[c("lon", "lat")],
                                        data        = BIAS.DATA.DF,
                                        proj4string = CRS.WGS.84)


BG.DATA.SP   = SpatialPointsDataFrame(coords      = BG.DATA.DF[c("lon", "lat")],
                                      data        = BG.DATA.DF,
                                      proj4string = CRS.WGS.84)


## All predictors
sdm.predictors <- c("Annual_mean_temp",    "Mean_diurnal_range",  "Isothermality",      "Temp_seasonality",  
                    "Max_temp_warm_month", "Min_temp_cold_month", "Temp_annual_range",  
                    "Mean_temp_wet_qu",    "Mean_temp_dry_qu",    
                    "Mean_temp_warm_qu",   "Mean_temp_cold_qu",   
                    
                    "Annual_precip",       "Precip_wet_month",   "Precip_dry_month",    "Precip_seasonality",  
                    "Precip_wet_qu",       "Precip_dry_qu",      
                    "Precip_warm_qu",      "Precip_col_qu")


## A-priori worldclim predictors
sdm.select     <- c("Annual_mean_temp", "Temp_seasonality",    "Max_temp_warm_month", "Min_temp_cold_month",
                    "Annual_precip",    "Precip_seasonality",  
                    "Precip_wet_month", "Precip_dry_month")  
 

#########################################################################################################################
## Check data subsets - subset the big dataframe to just the biased species
BIAS.DF = BIAS.DATA.SP#[BIAS.DATA.SP$searchTaxon %in% SPP.BIAS, ]
length(unique(BIAS.DF$searchTaxon))
dropList <- setdiff(sdm.predictors, sdm.select)
BIAS.DF  <- BIAS.DF[, !names(BIAS.DF) %in% dropList]
names(BIAS.DF)


## Now rename columns for ArcMap
BIAS.DF     = dplyr::rename(as.data.frame(BIAS.DF),
                            BIO1      = Annual_mean_temp,
                            BIO4      = Temp_seasonality,
                            BIO5      = Max_temp_warm_month,
                            BIO6      = Min_temp_cold_month,
                            BIO12     = Annual_precip,
                            BIO16     = Precip_wet_month,
                            BIO14     = Precip_dry_month,
                            BIO15     = Precip_seasonality,
                            LAT       = lat,
                            LON       = lon)


## Then create SPDF
BIAS.DF   = SpatialPointsDataFrame(coords      = BIAS.DF[c("LON", "LAT")],
                                   data        = BIAS.DF,
                                   proj4string = CRS.WGS.84)


## Project the SDM data into WGS
BIAS.DF <- spTransform(BIAS.DF, CRS("+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
projection(BIAS.DF)


## Save the shapefile, to be subsampled in ArcMap
names(BIAS.DF);head(BIAS.DF)
writeOGR(obj = BIAS.DF,    dsn = "./data/base/HIA_LIST/COMBO", layer = "SPP_ALL_DF_TREE_INV", driver = "ESRI Shapefile")
writeOGR(obj = BG.DATA.SP, dsn = "./data/base/HIA_LIST/COMBO", layer = "BG_POINTS", driver = "ESRI Shapefile")




#########################################################################################################################
## 2). PLOT BIASED SPECIES AND SAVE TABLE FOR MAXENT
#########################################################################################################################


#########################################################################################################################
## Loop over the species list and plot the occurrence data for each to check the data bias
for (i in 1:length(SPP.BIAS)) {
  
  ## Create points for each species
  spp.points <- BIAS.DATA.ALL[BIAS.DATA.ALL$searchTaxon == SPP.BIAS[i], ] %>%
    spTransform(CRS('+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'))
  
  aus.mol <- aus %>%
    spTransform(CRS('+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'))
  
  ## Print to file
  save_name = gsub(' ', '_', SPP.BIAS[i])
  save_dir  = "data/base/HIA_LIST/COMBO/BIASED"
  png(sprintf('%s/%s_%s.png', save_dir,
              save_name, "Australian_points"),
      3236, 2000, units = 'px', res = 300)

  ## set margins
  par(mar   = c(3, 3, 5, 3),  ## b, l, t, r
      #mgp   = c(9.8, 2.5, 0),
      oma   = c(1.5, 1.5, 1.5, 1.5))

  ## Plot just the Australian points
  plot(aus.mol, main = SPP.BIAS[i])
  points(spp.points, col = "red", cex = .15, pch = 19)
  
  # Finish the device
  dev.off()
  
}


#########################################################################################################################
## How many records do these species have?
names(TOT.GROW);names(TREE.EVERGREEN)
COUNTS = COMBO.NICHE.CONTEXT[c("searchTaxon", "COMBO.count", "AUS_RECORDS")]
TREES.BIAS  = join(TREE.EVERGREEN, TOT.GROW)
TREES.BIAS  = join(TREES.BIAS, COUNTS)
TREES.BIAS  = TREES.BIAS[c("searchTaxon", "Origin", "Plantings", "COMBO.count",   "AUS_RECORDS", "Number.of.States")]
#TREES.BIAS  = TREES.BIAS [TREES.BIAS$searchTaxon %in% SPP.BIAS, ]

dim(TREES.BIAS);head(TREES.BIAS, 40)





#########################################################################################################################
## Read the bias file back in :
## H:\green_cities_sdm\data\base\HIA_LIST\COMBO\RAREFY_SPP\SPAT_OUT\ALL_TREE_INV_2KM_spatially_rarified_locs.shp
BIAS.RAREFY = readOGR("H:/green_cities_sdm/data/base/HIA_LIST/COMBO/RAREFY_SPP/SPAT_OUT/ALL_TREE_INV_2KM_spatially_rarified_locs.shp",
                      layer = "ALL_TREE_INV_2KM_spatially_rarified_locs")
# saveRDS(BIAS.RAREFY,    'data/base/HIA_LIST/COMBO/RAREFY_SPP/SPAT_OUT/BIAS_RAREFY.rds')
# BIAS.RAREFY       = readRDS('./data/base/HIA_LIST/COMBO/RAREFY_SPP/SPAT_OUT/BIAS_RAREFY.rds')
class(BIAS.RAREFY)
names(BIAS.RAREFY)


## Rename
BIAS.RAREFY     = dplyr::rename(as.data.frame(BIAS.RAREFY),
                                searchTaxon         = srchTxn,
                                Annual_mean_temp    = BIO1,
                                Temp_seasonality    = BIO4,
                                Max_temp_warm_month = BIO5,
                                Min_temp_cold_month = BIO6,
                                Annual_precip       = BIO12,
                                Precip_wet_month    = BIO16,
                                Precip_dry_month    = BIO14,
                                Precip_seasonality  = BIO15,
                                lon = LON,
                                lat = LAT)

dropList <- c("lon_1", "lat_1", "RASTERVALU", "coords.x1", "coords.x2" )
BIAS.RAREFY  <- BIAS.RAREFY[, !names(BIAS.RAREFY) %in% dropList]
names(BIAS.RAREFY)


## Re-project
BIAS.RAREFY      = SpatialPointsDataFrame(coords      = BIAS.RAREFY[c("lon", "lat")],
                                          data        = BIAS.RAREFY,
                                          proj4string = CRS.WGS.84)
BIAS.RAREFY       = spTransform(BIAS.RAREFY , CRS('+init=ESRI:54009'))
projection(BIAS.RAREFY)
length(unique(BIAS.RAREFY$searchTaxon))
names(BIAS.RAREFY)



#########################################################################################################################
## Loop over the species list and plot the occurrence data for each to check the data bias
for (i in 1:length(SPP.BIAS)) {
  
  ## Create points for each species
  spp.points <- BIAS.RAREFY[BIAS.RAREFY$searchTaxon == SPP.BIAS[i], ] %>%
    spTransform(CRS('+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'))
  
  aus.mol <- aus %>%
    spTransform(CRS('+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'))
  
  ## Print to file
  save_name = gsub(' ', '_', SPP.BIAS[i])
  save_dir  = "data/base/HIA_LIST/COMBO/RAREFY"
  png(sprintf('%s/%s_%s.png', save_dir,
              save_name, "Australian_points"),
      3236, 2000, units = 'px', res = 300)
  
  ## set margins
  par(mar   = c(3, 3, 5, 3),  ## b, l, t, r
      #mgp   = c(9.8, 2.5, 0),
      oma   = c(1.5, 1.5, 1.5, 1.5))
  
  ## Plot just the Australian points
  plot(aus.mol, main = SPP.BIAS[i])
  points(spp.points, col = "red", cex = .15, pch = 19)
  
  # Finish the device
  dev.off()
  
}


## save as 
saveRDS(BIAS.RAREFY, 'data/base/HIA_LIST/COMBO/RAREFY_SPP/SPAT_OUT/SDM_TREES_RAREFY.rds')





#########################################################################################################################
#####################################################  TBC ##############################################################
#########################################################################################################################