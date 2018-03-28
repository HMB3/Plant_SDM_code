#########################################################################################################################
############################################# ESTIMATE SPECIES NICHES ################################################### 
#########################################################################################################################


#########################################################################################################################
## This code combines the GBIF records for all species with the ALA data into a single table, extracts environmental 
## values and final adds contextual info for each record (taxonomic and horticultural) 


## It creates two tables:

## 1). A large table with one row for each species record
## 2). A smaller table with one row for each species, including contextual data and species attributes (niches, traits, etc.)
 
## These tables are subsequently used to estimate the current global realised niche/climatic tolerance using the best 
## available data, and susequently model the niches using the maxent algorithm.  


#########################################################################################################################
## To save time, load in previous data
source('./R/HIA_LIST_MATCHING.R')
load("./data/base/HIA_LIST/GBIF/GBIF_LAND_POINTS.RData")
load("./data/base/HIA_LIST/ALA/ALA_LAND_POINTS.RData")


## 196 species from the risk list are missing...presumably due to inusfficient data
length(setdiff(RISK.BINOMIAL.CLEAN$Plant_name, GBIF.LAND$searchTaxon))


#########################################################################################################################
## 1). MERGE GBIF DATA WITH CLEANED AVH/ALA RECORDS
#########################################################################################################################


#########################################################################################################################
## Merge on the ALA data. Consider that GBIF has data for both sources. We are topping up the native ranges with the AVH. 
## So there could be duplicates between both sources
length(unique(GBIF.LAND$searchTaxon)) 
names(GBIF.LAND)
names(ALA.LAND)
setdiff(names(GBIF.LAND), names(ALA.LAND))


## Test if the new experimental species are there
'Swainsona formosa'  %in% GBIF.LAND$searchTaxon


#########################################################################################################################
## Rename a few fields: Add rename of cultivated field, and make the values the same as mine "CULT" or "UNKNOWN"
## X = CULTIVATED
ALA.LAND     = dplyr::rename(ALA.LAND, 
                             searchTaxon                   = scientificname,
                             coordinateUncertaintyInMeters = uncertainty_m)


## Rename the field values: cultivated and uncultivated?
#ALA.LAND$CULT  = gsub("x", "CULTIVATED",   ALA.LAND$CULT)
#ALA.LAND$CULT  = gsub("x", "UNCULTIVATED", ALA.LAND$CULT)


## Restrict ALA data to just those species on the big list
library(plyr)
HIA.SPP.JOIN     = HIA.SPP
HIA.SPP.JOIN     = dplyr::rename(HIA.SPP.JOIN, searchTaxon = Binomial)


## Set NA to blank, then sort by no. of growers to get them to the top
HIA.SPP.JOIN[is.na(HIA.SPP.JOIN)] <- 0
HIA.SPP.JOIN = HIA.SPP.JOIN[with(HIA.SPP.JOIN, rev(order(Number.of.growers))), ]
head(HIA.SPP.JOIN[, c("searchTaxon", "Number.of.growers")])
View(HIA.SPP.JOIN)


## Get just those ALA species which are on the combined list of HIA and planted/growing
## ALA.LAND.HIA  = ALA.LAND[ALA.LAND$searchTaxon %in% HIA.SPP.JOIN$searchTaxon, ] 
ALA.LAND.HIA  = ALA.LAND[ALA.LAND$searchTaxon %in% unique(GBIF.LAND$searchTaxon), ] ## OR, %in% all.taxa. Old data includes Paul's extra species
str(unique(ALA.LAND$searchTaxon))       ##
str(unique(ALA.LAND.HIA$searchTaxon))   ## Reduced from 30k to 4K


#########################################################################################################################
## Bind the rows together?
GBIF.ALA.COMBO.LAND = bind_rows(GBIF.LAND, ALA.LAND.HIA)
names(GBIF.ALA.COMBO.LAND)
identical((dim(GBIF.LAND)[1]+dim(ALA.LAND)[1]),dim(GBIF.ALA.COMBO.LAND)[1])   ## only adding the ovelap...
head(GBIF.ALA.COMBO.LAND)


## Check the new data frame has just the species on the HIA list
str(unique(GBIF.ALA.COMBO.LAND$searchTaxon))                                  ## ok


## Test if a particular species is there
'Swainsona formosa'  %in% GBIF.ALA.COMBO.LAND$searchTaxon


## What species are unique to each dataset?
# length(setdiff(unique(GBIF.LAND$searchTaxon),  unique(ALA.LAND$searchTaxon)))
# length(setdiff(unique(ALA.LAND$searchTaxon),   unique(GBIF.LAND$searchTaxon)))
# length(intersect(unique(ALA.LAND$searchTaxon), unique(GBIF.LAND$searchTaxon)))


#########################################################################################################################
## Now crunch the big dataset down to just the species on the 25 growers or more list: the extras are just overkill......
GBIF.ALA.COMBO.HIA  = GBIF.ALA.COMBO.LAND[GBIF.ALA.COMBO.LAND$searchTaxon %in% unique(c(test.spp, HIA.SPP$Binomial)), ]
'Swainsona formosa'  %in% GBIF.ALA.COMBO.HIA $searchTaxon


## This unique ID column can be applied across the project  
GBIF.ALA.COMBO.HIA$OBS <- 1:nrow(GBIF.ALA.COMBO.HIA)
dim(GBIF.ALA.COMBO.HIA)[1];length(GBIF.ALA.COMBO.HIA$OBS)  
names(GBIF.ALA.COMBO.HIA)


## Now move this out to a separet file, so the cleaning code can be run separately
saveRDS(GBIF.ALA.COMBO.HIA, file = paste("./data/base/HIA_LIST/COMBO/GBIF_ALA_COMBO_PRE_CLEAN.RData"))


#########################################################################################################################
## DON'T REMOVE SPATIAL OUTLIERS FROM NICHE CALCULATION
#########################################################################################################################


#########################################################################################################################
## Create points: consider changing the coordinate system here to a global projected system?
COMBO.POINTS   = SpatialPointsDataFrame(coords      = GBIF.ALA.COMBO.HIA[c("lon", "lat")], 
                                        data        = GBIF.ALA.COMBO.HIA[c("lon", "lat")],
                                        #proj4string = CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")
                                        proj4string = CRS('+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'))


## Check
dim(COMBO.POINTS)


#########################################################################################################################
## CREATE NEW COLUMNS FOR ALA CULTIVATED/NOT 
#########################################################################################################################


## Need the ALA data to have the columns for cultivated/not: Need the original data frame. Also we don't need to worry
## that much about the duplicate recrods between GBIF and ALA, given we will just take one records per grid cell. Having
## duplicate records when estimating the niche won't make much difference either, assuming the are close together in space.

## Multiple records of the same specimen from different herbaria could be a problem though. Rachel's criteria of same month,
## year, lat/long could help here though. Ask Stu for this code.





#########################################################################################################################
## 2). PROJECT RASTERS AND EXTRACT ALL WORLDCLIM DATA FOR SPECIES RECORDS
#########################################################################################################################


#########################################################################################################################
## Ignore edaphic variables


# BIO1  = Annual Mean Temperature                                     ## 
# BIO2  = Mean Diurnal Range (Mean of monthly (max temp - min temp))  ##
# BIO3  = Isothermality (BIO2/BIO7) (* 100)
# BIO4  = Temperature Seasonality (standard deviation *100)           ##
# BIO5  = Max Temperature of Warmest Month                            ## 
# BIO6  = Min Temperature of Coldest Month                            ## 
# BIO7  = Temperature Annual Range (BIO5-BIO6)
# BIO8  = Mean Temperature of Wettest Quarter
# BIO9  = Mean Temperature of Driest Quarter
# BIO10 = Mean Temperature of Warmest Quarter
# BIO11 = Mean Temperature of Coldest Quarter
# BIO12 = Annual Precipitation                                        ##
# BIO13 = Precipitation of Wettest Month                              ##
# BIO14 = Precipitation of Driest Month                               ##
# BIO15 = Precipitation Seasonality (Coefficient of Variation)        ##
# BIO16 = Precipitation of Wettest Quarter
# BIO17 = Precipitation of Driest Quarter
# BIO18 = Precipitation of Warmest Quarter
# BIO19 = Precipitation of Coldest Quarter


#########################################################################################################################
## Create a stack of rasters to sample: get all the Worldclim variables just for good measure
# env.grids.current = stack(
#   file.path('./data/base/worldclim/world/0.5/bio/current',
#             sprintf('bio_%02d', 1:19)))


## Change the directories in the maxent mapping function to 1km...........................................................


## These are the variables we will use for the SDMs
sdm.select <- c("Annual_mean_temp",   "Temp_seasonality",   "Max_temp_warm_month", "Min_temp_cold_month",
                "Annual_precip",      "Precip_seasonality", "Precip_wet_month",    "Precip_dry_month")


## Create names for all 19 predictors
pred_names <- c(
  'Annual_mean_temp',   ## To select a column with the cursor, hold ctrl+alt and use up or down arrow
  'Mean_diurnal_range',
  'Isothermality',
  'Temp_seasonality', 
  'Max_temp_warm_month', 
  'Min_temp_cold_month', 
  'Temp_annual_range',
  'Mean_temp_wet_qu', 
  'Mean_temp_dry_qu', 
  'Mean_temp_warm_qu',
  'Mean_temp_cold_qu',
  'Annual_precip', 
  'Precip_wet_month', 
  'Precip_dry_month', 
  'Precip_seasonality', 
  'Precip_wet_qu', 
  'Precip_dry_qu', 
  'Precip_warm_qu', 
  'Precip_col_qu') 


## Create an index for each name 
i  <- match(pred_names, pred_names)


## Create file paths using existing files, these will be projected
ff_current <- file.path('./data/base/worldclim/world/0.5/bio/current', sprintf('bio_%02d.tif', i))
ff_2030    <- file.path('./data/base/worldclim/aus/0.5/bio/2030',      sprintf('bio_%02d.tif', i))
ff_2050    <- file.path('./data/base/worldclim/aus/0.5/bio/2050',      sprintf('bio_%02d.tif', i))
ff_2070    <- file.path('./data/base/worldclim/aus/0.5/bio/2070',      sprintf('bio_%02d.tif', i))


## Create directories for the projected files :: 1km 
dir.create(dirname(sub('0.5', '1km', ff_current)[1]), recursive = TRUE)
dir.create(dirname(sub('0.5', '1km', ff_2030)[1]),    recursive = TRUE)
dir.create(dirname(sub('0.5', '1km', ff_2050)[1]),    recursive = TRUE)
dir.create(dirname(sub('0.5', '1km', ff_2070)[1]),    recursive = TRUE)


## Run a loop to warp the worldclim variables into the World Mollweide projected system, saving in the 1km folder
lapply(ff_current, function(f) {
  
  message(f)
  gdalwarp(f, sub('0.5', '1km', f), tr = c(1000, 1000),
           t_srs = '+init=esri:54009', r = 'bilinear', 
           multi = TRUE)
  
})


## Run a loop to warp the worldclim variables into the World Mollweide projected system, saving in the 1km folder
lapply(ff_current, function(f) {
  
  message(f)
  gdalwarp(f, sub('0.5', '1km', f), tr = c(1000, 1000),
           t_srs = '+init=esri:54009', r = 'bilinear', 
           multi = TRUE)
  
})


## Run a loop to warp the worldclim variables into the World Mollweide projected system, saving in the 1km folder
lapply(ff_2030, function(f) {
  
  message(f)
  gdalwarp(f, sub('0.5', '1km', f), tr = c(1000, 1000),
           t_srs = '+init=esri:54009', r = 'bilinear', 
           multi = TRUE)
  
})


## Run a loop to warp the worldclim variables into the World Mollweide projected system, saving in the 1km folder
lapply(ff_2050, function(f) {
  
  message(f)
  gdalwarp(f, sub('0.5', '1km', f), tr = c(1000, 1000),
           t_srs = '+init=esri:54009', r = 'bilinear', 
           multi = TRUE)
  
})


## Run a loop to warp the worldclim variables into the World Mollweide projected system, saving in the 1km folder
lapply(ff_2070, function(f) {
  
  message(f)
  gdalwarp(f, sub('0.5', '1km', f), tr = c(1000, 1000),
           t_srs = '+init=esri:54009', r = 'bilinear', 
           multi = TRUE)
  
})


## Create a raster stack of the projected grids
env.grids.current = stack(sub('0.5', '1km', ff))
names(env.grids.current) <- pred_names[i]


## Here, we should project all datasets into the final system used by maxent (World Mollweide https://epsg.io/54009)
## Using step 6...

# ## Can we resample the rasters to 10km?
# BIOW_10km = raster("./data/base/worldclim/world/0.5/bio/current/bio_01_10km.tif")
# BIOA_10km = raster("./data/base/worldclim/aus/0.5/bio/current/bio_aus_01_10km.tif")
# 
# 
# ## Use bilinear?
# test          = resample(env.grids.current, BIOW_10km, method = 'bilinear')
# raster_path   = "./data/base/worldclim/aus/0.5/bio/current/"
# current_10km  <- sprintf('%scurrent_env_10km.tif', raster_path)
# writeRaster(test, current_10km , overwrite = TRUE)


#########################################################################################################################
## Then use the extract function for all the rasters, and finaly bind on the COMBO data to the left of the raster values
## Can we use a cluster to speed this up?

## Best option to speed this up is to use only the unique cells
# COMBO.XY <- cellFromXY(world.temp, COMBO.POINTS) %>% 
#   
#   ## get the unique raster cells
#   unique %>% 
#   
#   ## Get coordinates of the center of raster cells for a row, column, or cell number of WORLDCLIM raster
#   xyFromCell(world.temp, .)


#########################################################################################################################
## Is there a way to speed this up?
## ## Test if the new experimental species are there : 'Swainsona formosa'  %in% GBIF.ALA.COMBO.HIA$searchTaxon
COMBO.RASTER <- extract(env.grids.current, COMBO.POINTS) %>% 
  cbind(GBIF.ALA.COMBO.HIA, .)


#'Swainsona formosa'  %in% COMBO.RASTER$searchTaxon


## Multiple rename using dplyr
# COMBO.RASTER = dplyr::rename(COMBO.RASTER,
#                              Annual_mean_temp     = bio_01, 
#                              Mean_diurnal_range   = bio_02,
#                              Isothermality        = bio_03,
#                              Temp_seasonality     = bio_04, 
#                              Max_temp_warm_month  = bio_05, 
#                              Min_temp_cold_month  = bio_06, 
#                              Temp_annual_range    = bio_07,
#                              Mean_temp_wet_qu     = bio_08, 
#                              Mean_temp_dry_qu     = bio_09, 
#                              Mean_temp_warm_qu    = bio_10,
#                              Mean_temp_cold_qu    = bio_11,
#                              
#                              Annual_precip        = bio_12, 
#                              Precip_wet_month     = bio_13, 
#                              Precip_dry_month     = bio_14, 
#                              Precip_seasonality   = bio_15, 
#                              Precip_wet_qu        = bio_16, 
#                              Precip_dry_qu        = bio_17, 
#                              Precip_warm_qu       = bio_18, 
#                              Precip_col_qu        = bio_19) 


## Save/load
saveRDS(COMBO.RASTER, file = paste("./data/base/HIA_LIST/GBIF/COMBO_GBIF_ALA_RASTER.rds"))
#load("./data/base/HIA_LIST/GBIF/COMBO_GBIF_ALA_RASTER.rds")
#COMBO.RASTER  = COMBO.RASTER[COMBO.RASTER$searchTaxon %in% HIA.SPP$Binomial, ]


## Check
dim(COMBO.RASTER)
names(COMBO.RASTER)





#########################################################################################################################
## 3). INTERSECT SPECIES RECORDS WITH LOCAL GOV AREAS AND SIGNIFICANT URBAN AREAS
#########################################################################################################################


## We want to know the count of species that occur in 'n' LGAs, across a range of climates. Read in LGA and SUA
COMBO.RASTER.SP   = SpatialPointsDataFrame(coords      = COMBO.RASTER[c("lon", "lat")], 
                                           data        = COMBO.RASTER,
                                           #proj4string = CRS("+init=epsg:4326"),
                                           proj4string = CRS('+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'))

SUA      = readOGR("F:/green_cities_sdm/data/base/CONTEXTUAL/SUA_2011_AUST.shp", layer = "SUA_2011_AUST")
LGA      = readOGR("F:/green_cities_sdm/data/base/CONTEXTUAL/LGA_2016_AUST.shp", layer = "LGA_2016_AUST")

names(SUA)
names(LGA)


## Project using a projected rather than geographic coordinate system
#CRS.new  <- CRS("+init=epsg:4326") # EPSG:3577
CRS.new  <- CRS('+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')
LGA.WGS  = spTransform(LGA, CRS.new)
SUA.WGS  = spTransform(SUA, CRS.new)


## Double check they are the same
projection(COMBO.RASTER.SP);projection(LGA.WGS);projection(SUA.WGS)


## Remove the LGA columns we don't need
LGA.WGS = LGA.WGS[, c("LGA_CODE16", "LGA_NAME16")] 


## Now create a shapefile which is just the SUA: in or out...
# IN.SUA <- SUA[ !(SUA$SUA_NAME11 %in% c("Not in any Significant Urban Area (NSW)", 
#                                        "Not in any Significant Urban Area (Vic.)",
#                                        "Not in any Significant Urban Area (Qld)",
#                                        "Not in any Significant Urban Area (SA)",
#                                        "Not in any Significant Urban Area (WA)",
#                                        "Not in any Significant Urban Area (Tas.)",
#                                        "Not in any Significant Urban Area (NT)",
#                                        "Not in any Significant Urban Area (OT)",
#                                        "Not in any Significant Urban Area (ACT)")), ]
# 
# plot(IN.SUA)
# writeOGR(obj = IN.SUA, dsn = "./data/base/CONTEXTUAL", layer = "IN.SUA", driver = "ESRI Shapefile")


## Then, we want to create a layer which is just in the urban area, or not. This would need to combine the above fields into one
IN.SUA   = readOGR("F:/green_cities_sdm/data/base/CONTEXTUAL/INSIDE_AUS_SUA.shp", layer = "INSIDE_AUS_SUA")
IN.SUA   = IN.SUA[,-(1)];
IN.SUA   = IN.SUA[,-(2)]
names(IN.SUA)
IN.SUA  = spTransform(IN.SUA, CRS.new)
plot(IN.SUA)


#########################################################################################################################
## Run join between species records and LGAs/SUAs
SUA.JOIN      = over(COMBO.RASTER.SP, IN.SUA)  
LGA.JOIN      = over(COMBO.RASTER.SP, LGA.WGS)   
COMBO.SUA.LGA = cbind.data.frame(COMBO.RASTER.SP, SUA.JOIN, LGA.JOIN) ## Include below, just get the LGA columns?  

head(SUA.JOIN)
head(LGA.JOIN)


#########################################################################################################################
## AGGREGATE THE NUMBER OF LGAs EACH SPECIES IS FOUND IN 
LGA.AGG   = tapply(COMBO.SUA.LGA$LGA_NAME16, COMBO.SUA.LGA$searchTaxon, function(x) length(unique(x))) ## group LGA by species name
LGA.AGG   = as.data.frame(LGA.AGG)
head(LGA.AGG)


## Save
saveRDS(COMBO.SUA.LGA, file = paste("./data/base/HIA_LIST/GBIF/COMBO_SUA_LGA.rds"))
saveRDS(LGA.AGG,   file = paste("./data/base/HIA_LIST/GBIF/LGA_AGG.rds"))


## 
str(COMBO.SUA.LGA)
head(COMBO.SUA.LGA)
names(COMBO.SUA.LGA)
COMBO.SUA.LGA = subset(COMBO.SUA.LGA, select = -c(lon.1, lat.1))




#########################################################################################################################
## 4). CREATE NICHES FOR SELECTED TAXA
#########################################################################################################################


#########################################################################################################################
## Now summarise the niches. But figure out a cleaner way of doing this
# env.variables = c("Annual_mean_temp",     
#                   "Mean_diurnal_range",   
#                   "Isothermality",        
#                   "Temp_seasonality",     
#                   "Max_temp_warm_month",  
#                   "Min_temp_cold_month",  
#                   "Temp_annual_range",    
#                   "Mean_temp_wet_qu",     
#                   "Mean_temp_dry_qu",     
#                   "Mean_temp_warm_qu",    
#                   "Mean_temp_cold_qu",
#                   
#                   "Annual_precip",        
#                   "Precip_wet_month",     
#                   "Precip_dry_month",     
#                   "Precip_seasonality",   
#                   "Precip_wet_qu",        
#                   "Precip_dry_qu",        
#                   "Precip_warm_qu",       
#                   "Precip_col_qu")


#########################################################################################################################
## Change the raster values here: See http://worldclim.org/formats1 for description of the interger conversion. 
## All temperature variables wer multiplied by 10, so divide by 10 to reverse it.
COMBO.RASTER.CONVERT = as.data.table(COMBO.SUA.LGA)                           ## Check this works, also inefficient
COMBO.RASTER.CONVERT[, (pred_names [c(1:11)]) := lapply(.SD, function(x) 
  x / 10 ), .SDcols = pred_names [c(1:11)]]
COMBO.RASTER.CONVERT = as.data.frame(COMBO.RASTER.CONVERT)                   ## Find another method without using data.table


## Check Looks ok?
summary(COMBO.RASTER.CONVERT$Annual_mean_temp)  ## -23 looks too low/. Check where these are ok
summary(COMBO.RASTER$Annual_mean_temp)

summary(COMBO.RASTER.CONVERT$Isothermality)
summary(COMBO.RASTER$Isothermality)


## Plot a few points to see
plot(LAND, col = 'grey', bg = 'sky blue')
points(COMBO.RASTER.CONVERT[ which(COMBO.RASTER.CONVERT$Annual_mean_temp < -5), ][, c("lon", "lat")], 
       pch = ".", col = "red", cex = 3, asp = 1)



#########################################################################################################################
## Create niche summaries for each environmental condition like this...
## Here's what the function will produce :
library(plyr)
head(niche_estimate (DF = COMBO.RASTER.CONVERT, colname = "Annual_mean_temp"))  ## including the q05 and q95


## So lets use lapply on the "Search Taxon". Note additonal flags are needed, and the taxonomic lists need to be managed better...
## test = run_function_concatenate(list, DF, "DF, colname = x") 
COMBO.NICHE <- pred_names[c(1:length(pred_names))] %>% 
  
  ## Pipe the list into lapply
  lapply(function(x) {
    
    ## Now use the niche width function on each colname (so 8 environmental variables)
    ## Also, need to figure out how to make the aggregating column generic (species, genus, etc.)
    ## currently it only works hard-wired
    niche_estimate (DF = COMBO.RASTER.CONVERT, colname = x)
    
    ## would be good to remove the duplicate columns here
    
  }) %>% 
  
  ## finally, create one dataframe for all niches
  as.data.frame


## Remove duplicate Taxon columns and check the output :: would be great to skip these columns when running the function
names(COMBO.NICHE)
COMBO.NICHE = subset(COMBO.NICHE, select = -c(searchTaxon.1,  searchTaxon.2,  searchTaxon.3,  searchTaxon.4,
                                              searchTaxon.5,  searchTaxon.6,  searchTaxon.7,  searchTaxon.7,
                                              searchTaxon.8,  searchTaxon.9,  searchTaxon.10, searchTaxon.11,
                                              searchTaxon.12, searchTaxon.13, searchTaxon.14, searchTaxon.15,
                                              searchTaxon.16, searchTaxon.17, searchTaxon.18))


## Add counts for each species, and record the total number of taxa processed
COMBO.count = as.data.frame(table(COMBO.RASTER.CONVERT$searchTaxon))$Freq
Total.taxa.processed = dim(COMBO.NICHE)[1]
COMBO.NICHE  = cbind(COMBO.count, COMBO.NICHE)
names(COMBO.NICHE)





#########################################################################################################################
## 5). CALCULATE AREA OF OCCUPANCY RANGES 
#########################################################################################################################


## These numbers don't look that accurate: Try convert into a sp data frame and projecting into a projected system?
## Create a species list to estimate the ranges for
spp.geo = as.character(unique(COMBO.RASTER$searchTaxon)) 
data    = COMBO.RASTER


#########################################################################################################################
## AREA OF OCCUPANCY (AOO)
## For every species in the list: calculate the AOO
GBIF.AOO <- spp.geo[c(1:length(spp.geo))] %>% 
  
  ## Pipe the list into lapply
  lapply(function(x) {
    
    ## Subset the the data frame 
    DF      = subset(data, searchTaxon == x)[, c("lon", "lat")]
    
    ## Calculate area of occupancy according the the "red" package
    aoo (DF)
    
    ## Warning messages: Ask John if this is a problem
    ## In rgdal::project(longlat, paste("+proj=utm +zone=", zone,  ... :
    ## 3644 projected point(s) not finite
    
  }) %>% 
  
  ## Finally, create one dataframe for all niches
  as.data.frame


#########################################################################################################################
## EXTENT OF OCCURRENCE
## For every species in the list: calculate the EOO
# GBIF.EOO <- spp.geo[c(1:length(spp.geo))] %>% 
#   
#   ## Pipe the list into lapply
#   lapply(function(x) {
#     
#     ## Subset the the data frame 
#     DF      = subset(data, searchTaxon == x)[, c("searchTaxon", "lon", "lat")]
#     DF.GEO  = dplyr::rename(DF, 
#                             identifier = searchTaxon,
#                             XCOOR      = lon,
#                             YCOOR      = lat)
#     
#     ## Calculate area of occupancy according the the "red" package
#     CalcRange (DF.GEO)
#     
#     ## Warning messages: Ask John if this is a problem
#     ## In rgdal::project(longlat, paste("+proj=utm +zone=", zone,  ... :
#     ## 3644 projected point(s) not finite
#     
#   }) %>% 
#   
#   ## Finally, create one dataframe for all niches
#   as.data.frame


#########################################################################################################################
## Clean it up :: the order of species should be preserved
GBIF.AOO = gather(GBIF.AOO)
str(GBIF.AOO)
str(unique(GBIF.ALA.COMBO.HIA$searchTaxon))                     ## same number of species...


## Now join on the GEOGRAPHIC RANGE
COMBO.NICHE$AREA_OCCUPANCY = GBIF.AOO$value                     ## vectors same length so don't need to match


## AOO is calculated as the area of all known or predicted cells for the species. The resolution will be 2x2km as 
## required by IUCN. A single value in km2.
## Add the counts of LGAs for each species in here:
COMBO.LGA = cbind.data.frame(COMBO.NICHE, LGA.AGG)              ## The tapply needs to go where the niche summaries are
names(COMBO.LGA)





#########################################################################################################################
## 6). JOIN ON CONTEXTUAL DATA
#########################################################################################################################


#########################################################################################################################
## Now join the horticultural contextual data onto one or both tables ()
COMBO.RASTER.CONTEXT = join(COMBO.RASTER.CONVERT, HIA.SPP.JOIN, 
                            by = "searchTaxon", type = "left", match = "all")


## Now join hort context to all the niche
COMBO.NICHE.CONTEXT = join(COMBO.LGA, HIA.SPP.JOIN, 
                           by = "searchTaxon", type = "left", match = "all")


## Check taxa again
'Kennedia beckxiana' %in% COMBO.NICHE.CONTEXT$searchTaxon  


## Changing the order is not needed for the raster data
#COMBO.RASTER.CONTEXT = COMBO.RASTER.CONTEXT[, c(1:5,  45:57, 6:43)]                                      
#COMBO.NICHE.CONTEXT  = COMBO.NICHE.CONTEXT[,  c(176, 2, 1, 174:175, 177:189, 3:173)]                ## REDO


## Set NA to blank, then sort by no. of growers. The extra species are ones I dowloaded accidentally
COMBO.NICHE.CONTEXT$Number.of.growers[is.na(COMBO.NICHE.CONTEXT$Number.of.growers)] <- 0
COMBO.NICHE.CONTEXT = COMBO.NICHE.CONTEXT[with(COMBO.NICHE.CONTEXT, rev(order(Number.of.growers))), ]


## View the data
names(COMBO.RASTER.CONTEXT)
names(COMBO.NICHE.CONTEXT)
dim(COMBO.RASTER.CONTEXT)
dim(COMBO.NICHE.CONTEXT)








#########################################################################################################################
## quickly check how many species match from the original 605. Only 553 are currently there.


## Which species from those with > 25 growers are missing?
missing.25     = setdiff(unique(HIA.SPP.JOIN[ which(HIA.SPP.JOIN$Number.of.growers >= 25), ][["searchTaxon"]]),
                         unique(COMBO.NICHE.CONTEXT[ which(COMBO.NICHE.CONTEXT$Number.of.growers >= 25), ][["searchTaxon"]]))


## Same as for
missing.all    = setdiff(unique(HIA.SPP.JOIN[["searchTaxon"]]),
                         unique(COMBO.NICHE.CONTEXT[["searchTaxon"]]))

## Which species from the top 200 are missing?
missing.200    = setdiff(unique(spp.200$Binomial),
                         unique(COMBO.NICHE.CONTEXT[ which(COMBO.NICHE.CONTEXT$Number.of.growers >= 25), ][["searchTaxon"]]))


## Which species from the experimental list are missing?
missing.renee  = setdiff(unique(test.spp),
                         unique(COMBO.NICHE.CONTEXT[["searchTaxon"]]))


## Overall list
missing.taxa   = sort(unique(c(missing.25, missing.200, missing.renee)))
missing.taxa   = as.data.frame(missing.taxa)
missing.taxa


## The missing species are due to too few records, too many, or taxonomy problems. EG some of the species are varieties, so they 
## only match to the genus. So 605 - 30 = 575.




#########################################################################################################################
## Combine the existing data with the new species :: use this approach to just process a subset 
## COMBO.RASTER.CONTEXT.UPDATE = COMBO.RASTER.CONTEXT
## saveRDS(COMBO.RASTER.CONTEXT.UPDATE, file = paste("./data/base/HIA_LIST/COMBO/COMBO_RASTER_CONTEXT_UPDATE.rds", sep = ""))
## load("./data/base/HIA_LIST/COMBO/COMBO_RASTER_CONTEXT.rds")
## load("./data/base/HIA_LIST/COMBO/COMBO_RASTER_CONTEXT_UPDATE.rds")
# names(COMBO.RASTER.CONTEXT);names(COMBO.RASTER.CONTEXT.UPDATE)
# COMBO.RASTER.CONTEXT.UPDTATE = bind_rows(COMBO.RASTER.CONTEXT, COMBO.RASTER.CONTEXT.UPDATE)
       

## Save the summary datasets
saveRDS(missing.taxa,         file = paste("./data/base/HIA_LIST/COMBO/MISSING_TAXA.rds",                   sep = ""))
saveRDS(COMBO.RASTER.CONTEXT, file = paste("./data/base/HIA_LIST/COMBO/COMBO_RASTER_CONTEXT_2703_2018.rds", sep = ""))
saveRDS(COMBO.NICHE.CONTEXT,  file = paste("./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT_2703_2018.rds",  sep = ""))
write.csv(COMBO.NICHE.CONTEXT, "./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT_2703_2018.csv",       row.names = FALSE)


## Now save .RData file for the next session...
save.image("STEP_4_NICHES.RData")





#########################################################################################################################
## OUTSTANDING NICHE TASKS:
#########################################################################################################################


## Check geographic range: doesn't look right for some species        - Keep a spreadsheet of all species...

## Estimate native/naturalised ranges as a separate colum             - APC data + Ubran polygons for AUS, USA, EU: only a subset

## GBIF taxonomic errors                                              - Use taxonstand

## Keep cultivated records as a separate column/file                  - Get cultivated column from ALA data...

## Duplicates between GBIF and ALA                                    - See email from CSIRO - only a problem for niches...



#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################