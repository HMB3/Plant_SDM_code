#########################################################################################################################
############################################# ESTIMATE SPECIES NICHES ################################################### 
#########################################################################################################################


#########################################################################################################################
## This code combines the GBIF records for all species with the ALA data into a single table, and adds contextual info  


## It creates two tables:

## 1). A table with one row for each species record
## 2). A table with one row for each species, including contextual data and species attributes (niches, traits, etc.)

## These tables are used to estimate the current global realised niche/climatic tolerance using the best available data 


#########################################################################################################################
## To save time, load in previous data
source('./R/HIA_LIST_MATCHING.R')
load("./data/base/HIA_LIST/GBIF/GBIF_LAND_POINTS.RData")
load("./data/base/HIA_LIST/ALA/ALA_LAND_POINTS.RData")


#########################################################################################################################
## 1). MERGE GBIF DATA WITH CLEANED AVH/ALA RECORDS
#########################################################################################################################


#########################################################################################################################
## Merge on the ALA data. Consider that GBIF has data for both sources. We are topping up the native ranges with the AVH. 
## So there could be duplicates between both sources
length(unique(GBIF.LAND$searchTaxon)) 


## Rename a few fields:
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
ALA.LAND.HIA  = ALA.LAND[ALA.LAND$searchTaxon %in% kop.spp, ] ## OR, %in% all.taxa. Old data includes Paul's extra species
str(unique(ALA.LAND$searchTaxon))       ##
str(unique(ALA.LAND.HIA$searchTaxon))   ## 


#########################################################################################################################
## Bind the rows together?
GBIF.ALA.COMBO.LAND = bind_rows(GBIF.LAND, ALA.LAND.HIA)


## Check the new data frame has just the species on the HIA list
str(unique(GBIF.ALA.COMBO.LAND$searchTaxon))                                  ## ok


#########################################################################################################################
## Now crunch the big dataset down to just the species on the 25 growers or more list: the extras are just overkill......
GBIF.ALA.COMBO.HIA  = GBIF.ALA.COMBO.LAND[GBIF.ALA.COMBO.LAND$searchTaxon %in% kop.spp, ]
dim(GBIF.ALA.COMBO.HIA)


# ## Now move this out to a separet file, so the cleaning code can be run separately
# save(GBIF.ALA.COMBO.HIA, file = paste("./data/base/HIA_LIST/COMBO/GBIF_ALA_COMBO_PRE_CLEAN.RData"))


#########################################################################################################################
## Create points: consider changing the coordinate system here to a global projected system?
COMBO.POINTS   = SpatialPointsDataFrame(coords      = GBIF.ALA.COMBO.HIA[c("lon", "lat")], 
                                        data        = GBIF.ALA.COMBO.HIA[c("lon", "lat")],
                                        proj4string = CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))
dim(COMBO.POINTS)





#########################################################################################################################
## 2). EXTRACT ALL WORLDCLIM DATA FOR SPECIES RECORDS
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
## Create a stack of rasters to sample: get all the World clim variables just for good measure
env.grids.current = stack(
  file.path('./data/base/worldclim/world/0.5/bio/current', #/10km
            sprintf('bio_%02d', 1:19)))


## Can we resample the rasters to 10km?
# BIOW_10km = raster("./data/base/worldclim/world/0.5/bio/current/bio_01_10km.tif")
# BIOA_10km = raster("./data/base/worldclim/aus/0.5/bio/current/bio_aus_01_10km.tif")
# 
# 
# ## Use bilinear?
# ten.km.grids.current = resample(env.grids.current, BIOW_10km, method = 'bilinear')
# writeRaster(ten.km.grids.current, 
#             filename = "./data/base/worldclim/aus/0.5/bio/current/10km_current_grids.current.tif", 
#             options = "INTERLEAVE = BAND", overwrite = TRUE)
# 
# 
# # raster_path   = "./data/base/worldclim/aus/0.5/bio/current/"
# # current_10km  <- sprintf('%scurrent_env_10km.tif', raster_path)
# # writeRaster(test, current_10km , overwrite = TRUE)
# xres(ten.km.grids.current)
# yres(ten.km.grids.current)


#########################################################################################################################
## Is there a way to speed this up?
COMBO.RASTER <- extract(env.grids.current, COMBO.POINTS) %>% 
  cbind(GBIF.ALA.COMBO.HIA, .)


## Multiple rename using dplyr
COMBO.RASTER = dplyr::rename(COMBO.RASTER,
                             Annual_mean_temp     = bio_01, 
                             Mean_diurnal_range   = bio_02,
                             Isothermality        = bio_03,
                             Temp_seasonality     = bio_04, 
                             Max_temp_warm_month  = bio_05, 
                             Min_temp_cold_month  = bio_06, 
                             Temp_annual_range    = bio_07,
                             Mean_temp_wet_qu     = bio_08, ## remove
                             Mean_temp_dry_qu     = bio_09, ## remove
                             Mean_temp_warm_qu    = bio_10,
                             Mean_temp_cold_qu    = bio_11,
                             
                             Annual_precip        = bio_12, 
                             Precip_wet_month     = bio_13, 
                             Precip_dry_month     = bio_14, 
                             Precip_seasonality   = bio_15, 
                             Precip_wet_qu        = bio_16, 
                             Precip_dry_qu        = bio_17, 
                             Precip_warm_qu       = bio_18, ## remove
                             Precip_col_qu        = bio_19) ## remove


## Check
dim(COMBO.RASTER)
names(COMBO.RASTER)


#########################################################################################################################
## Change the raster values here: See http://worldclim.org/formats1 for description of the interger conversion. 
## All temperature variables wer multiplied by 10, so divide by 10 to reverse it.
COMBO.RASTER.CONVERT = as.data.table(COMBO.RASTER)                           ## Check this works, also inefficient
COMBO.RASTER.CONVERT[, (env.variables[c(1:11)]) := lapply(.SD, function(x) 
  x / 10 ), .SDcols = env.variables[c(1:11)]]
COMBO.RASTER.CONVERT = as.data.frame(COMBO.RASTER.CONVERT)                    ## Find another method without using data.table


## Check Looks ok?
summary(COMBO.RASTER.CONVERT$Annual_mean_temp)  ## -23 looks too low/. Check where these are ok
summary(COMBO.RASTER$Annual_mean_temp)

summary(COMBO.RASTER.CONVERT$Isothermality)
summary(COMBO.RASTER$Isothermality)




#########################################################################################################################
## 3). FLAG OUTLIERS
#########################################################################################################################


##






## Save the summary datasets
save(COMBO.RASTER.CONTEXT, file = paste("./data/base/HIA_LIST/COMBO/COMBO_RASTER_CONTEXT_10km.RData", sep = ""))
save.image("STEP_4_NICHES_TEST_SPP.RData")





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