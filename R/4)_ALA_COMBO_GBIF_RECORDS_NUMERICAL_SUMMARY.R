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
#source('./R/HIA_LIST_MATCHING.R')
#rasterOptions(tmpdir = file.path("'H:/green_cities_sdm/RTEMP")) 


## GBIF and ALA data
#GBIF.LAND = readRDS("./data/base/HIA_LIST/GBIF/GBIF_LAND_POINTS.rds")
load("./data/base/HIA_LIST/ALA/ALA_LAND_POINTS.RData")


## Alessandros urban inventory data......................................................................................
# AUS.URBAN   = readRDS("./data/base/HIA_LIST/URBAN/AUS_URBAN_POINTS.rds")
# AUS.NURSE   = readRDS("./data/base/HIA_LIST/URBAN/BUSH_NURSERY_POINTS.rds")
# AUS.INAT    = readRDS("./data/base/HIA_LIST/URBAN/INAT_POINTS.rds")





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


#########################################################################################################################
## Rename a few fields: Add rename of cultivated field, and make the values the same as mine "CULT" or "UNKNOWN"
## X = CULTIVATED
ALA.LAND     = dplyr::rename(ALA.LAND, 
                             searchTaxon                   = scientificname,
                             coordinateUncertaintyInMeters = uncertainty_m)


## Restrict ALA data to just those species on the big list
HIA.SPP.JOIN     = CLEAN.SPP
HIA.SPP.JOIN     = dplyr::rename(HIA.SPP.JOIN, searchTaxon = Binomial)


## Set NA to blank, then sort by no. of growers to get them to the top
HIA.SPP.JOIN[is.na(HIA.SPP.JOIN)] <- 0
HIA.SPP.JOIN = HIA.SPP.JOIN[with(HIA.SPP.JOIN, rev(order(Number.of.growers))), ]
head(HIA.SPP.JOIN[, c("searchTaxon", "Number.of.growers")])
View(HIA.SPP.JOIN)


## Get just those ALA species which are on the combined list of HIA and planted/growing
## ALA.LAND.HIA  = ALA.LAND[ALA.LAND$searchTaxon %in% HIA.SPP.JOIN$searchTaxon, ] 
ALA.LAND.HIA  = ALA.LAND[ALA.LAND$searchTaxon %in% unique(GBIF.LAND$searchTaxon), ] ## 
str(unique(ALA.LAND$searchTaxon))       ##
str(unique(ALA.LAND.HIA$searchTaxon))   ## Reduced from 30k to 4K


#########################################################################################################################
## Bind the rows together?
GBIF.ALA.COMBO.LAND = bind_rows(GBIF.LAND, ALA.LAND.HIA)
names(GBIF.ALA.COMBO.LAND)
identical((dim(GBIF.LAND)[1]+dim(ALA.LAND)[1]),dim(GBIF.ALA.COMBO.LAND)[1])      ## Only adding the ovelap...
head(GBIF.ALA.COMBO.LAND)


## Check the new data frame has just the species on the HIA list
length(unique(GBIF.ALA.COMBO.LAND$searchTaxon))                                  ## ok


## How many risky species are there?
summary(RISK.BINOMIAL  %in% GBIF.LAND$searchTaxon)


#########################################################################################################################
## Now crunch the big dataset down to just the species on the 25 growers or more list: the extras are just overkill......
#GBIF.ALA.COMBO.HIA  = GBIF.ALA.COMBO.LAND[GBIF.ALA.COMBO.LAND$searchTaxon %in% unique(c(test.spp, HIA.SPP$Binomial)), ]
GBIF.ALA.COMBO.HIA  = GBIF.ALA.COMBO.LAND[GBIF.ALA.COMBO.LAND$searchTaxon %in% unique(GBIF.LAND$searchTaxon), ]
length(unique(GBIF.ALA.COMBO.HIA$searchTaxon))


## This unique ID column can be applied across the project  
GBIF.ALA.COMBO.HIA$OBS <- 1:nrow(GBIF.ALA.COMBO.HIA)
dim(GBIF.ALA.COMBO.HIA)[1];length(GBIF.ALA.COMBO.HIA$OBS)  
names(GBIF.ALA.COMBO.HIA)


#########################################################################################################################
## Create points: consider changing the coordinate system here to a global projected system?
CRS.MOL      <- CRS('+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs')
CRS.WGS.84   <- CRS("+init=epsg:4326")
CRS.AUS.ALB  <- CRS("+init=EPSG:3577")


## The points table can be in global mollweide :: how can we 
COMBO.POINTS   = SpatialPointsDataFrame(coords      = GBIF.ALA.COMBO.HIA[c("lon", "lat")], 
                                        data        = GBIF.ALA.COMBO.HIA[c("lon", "lat")],
                                        proj4string = CRS.WGS.84)


## Check
dim(COMBO.POINTS)
projection(COMBO.POINTS)
names(COMBO.POINTS)





#########################################################################################################################
## 2). MERGE OCCURRENCE DATA WITH COUNCIL AND OTHER URBAN DATA
#########################################################################################################################


# ## what does the data look like?
# names(COMBO.POINTS);names(AUS.URBAN);names(AUS.NURSE);names(AUS.INAT)
# 
# 
# ## Bind the rows together?
# GBIF.ALA.COMBO.LAND = bind_rows(COMBO.POINTS, AUS.URBAN, AUS.NURSE, AUS.INAT)
# names(GBIF.ALA.COMBO.LAND)





#########################################################################################################################
## 3). PROJECT RASTERS AND EXTRACT ALL WORLDCLIM DATA FOR SPECIES RECORDS
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
## Use the Mollweide projection for the points and rasters 
env.grids.current = stack(
  file.path('./data/base/worldclim/world/0.5/bio/current',
            sprintf('bio_%02d', 1:19))) 

## Also get the PET raster
PET               = raster("./data/base/worldclim/world/1km/pet_he_yr1.tif")


#########################################################################################################################
## Is there a way to speed this up?
projection(COMBO.POINTS);projection(env.grids.current)
COMBO.RASTER <- extract(env.grids.current, COMBO.POINTS) %>% 
  cbind(GBIF.ALA.COMBO.HIA, .)


## Multiple rename using dplyr
COMBO.RASTER = dplyr::rename(COMBO.RASTER,
                             
                             ## Temperature
                             Annual_mean_temp     = bio_01,
                             Mean_diurnal_range   = bio_02,
                             Isothermality        = bio_03,
                             Temp_seasonality     = bio_04,
                             Max_temp_warm_month  = bio_05,
                             Min_temp_cold_month  = bio_06,
                             Temp_annual_range    = bio_07,
                             Mean_temp_wet_qu     = bio_08,
                             Mean_temp_dry_qu     = bio_09,
                             Mean_temp_warm_qu    = bio_10,
                             Mean_temp_cold_qu    = bio_11,
                             
                             ## Rainfall
                             Annual_precip        = bio_12,
                             Precip_wet_month     = bio_13,
                             Precip_dry_month     = bio_14,
                             Precip_seasonality   = bio_15,
                             Precip_wet_qu        = bio_16,
                             Precip_dry_qu        = bio_17,
                             Precip_warm_qu       = bio_18,
                             Precip_col_qu        = bio_19)


## Save/load
#summary(COMBO.RASTER)
#saveRDS(COMBO.RASTER, file = paste("./data/base/HIA_LIST/GBIF/COMBO_GBIF_ALA_RASTER.rds"))
gc();gc()


#########################################################################################################################
## Extract the raster data for PET
COMBO.POINTS   = SpatialPointsDataFrame(coords      = COMBO.RASTER[c("lon", "lat")], 
                                        data        = COMBO.RASTER[c("lon", "lat")],
                                        proj4string = CRS.WGS.84)


projection(COMBO.POINTS);projection(PET)
POINTS.PET <- extract(PET, COMBO.POINTS) %>% 
  cbind(COMBO.RASTER, .)
#saveRDS(POINTS.PET, file = paste("./data/base/HIA_LIST/GBIF/COMBO_GBIF_ALA_RASTER_PET.rds"))
COMBO.RASTER = POINTS.PET
colnames(COMBO.RASTER)[62] <- "PET"


## Check 
dim(COMBO.RASTER)
names(COMBO.RASTER)
projection(COMBO.RASTER)





#########################################################################################################################
## 4). CREATE NICHES FOR SELECTED TAXA
#########################################################################################################################


#########################################################################################################################
## Now summarise the niches. But figure out a cleaner way of doing this
env.variables = c("Annual_mean_temp",
                  "Mean_diurnal_range",
                  "Isothermality",
                  "Temp_seasonality",
                  "Max_temp_warm_month",
                  "Min_temp_cold_month",
                  "Temp_annual_range",
                  "Mean_temp_wet_qu",
                  "Mean_temp_dry_qu",
                  "Mean_temp_warm_qu",
                  "Mean_temp_cold_qu",

                  "Annual_precip",
                  "Precip_wet_month",
                  "Precip_dry_month",
                  "Precip_seasonality",
                  "Precip_wet_qu",
                  "Precip_dry_qu",
                  "Precip_warm_qu",
                  "Precip_col_qu",
                  "PET")


#########################################################################################################################
## Change the raster values here: See http://worldclim.org/formats1 for description of the interger conversion. 
## All temperature variables wer multiplied by 10, so divide by 10 to reverse it.
COMBO.RASTER.CONVERT = as.data.table(COMBO.RASTER)                           ## Check this works, also inefficient
COMBO.RASTER.CONVERT[, (env.variables [c(1:11)]) := lapply(.SD, function(x) 
  x / 10 ), .SDcols = env.variables [c(1:11)]]
COMBO.RASTER.CONVERT = as.data.frame(COMBO.RASTER.CONVERT)                   ## Find another method without using data.table


## Check Looks ok?
summary(COMBO.RASTER.CONVERT$Annual_mean_temp)  ## -23 looks too low/. Check where these are ok
summary(COMBO.RASTER$Annual_mean_temp)

summary(COMBO.RASTER.CONVERT$Isothermality)
summary(COMBO.RASTER$Isothermality)


## Plot a few points to see :: do those look reasonable?
# plot(LAND, col = 'grey', bg = 'sky blue')
# points(COMBO.RASTER.CONVERT[ which(COMBO.RASTER.CONVERT$Annual_mean_temp < -5), ][, c("lon", "lat")], 
#        pch = ".", col = "red", cex = 3, asp = 1, main = "temp records < -5")




#########################################################################################################################
## Combine the existing data with the new species :: use this approach to just process a subset 
## COMBO.RASTER.CONTEXT.UPDATE = COMBO.RASTER.CONTEXT
## saveRDS(COMBO.RASTER.CONTEXT.UPDATE, file = paste("./data/base/HIA_LIST/COMBO/COMBO_RASTER_CONTEXT_UPDATE.rds", sep = ""))
## load("./data/base/HIA_LIST/COMBO/COMBO_RASTER_CONTEXT.rds")
## load("./data/base/HIA_LIST/COMBO/COMBO_RASTER_CONTEXT_UPDATE.rds")
# names(COMBO.RASTER.CONTEXT);names(COMBO.RASTER.CONTEXT.UPDATE)
# COMBO.RASTER.CONTEXT.UPDTATE = bind_rows(COMBO.RASTER.CONTEXT, COMBO.RASTER.CONTEXT.UPDATE)
       

## Save the summary datasets
#saveRDS(COMBO.RASTER.CONVERT, file = paste("./data/base/HIA_LIST/COMBO/COMBO_RASTER_CONVERT_APRIL_2018.rds", sep = ""))
#write.csv(COMBO.NICHE.CONTEXT, "./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT_APRIL_2018_STANDARD_CLEAN.csv",        row.names = FALSE)



## Now save .RData file for the next session...
#save.image("STEP_4_NICHES.RData")





#########################################################################################################################
## OUTSTANDING NICHE TASKS:
#########################################################################################################################


## Integrate Ale's dataset

## Integrate Manuel's dataset

## Integrate Bush's data, etc.





#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################