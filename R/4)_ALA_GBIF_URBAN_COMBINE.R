#########################################################################################################################
############################################# ESTIMATE SPECIES NICHES ################################################### 
#########################################################################################################################


#########################################################################################################################
## This code combines the GBIF records with the ALA records, then extracts environmental values. 
## It creates :

## 1). A large table with one row for each species record


## This table is subsequently used to estimate the current global realised niche/climatic tolerance, snd susequently model 
## the niches using the maxent algorithm.  


#########################################################################################################################
## Read in all data to run the SDM code :: species lists, shapefile, rasters & tables
#source('./R/HIA_LIST_MATCHING.R')
rasterTmpFile()


## Print the species run to the screen
message('Extracting Worldclim data for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")


if(read_data == "TRUE") {
  
  ## Load GBIF and ALA data
  GBIF.LAND      = readRDS(paste0('data/base/HIA_LIST/GBIF/GBIF_TREES_LAND_', save_run, '.rds'))
  ALA.TREES.LAND = readRDS(paste0('data/base/HIA_LIST/ALA/ALA_TREES_LAND_',   save_run, '.rds'))
  
} else {
  
  message(' skip file reading, not many species analysed')   ##
  
}




#########################################################################################################################
## 1). MERGE GBIF DATA WITH CLEANED AVH/ALA RECORDS
#########################################################################################################################


#########################################################################################################################
## Merge on the ALA data. Consider that GBIF has data for both sources. We are topping up the native ranges with the AVH. 
## So there could be duplicates between both sources.
length(unique(GBIF.LAND$searchTaxon)) 
length(unique(ALA.TREES.LAND$searchTaxon)) 
unique(GBIF.LAND$SOURCE)
unique(ALA.TREES.LAND$SOURCE)
sort(intersect(sort(names(GBIF.LAND)), sort(names(ALA.TREES.LAND))))


dim(GBIF.LAND);dim(ALA.TREES.LAND)
setdiff(names(GBIF.LAND), names(ALA.TREES.LAND))


## Rename a few fields
# GBIF.LAND     = dplyr::rename(GBIF.LAND, 
#                               #coordinateUncertaintyInMetres = coordinateUncertaintyInMeters,
#                               rank = taxonRank)


#########################################################################################################################
## Bind the rows together
common.cols = intersect(names(GBIF.LAND), names(ALA.TREES.LAND))
GBIF.ALA.COMBO = bind_rows(GBIF.LAND, ALA.TREES.LAND)
GBIF.ALA.COMBO = GBIF.ALA.COMBO %>% 
  select(one_of(common.cols))


names(GBIF.ALA.COMBO)
unique(GBIF.ALA.COMBO$SOURCE)
identical((dim(GBIF.LAND)[1]+dim(ALA.TREES.LAND)[1]),dim(GBIF.ALA.COMBO)[1])      ## Only adding the ovelap...

head(GBIF.ALA.COMBO)
length(unique(GBIF.ALA.COMBO$searchTaxon))
length(unique(GBIF.ALA.COMBO$scientificName)) 


#########################################################################################################################
## Check these GBIF/ALA data :: No na year or lat/long, all years > 1950
summary(GBIF.ALA.COMBO$year)
summary(GBIF.ALA.COMBO$lat)
summary(GBIF.ALA.COMBO$lon)


## Check NAs again
(sum(is.na(GBIF.ALA.COMBO$scientificName)) + dim(subset(GBIF.ALA.COMBO, scientificName == ""))[1])/dim(GBIF.ALA.COMBO)[1]*100


#########################################################################################################################
## Now create table of species counts
COMBO.LUT = as.data.frame(table(GBIF.ALA.COMBO$scientificName))
names(COMBO.LUT) = c("scientificName", "FREQUENCY")
COMBO.LUT = COMBO.LUT[with(COMBO.LUT, rev(order(FREQUENCY))), ] 
head(COMBO.LUT);dim(COMBO.LUT)


#########################################################################################################################
## save data
if(save_data == "TRUE") {
  
  ## save .rds file for the next session
  write.csv(COMBO.LUT, "./data/base/HIA_LIST/COMBO/SUA_TREES_GBIF_ALA_LUT.csv", row.names = FALSE)
  
} else {
  
  message(' skip file saving, not many species analysed')   ##
  
}


#########################################################################################################################
## Create points: the 'over' function seems to need geographic coordinates for this data...
COMBO.POINTS   = SpatialPointsDataFrame(coords      = GBIF.ALA.COMBO[c("lon", "lat")],
                                        data        = GBIF.ALA.COMBO[c("lon", "lat")],
                                        proj4string = CRS.WGS.84)


## Check the logic of if the niche is different when calculate this way, insted of using all data.
## Would the median, 95-5, etc, be the same using only unique values?...................................................


## Now get the XY centroids of the unique 1km * 1km WORLDCLIM blocks where ALA records are found
## Get cell number(s) of WORLDCLIM raster from row and/or column numbers. Cell numbers start at 1 in the upper left corner, 
## and increase from left to right, and then from top to bottom. The last cell number equals the number of raster cells 
# COMBO.POINTS <- cellFromXY(world.grids.current, GBIF.ALA.COMBO[c("lon", "lat")]) %>% 
#   
#   ## get the unique raster cells
#   unique %>% 
#   
#   ## Get coordinates of the center of raster cells for a row, column, or cell number of WORLDCLIM raster
#   xyFromCell(world.grids.current, .) %>%
#   
#   as.data.frame() %>%
#   
#   SpatialPointsDataFrame(coords = ., data = .,
#                          proj4string = CRS.WGS.84)


## Check
dim(COMBO.POINTS)
projection(COMBO.POINTS)
names(COMBO.POINTS)





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
## This step is a bottleneck - create niche summaries using both unique cells, and also all cells, and plot them (LM).
## If they are not different, then consider using the unique version. Email john, and include both code and output.
## Then we can decide what makes sense


## Also get the PET raster
PET               = raster("./data/base/worldclim/world/1km/pet_he_yr1.tif")


#########################################################################################################################
## Extract worldclim data
## This step is a bottle neck ...........................................................................................
message('Extracting raster values for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")
projection(COMBO.POINTS);projection(world.grids.current)
dim(COMBO.POINTS);dim(GBIF.ALA.COMBO)

COMBO.RASTER <- raster::extract(world.grids.current, COMBO.POINTS) %>% 
  
  cbind(GBIF.ALA.COMBO, .) %>% 
  
  dplyr::rename(
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


## Free some memory
gc();gc()


#########################################################################################################################
## Extract the raster data for PET
projection(COMBO.POINTS);projection(PET)
dim(COMBO.POINTS)

POINTS.PET <- raster::extract(PET, COMBO.POINTS) %>% 
  cbind(COMBO.RASTER, .)
COMBO.RASTER = POINTS.PET
names(COMBO.RASTER)[names(COMBO.RASTER) == "."] <- 'PET'


## Check 
dim(COMBO.RASTER)
names(COMBO.RASTER)
summary(COMBO.RASTER$Annual_mean_temp)
summary(COMBO.RASTER$PET)





#########################################################################################################################
## 3). CONVERT RASTER VALUES
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
## All temperature variables were multiplied by 10, so divide by 10 to reverse it.
COMBO.RASTER.CONVERT = as.data.table(COMBO.RASTER)                           ## Check this works, also inefficient
COMBO.RASTER.CONVERT[, (env.variables [c(1:11)]) := lapply(.SD, function(x) 
  x / 10 ), .SDcols = env.variables [c(1:11)]]
COMBO.RASTER.CONVERT = as.data.frame(COMBO.RASTER.CONVERT)                   ## Find another method without using data.table


## Check Looks ok?
summary(COMBO.RASTER.CONVERT$Annual_mean_temp)  ## -23 looks too low/. Check where these are ok
summary(COMBO.RASTER$Annual_mean_temp)

summary(COMBO.RASTER.CONVERT$Isothermality)
summary(COMBO.RASTER$Isothermality)

summary(COMBO.RASTER.CONVERT$PET)


## Print the dataframe dimensions to screen :: format to recognise millions, hundreds of thousands, etc.
names(COMBO.RASTER.CONVERT)
dim(COMBO.RASTER.CONVERT)
formatC(dim(COMBO.RASTER.CONVERT)[1], format = "e", digits = 2)
length(unique(COMBO.RASTER.CONVERT$searchTaxon));length(GBIF.spp)


## Plot a few points to see :: do those look reasonable?
# plot(LAND, col = 'grey', bg = 'sky blue')
# points(COMBO.RASTER.CONVERT[ which(COMBO.RASTER.CONVERT$Annual_mean_temp < -5), ][, c("lon", "lat")], 
#        pch = ".", col = "red", cex = 3, asp = 1, main = "temp records < -5")



#########################################################################################################################
## save data
if(save_data == "TRUE") {
  
  ## save .rds file for the next session
  saveRDS(COMBO.RASTER.CONVERT, paste0('data/base/HIA_LIST/COMBO/COMBO_RASTER_CONVERT_', save_run, '.rds'))
  
} else {
  
  message(' skip file saving, not many species analysed')   ##
  
}

## get rid of some memory
gc()





#########################################################################################################################
## OUTSTANDING NICHE TASKS:
#########################################################################################################################


## Use better ALA data when available



#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################