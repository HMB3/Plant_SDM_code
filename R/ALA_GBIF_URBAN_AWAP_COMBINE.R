#########################################################################################################################
############################################# ESTIMATE SPECIES NICHES ################################################### 
#########################################################################################################################


#########################################################################################################################
## This code combines the GBIF records with the ALA records, then extracts environmental values. 


## It creates :

## 1). A large table with one row for each species record


## This tables are subsequently used to estimate the current global realised niche/climatic tolerance 
## and susequently model the niches using the maxent algorithm.  


#########################################################################################################################
## Read in all data to run the SDM code :: species lists, shapefile, rasters & tables
#source('./R/HIA_LIST_MATCHING.R')
rasterTmpFile()

GBIF.LAND = readRDS(paste0('data/base/HIA_LIST/GBIF/GBIF_TREES_LAND_', save_run, '.rds'))
ALA.TREES.LAND = readRDS(paste0('data/base/HIA_LIST/ALA/ALA_TREES_LAND_', save_run, '.rds'))



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
GBIF.LAND     = dplyr::rename(GBIF.LAND, 
                              coordinateUncertaintyInMetres = coordinateUncertaintyInMeters,
                              rank = taxonRank)


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


## Write out the table
#write.csv(COMBO.LUT, "./data/base/HIA_LIST/COMBO/SUA_TREES_GBIF_ALA_LUT.csv", row.names = FALSE)


#########################################################################################################################
## Create points: the 'over' function seems to need geographic coordinates for this data...
COMBO.POINTS   = SpatialPointsDataFrame(coords      = GBIF.ALA.COMBO[c("lon", "lat")], 
                                        data        = GBIF.ALA.COMBO[c("lon", "lat")],
                                        proj4string = CRS.WGS.84)


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
## Create a stack of rasters to sample: get all the Worldclim variables just for good measure
## Use the Mollweide projection for the points and rasters 
env.grids.current = stack(
  file.path('./data/base/worldclim/world/0.5/bio/current',
            sprintf('bio_%02d', 1:19)))


## Also get the PET raster
PET               = raster("./data/base/worldclim/world/1km/pet_he_yr1.tif")


#########################################################################################################################
## Also extract the netCDF AWAP drought rasters
awap.extreme  = stack(list.files(as.character('./data/base/AWAP'), pattern = 'extreme_droughts.nc$',  full.names = TRUE))
awap.severe   = stack(list.files(as.character('./data/base/AWAP'), pattern = 'severe_droughts.nc$',   full.names = TRUE))
awap.moderate = stack(list.files(as.character('./data/base/AWAP'), pattern = 'moderate_droughts.nc$', full.names = TRUE))

names(awap.extreme)
names(awap.severe)
names(awap.moderate)

names(awap.extreme) = c("Drought_freq_extr", "Drought_max_dur_extr", "Drought_max_int_extr", "Drought_max_rel_int_extr",
                        "Drought_mean_dur_extr", "Drought_mean_int_extr", "Drought_mean_rel_int_extr")

names(awap.severe) = c("Drought_freq_extr", "Drought_max_dur_extr", "Drought_max_int_extr", "Drought_max_rel_int_extr",
                       "Drought_mean_dur_extr", "Drought_mean_int_extr", "Drought_mean_rel_int_extr")

names(awap.moderate) = c("Drought_freq_extr", "Drought_max_dur_extr", "Drought_max_int_extr", "Drought_max_rel_int_extr",
                         "Drought_mean_dur_extr", "Drought_mean_int_extr", "Drought_mean_rel_int_extr")



names(awap.extreme)
projection(awap.extreme)
plot(awap.extreme[["Drought_max_dur_extr"]])
plot(awap.extreme[["Drought_max_int_extr"]])
plot(awap.extreme[["Drought_mean_rel_int_extr"]])


#########################################################################################################################
## Also extract the netCDF BOM Heatwave rasters
## read in each file using nc_open
#awap.heatwave  = stack(list.files(as.character('./data/base/AWAP'), pattern = 'AWAP_HW',  full.names = TRUE))

HWF = raster('./data/base/AWAP/AWAP_HWF.nc')  #number of heatwave days in a season
names(HWF) = "HWF"

HWA = raster('./data/base/AWAP/AWAP_HWA.nc')  #hottest day of hottest heatwave
names(HWA) = "HWA"

HWM =  raster('./data/base/AWAP/AWAP_HWM.nc')  #average heatwave intensity
names(HWM) = "HWM"

HWD =  raster('./data/base/AWAP/AWAP_HWD.nc') #length of longest event
names(HWD) = "HWD"

HWN = raster('./data/base/AWAP/AWAP_HWN.nc')  #number of individual heatwaves
names(HWN) = "HWN"

cum_all =  raster('./data/base/AWAP/AWAP_HW_cum_all.nc') #the total temperature across all heatwaves in the season
names(cum_all) = "HW_CUM_ALL"

cum_av = raster('./data/base/AWAP/AWAP_HW_cum_av.nc') #the average total temperature across all heatwaves in the season
names(cum_av) = "HW_CUM_AV"

cum_hot = raster('./data/base/AWAP/AWAP_HW_cum_hottest.nc') #the cumulative heat during the hottest heatwave.
names(cum_hot) = "HW_CUM_HOT"


## rename
awap.heatwave = stack(HWF, HWA, HWM, HWD, HWN, cum_all, cum_av, cum_hot)
names(AWAP.heatwave)
plot(HWF)
plot(cum_hot)


#########################################################################################################################
## Extract worldclim data
projection(COMBO.POINTS);projection(env.grids.current)
dim(COMBO.POINTS);dim(GBIF.ALA.COMBO)

COMBO.RASTER <- raster::extract(env.grids.current, COMBO.POINTS) %>% 
  
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


## Save/load
#summary(COMBO.RASTER)
#saveRDS(COMBO.RASTER, file = paste("./data/base/HIA_LIST/GBIF/COMBO_GBIF_ALA_RASTER.rds"))
gc();gc()


#########################################################################################################################
## Extract the raster data for PET
projection(COMBO.POINTS);projection(PET)
dim(COMBO.POINTS)

POINTS.PET <- raster::extract(PET, COMBO.POINTS) %>% 
  cbind(COMBO.RASTER, .)
COMBO.RASTER = POINTS.PET
names(COMBO.RASTER)[names(COMBO.RASTER) == "."] <- 'PET'


#########################################################################################################################
## Extract AWAP drought and heatwave data
projection(COMBO.POINTS);projection(awap.extreme);projection(awap.heatwave)

COMBO.DROUGHT <- raster::extract(awap.extreme, COMBO.POINTS)
COMBO.HEAT    <- raster::extract(awap.heatwave, COMBO.POINTS) 
COMBO.AWAP = cbind(COMBO.RASTER, COMBO.DROUGHT, COMBO.HEAT)



## Check 
dim(COMBO.AWAP)
names(COMBO.AWAP)

summary(COMBO.AWAP$Annual_mean_temp)
summary(COMBO.AWAP$PET)

summary(COMBO.AWAP$Drought_mean_rel_int_extr)  ## 8365 NA values for drought
summary(COMBO.AWAP$Drought_max_int_ext)        ## 45089 NA values for heatwave

summary(COMBO.AWAP$HWA)
summary(COMBO.AWAP$HW_CUM_HOT)



#########################################################################################################################
## 3). CONVERT AWAP VALUES
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
COMBO.AWAP.CONVERT = as.data.table(COMBO.AWAP)                           ## Check this works, also inefficient
COMBO.AWAP.CONVERT[, (env.variables [c(1:11)]) := lapply(.SD, function(x) 
  x / 10 ), .SDcols = env.variables [c(1:11)]]
COMBO.AWAP.CONVERT = as.data.frame(COMBO.AWAP.CONVERT)                   ## Find another method without using data.table


## Check Looks ok?
summary(COMBO.AWAP.CONVERT$Annual_mean_temp)  ## -23 looks too low/. Check where these are ok
summary(COMBO.AWAP$Annual_mean_temp)

summary(COMBO.AWAP.CONVERT$Isothermality)
summary(COMBO.AWAP$Isothermality)

summary(COMBO.AWAP.CONVERT$PET)


## Print the dataframe dimensions to screen :: format to recognise millions, hundreds of thousands, etc.
names(COMBO.AWAP.CONVERT)
dim(COMBO.AWAP.CONVERT)
formatC(dim(COMBO.AWAP.CONVERT)[1], format = "e", digits = 2)
length(unique(COMBO.AWAP.CONVERT$searchTaxon));length(GBIF.spp)


## Plot a few points to see :: do those look reasonable?
# plot(LAND, col = 'grey', bg = 'sky blue')
# points(COMBO.AWAP.CONVERT[ which(COMBO.AWAP.CONVERT$Annual_mean_temp < -5), ][, c("lon", "lat")], 
#        pch = ".", col = "red", cex = 3, asp = 1, main = "temp records < -5")


#########################################################################################################################
## Save the summary datasets
saveRDS(COMBO.AWAP.CONVERT, paste0('data/base/HIA_LIST/COMBO/COMBO_AWAP_CONVERT_', save_run, '.rds'))


## Now save .RData file for the next session...
## save.image("TRAITS_RASTER.RData")





#########################################################################################################################
## OUTSTANDING NICHE TASKS:
#########################################################################################################################


## Use better ALA data when available



#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################