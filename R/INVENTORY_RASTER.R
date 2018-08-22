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
## Read in all data to run the SDM code :: species lists, shapefile, rasters & tables
#source('./R/HIA_LIST_MATCHING.R')
rasterTmpFile()





#########################################################################################################################
## 1). CREATE SPDF FOR INVENTORY DATA
#########################################################################################################################


#########################################################################################################################
## Restrict the inventory species to just the analysed species
TI.XY.SPP = TI.XY[TI.XY$searchTaxon %in% GBIF.spp, ]
 

#########################################################################################################################
## Create points: the over function seems to need geographic coordinates for this data
TI.POINTS   = SpatialPointsDataFrame(coords      = TI.XY.SPP[c("lon", "lat")], 
                                     data        = TI.XY.SPP[c("lon", "lat")],
                                     proj4string = CRS.WGS.84)


## Check
dim(TI.POINTS)
projection(TI.POINTS)
names(TI.POINTS)





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
aus.grids.current = stack(
  file.path('./data/base/worldclim/aus/1km/bio/current', 
            sprintf('bio_%02d.tif', 1:19))) 


## Also get the PET raster
PET               = raster("./data/base/worldclim/world/1km/pet_he_yr1.tif")
PET <- PET %>%
  projectRaster(crs = CRS.WGS.84)


#########################################################################################################################
## Check the projection and raster extents for worldclim vs aus data
aus.grids.current <- aus.grids.current %>%
  projectRaster(crs = CRS.WGS.84)
projection(TI.POINTS);projection(aus.grids.current)
saveRDS("./data/base/worldclim/aus/1km/bio/current/aus_grids_current.rds")

TI.RASTER <- extract(aus.grids.current, TI.POINTS) %>% 
  cbind(TI.XY.SPP, .)
summary(TI.RASTER)
class(TI.RASTER)



#########################################################################################################################
## Now try to extract the nearest values for the NA rasters. This is just too slow to bother for only 112 points
# TI.NA  = TI.RASTER[is.na(TI.RASTER$bio_01),]
# NA.XY  = TI.NA[c("lon", "lat")]
# 
# 
# ##
# TI.NA <- names(aus.grids.current) %>%
# 
#   ## pipe the list of raster names into lapply
#   lapply(function(x) {
# 
#     ## Create the raster values
#     r = aus.grids.current[[x]]
#     r@data@values = values(r)
#     str(r@data@values)
#     xy = NA.XY
#     dim(xy)
# 
#     ## The sample the nearest non-NA raster values
#     message("sampling ", dim(xy)[1], "points from raster") 
#     sampled = apply(X = xy, MARGIN = 1, FUN = function(xy) r@data@values[which.min(replace(distanceFromPoints(r, xy), is.na(r), NA))])
# 
#     ## Check the data
#     head(sampled)
# 
#   }) %>%
# 
#   ## Finally, bind all the rows together
#   cbind
# 
# 
# ## Save data
# saveRDS(TI.NA, 'data/base/HIA_LIST/COMBO/SPAT_OUT/TREE_INV_NA_RASTER.rds')
# TI.ALL = cbind(TI.XY.SPP, TI.RASTER, TI.NA)
# indentical(dim(TI.ALL)[1], dim(TI.RASTER)[1])
# 
# 
# saveRDS(TI.ALL, 'data/base/HIA_LIST/COMBO/SPAT_OUT/TREE_INV_ALL_POINTS.rds')



#########################################################################################################################
## Multiple rename using dplyr
TI.RASTER = dplyr::rename(TI.RASTER,
                          
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
#summary(TI.RASTER)
#saveRDS(TI.RASTER, file = paste("./data/base/HIA_LIST/GBIF/TI_GBIF_ALA_RASTER.rds"))
gc();gc()


#########################################################################################################################
## Extract the raster data for PET
TI.POINTS   = SpatialPointsDataFrame(coords         = TI.RASTER[c("lon", "lat")], 
                                     data           = TI.RASTER[c("lon", "lat")],
                                     proj4string    = CRS.WGS.84)


projection(TI.POINTS);projection(PET)
POINTS.PET <- extract(PET, TI.POINTS) %>% 
  cbind(TI.RASTER, .)
TI.RASTER = POINTS.PET
names(TI.RASTER)[names(TI.RASTER) == "."] <- 'PET'


## Check 
dim(TI.RASTER)
names(TI.RASTER)
projection(TI.RASTER)


## Check the raster values here..........................................................................................
summary(TI.RASTER$Annual_mean_temp)





#########################################################################################################################
## 4). CONVERT RASTER VALUES
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
TI.RASTER.CONVERT = as.data.table(TI.RASTER)                           ## Check this works, also inefficient
TI.RASTER.CONVERT[, (env.variables [c(1:11)]) := lapply(.SD, function(x) 
  x / 10 ), .SDcols = env.variables [c(1:11)]]
TI.RASTER.CONVERT = as.data.frame(TI.RASTER.CONVERT)                   ## Find another method without using data.table


## Check Looks ok?
summary(TI.RASTER.CONVERT$Annual_mean_temp)  ## -23 looks too low/. Check where these are ok
summary(TI.RASTER$Annual_mean_temp)

summary(TI.RASTER.CONVERT$Isothermality)
summary(TI.RASTER$Isothermality)


#########################################################################################################################
## Which records are NA?
TI.NA      = TI.RASTER.CONVERT[is.na(TI.RASTER.CONVERT$Annual_mean_temp),]
TI.NA.DF   = SpatialPointsDataFrame(coords      = TI.NA[c("lon", "lat")],
                                    data        = TI.NA,
                                    proj4string = CRS.WGS.84)

aus.mol <- aus %>%
  spTransform(CRS.WGS.84)

plot(aus.mol, main = "NA Tree inventory records")
points(TI.NA.DF, col = "red", cex = .15, pch = 19)


## Save the shapefile, to be subsampled in ArcMap
# names(TI.NA.DF)
# writeOGR(obj = TI.NA.DF, dsn = "./data/base/HIA_LIST/COMBO", layer = "TI_NA", driver = "ESRI Shapefile")


########################################################################################################################
## How big is the dataset?
message(round(dim(TI.NA)[1]/dim(TI.RASTER.CONVERT)[1]*100, 2), "% of points are outside worldclim extent")
dim(TI.RASTER.CONVERT)
identical(length(unique(TI.RASTER.CONVERT$searchTaxon)), length(GBIF.spp))


## Plot a few points to see :: do those look reasonable?
# plot(LAND, col = 'grey', bg = 'sky blue')
# points(TI.RASTER.CONVERT[ which(TI.RASTER.CONVERT$Annual_mean_temp < -5), ][, c("lon", "lat")], 
#        pch = ".", col = "red", cex = 3, asp = 1, main = "temp records < -5")


#########################################################################################################################
## Save the raster datasets
TI.RASTER.CONVERT = na.omit(TI.RASTER.CONVERT)
saveRDS(TI.RASTER.CONVERT, file = paste("./data/base/HIA_LIST/COMBO/TI_RASTER_CONVERT.rds", sep = ""))





#########################################################################################################################
## OUTSTANDING NICHE TASKS:
#########################################################################################################################


## Clean species data with Alessandro's method 





#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################