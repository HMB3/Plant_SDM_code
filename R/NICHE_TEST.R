#########################################################################################################################
####################################################  NICHE TEST ######################################################## 
#########################################################################################################################


## See how different the niches are, using all data vs. unique cells.
GBIF.ALA.COMBO$OBS = 1:nrow(GBIF.ALA.COMBO)


#########################################################################################################################
## 1). CREATE POINT DATA
#########################################################################################################################


#########################################################################################################################
## Create all points: the 'over' function seems to need geographic coordinates for this data...
COMBO.POINTS   = SpatialPointsDataFrame(coords      = GBIF.ALA.COMBO[c("lon", "lat")],
                                        data        = GBIF.ALA.COMBO[c("lon", "lat")],
                                        proj4string = CRS.WGS.84)


## Check the logic of if the niche is different when calculate this way, insted of using all data.
## Would the median, 95-5, etc, be the same using only unique values?...................................................


## Now get the XY centroids of the unique 1km * 1km WORLDCLIM blocks where ALA records are found
## Get cell number(s) of WORLDCLIM raster from row and/or column numbers. Cell numbers start at 1 in the upper left corner, 
## and increase from left to right, and then from top to bottom. The last cell number equals the number of raster cells 
COMBO.UNIQUE <- cellFromXY(world.grids.current, COMBO.POINTS) %>%

  ## get the unique raster cells
  unique %>%

  ## Get coordinates of the center of raster cells for a row, column, or cell number of WORLDCLIM raster
  xyFromCell(world.grids.current, .) %>%

  as.data.frame() %>%

  SpatialPointsDataFrame(coords = ., data = .,
                         proj4string = CRS.WGS.84)
  
names(COMBO.UNIQUE) <- c("lon", "lat")


## Check
message(round(dim(COMBO.UNIQUE)[1]/dim(COMBO.POINTS)[1]*100, 2), "% of all records utilised for unique cells")
head(COMBO.POINTS);head(COMBO.UNIQUE)

dim(COMBO.POINTS)
projection(COMBO.POINTS)
names(COMBO.POINTS)





#########################################################################################################################
## 2). EXTRACT ALL WORLDCLIM DATA FOR SPECIES RECORDS
#########################################################################################################################


#########################################################################################################################
## This step is a bottleneck - create niche summaries using both unique cells, and also all cells, and plot them (LM).
## If they are not different, then consider using the unique version. Email john, and include both code and output.
## Then we can decide what makes sense


#########################################################################################################################
## Extract worldclim data for all records
COMBO.RASTER <- raster::extract(world.grids.current, COMBO.POINTS) %>% 
  
  ## Join on the species column, using the unique identifier "OBS"
  cbind(GBIF.ALA.COMBO[c("searchTaxon", "lon", "lat", "OBS")], .) %>% 
  
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
head(COMBO.RASTER);dim(COMBO.RASTER)
gc();gc()


#########################################################################################################################
## Extract worldclim data for all unique records
COMBO.RASTER.UNIQUE <- raster::extract(world.grids.current, COMBO.UNIQUE) %>% 
  
  ## cbind
  cbind(COMBO.UNIQUE, .) %>% 
  
  as.data.frame() %>% 
  
  ## Join on the species column, using the unique identifier "OBS"
  join(COMBO.RASTER[c("lon", "lat", "searchTaxon")], ., type = "right") %>%
  #join(GBIF.ALA.COMBO[c("lon", "lat", "searchTaxon")], ., type = "right")
  
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
head(COMBO.RASTER);head(COMBO.RASTER.UNIQUE)
gc();gc()





#########################################################################################################################
## 3). ESTIMATE NICHES
#########################################################################################################################

## Estimate niches using all records
message('Estimating niche differences for ', length(GBIF.spp), ' species across ', length(env.variables), ' climate variables')
ALL.NICHE.DF    = completeFun(COMBO.RASTER,        "Annual_mean_temp")
UNIQUE.NICHE.DF = completeFun(COMBO.RASTER.UNIQUE, "Annual_mean_temp")


#########################################################################################################################
## Pipe the environmental variables into niche estimate function
ALL.NICHE <- env.variables %>% 
  
  ## Pipe the list into lapply
  lapply(function(x) {
    
    ## Now use the niche width function on each colname (so 8 environmental variables)
    ## Also, need to figure out how to make the aggregating column generic (species, genus, etc.)
    ## currently it only works hard-wired
    niche_estimate (DF = ALL.NICHE.DF, colname = x)
    
    ## would be good to remove the duplicate columns here
    
  }) %>% 
  
  ## finally, create one dataframe for all niches
  as.data.frame


#########################################################################################################################
## Pipe the environmental variables into niche estimate function
UNIQUE.NICHE <- env.variables %>% 
  
  ## Pipe the list into lapply
  lapply(function(x) {
    
    ## Now use the niche width function on each colname (so 8 environmental variables)
    ## Also, need to figure out how to make the aggregating column generic (species, genus, etc.)
    ## currently it only works hard-wired
    niche_estimate (DF = UNIQUE.NICHE.DF, colname = x)
    
    ## would be good to remove the duplicate columns here
    
  }) %>% 
  
  ## finally, create one dataframe for all niches
  as.data.frame


#########################################################################################################################
## Now show niches and boxplots
View(ALL.NICHE)
View(UNIQUE.NICHE)


## Updte from here



