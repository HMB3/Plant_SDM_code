#########################################################################################################################
############################################# ESTIMATE SPECIES NICHES ################################################### 
#########################################################################################################################


#########################################################################################################################
## This code take the table of inventory data for the species list, extracts environmental values 


## 1). This code takes the table of inventory data for the species list, extracts environmental values and creates 
## a large table with one row for each species record.


## These tables are subsequently used to estimate the current global realised niche/climatic tolerance using the best 
## available data, and susequently model the niches using the maxent algorithm. 


## Could speed up the extract step, by creating an index of all possible raster cells, then only using the unique cells 
## where occurrence data is found..................................................................................... 


#########################################################################################################################
## Read in all data to run the SDM code :: species lists, shapefile, rasters & tables
#source('./R/HIA_LIST_MATCHING.R')
rasterTmpFile()


#########################################################################################################################
## 1). PROJECT RASTERS AND EXTRACT ALL WORLDCLIM DATA FOR SPECIES RECORDS
#########################################################################################################################


#########################################################################################################################
## Restrict the inventory species to just the analysed species
TI.XY.SPP = TI.XY[TI.XY$searchTaxon %in% GBIF.spp, ]
length(unique(TI.XY.SPP$searchTaxon))


message('Extracting raster values for ', 
        length(unique(TI.XY.SPP$searchTaxon)), ' urban species across ',
        length(unique(TI.XY.SPP$INVENTORY)),   ' Councils ')


#########################################################################################################################
## Check the inventory table works
if(nrow(TI.XY.SPP) > 0) {
  
  
  #########################################################################################################################
  ## Create points: the over function seems to need geographic coordinates for this data
  TI.POINTS   = SpatialPointsDataFrame(coords      = TI.XY.SPP[c("lon", "lat")],
                                       data        = TI.XY.SPP[c("lon", "lat")],
                                       proj4string = CRS(projection(world.grids.current)))

  
  ## Check
  dim(TI.POINTS)
  projection(TI.POINTS)
  names(TI.POINTS)
  
  
  ## Now get the XY centroids of the unique 1km * 1km WORLDCLIM blocks where ALA records are found
  ## Get cell number(s) of WORLDCLIM raster from row and/or column numbers. Cell numbers start at 1 in the upper left corner, 
  ## and increase from left to right, and then from top to bottom. The last cell number equals the number of raster cells 
  # TI.POINTS <- cellFromXY(world.grids.current, TI.XY.SPP[c("lon", "lat")]) %>% 
  #   
  #   ## get the unique raster cells
  #   unique %>% 
  #   
  #   ## Get coordinates of the center of raster cells for a row, column, or cell number of WORLDCLIM raster
  #   xyFromCell(world.temp, .) %>%
  #   
  #   as.data.frame() %>%
  #   
  #   SpatialPointsDataFrame(coords = ., data = .,
  #                          proj4string = CRS.WGS.84)
  
  
  ## Check
  dim(TI.POINTS)
  projection(TI.POINTS)
  names(TI.POINTS)
  
  #########################################################################################################################
  ## Create a stack of rasters to sample: get all the Worldclim variables just for good measure
  ## Use the Mollweide projection for the points and rasters 
  

  ## Also get the PET raster
  projection(PET);projection(world.grids.current)
  
  #########################################################################################################################
  ## Check the projection and raster extents for worldclim vs aus data
  # world.grids.current <- world.grids.current %>%
  #   projectRaster(crs = CRS.WGS.84)
  # saveRDS(world.grids.current, "./data/base/worldclim/aus/1km/bio/current/aus_grids_current.rds")
  # projection(TI.POINTS);projection(world.grids.current)
  # 
  # 
  # ## Save projected files
  # saveRDS(world.grids.current, "./data/base/worldclim/aus/1km/bio/current/aus_grids_current.rds")
  # world.grids.current = readRDS("./data/base/worldclim/aus/1km/bio/current/aus_grids_current.rds")
  
  
  #########################################################################################################################
  ## Now extract data
  ## > 200k points are coming out as NA. What could be the reason? Outside the raster extent 
  projection(TI.POINTS);projection(world.grids.current);projection(PET)
  TI.RASTER <- raster::extract(world.grids.current, TI.POINTS) %>% 
    cbind(TI.XY.SPP, .)
  
  
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
  
  gc()
  
  
  #########################################################################################################################
  ## Extract the raster data for PET
  projection(TI.POINTS);projection(PET)
  POINTS.PET <- raster::extract(PET, TI.POINTS) %>%
    cbind(TI.RASTER, .)
  TI.RASTER = POINTS.PET
  names(TI.RASTER)[names(TI.RASTER) == "."] <- 'PET'
  
  ## Check 
  dim(TI.RASTER)
  dim(TI.POINTS)
  names(TI.RASTER)
  
  
  ## Check the raster values here..........................................................................................
  ## Need a measure of how many points are outside the raster exent
  summary(TI.RASTER$Annual_mean_temp)
  
  
  
  #########################################################################################################################
  ## 2). CONVERT RASTER VALUES
  #########################################################################################################################
  
  
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
  ## Save the raster datasets
  TI.RASTER.CONVERT = completeFun(TI.RASTER.CONVERT, "PET")
  summary(TI.RASTER.CONVERT)
  
  
  ## Need a measure of how many points are outside the raster exent
  message(round(nrow(TI.RASTER.CONVERT)/nrow(TI.RASTER)*100, 2), " % records retained after worldclim extract") 
  
  
  names(TI.RASTER.CONVERT)
  unique(TI.RASTER.CONVERT$INVENTORY)
  unique(TI.RASTER.CONVERT$searchTaxon)
  saveRDS(TI.RASTER.CONVERT, paste0(DATA_path, 'TI_RASTER_CONVERT_', save_run, '.rds'))
  
  
} else {
  
  message(' skipped urban raster extraction, these species have no inventory records')   ##
  
}





#########################################################################################################################
## OUTSTANDING RASTER TASKS:
#########################################################################################################################


## Figure out how to extract points which are NA with a buffer, in either R or ArcMap



#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################