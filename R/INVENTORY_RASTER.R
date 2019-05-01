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
unique(TI.XY.SPP$searchTaxon)


message('Extracting raster values for ', 
        length(unique(TI.XY.SPP$searchTaxon)), ' urban species across ',
        length(unique(TI.XY.SPP$INVENTORY)),   ' Councils ')


#########################################################################################################################
## Check the inventory table works
if(nrow(TI.XY.SPP) > 0) {
  
  
  #########################################################################################################################
  ## Create points: the 'over' function seems to need geographic coordinates for this data...
  TI.XY.84   = SpatialPointsDataFrame(coords      = TI.XY.SPP[c("lon", "lat")],
                                      data        = TI.XY.SPP,
                                      proj4string = CRS.WGS.84)
  
  
  ## Now split using the data using the species column, and get the unique occurrence cells
  ## Create a template raster in WGS84 projection
  template.raster.1km.84 = raster("./data/world_koppen/template_1km_WGS84.tif")
  
  TI.XY.84.SPLIT.ALL <- split(TI.XY.84, as.character(unique(TI.XY.84$searchTaxon)))
  inventory_cells_all  <- lapply(TI.XY.84.SPLIT.ALL, function(x) cellFromXY(template.raster.1km.84, x))
  length(inventory_cells_all)   ## this is a list of dataframes, where the number of rows for each being the species table
  
  
  #########################################################################################################################
  ## Now get just one record within each 10*10km cell.
  TI.XY.84.1KM <- mapply(function(x, cells) {
    x[!duplicated(cells), ]
  }, TI.XY.84.SPLIT.ALL, inventory_cells_all, SIMPLIFY = FALSE) %>% do.call(rbind, .)
  
  
  ## Check to see we have 19 variables + the species for the standard predictors, and 19 for all predictors
  ## create two more template.raster files: 5km amnd 10km
  ## Check data :: template, data table and species
  message(round(nrow(TI.XY.84.1KM)/nrow(TI.XY.84)*100, 2), " % records retained at 1km resolution") 
  TI.POINTS   = TI.XY.84.1KM[c("lon", "lat")]

  
  ## Check
  dim(TI.POINTS);dim(TI.XY.84.1KM)
  projection(TI.POINTS)
  names(TI.POINTS)
  
  
  ## Also get the PET raster
  projection(PET);projection(world.grids.current)

  
  
  #########################################################################################################################
  ## Now extract data
  projection(TI.POINTS);projection(world.grids.current);projection(PET)
  TI.RASTER <- raster::extract(world.grids.current, TI.POINTS) %>% 
    cbind(as.data.frame(TI.XY.84.1KM), .)
  
  
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
  
  
  ## How many species were processed?
  length(unique(TI.RASTER.CONVERT$INVENTORY))
  message(length(unique(TI.RASTER.CONVERT$searchTaxon)), 
          ' species processed in ', length(unique(TI.RASTER.CONVERT$INVENTORY)), ' councils')
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