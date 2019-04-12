#########################################################################################################################
################################################# FLAG SPATIAL OUTLIERS #################################################
#########################################################################################################################


#########################################################################################################################
## This code uses the CoordinateCleaner package to Automatcially screen out the dodgy spatial records. See ::
## https://github.com/azizka/CoordinateCleaner


## It then creates a table of the species niche :

## 1). A table with one row for each species, including contextual data and species attributes (niches, traits, etc.)

## These tables are subsequently used to estimate the current global realised niche/climatic tolerance using the best 
## available data, and susequently model the niches using the maxent algorithm.


## Update line 548 RE niches.............................................................................................


#########################################################################################################################
## Michelle's advice :: use all the global data.
## Does including the global data improve the bad species?
## Re-import Alessandro's inventory data and re-run the models...........................................................


## Print the species run to the screen
message('Cleaning geogrphic outliers for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")


#########################################################################################################################
## Read in the global data
if(read_data == "TRUE") {
  
  ## Read in RDS files from previous step
  GBIF.ALA.MATCH = readRDS(paste0(DATA_path, 'GBIF_ALA_MATCH_',  save_run, '.rds'))
  
} else {
  
  message(' skip file reading, not many species analysed')   ##
  
}


#########################################################################################################################
## Restrict the tree inventory data to just the species analysed
TI.XY.SPP = TI.XY[TI.XY$searchTaxon %in% GBIF.spp, ]


## Restrict the tree inventory data to just the species analysed
message('Using tree inventory data for ', 
        length(unique(TI.XY.SPP$searchTaxon)), ' urban species across ',
        length(unique(TI.XY.SPP$INVENTORY)),   ' Councils ')

## Check dimensions of the occurrence and inventory data tables.
length(unique(GBIF.ALA.MATCH$searchTaxon))
formatC(nrow(GBIF.ALA.MATCH), format = "e", digits = 2)
names(GBIF.ALA.MATCH)


#########################################################################################################################
## Create a unique identifier. This is used for automated cleaing of the records, and also saving shapefiles
## But this will not be run for all species linearly. So, it probably needs to be a combination of species and number
GBIF.ALA.MATCH$CC.OBS <- 1:nrow(GBIF.ALA.MATCH)
GBIF.ALA.MATCH$CC.OBS <- paste0(GBIF.ALA.MATCH$CC.OBS, "_CC_", GBIF.ALA.MATCH$searchTaxon)
GBIF.ALA.MATCH$CC.OBS <- gsub(" ",     "_",  GBIF.ALA.MATCH$CC.OBS, perl = TRUE)
length(GBIF.ALA.MATCH$CC.OBS);length(unique(GBIF.ALA.MATCH$CC.OBS))





#########################################################################################################################
## 1). FLAG GBIF OUTLIERS
#########################################################################################################################


#########################################################################################################################
## Rename the columns to fit the CleanCoordinates format and create a tibble. 
TIB.GBIF <- GBIF.ALA.MATCH %>% dplyr::rename(coord_spp        = searchTaxon,
                                             decimallongitude = lon, 
                                             decimallatitude  = lat) %>%
  
  ## The create a tibble for running the spatial outlier cleaning
  timetk::tk_tbl() %>% 
  
  ## Consider the arguments. We've already stripped out the records that fall outside
  ## the worldclim raster boundaries, so the sea test is probably not the most important
  ## Study area is the globe, but we are only projecting models onto Australia
  
  ## The geographic outlier detection is not working here. Try it in setp 6
  ## Also, the duplicates step is not working, flagging too many records
  clean_coordinates(.,
                    verbose         = TRUE,
                    tests = c("capitals",     "centroids", "equal", "gbif", 
                              "institutions", "zeros"), ## duplicates flagged too many
                    
                    capitals_rad    = 10000,  ## remove records within 10km  of capitals
                    centroids_rad   = 5000    ## remove records within 5km of country centroids
                    # outliers_method = "distance", ## Save the other checks for step 6
                    # outliers_td     = 800,
                    # outliers_mtp    = 5
  ) %>%
  
  ## The select the relevant columns and rename
  select(., coord_spp, CC.OBS, .val,  .equ, .zer,  .cap,   
         .cen,    .gbf,   .inst, .summary) 

## Then rename
summary(TIB.GBIF)
names(TIB.GBIF) = c("coord_spp", "CC.OBS",    "coord_val",  "coord_equ",  "coord_zer",  "coord_cap", 
                    "coord_cen", "coord_gbf", "coord_inst", "coord_summary")


## Flagging ~ x%, excluding the spatial outliers. Seems reasonable?
message(round(with(TIB.GBIF, table(coord_summary)/sum(table(coord_summary))*100), 2), " % records removed")


#########################################################################################################################
## Check the order still matches
identical(GBIF.ALA.MATCH$CC.OBS, TIB.GBIF$CC.OBS)                                               ## order matches
identical(GBIF.ALA.MATCH$searchTaxon, TIB.GBIF$coord_spp)                                       ## order matches


## Is the species column the same as the searchTaxon column?
TEST.GEO   = join(GBIF.ALA.MATCH, TIB.GBIF)
identical(GBIF.ALA.MATCH$searchTaxon, TEST.GEO$coord_spp)                                       ## order matches 
identical(GBIF.ALA.MATCH$CC.OBS,      TEST.GEO$CC.OBS)                                          ## order matches


## Now subset to records that are flagged as outliers
CLEAN.TRUE  = subset(TEST.GEO, coord_summary == "TRUE")
CLEAN.FALSE = subset(TEST.GEO, coord_summary == "FALSE")
table(CLEAN.TRUE$coord_summary) 


## What percentage of records are retained?
length(unique(CLEAN.TRUE$searchTaxon))
message(round(nrow(CLEAN.TRUE)/nrow(TEST.GEO)*100, 2), " % records retained")                                               





#########################################################################################################################
## 2). JOIN ALA, GBIF & INVENTORY DATA
#########################################################################################################################


#########################################################################################################################
## Now bind on the urban tree inventory data. We are assuming this data is clean, after we manually fix the taxonomy
if(nrow(TI.XY.SPP) > 0) {
  
  message('Combining Australian inventory data with occurrence data') 
  intersect(names(TI.XY.SPP), names(CLEAN.TRUE))
  CLEAN.INV = bind_rows(CLEAN.TRUE, TI.XY.SPP)
  
} else {
  
  ## Update with global data
  message('No Australian inventory data for these species')   ##
  CLEAN.INV = CLEAN.TRUE
  
}


## How many local municipalities have data for the species analysed
## Remove duplicate coordinates
drops <- c("lon.1", "lat.1")
CLEAN.INV <- CLEAN.INV[ , !(names(CLEAN.INV) %in% drops)]
unique(CLEAN.INV$SOURCE) 
length(unique(CLEAN.INV$INVENTORY))


length(unique(CLEAN.INV$searchTaxon))
summary(CLEAN.INV$Annual_mean_temp)
summary(CLEAN.INV$Annual_mean_temp)


## By how many % does including tree inventories increase the overal number of records?
message("Tree inventory data increases records by ", round(nrow(CLEAN.INV)/nrow(TEST.GEO)*100, 2), " % ")


#########################################################################################################################
## Plot GBIF outliers to check. This might be overkill, but useful to interrogate why species didn't work
#source('./R/CC_CLEAN_TEST.R', echo = TRUE)


#########################################################################################################################
## save data
if(save_data == "TRUE") {
  
  ## save .rds file for the next session
  saveRDS(CLEAN.INV, paste0(DATA_path, 'CLEAN_INV_', save_run, '.rds'))
  
} else {
  
  message(' skip file saving, not many species analysed')   ##
  
}


#########################################################################################################################
## Plot GBIF outliers to check. This might be overkill, but useful to interrogate why species didn't work
if(check_maps == "TRUE") {
  
  message('Writing shapefiles and maps to checking directory') 
  source('./R/CC_CLEAN_TEST.R')
  
} else {
  
  ## Update with global data
  message('Dont create maps and shapefile of points')   ##
  
}





#########################################################################################################################
#####################################################  TBC ##############################################################
#########################################################################################################################