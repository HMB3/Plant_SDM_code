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
message('Cleaning outliers and creating niches for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")


#########################################################################################################################
## Read in the global data
if(read_data == "TRUE") {
  
  ## Read in RDS files from previous step
  COMBO.RASTER.CONVERT = readRDS(paste0(DATA_path, 'COMBO_RASTER_CONVERT_',  save_run, '.rds'))
  
} else {
  
  message(' skip file reading, not many species analysed')   ##
  
}


#########################################################################################################################
## Read in the global data
if(dim(TI.XY.SPP)[1] > 0 & read_data == "TRUE") {
  
  ## Read in RDS files from previous step
  TI.RASTER.CONVERT = readRDS(paste0(DATA_path, 'TI_RASTER_CONVERT_', save_run, '.rds'))
  
} else {
  
  message(' skip file reading, no urban records')   ##
  
}


## Check dimensions of the occurrence and inventory data tables.
length(unique(COMBO.RASTER.CONVERT$searchTaxon))
formatC(dim(COMBO.RASTER.CONVERT)[1], format = "e", digits = 2)
names(COMBO.RASTER.CONVERT)


#########################################################################################################################
## Create a unique identifier
COMBO.RASTER.CONVERT$CC.OBS <- 1:nrow(COMBO.RASTER.CONVERT)





#########################################################################################################################
## 1). FLAG GBIF OUTLIERS
#########################################################################################################################


#########################################################################################################################
## Rename the columns to fit the CleanCoordinates format and create a tibble. 
## A Tibble is needed for running the spatial outlier cleaning
TIB.GBIF <- COMBO.RASTER.CONVERT %>% dplyr::rename(species          = searchTaxon,
                                                   decimallongitude = lon, 
                                                   decimallatitude  = lat) %>%
  timetk::tk_tbl()


## Add a column for unique observation so we can check the records match up after joining
#TIB.GBIF$CC.OBS <- 1:nrow(TIB.GBIF)
identical(length(TIB.GBIF$CC.OBS), dim(TIB.GBIF)[1])


## We've already stripped out the records that fall outside
## the worldclim raster boundaries, so the sea test is probably not the most important
## Study area is the globe, but we are only projecting models onto Australia


#########################################################################################################################
## Don't run the outliers test here, it is slower. Also, can't run cleaning on the urban tree inventory data, because this
## removes all the records near capital cities 
message('Flagging GBIF outliers for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")
FLAGS  <- CleanCoordinates(TIB.GBIF,
                           #countries        = "country",    ## too many flagged here...
                           capitals.rad     = 0.12,
                           countrycheck     = TRUE,
                           duplicates       = TRUE,
                           seas             = FALSE,
                           verbose          = FALSE)


## save/load the flags
identical(dim(FLAGS)[1], dim(TIB.GBIF)[1])


## Flagging ~ 1.64%, excluding the spatial outliers. Seems reasonable?
summary(FLAGS)
FLAGS = FLAGS[ ,!(colnames(FLAGS) == "decimallongitude" | colnames(FLAGS) == "decimallatitude")]
message(round(summary(FLAGS)[8]/dim(FLAGS)[1]*100, 2), " % records removed")





#########################################################################################################################
## 2). FILTER DATA TO THE CLEAN RECORDS
#########################################################################################################################


#########################################################################################################################
## Check the order still matches
identical(dim(TIB.GBIF)[1],dim(FLAGS)[1])
names(FLAGS)[1]    = c("coord_spp")
identical(COMBO.RASTER.CONVERT$searchTaxon, FLAGS$coord_spp)                                            ## order matches


## Is the species column the same as the searchTaxon column?
TEST.GEO   = cbind(COMBO.RASTER.CONVERT, FLAGS)
identical(TEST.GEO$searchTaxon, TEST.GEO$coord_spp)                                                     ## order matches


## Not sure why the inverse did not work :: get only the records which were not flagged as being dodgy.
dim(subset(TEST.GEO, summary == "TRUE"))
CLEAN.TRUE = subset(TEST.GEO, summary == "TRUE")
unique(CLEAN.TRUE$summary)   


## What percentage of records are retained?
identical(dim(CLEAN.TRUE)[1], (dim(COMBO.RASTER.CONVERT)[1] - dim(subset(TEST.GEO, summary == "FALSE"))[1]))
length(unique(CLEAN.TRUE$searchTaxon))
message(round(dim(CLEAN.TRUE)[1]/dim(TEST.GEO)[1]*100, 2), " % records retained")                                               


#########################################################################################################################
## Now bind on the urban tree inventory data. We are assuming this data is clean, after we manually fix the taxonomy
if(dim(TI.XY.SPP)[1] > 0) {
  
  message('Combining Australian inventory data with occurrence data') 
  intersect(names(TI.RASTER.CONVERT), names(CLEAN.TRUE))
  CLEAN.TRUE = bind_rows(CLEAN.TRUE, TI.RASTER.CONVERT)
  
} else {
  
  ## Update with global data
  message('No Australian inventory data for these species')   ##
  CLEAN.TRUE = CLEAN.TRUE
  
}


## How many local municipalities have data for the species analysed
names(CLEAN.TRUE)
unique(CLEAN.TRUE$SOURCE) 
unique(CLEAN.TRUE$INVENTORY)
length(unique(CLEAN.TRUE$INVENTORY))


length(unique(CLEAN.TRUE$searchTaxon))
summary(CLEAN.TRUE$Annual_mean_temp)
summary(CLEAN.TRUE$Annual_mean_temp)


## Then create a unique ID column which can be used to identify outlier records 
CLEAN.TRUE$OBS <- 1:nrow(CLEAN.TRUE)
dim(CLEAN.TRUE)[1];length(CLEAN.TRUE$OBS)  
identical(dim(CLEAN.TRUE)[1], length(CLEAN.TRUE$OBS))


## By how many % does including tree inventories increase the overal number of records?
message("Tree inventory data increases records by ", round(dim(CLEAN.TRUE)[1]/dim(TEST.GEO)[1]*100, 2), " % ")   


#########################################################################################################################
## save data
if(save_data == "TRUE") {
  
  ## save .rds file for the next session
  saveRDS(CLEAN.TRUE, paste0(DATA_path, 'CLEAN_TRUE_', save_run, '.rds'))
  
} else {
  
  message(' skip file saving, not many species analysed')   ##
  
}



#########################################################################################################################
## 4). CREATE SHAPEFILES TO CHECK OUTLIERS ARCMAP
#########################################################################################################################


## This code could be used to manually clean outlying records, using the OBS column......................................
## Select columns
GBIF.ALA.CHECK  = select(CLEAN.TRUE,     OBS, searchTaxon, scientificName, lat, lon, SOURCE, INVENTORY, year, 
                         country, locality, basisOfRecord, institutionCode, 
                         Taxonomic.status, New.Taxonomic.status)


## Rename the fields so that ArcMap can handle them
GBIF.ALA.CHECK  = dplyr::rename(GBIF.ALA.CHECK, 
                                TAXON     = searchTaxon,
                                LAT       = lat,
                                LON       = lon,
                                SC_NAME   = scientificName,
                                BASIS     = basisOfRecord,                
                                LOCAL     = locality,                      
                                INSTIT    = institutionCode,                
                                COUNTRY   = country,                
                                YEAR      = year,
                                TAX_STAT  = Taxonomic.status,
                                NEW_STA   = New.Taxonomic.status)
names(GBIF.ALA.CHECK)


## Then create a SPDF
GBIF.ALA.SPDF    = SpatialPointsDataFrame(coords      = GBIF.ALA.CHECK[c("LON", "LAT")],
                                          data        = GBIF.ALA.CHECK,
                                          proj4string = CRS.WGS.84)


## Create list of species 
TAXA <- unique(GBIF.ALA.SPDF$TAXON)


#########################################################################################################################
## Then, loop over the species list and create a shapefile for each 
## For lots of taxa, this would be needed to 
# for (i in 1:length(TAXA)) {
#   
#   ## Need to check the OBS column matches up - or do we not need this again?
#   message("Writing shapefile for ", TAXA[i])
#   tmp <- GBIF.ALA.SPDF[GBIF.ALA.SPDF$TAXON == TAXA[i], ] 
#   writeOGR(tmp, dsn = "./data/ANALYSIS/CLEAN_GBIF", TAXA[i], driver = "ESRI Shapefile", overwrite_layer = TRUE)
#   
# }


## Then write the shapefile out just in case
# writeOGR(obj = GBIF.ALA.SPDF, dsn = "./data/ANALYSIS/CLEAN_GBIF", 
#          layer = paste0("CLEAN_SPDF_", save_run), driver = "ESRI Shapefile")





#########################################################################################################################
## 5). INTERSECT SPECIES RECORDS WITH LOCAL GOV AREAS AND SIGNIFICANT URBAN AREAS
#########################################################################################################################


#########################################################################################################################
## Only run this code if we are saving niches 
if(calc_niche == "TRUE") {
  
  ## We want to know the count of species that occur in 'n' LGAs, across a range of climates. Read in LGA and SUA
  length(unique(CLEAN.TRUE$searchTaxon))
  projection(LGA);projection(AUS);projection(SUA_2016 )
  
  
  ## Convert the raster data back into a spdf
  COMBO.RASTER.SP   = SpatialPointsDataFrame(coords      = CLEAN.TRUE[c("lon", "lat")], 
                                             data        = CLEAN.TRUE,
                                             proj4string = CRS.WGS.84)
  
  
  ## Project using a projected rather than geographic coordinate system
  ## Not sure why, but this is needed for the extraction step
  LGA.WGS  = spTransform(LGA,       CRS.WGS.84)
  SUA.WGS  = spTransform(SUA_2016 , CRS.WGS.84)
  AUS.WGS  = spTransform(AUS,       CRS.WGS.84)
  
  
  ## Remove the columns we don't need
  LGA.WGS = LGA.WGS[, c("LGA_CODE16", "LGA_NAME16")] 
  head(SUA.WGS)
  
  
  #########################################################################################################################
  ## Run join between species records and LGAs/SUAs :: Double check they are the same
  ## See the ABS for details :: there are 563 LGAs
  ## http://www.abs.gov.au/ausstats/abs@.nsf/Lookup/by%20Subject/1270.0.55.003~July%202016~Main%20Features~Local%20Government%20Areas%20(LGA)~7
  message('Joining occurence data to SUAs for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")
  
  
  projection(COMBO.RASTER.SP);projection(LGA.WGS);projection(SUA.WGS);projection(AUS.WGS)
  SUA.JOIN      = over(COMBO.RASTER.SP,   SUA.WGS)
  COMBO.SUA.LGA = cbind.data.frame(COMBO.RASTER.SP, SUA.JOIN) 
  
  ## save .rds file for the next session
  #saveRDS(COMBO.SUA.LGA, file = paste0(DATA_path, 'COMBO_SUA_OVER_', save_run, '.rds'))
  
  
  #########################################################################################################################
  ## AGGREGATE THE NUMBER OF SUAs EACH SPECIES IS FOUND IN. NA LGAs ARE OUTSIDE AUS
  SUA.AGG   = tapply(COMBO.SUA.LGA$SUA_NAME16, COMBO.SUA.LGA$searchTaxon, function(x) length(unique(x))) ## group SUA by species name
  SUA.AGG   = as.data.frame(SUA.AGG)
  AUS.AGG   = aggregate(SUA_NAME16 ~ searchTaxon, data = COMBO.SUA.LGA,   function(x) {sum(!is.na(x))}, na.action = NULL)
  SUA.AGG   = cbind.data.frame(AUS.AGG, SUA.AGG)
  names(SUA.AGG) = c("searchTaxon", "AUS_RECORDS", "SUA_COUNT")
  
  
  ## Now create a table of all the SUA's that each species occurrs
  SUA.SPP.COUNT = as.data.frame(table(COMBO.SUA.LGA[["SUA_NAME16"]], COMBO.SUA.LGA[["searchTaxon"]]))
  names(SUA.SPP.COUNT) = c("SUA", "SPECIES", "SUA_RECORDS")
  
  
  ## Save .rds file for the next session
  saveRDS(SUA.SPP.COUNT, paste0(DATA_path, 'SUA_SPP_COUNT_', save_run, '.rds'))
  
  
  ## Check : That's ok, but we want a table of which SUA each species is actually in.
  dim(SUA.AGG)
  head(SUA.AGG)
  
  
  ## Remove duplicate coordinates
  COMBO.SUA.LGA = subset(COMBO.SUA.LGA, select = -c(lon.1, lat.1))
  names(COMBO.SUA.LGA)
  dim(COMBO.SUA.LGA)
  str(unique(COMBO.SUA.LGA$searchTaxon))
  unique(COMBO.SUA.LGA$SOURCE)
  
  
  
  
  
  #########################################################################################################################
  ## 6). CREATE NICHES FOR SELECTED TAXA
  #########################################################################################################################
  
  
  #########################################################################################################################
  ## Now summarise the niches. But figure out a cleaner way of doing this
  message('Estimating global niches for ', length(GBIF.spp), ' species across ', length(env.variables), ' climate variables')
  
  
  if(OCC_SOURCE == "ALA") {
    
    ## save .rds file for the next session
    message('Create niches from ', OCC_SOURCE, ' records')
    COMBO.SUA.LGA = COMBO.SUA.LGA[COMBO.SUA.LGA$SOURCE %in% OCC_SOURCE , ]
    message(unique(COMBO.SUA.LGA$SOURCE))
    
  } else {
    
    message('Create niches from all sources')
    message(unique(COMBO.SUA.LGA$SOURCE))
    
  }
  
  
  #########################################################################################################################
  ## Create niche summaries for each environmental condition like this...
  ## Here's what the function will produce :
  NICHE.DF = completeFun(COMBO.SUA.LGA, "PET")
  dim(NICHE.DF)
  head(niche_estimate (DF = NICHE.DF, colname = "Annual_mean_temp"))  ## including the q05 and q95
  
  
  ## So lets use lapply on the "SearchTaxon"
  ## test = run_function_concatenate(list, DF, "DF, colname = x") 
  message('Estimating global niches for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")
  COMBO.NICHE <- env.variables %>% 
    
    ## Pipe the list into lapply
    lapply(function(x) {
      
      ## Now, use the niche width function on each colname
      ## Also, need to figure out how to make the aggregating column generic (species, genus, etc.)
      ## currently it only works hard-wired..............................................................................
      niche_estimate (DF = NICHE.DF, colname = x)
      
      ## would be good to remove the duplicate columns here..............................................................
      
    }) %>% 
    
    ## finally, create one dataframe for all niches
    as.data.frame
  
  
  ## Remove duplicate Taxon columns and check the output :: would be great to skip these columns when running the function
  names(COMBO.NICHE)
  COMBO.NICHE = subset(COMBO.NICHE, select = -c(searchTaxon.1,  searchTaxon.2,  searchTaxon.3,  searchTaxon.4,
                                                searchTaxon.5,  searchTaxon.6,  searchTaxon.7,  searchTaxon.7,
                                                searchTaxon.8,  searchTaxon.9,  searchTaxon.10, searchTaxon.11,
                                                searchTaxon.12, searchTaxon.13, searchTaxon.14, searchTaxon.15,
                                                searchTaxon.16, searchTaxon.17, searchTaxon.18, searchTaxon.19))
  
  
  #########################################################################################################################
  ## Add counts for each species, and record the total number of taxa processed
  ## dim(COMBO.RASTER.CONVERT);dim(CLEAN.TRUE)
  GLOBAL_RECORDS = as.data.frame(table(COMBO.SUA.LGA$searchTaxon))
  names(GLOBAL_RECORDS) = c("searchTaxon", "GLOBAL_RECORDS")
  identical(dim(GLOBAL_RECORDS)[1], dim(COMBO.NICHE)[1])
  
  
  #########################################################################################################################
  ## Add the counts of Australian records for each species to the niche database
  Total.taxa.processed = dim(COMBO.NICHE)[1]
  COMBO.NICHE  = join(GLOBAL_RECORDS, COMBO.NICHE, type = "right")
  COMBO.NICHE  = join(SUA.AGG, COMBO.NICHE, type = "right")
  
  head(COMBO.NICHE$AUS_RECORDS)
  head(COMBO.NICHE$SUA_COUNT)
  
  names(COMBO.NICHE)
  dim(COMBO.NICHE)
  
  
  #########################################################################################################################
  ## Add the counts of Australian records for each species to the niche database
  # names(COMBO.NICHE)
  # names(SUA.AGG)
  # dim(COMBO.NICHE)
  # dim(SUA.AGG)
  # 
  # COMBO.LGA = join(COMBO.NICHE, SUA.AGG)
  # names(COMBO.LGA)
  # 
  # dim(COMBO.LGA)
  # head(COMBO.LGA$AUS_RECORDS)
  # head(COMBO.LGA$SUA_COUNT)
  
  
  #########################################################################################################################
  ## 7). CALCULATE AREA OF OCCUPANCY
  #########################################################################################################################
  
  
  ## Given a dataframe of georeferenced occurrences of one, or more, taxa, this function provide statistics values 
  ## (Extent of Occurrence, Area of Occupancy, number of locations, number of subpopulations) and provide a preliminary 
  ## conservation status following Criterion B of IUCN. A graphical map output is also available.
  
  
  ## Create a species list to estimate the ranges for. Currently we can't estimate ranges for widely distributed species
  ## A workaround for the trait species is to just use the ALA records. But for the full list of species, we will need to 
  ## change this. Otherwise, we could report the global niche, but the Australian range only. 
  ## Change when Gilles Dauby updates......................................................................................
  spp.geo = as.character(unique(COMBO.SUA.LGA$searchTaxon))
  AOO.DAT = COMBO.SUA.LGA %>%
    subset(., SOURCE == "ALA")
  
  
  #########################################################################################################################
  ## AREA OF OCCUPANCY (AOO)
  ## For every species in the list: calculate the AOO
  ## The extent of occurrence function is causing problems.................................................................
  ## Change this once Giles Dauby updates..................................................................................
  ## x = spp.geo[1]
  GBIF.AOO <- spp.geo %>%
    
    ## Pipe the list into lapply
    lapply(function(x) {
      
      ## Subset the the data frame to calculate area of occupancy according the IUCN.eval
      DF = subset(AOO.DAT, searchTaxon == x)[, c("lat", "lon", "searchTaxon")]
      IUCN.eval(DF, Cell_size_AOO = 10, DrawMap = FALSE, country_map = land, SubPop = FALSE) 
      
    }) %>%
    
    ## Finally, create one dataframe for all niches
    bind_rows
  
  
  names(GBIF.AOO)[names(GBIF.AOO) == 'taxa'] <- 'searchTaxon'
  GBIF.AOO             = select(GBIF.AOO, searchTaxon, EOO, AOO, Category_CriteriaB)
  names(GBIF.AOO)      = c("searchTaxon",  "EOO", "AOO", "ICUN_catb")
  head(GBIF.AOO)
  
  
  #########################################################################################################################
  ## Now join on the geographic range and glasshouse data
  identical(length(GBIF.AOO$AOO), length(GBIF.spp))
  COMBO.NICHE = join(GBIF.AOO, COMBO.NICHE, type = "right")

  
  
  #########################################################################################################################
  ## 8). JOIN ON CONTEXTUAL DATA
  #########################################################################################################################
  
  
  #########################################################################################################################
  ## Now join the horticultural contextual data onto one or both tables ()
  message('Joining contextual data for raster and niche files', length(GBIF.spp), ' species in the set ', "'", save_run, "'")
  COMBO.RASTER.CONTEXT = CLEAN.TRUE
  names(COMBO.RASTER.CONTEXT)
  
  
  #########################################################################################################################
  ## Now join the hort context data and the tree inventory plantings to the niche
  COMBO.NICHE.CONTEXT = join(CLEAN.GROW, COMBO.NICHE, type = "right") 
  COMBO.NICHE.CONTEXT = join(TI.LIST, COMBO.NICHE.CONTEXT, type = "right")
  head(COMBO.NICHE.CONTEXT$AUS_RECORDS)
  head(COMBO.NICHE.CONTEXT$SUA_COUNT)
  
  
  #########################################################################################################################
  ## Now combine the SDM output with the niche context data 
  ## CLEAN.GROW needs to be put through GBIF and TPL........................................................................
  COMBO.NICHE.CONTEXT = COMBO.NICHE.CONTEXT [order(COMBO.NICHE.CONTEXT $searchTaxon),]
  
  
  ## Set NA to blank, then sort by no. of growers
  COMBO.NICHE.CONTEXT$Total_growers[is.na(COMBO.NICHE.CONTEXT$Total.growers)] <- 0
  COMBO.NICHE.CONTEXT = COMBO.NICHE.CONTEXT[with(COMBO.NICHE.CONTEXT, rev(order(Total_growers))), ]
  
  
  ## View the data
  dim(COMBO.RASTER.CONTEXT)
  dim(COMBO.NICHE.CONTEXT)
  #View(COMBO.NICHE.CONTEXT)
  
  
  ## Print the dataframe dimensions to screen
  dim(CLEAN.TRUE)
  dim(COMBO.NICHE.CONTEXT)
  dim(COMBO.RASTER.CONTEXT)
  
  length(unique(CLEAN.TRUE$searchTaxon))
  length(COMBO.NICHE.CONTEXT$searchTaxon)
  length(unique(COMBO.RASTER.CONTEXT$searchTaxon))
  
  unique(COMBO.RASTER.CONTEXT$SOURCE)
  
  
  #########################################################################################################################
  ## save .rds file for the next session
  message('Writing niche and raster data for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")
  saveRDS(COMBO.NICHE.CONTEXT,    paste0(DATA_path, 'COMBO_NICHE_CONTEXT_',  OCC_SOURCE, '_RECORDS_', save_run, '.rds'))
  saveRDS(COMBO.RASTER.CONTEXT,   paste0(DATA_path, 'COMBO_RASTER_CONTEXT_', OCC_SOURCE, '_RECORDS_', save_run, '.rds'))
  write.csv(COMBO.NICHE.CONTEXT,  paste0(DATA_path, 'COMBO_NICHE_CONTEXT_',  OCC_SOURCE, '_RECORDS_', save_run, '.csv'), row.names = FALSE)
  
  
} else {
  
  message(' skip file saving, ', length(GBIF.spp), ' species analysed')   ##
  
}


#########################################################################################################################
## OUTSTANDING CLEANING TASKS:
#########################################################################################################################


## 1). Check the settings of coordinate cleaner..........................................................................


## 2). Check the automatic cleaning results against manual cleaning results..............................................


## 3). Add a calculation of geographic ranges back in ..............................................................................



#########################################################################################################################
#####################################################  TBC ##############################################################
#########################################################################################################################