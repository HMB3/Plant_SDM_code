#########################################################################################################################
################################################# FLAG SPATIAL OUTLIERS #################################################
#########################################################################################################################


## This code uses the CoordinateCleaner package to Automatcially screen out the dodgy spatial records. See ::
## https://github.com/azizka/CoordinateCleaner


## It then creates a table of the species niche :

## 1). A table with one row for each species, including contextual data and species attributes (niches, traits, etc.)

## These tables are subsequently used to estimate the current global realised niche/climatic tolerance using the best 
## available data, and susequently model the niches using the maxent algorithm.


## Print the species run to the screen
message('Cleaning outliers and creating niches for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")


#########################################################################################################################
## Read in the two data tables
if(read_data == "TRUE") {
  
  ## read in RDS files from previous step
  TI.RASTER.CONVERT = readRDS(paste0('data/base/HIA_LIST/COMBO/TI_RASTER_CONVERT_', save_run, '.rds'))
  COMBO.RASTER.CONVERT = readRDS(paste0('data/base/HIA_LIST/COMBO/COMBO_RASTER_CONVERT_', save_run, '.rds'))
  
} else {
  
  message(' skip file reading, not many species analysed')   ##
  
}


## Check dimensions of the occurrence and inventory data tables.
length(unique(COMBO.RASTER.CONVERT$searchTaxon))
formatC(dim(COMBO.RASTER.CONVERT)[1], format = "e", digits = 2)
names(COMBO.RASTER.CONVERT)


length(unique(TI.RASTER.CONVERT$searchTaxon))
formatC(dim(TI.RASTER.CONVERT)[1], format = "e", digits = 2)
names(TI.RASTER.CONVERT)


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
                           seas             = FALSE)


## save/load the flags
identical(dim(FLAGS)[1], dim(TIB.GBIF)[1])
#saveRDS(FLAGS, paste0('data/base/HIA_LIST/COMBO/ALA_GBIF_FLAGS_', save_run, '.rds'))
#FLAGS = readRDS(paste0('data/base/HIA_LIST/COMBO/ALA_GBIF_FLAGS_', save_run, '.rds'))


## Flagging ~ 1.64%, excluding the spatial outliers. Seems reasonable?
summary(FLAGS)
FLAGS = FLAGS[ ,!(colnames(FLAGS) == "decimallongitude" | colnames(FLAGS) == "decimallatitude")]
message(round(summary(FLAGS)[8]/dim(FLAGS)[1]*100, 2), " % records removed")



#########################################################################################################################
## 2). FLAG SPATIAL OUTLIERS
#########################################################################################################################


## Alex Zizka ::
## The outlier function is limited in the amount of records it can process. It uses a distance matrix of all records per 
## species, which means that a species with 200k records will result in a 200,000x200,000 cells matrix, which will 
## probably choke your computer. The latest version of cc_outl includes a subsampling heuristic to address this problem. 
## I think this will work for you case, but it might run for a while


## Check the frequency table first, to see if any species are likely to hit this threshold ::
## So there are a few species with +100k records, this will be hard for the computer 
## Could split them up into species under 200k or not?
# COMBO.LUT = as.data.frame(table(TIB.GBIF$species))
# names(COMBO.LUT) = c("species", "FREQUENCY")
# COMBO.LUT = COMBO.LUT[with(COMBO.LUT, rev(order(FREQUENCY))), ] 
# head(COMBO.LUT);summary(COMBO.LUT$FREQUENCY)  ## Quercus robur, 214, and Fraxinus excelsior, 156
# LUT.100K = as.character(subset(COMBO.LUT, FREQUENCY < 100000)$species)
# LUT.100K = LUT.100K [order(LUT.100K)]


## Unfortunately, the cc_outl function can't handle vectors of a certain size - over 40 GB at least.
## So we are problably better off moving this code to after the data has been thinned by the SDM
## step 


# ## Create a data frame of species name and spatial outlier
# SPAT.OUT <- unique(TIB.GBIF$species) %>%  ## LUT.100K
#   
#   ## pipe the list of species into lapply
#   lapply(function(x) {
#     
#     ## Create the species df by subsetting by species
#     f <- subset(TIB.GBIF, species == x)
#     
#     ## Run the spatial outlier detection
#     message("Running spatial outlier detection for ", x)
#     message(dim(f)[1], " records for ", x)
#     sp.flag <- cc_outl(f,
#                        lon     = "decimallongitude",
#                        lat     = "decimallatitude",
#                        species = "species",
#                        method  = "distance",
#                        #mltpl   = 5,
#                        tdi     = 300,
#                        value   = "flags")
#     
#     ## Now add attache column for species, and the flag for each record
#     #d = data.frame(searchTaxon = x, SPAT_OUT = sp.flag)
#     d = cbind(species = x,
#               SPAT_OUT = sp.flag, f)[c("species", "SPAT_OUT", "CC.OBS")]
#     
#   }) %>%
#   
#   ## Finally, bind all the rows together
#   bind_rows


## Save data
# identical(dim(SPAT.OUT)[1], dim(TIB.GBIF)[1])
# unique(SPAT.OUT$species)
# summary(SPAT.OUT$SPAT_OUT)
# head(SPAT.OUT)
# saveRDS(SPAT.OUT, paste0('data/base/HIA_LIST/COMBO/ALA_GBIF_SPAT_OUT_', save_run, '.rds'))





#########################################################################################################################
## 3). FILTER DATA TO THE CLEAN RECORDS
#########################################################################################################################


#########################################################################################################################
## Join data :: exclude the decimal lat/long, check the length 
identical(dim(TIB.GBIF)[1],dim(FLAGS)[1])#;length(GBIF.SPAT.OUT)
names(FLAGS)[1]    = c("coord_spp")
identical(COMBO.RASTER.CONVERT$searchTaxon, FLAGS$coord_spp)                                            ## order matches


#########################################################################################################################
## Is the species column the same as the searchTaxon column?
TEST.GEO   = cbind(COMBO.RASTER.CONVERT, FLAGS)
identical(TEST.GEO$searchTaxon, TEST.GEO$coord_spp)                                                     ## order matches


## Not sure why the inverse did not work :: get only the records which were not flagged as being dodgy.
dim(subset(TEST.GEO, summary == "TRUE"))
CLEAN.TRUE = subset(TEST.GEO, summary == "TRUE")
unique(CLEAN.TRUE$summary)   
#unique(CLEAN.TRUE$SPAT_OUT)   


## What percentage of records are retained?
identical(dim(CLEAN.TRUE)[1], (dim(COMBO.RASTER.CONVERT)[1] - dim(subset(TEST.GEO, summary == "FALSE"))[1]))
length(unique(CLEAN.TRUE$searchTaxon))
message(round(dim(CLEAN.TRUE)[1]/dim(TEST.GEO)[1]*100, 2), " % records retained")                                               


#########################################################################################################################
## Record how many points were removed, just from the analysis species
## GBIF, ALA 1144726 + 1257152
length(native.sua)
CLEAN.NATIVE = CLEAN.TRUE[CLEAN.TRUE$searchTaxon %in% native.sua, ]
length(unique(CLEAN.NATIVE$searchTaxon))
message(round(dim(CLEAN.NATIVE)[1]/2401878*100, 2), " % records retained")


#########################################################################################################################
## Now bind on the urban tree inventory data. We are assuming this data is clean, after we manually fix the taxonomy
## Check the NAs
intersect(names(TI.RASTER.CONVERT), names(CLEAN.TRUE))
CLEAN.TRUE = bind_rows(CLEAN.TRUE, TI.RASTER.CONVERT)


names(CLEAN.TRUE)
unique(CLEAN.TRUE$SOURCE) 
unique(CLEAN.TRUE$INVENTORY) 
length(unique(CLEAN.TRUE$searchTaxon))
summary(CLEAN.TRUE$Annual_mean_temp)


## Then create a unique ID column which can be used to identify outlier records 
CLEAN.TRUE$OBS <- 1:nrow(CLEAN.TRUE)
dim(CLEAN.TRUE)[1];length(CLEAN.TRUE$OBS)  
identical(dim(CLEAN.TRUE)[1], length(CLEAN.TRUE$OBS))


## How many records are added by including the tree inventories?
message("Tree inventory data increases records by ", round(dim(CLEAN.TRUE)[1]/dim(TEST.GEO)[1]*100, 2), " % ")   


#########################################################################################################################
## save data
if(save_data == "TRUE") {
  
  ## save .rds file for the next session
  saveRDS(CLEAN.TRUE, paste0('data/base/HIA_LIST/COMBO/CLEAN_TRUE_', save_run, '.rds'))
  
} else {
  
  message(' skip file saving, not many species analysed')   ##
  
}



#########################################################################################################################
## 4). CREATE SHAPEFILES TO CHECK OUTLIERS ARCMAP
#########################################################################################################################


## Select columns
GBIF.ALA.CHECK  = select(CLEAN.TRUE,     OBS, searchTaxon, scientificName, lat, lon, SOURCE, INVENTORY, year, 
                         #coordinateUncertaintyInMetres, #geodeticDatum,  
                         country, locality, basisOfRecord, institutionCode, 
                         #rank, 
                         Taxonomic.status, New.Taxonomic.status)


## Rename the fields so that ArcMap can handle them
GBIF.ALA.CHECK     = dplyr::rename(GBIF.ALA.CHECK, 
                                   TAXON     = searchTaxon,
                                   LAT       = lat,
                                   LON       = lon,
                                   SC_NAME   = scientificName,
                                   #RANK      = rank,
                                   BASIS     = basisOfRecord,                
                                   LOCAL     = locality,                      
                                   INSTIT    = institutionCode,                
                                   COUNTRY   = country,                
                                   #COORD_UN  = coordinateUncertaintyInMetres,
                                   #DATUM     = geodeticDatum,                 
                                   YEAR      = year,
                                   TAX_STAT  = Taxonomic.status,
                                   NEW_STA   = New.Taxonomic.status)
names(GBIF.ALA.CHECK)


#########################################################################################################################
## Then create SPDF
GBIF.ALA.SPDF    = SpatialPointsDataFrame(coords      = GBIF.ALA.CHECK[c("LON", "LAT")],
                                          data        = GBIF.ALA.CHECK,
                                          proj4string = CRS.WGS.84)


## Create list of species 
TAXA <- unique(GBIF.ALA.SPDF$TAXON)


#########################################################################################################################
## Then, loop over the species list and create a shapefile for each 
# for (i in 1:length(TAXA)) {
#   
#   ## Need to check the OBS column matches up - or do we not need this again?
#   message("Writing shapefile for ", TAXA[i])
#   tmp <- GBIF.ALA.SPDF[GBIF.ALA.SPDF$TAXON == TAXA[i], ] 
#   writeOGR(tmp, dsn = "./data/base/HIA_LIST/COMBO/CLEAN_GBIF", TAXA[i], driver = "ESRI Shapefile", overwrite_layer = TRUE)
#   
# }


## Write the shapefile out just in case
# writeOGR(obj = GBIF.ALA.SPDF, dsn = "./data/base/HIA_LIST/COMBO", layer = paste0("CLEAN_SPDF_", save_run), driver = "ESRI Shapefile")





#########################################################################################################################
## 5). INTERSECT SPECIES RECORDS WITH LOCAL GOV AREAS AND SIGNIFICANT URBAN AREAS
#########################################################################################################################


#########################################################################################################################
## Only run this code if we are saving niches 
if(save_data == "TRUE") {
  
  ## We want to know the count of species that occur in 'n' LGAs, across a range of climates. Read in LGA and SUA
  length(unique(CLEAN.TRUE$searchTaxon))
  projection(LGA);projection(AUS);projection(SUA.16)
  
  
  ## Convert the raster data back into a spdf
  COMBO.RASTER.SP   = SpatialPointsDataFrame(coords      = CLEAN.TRUE[c("lon", "lat")], 
                                             data        = CLEAN.TRUE,
                                             proj4string = CRS.WGS.84)
  
  
  ## Project using a projected rather than geographic coordinate system
  LGA.WGS  = spTransform(LGA,    CRS.WGS.84)
  SUA.WGS  = spTransform(SUA.16, CRS.WGS.84)
  AUS.WGS  = spTransform(AUS,    CRS.WGS.84)
  
  
  ## Remove the columns we don't need
  LGA.WGS = LGA.WGS[, c("LGA_CODE16", "LGA_NAME16")] 
  head(SUA.WGS)
  
  
  #########################################################################################################################
  ## Run join between species records and LGAs/SUAs :: Double check they are the same
  ## See the ABS for details :: there are 563 LGAs
  ## http://www.abs.gov.au/ausstats/abs@.nsf/Lookup/by%20Subject/1270.0.55.003~July%202016~Main%20Features~Local%20Government%20Areas%20(LGA)~7
  message('Joining occurence data to SUAs for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")
  
  
  projection(COMBO.RASTER.SP);projection(LGA.WGS);projection(SUA.WGS);projection(AUS.WGS)
  SUA.JOIN      = over(COMBO.RASTER.SP, SUA.WGS)              
  COMBO.SUA.LGA = cbind.data.frame(COMBO.RASTER.SP, SUA.JOIN) 
  
  
  #########################################################################################################################
  ## save .rds file for the next session
  saveRDS(COMBO.SUA.LGA, file = paste0('data/base/HIA_LIST/COMBO/COMBO_SUA_OVER_', save_run, '.rds'))
  
  
  #########################################################################################################################
  ## AGGREGATE THE NUMBER OF SUAs EACH SPECIES IS FOUND IN. NA LGAs ARE OUTSIDE AUS
  SUA.AGG   = tapply(COMBO.SUA.LGA$SUA_NAME16, COMBO.SUA.LGA$searchTaxon, function(x) length(unique(x))) ## group SUA by species name
  SUA.AGG   = as.data.frame(SUA.AGG)
  AUS.AGG   = aggregate(SUA_NAME16 ~ searchTaxon, data = COMBO.SUA.LGA, function(x) {sum(!is.na(x))}, na.action = NULL)
  
  SUA.AGG   = cbind.data.frame(AUS.AGG, SUA.AGG)
  names(SUA.AGG) = c("searchTaxon", "AUS_RECORDS", "SUA_COUNT")
  
  
  ## Now create a table of all the SUA's that each species occurrs
  SUA.SPP.COUNT = as.data.frame(table(COMBO.SUA.LGA[["SUA_NAME16"]], COMBO.SUA.LGA[["searchTaxon"]]))
  names(SUA.SPP.COUNT) = c("SUA", "SPECIES", "SUA_RECORDS")
  
  
  #########################################################################################################################
  ## save data
  ## Save .rds file for the next session
  saveRDS(SUA.SPP.COUNT, paste0('data/base/HIA_LIST/COMBO/SUA_SPP_COUNT_', save_run, '.rds'))
  
  
  
  ## Check : That's ok, but we want a table of which SUA each species is actually in.
  dim(SUA.AGG)
  head(SUA.AGG)
  
  
  ## 
  names(COMBO.SUA.LGA)
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
      
      ## Now use the niche width function on each colname (so 8 environmental variables)
      ## Also, need to figure out how to make the aggregating column generic (species, genus, etc.)
      ## currently it only works hard-wired
      niche_estimate (DF = NICHE.DF, colname = x)
      
      ## would be good to remove the duplicate columns here
      
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
  GLOBAL_RECORDS = as.data.frame(table(COMBO.SUA.LGA$searchTaxon))$Freq
  identical(length(GLOBAL_RECORDS), dim(COMBO.NICHE)[1])
  
  Total.taxa.processed = dim(COMBO.NICHE)[1]
  COMBO.NICHE  = cbind(GLOBAL_RECORDS, COMBO.NICHE)
  names(COMBO.NICHE)
  dim(COMBO.NICHE)
  
  
  #########################################################################################################################
  ## Add the counts of Australian records for each species to the niche database
  names(COMBO.NICHE)
  names(SUA.AGG)
  dim(COMBO.NICHE)
  dim(SUA.AGG)
  
  COMBO.LGA = join(COMBO.NICHE, SUA.AGG)                            ## The tapply needs to go where the niche summaries are
  names(COMBO.LGA)
  
  dim(COMBO.LGA)
  head(COMBO.LGA$AUS_RECORDS)
  head(COMBO.LGA$SUA_COUNT)
  
  ## Pick up tidying from here.............................................................................................
  
  
  
  
  
  #########################################################################################################################
  ## 8). JOIN ON CONTEXTUAL DATA
  #########################################################################################################################
  
  
  #########################################################################################################################
  ## Now join the horticultural contextual data onto one or both tables ()
  message('Joining contextual data for raster and niche files', length(GBIF.spp), ' species in the set ', "'", save_run, "'")
  COMBO.RASTER.CONTEXT = CLEAN.TRUE
  names(COMBO.RASTER.CONTEXT)
  
  
  ## Now join hort context to all the niche
  names(CLEAN.SPP)
  COMBO.NICHE.CONTEXT = join(COMBO.LGA, TOT.GROW)
  #COMBO.NICHE.CONTEXT =  COMBO.NICHE.CONTEXT[, c(2, 185, 1, 183:184, 186:197, 3:182)]
  head(COMBO.NICHE.CONTEXT$AUS_RECORDS)
  head(COMBO.NICHE.CONTEXT$LGA_COUNT)
  
  
  #########################################################################################################################
  ## Now combine the SDM output with the niche context data 
  ## Get the number of aus records too ....................................................................................
  
  #NICHE.CONTEXT   = TOT.GROW[, c("searchTaxon", "Plant.type", "Origin", "Total.growers", "Number.of.States")]
  TREE.PLANTINGS            = TREE.EVERGREEN[, c("searchTaxon",
                                                 "Plantings")]
  COMBO.NICHE.CONTEXT       = join(COMBO.LGA, TOT.GROW,  type = "left")  ## join does not support the sorting
  COMBO.NICHE.CONTEXT       = join(COMBO.NICHE.CONTEXT, TREE.PLANTINGS, type = "left")  
  COMBO.NICHE.CONTEXT       = COMBO.NICHE.CONTEXT [order(COMBO.NICHE.CONTEXT $searchTaxon),]
  
  
  ## Re-order table
  #COMBO.NICHE.CONTEXT =  COMBO.NICHE.CONTEXT[, c(2, 185, 1, 183:184, 186:197, 3:182)]
  
  ## Set NA to blank, then sort by no. of growers
  COMBO.NICHE.CONTEXT$Number.of.growers[is.na(COMBO.NICHE.CONTEXT$Total.growers)] <- 0
  COMBO.NICHE.CONTEXT = COMBO.NICHE.CONTEXT[with(COMBO.NICHE.CONTEXT, rev(order(Total.growers))), ]
  
  
  ## View the data
  dim(COMBO.RASTER.CONTEXT)
  dim(COMBO.NICHE.CONTEXT)
  View(COMBO.NICHE.CONTEXT)
  
  
  ## Print the dataframe dimensions to screen
  dim(CLEAN.TRUE)
  dim(COMBO.NICHE.CONTEXT)
  length(unique(CLEAN.TRUE$searchTaxon))
  length(COMBO.NICHE.CONTEXT$searchTaxon)
  
  
  #########################################################################################################################
  ## save .rds file for the next session
  saveRDS(COMBO.NICHE.CONTEXT,   paste0('data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT_',  save_run, '.rds'))
  saveRDS(COMBO.RASTER.CONTEXT,  paste0('data/base/HIA_LIST/COMBO/COMBO_RASTER_CONTEXT_', save_run, '.rds'))
  
  
} else {
  
  message(' skip file reading, not many species analysed')   ##
  
}


#########################################################################################################################
## OUTSTANDING CLEANING TASKS:
#########################################################################################################################


## Check the automatic cleaning results against manual cleaning results




#########################################################################################################################
#####################################################  TBC ##############################################################
#########################################################################################################################