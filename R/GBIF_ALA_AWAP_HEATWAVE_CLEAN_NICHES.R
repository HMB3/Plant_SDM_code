#########################################################################################################################
################################################# FLAG SPATIAL OUTLIERS #################################################
#########################################################################################################################


## This code uses the CoordinateCleaner package to Automatcially screen out the dodgy spatial records. See ::
## https://github.com/azizka/CoordinateCleaner


## It then creates a table of the species niche :

## 2). A table with one row for each species, including contextual data and species attributes (niches, traits, etc.)

## These tables are subsequently used to estimate the current global realised niche/climatic tolerance using the best 
## available data, and susequently model the niches using the maxent algorithm.


#########################################################################################################################
## Read in the three data tables
if(read_data == "TRUE") {
  
  ## Load GBIF and ALA data
  COMBO.AWAP.CONVERT = readRDS( paste0(DATA_path, 'COMBO_RASTER_CONVERT_',  save_run, '.rds'))
  
  
} else {
  
  message(' skip file reading, not many species analysed')   ##
  
}

## Check dimensions of the occurrence and inventory data tables.
length(unique(COMBO.AWAP.CONVERT$searchTaxon))
formatC(dim(COMBO.AWAP.CONVERT)[1], format = "e", digits = 2)
names(COMBO.AWAP.CONVERT)


#########################################################################################################################
## Create a unique identifier
COMBO.AWAP.CONVERT$CC.OBS <- 1:nrow(COMBO.AWAP.CONVERT)





#########################################################################################################################
## 1). FLAG GBIF OUTLIERS
#########################################################################################################################


#########################################################################################################################
## Rename the columns to fit the CleanCoordinates format and create a tibble. 
## A Tibble is needed for running the spatial outlier cleaning
TIB.GBIF = COMBO.AWAP.CONVERT[, !duplicated(names(COMBO.AWAP.CONVERT))]
TIB.GBIF <- TIB.GBIF %>% dplyr::rename(species          = searchTaxon,
                                       decimallongitude = lon, 
                                       decimallatitude  = lat) %>%
  timetk::tk_tbl()


## We've already stripped out the records that fall outside
## the worldclim raster boundaries, so the sea test is probably not the most important
## Study area is the globe, but we are only projecting models onto Australia


#########################################################################################################################
## Don't run the outliers test here, it is slower. Also, can't run cleaning on the urban tree inventory data, because this
## removes all the records near capital cities 
FLAGS  <- CleanCoordinates(TIB.GBIF,
                           #countries        = "country",    ## too many flagged here...
                           capitals.rad     = 0.12,
                           countrycheck     = TRUE,
                           duplicates       = TRUE,
                           seas             = FALSE)


## save/load the flags
identical(dim(FLAGS)[1], dim(TIB.GBIF)[1])
#saveRDS(FLAGS, paste0('data/base/HIA_LIST/COMBO/ALA_GBIF_FLAGS_', save_run, '.rds'))
#FLAGS = readRDS('data/base/HIA_LIST/COMBO/ALA_GBIF_FLAGS_OLD_ALA.rds')


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
COMBO.LUT = as.data.frame(table(TIB.GBIF$species))
names(COMBO.LUT) = c("species", "FREQUENCY")
COMBO.LUT = COMBO.LUT[with(COMBO.LUT, rev(order(FREQUENCY))), ] 
head(COMBO.LUT);summary(COMBO.LUT$FREQUENCY)  ## Quercus robur, 214, and Fraxinus excelsior, 156
LUT.100K = as.character(subset(COMBO.LUT, FREQUENCY < 100000)$species)
LUT.100K = LUT.100K [order(LUT.100K)]


## Unfortunately, the cc_outl function can't handle vectors of a certain size - over 40 GB at least.
## So we are problably better off moving this code to after the data has been thinned by the SDM
## step 


## Create a data frame of species name and spatial outlier
SPAT.OUT <- unique(TIB.GBIF$species) %>%  ## LUT.100K
  
  ## pipe the list of species into lapply
  lapply(function(x) {
    
    ## Create the species df by subsetting by species
    f <- subset(TIB.GBIF, species == x)
    
    ## Run the spatial outlier detection
    message("Running spatial outlier detection for ", x)
    message(dim(f)[1], " records for ", x)
    sp.flag <- cc_outl(f,
                       lon     = "decimallongitude",
                       lat     = "decimallatitude",
                       species = "species",
                       method  = "distance",
                       #mltpl   = 5,
                       tdi     = 300,
                       value   = "flags")
    
    ## Now add attache column for species, and the flag for each record
    #d = data.frame(searchTaxon = x, SPAT_OUT = sp.flag)
    d = cbind(species = x,
              SPAT_OUT = sp.flag, f)[c("species", "SPAT_OUT", "CC.OBS")]
    
  }) %>%
  
  ## Finally, bind all the rows together
  bind_rows


## Check spatial data
identical(dim(SPAT.OUT)[1], dim(TIB.GBIF)[1])
unique(SPAT.OUT$species)
summary(SPAT.OUT$SPAT_OUT)
head(SPAT.OUT)




#########################################################################################################################
## 3). FILTER DATA TO THE CLEAN RECORDS
#########################################################################################################################


#########################################################################################################################
## Join data :: exclude the decimal lat/long, check the length 
identical(dim(TIB.GBIF)[1],dim(FLAGS)[1])#;length(GBIF.SPAT.OUT)
names(FLAGS)[1]    = c("coord_spp")
identical(COMBO.AWAP.CONVERT$searchTaxon, FLAGS$coord_spp)                                            ## order matches


#########################################################################################################################
## Is the species column the same as the searchTaxon column?
TEST.GEO   = cbind(COMBO.AWAP.CONVERT, FLAGS)
identical(TEST.GEO$searchTaxon, TEST.GEO$coord_spp)                                                     ## order matches


## Keep records which passed the GBIF and spatial test
dim(subset(TEST.GEO, summary == "TRUE")) 
CLEAN.TRUE = subset(TEST.GEO, summary == "TRUE")
unique(CLEAN.TRUE$summary)

                                    
## What percentage of records are retained?
identical(dim(CLEAN.TRUE)[1], (dim(COMBO.AWAP.CONVERT)[1] - dim(subset(TEST.GEO, summary == "FALSE"))[1]))
length(unique(CLEAN.TRUE$searchTaxon))
message(round(dim(CLEAN.TRUE)[1]/dim(TEST.GEO)[1]*100, 2), " % records retained")                                               





#########################################################################################################################
## 4). CREATE SHAPEFILES TO CHECK OUTLIERS ARCMAP
#########################################################################################################################


## Select columns
GBIF.ALA.CHECK  = select(CLEAN.TRUE, searchTaxon, scientificName, lat, lon, SOURCE, year,
                         country, locality, basisOfRecord, institutionCode)


## Rename the fields so that ArcMap can handle them
GBIF.ALA.CHECK     = dplyr::rename(GBIF.ALA.CHECK,
                                   TAXON     = searchTaxon,
                                   LAT       = lat,
                                   LON       = lon,
                                   SC_NAME   = scientificName,
                                   BASIS     = basisOfRecord,
                                   LOCAL     = locality,
                                   INSTIT    = institutionCode,
                                   COUNTRY   = country,
                                   YEAR      = year)
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
#   writeOGR(tmp, dsn = "./data/ANALYSIS/CLEAN_GBIF", TAXA[i], driver = "ESRI Shapefile", overwrite_layer = TRUE)
# 
# }


## Write the whole shapefile out just in case
writeOGR(obj = GBIF.ALA.SPDF, dsn = "./data/ANALYSIS/CLEAN_GBIF", 
         layer = paste0("CLEAN_SPDF_", save_run), driver = "ESRI Shapefile")





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
  
  
  #########################################################################################################################
  ## save .rds file for the next session
  saveRDS(COMBO.SUA.LGA, file = paste0(DATA_path, 'COMBO_SUA_OVER_', save_run, '.rds'))
  
  
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
  
  
  #########################################################################################################################
  ## Save .rds file for the next session
  saveRDS(SUA.SPP.COUNT, paste0(DATA_path, 'SUA_SPP_COUNT_', save_run, '.rds'))
  
  
  ## Check : That's ok, but we want a table of which SUA each species is actually in.
  dim(SUA.AGG)
  head(SUA.AGG)
  
  
  ## Remove coordinates
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
  ## Using only complete cases for the Ukola and Fitzpatrick data will make the dataset smaller
  ## Consider adding NA condition to the niche_estimate function...........................................................
  NICHE.DF = completeFun(COMBO.SUA.LGA, c("PET", "Drought_freq_extr", "Drought_mean_int_extr", "HWM", "HW_CUM_HOT"))
  dim(NICHE.DF)
  head(niche_estimate (DF = NICHE.DF, colname = "Annual_mean_temp"))  ## including the q05 and q95
  
  
  ## So lets use lapply on the "SearchTaxon"
  ## test = run_function_concatenate(list, DF, "DF, colname = x") 
  message('Estimating global niches for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")
  COMBO.NICHE <- awap.variables %>% 
    
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
                                                searchTaxon.16, searchTaxon.17, searchTaxon.18, searchTaxon.19,
                                                searchTaxon.20,  searchTaxon.21,  searchTaxon.22,  searchTaxon.23,
                                                searchTaxon.24,  searchTaxon.25,  searchTaxon.26,  searchTaxon.27,
                                                searchTaxon.28,  searchTaxon.29,  searchTaxon.30, searchTaxon.31,
                                                searchTaxon.32))
  
  
  #########################################################################################################################
  ## Add counts for each species, and record the total number of taxa processed
  ## dim(COMBO.RASTER.CONVERT);dim(CLEAN.TRUE)
  GLOBAL_RECORDS = as.data.frame(table(COMBO.SUA.LGA$searchTaxon))
  names(GLOBAL_RECORDS) = c("searchTaxon", "GLOBAL_RECORDS")
  identical(dim(GLOBAL_RECORDS)[1], dim(COMBO.NICHE)[1])
  
  Total.taxa.processed = dim(COMBO.NICHE)[1]
  COMBO.NICHE  = join(GLOBAL_RECORDS, COMBO.NICHE)
  names(COMBO.NICHE)
  dim(COMBO.NICHE)
  
  
  #########################################################################################################################
  ## Add the counts of Australian records for each species to the niche database
  names(COMBO.NICHE)
  names(SUA.AGG)
  dim(COMBO.NICHE)
  dim(SUA.AGG)
  
  COMBO.LGA = join(COMBO.NICHE, SUA.AGG)
  names(COMBO.LGA)
  
  dim(COMBO.LGA)
  head(COMBO.LGA$AUS_RECORDS)
  head(COMBO.LGA$SUA_COUNT)
  
  
  #########################################################################################################################
  ## 7). CALCULATE AREA OF OCCUPANCY
  #########################################################################################################################
  
  
  ## Given a dataframe of georeferenced occurrences of one, or more, taxa, this function provide statistics values 
  ## (Extent of Occurrence, Area of Occupancy, number of locations, number of subpopulations) and provide a preliminary 
  ## conservation status following Criterion B of IUCN. A graphical map output is also available.
  
  ## Create a species list to estimate the ranges for. Currently we can't estimate ranges for widely distributed species
  ## A workaround for the trait species is to just use the ALA records. But for the full list of species, we will need to 
  ## change this. Otherwise, we could report the global niche, but the Australian range only. This will mean that 
  spp.geo = as.character(unique(COMBO.SUA.LGA$searchTaxon))
  AOO.DAT = COMBO.SUA.LGA %>%
    subset(., SOURCE == "ALA")
  
  
  #########################################################################################################################
  ## AREA OF OCCUPANCY (AOO)
  ## For every species in the list: calculate the AOO
  ## The extent of occurrence function is causing problems.................................................................
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
  COMBO.NICHE       = join(GBIF.AOO, COMBO.NICHE)
  COMBO.NICHE.TRAIT = join(TCRIT, COMBO.NICHE)

  
  ## AOO is calculated as the area of all known or predicted cells for the species. The resolution will be 10x10km as 
  ## required by IUCN. A single value in km2.


  #########################################################################################################################
  ## Now join hort context to all the niche
  # TCRIT = select(DIANA.TC, Names, Tcrit, ster)
  # names(TCRIT)[names(TCRIT) == 'Names'] <- 'searchTaxon'


  #########################################################################################################################
  ## save .rds file for the next session
  message('Writing niche and raster data for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")
  saveRDS(COMBO.NICHE,    paste0(DATA_path, 'COMBO_NICHE_CONTEXT_',  OCC_SOURCE, '_RECORDS_', save_run, '.rds'))
  write.csv(COMBO.NICHE.TRAIT,  paste0(DATA_path, 'COMBO_NICHE_CONTEXT_',  OCC_SOURCE, '_RECORDS_', save_run, '.csv'), row.names = FALSE)
  
} else {
  
  message(' skip file saving, ', length(GBIF.spp), ' species analysed')   ##
  
}




#########################################################################################################################
## OUTSTANDING CLEANING TASKS:
#########################################################################################################################


## Check the automatic cleaning results against manual cleaning results




#########################################################################################################################
#####################################################  TBC ##############################################################
#########################################################################################################################