#########################################################################################################################
#############################################  CREATE GBIF CLEANING COLUMNS ############################################# 
#########################################################################################################################


#########################################################################################################################
## 1). CHECK SEARCHED AND RETURNED NAMES FROM GBIF, AND DUPLICATES
#########################################################################################################################


## create a dummy file to run the records through
TEST = GBIF.HIA.SPP.RECORDS.ALL


## check how many records match the search term?
head(GBIF.HIA.SPP.RECORDS.ALL$scientificName, 50)
head(GBIF.HIA.SPP.RECORDS.ALL$searchTaxon, 50)


## remove the varieties, etc, from the scientific name returned by GBIF.
GBIF.HIA.SPP.RECORDS.ALL$Returned.binomial <- unlist(lapply(GBIF.HIA.SPP.RECORDS.ALL$scientificName, string_fun_first_two_words))


## Check for duplicate record ids in occurrence download
## Not sure which field indicate duplicated records. I think it is gbifID?
length(GBIF.HIA.SPP.RECORDS.ALL$gbifID) - length(unique(GBIF.HIA.SPP.RECORDS.ALL$gbifID))
(GBIF.HIA.SPP.RECORDS.ALL$gbifID %>% unique %>% length) == (GBIF.HIA.SPP.RECORDS.ALL %>% nrow)





#########################################################################################################################
## 2). pre-clean flags based on occurrence data attributes & assertions
#########################################################################################################################


## which columns overlap between GBIF and ALA/AVH? The problem here is that we don't what each GBIF field means, and the
## documentation is a bit dodgy.

## dataProvider (ALA), GBIF ?
## basisOfRecord (ALA), basisOfRecord (GBIF)
## occCultivatedEscapee (ALA), establishmentMeans (GBIF)
## year = year


## check for invalid species name
TEST %>%
  
  filter(!grepl("^[A-Z][a-z]* [a-z\\-]+( (subsp\\.|var\\.|f\\.) [a-z\\-]+)?$", scientificName)) %>%
  
  nrow  ## 


## now create flags for removing data
TEST %<>%
  
  ## all the columns names on the right hand side must match GBIF
  mutate(#CLEAN_REMOVE_NOT_AVH = !(dataProvider == "WORLDralia's Virtual Herbarium" & !is.na(dataProvider)),
         
         #CLEAN_REMOVE_NOT_IN_TAXON_LIST = !scientificName %in% taxa$taxonName,
         
         CLEAN_REMOVE_BASISOFRECORD = basisOfRecord == "HUMAN_OBSERVATION",
         
         CLEAN_REMOVE_MISSING_LONLAT = is.na(lon) | is.na(lat),
         
         #CLEAN_REMOVE_COORD = coordinateUncertaintyInMeters <= 1000,
         
         #CLEAN_REMOVE_OUTLIER = !is.na(outlierForLayer),
         
         #CLEAN_REMOVE_CULTIVATED_FLAG = occCultivatedEscapee & !is.na(occCultivatedEscapee),
         
         #CLEAN_REMOVE_CULTIVATED_LOCALITY = NA,
         
         CLEAN_REMOVE_PRE_1950 = year < 1950 & !is.na(year),
         
         FLAG_MISSING_YEAR = is.na(year),
         
         #ASSERT_TAXON_ID_ISSUE = !taxonIdentificationIssue %in% c("noIssue", "matchedToMisappliedName"),
         
         # ASSERT_COORD_ISSUE = ifelse(countryCoordinateMismatch | stateCoordinateMismatch | coordinatesCentreOfStateProvince == "TRUE" |
         #                             coordinatesCentreOfCountry == "TRUE" | negatedLongitude == "TRUE" | negatedLatitude == "TRUE" |
         #                             unrecognizedGeodeticDatum, TRUE, FALSE),
         
         #ASSERT_COORD_ISSUE = ifelse(is.na(ASSERT_COORD_ISSUE), FALSE, ASSERT_COORD_ISSUE),
         
         #ASSERT_DATE_FIRST_CENTURY = ifelse(!is.na(firstOfCentury) & firstOfCentury & year == 1900, TRUE, FALSE)
         
  )





#########################################################################################################################
## 3). Spatial outliers?
#########################################################################################################################


## can write the GBIF data to shapefile for mapping?
GBIF.ALL.POINTS = GBIF.HIA.SPP.RECORDS.ALL[c("lon", "lat")]
GBIF.ALL.POINTS = na.omit(GBIF.ALL.POINTS)

coordinates(GBIF.ALL.POINTS) = c("lon", "lat")
proj4string(GBIF.ALL.POINTS) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")



## clip to coastline
WORLD <- readOGR("./data/base/URBAN/TM_WORLD_BORDERS-0.3.shp", layer = "TM_WORLD_BORDERS-0.3")
plot(WORLD)
points(GBIF.HIA.SPP.RECORDS.ALL[c("lon", "lat")], cex = 0.1, col = "blue", pch = 19)


## now clip the GBIF points to those on land
proj4string(GBIF.ALL.POINTS)
proj4string(WORLD)
GBIF.LAND.POINTS <- GBIF.ALL.POINTS[WORLD, ]



