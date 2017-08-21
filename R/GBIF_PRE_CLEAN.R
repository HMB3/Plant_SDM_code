#########################################################################################################################
#############################################  CREATE GBIF CLEANING COLUMNS ############################################# 
#########################################################################################################################


#########################################################################################################################
## 1). CHECK SEARCHED AND RETURNED NAMES FROM GBIF, AND DUPLICATES
#########################################################################################################################


## create a dummy file to run the records through
GBIF.ALL = GBIF.HIA.SPP.RECORDS.ALL


## check how many records match the search term?
head(GBIF.ALL$scientificName, 50)
head(GBIF.ALL$searchTaxon, 50)


## remove the varieties, etc, from the scientific name returned by GBIF.
GBIF.ALL$Returned.binomial <- unlist(lapply(GBIF.ALL$scientificName, string_fun_first_two_words))


## Check for duplicate record ids in occurrence download
## Not sure which field indicate duplicated records. I think it is gbifID?
length(GBIF.ALL$gbifID) - length(unique(GBIF.ALL$gbifID))
(GBIF.ALL$gbifID %>% unique %>% length) == (GBIF.ALL %>% nrow)





#########################################################################################################################
## 2). pre-clean flags based on occurrence data attributes & assertions
#########################################################################################################################


## which columns overlap between GBIF and ALA/AVH? The problem here is that we don't what each GBIF field means, and the
## documentation is a bit dodgy.

## dataProvider (ALA), GBIF ?
## basisOfRecord (ALA), basisOfRecord (GBIF)
## occCultivatedEscapee (ALA), establishmentMeans (GBIF)
## year = year


#########################################################################################################################
## create flags for removing data
#########################################################################################################################


## check the data quality columns for the combined GBIF data set
## unique identifier?
head(GBIF.ALL$occurrenceID)
head(GBIF.ALL$gbifID)


##
unique(GBIF.ALL$basisOfRecord)
unique(GBIF.ALL$establishmentMeans)


## coordinates
summary(GBIF.ALL$lon)
summary(GBIF.ALL$lat)
summary(GBIF.ALL$year)


## check data
summary(GBIF.ALL$lon)[7]/dim(GBIF.ALL)[1]*100   ## 11% of data has no lat/long
summary(GBIF.ALL$year)[7]/dim(GBIF.ALL)[1]*100  ## 23% of data has no year


#########################################################################################################################
## Now create columns on which to exclude records
GBIF.ALL %<>%
  
  ## all the columns names on the right hand side must match GBIF
  mutate(
    
    #CLEAN_REMOVE_NOT_AVH = !(dataProvider == "WORLDralia's Virtual Herbarium" & !is.na(dataProvider)),
    
    #CLEAN_REMOVE_NOT_IN_TAXON_LIST = !scientificName %in% taxa$taxonName,
    
    CLEAN_REMOVE_BASISOFRECORD = basisOfRecord == "HUMAN_OBSERVATION",  ## TRUE for HUMAN OBSERVATION
    
    CLEAN_REMOVE_MISSING_LONLAT = is.na(lon) | is.na(lat),              ## TRUE = na(lon/lat)
    
    #CLEAN_REMOVE_COORD = coordinateUncertaintyInMeters <= 1000,
    
    CLEAN_REMOVE_CULTIVATED_FLAG = establishmentMeans == "MANAGED",
    
    #CLEAN_REMOVE_CULTIVATED_LOCALITY = NA,
    
    CLEAN_REMOVE_PRE_1950 = year < 1950 & !is.na(year),                 ## TRUE = na(year)
    
    FLAG_MISSING_YEAR = is.na(year)
    
    #ASSERT_TAXON_ID_ISSUE = !taxonIdentificationIssue %in% c("noIssue", "matchedToMisappliedName"),
    
    # ASSERT_COORD_ISSUE = ifelse(countryCoordinateMismatch | stateCoordinateMismatch | coordinatesCentreOfStateProvince == "TRUE" |
    #                             coordinatesCentreOfCountry == "TRUE" | negatedLongitude == "TRUE" | negatedLatitude == "TRUE" |
    #                             unrecognizedGeodeticDatum, TRUE, FALSE),
    
    #ASSERT_COORD_ISSUE = ifelse(is.na(ASSERT_COORD_ISSUE), FALSE, ASSERT_COORD_ISSUE),
    
    #ASSERT_DATE_FIRST_CENTURY = ifelse(!is.na(firstOfCentury) & firstOfCentury & year == 1900, TRUE, FALSE)
    
  )


#########################################################################################################################
## Check what flags are doing?
table(GBIF.ALL$establishmentMeans)[3]
summary(GBIF.ALL$CLEAN_REMOVE_CULTIVATED_FLAG)[3]   ## TRUE = HUMAN OBSERVATION

table(GBIF.ALL$basisOfRecord)[2]
summary(GBIF.ALL$CLEAN_REMOVE_BASISOFRECORD)[3]     ## TRUE = MANAGED

summary(GBIF.ALL$lat)[7]                            
summary(GBIF.ALL$CLEAN_REMOVE_MISSING_LONLAT)[3]    ## TRUE = na(lon/lat)

summary(GBIF.ALL$year)[7]                           
summary(GBIF.ALL$FLAG_MISSING_YEAR)[3]              ## TRUE = na(year)

dim(GBIF.ALL[which(GBIF.ALL$year >= 1950), ])[1]                           
summary(GBIF.ALL$CLEAN_REMOVE_PRE_1950)[3]          ## TRUE = year > 1950 and not NA







#########################################################################################################################
## 3). Remove spatial outliers
#########################################################################################################################


## can write the GBIF data to shapefile for mapping?
GBIF.ALL.POINTS = GBIF.ALL[c("lon", "lat")]
GBIF.ALL.POINTS = na.omit(GBIF.ALL.POINTS)
coordinates(GBIF.ALL.POINTS) = c("lon", "lat")
proj4string(GBIF.ALL.POINTS) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")


## clip to coastline
WORLD <- readOGR("./data/base/URBAN/TM_WORLD_BORDERS-0.3.shp", layer = "TM_WORLD_BORDERS-0.3")
GBIF.LAND.POINTS <- GBIF.ALL.POINTS[WORLD, ]


## hard to tell if the points in the ocean are on islands?
plot(WORLD)
points(GBIF.LAND.POINTS, cex = 0.05, col = "blue", pch = 19)





#########################################################################################################################
## 4). CLEAN RECORDS
#########################################################################################################################


## Apply the filters created above
GBIF.CLEAN = GBIF.ALL


GBIF.CLEAN %<>%
  
  ## can't apply the filters that only apply to ALA data
  filter(#!CLEAN_REMOVE_NOT_IN_TAXON_LIST,
         #!CLEAN_REMOVE_NOT_OEH_AVH,
         !CLEAN_REMOVE_BASISOFRECORD,
         !CLEAN_REMOVE_MISSING_LONLAT,
         !CLEAN_REMOVE_CULTIVATED_FLAG,
         !CLEAN_REMOVE_CULTIVATED_LOCALITY,
         !CLEAN_REMOVE_PRE_1950,
         !FLAG_MISSING_YEAR,
         #!ASSERT_TAXON_ID_ISSUE,
         #!ASSERT_COORD_ISSUE,
         #!ASSERT_DATE_FIRST_CENTURY,
         #!is.na(clip_0.05)
         )


## How many records are knocked out by each filter?
GBIF.CLEAN %>% 
  
  ## 
  filter(CLEAN_REMOVE_BASISOFRECORD) %>% 
  
  ## how many rows?
  nrow


