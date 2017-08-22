#########################################################################################################################
#############################################  CREATE GBIF CLEANING COLUMNS ############################################# 
#########################################################################################################################


#########################################################################################################################
## 1). pre-clean flags based on occurrence data attributes & assertions
#########################################################################################################################


## Which columns overlap between GBIF and ALA/AVH? The problem here is that we don't what each GBIF field means, and the
## documentation is a bit dodgy.

## dataProvider             (ALA), GBIF ?
## basisOfRecord            (ALA), basisOfRecord (GBIF)
## occCultivatedEscapee     (ALA), establishmentMeans (GBIF)
## year
## taxonIdentificationIssue (ALA), GBIF no equivalent
## fatal assertions         (ALA), GBIF no equivalent


#########################################################################################################################
## create flags for removing data
#########################################################################################################################


# ## check the data quality columns for the combined GBIF data set
# ## unique identifier?
# head(GBIF.TRIM$occurrenceID)
# head(GBIF.TRIM$gbifID)
# 
# 
# ##
# unique(GBIF.TRIM$basisOfRecord)
# unique(GBIF.TRIM$establishmentMeans)
# 
# 
# ## coordinates
# summary(GBIF.TRIM$lon)
# summary(GBIF.TRIM$lat)
# summary(GBIF.TRIM$year)
# 
# 
# ## check data
# summary(GBIF.TRIM$lon)[7]/dim(GBIF.TRIM)[1]*100   ## 11% of data has no lat/long
# summary(GBIF.TRIM$year)[7]/dim(GBIF.TRIM)[1]*100  ## 23% of data has no year


# #########################################################################################################################
# ## Now create columns on which to exclude records
# GBIF.TRIM %<>%
#   
#   ## all the columns names on the right hand side must match GBIF
#   mutate(
#     
#     #CLEAN_REMOVE_NOT_AVH = !(dataProvider == "WORLDralia's Virtual Herbarium" & !is.na(dataProvider)),
#     
#     #CLEAN_REMOVE_NOT_IN_TAXON_LIST = !scientificName %in% taxa$taxonName,
#     
#     CLEAN_REMOVE_BASISOFRECORD = basisOfRecord == "HUMAN_OBSERVATION",  ## TRUE for HUMAN OBSERVATION
#     
#     CLEAN_REMOVE_MISSING_LONLAT = is.na(lon) | is.na(lat),              ## TRUE = na(lon/lat)
#     
#     #CLEAN_REMOVE_COORD = coordinateUncertaintyInMeters <= 1000,
#     
#     CLEAN_REMOVE_CULTIVATED_FLAG = establishmentMeans == "MANAGED",
#     
#     #CLEAN_REMOVE_CULTIVATED_LOCALITY = NA,
#     
#     CLEAN_REMOVE_PRE_1950 = year < 1950 & !is.na(year),                 ## TRUE = na(year)
#     
#     FLAG_MISSING_YEAR = is.na(year)
#     
#     #ASSERT_TAXON_ID_ISSUE = !taxonIdentificationIssue %in% c("noIssue", "matchedToMisappliedName"),
#     
#     # ASSERT_COORD_ISSUE = ifelse(countryCoordinateMismatch | stateCoordinateMismatch | coordinatesCentreOfStateProvince == "TRUE" |
#     #                             coordinatesCentreOfCountry == "TRUE" | negatedLongitude == "TRUE" | negatedLatitude == "TRUE" |
#     #                             unrecognizedGeodeticDatum, TRUE, FALSE),
#     
#     #ASSERT_COORD_ISSUE = ifelse(is.na(ASSERT_COORD_ISSUE), FALSE, ASSERT_COORD_ISSUE),
#     
#     #ASSERT_DATE_FIRST_CENTURY = ifelse(!is.na(firstOfCentury) & firstOfCentury & year == 1900, TRUE, FALSE)
#     
#   )





#########################################################################################################################
## 2). CLEAN RECORDS
#########################################################################################################################


## create copy and free memory
GBIF.CLEAN = GBIF.TRIM
gc()


# GBIF.CLEAN %<>%
#   
#   ## can't apply the filters that only apply to ALA data
#   filter(!CLEAN_REMOVE_MISSING_LONLAT)


## create a table which counts the number of records meeting a comination of all the criteria
## in this table, 
GBIF.PROBLEMS <- with(GBIF.TRIM,
                 
                 table(
                   
                   is.na(lon)|is.na(lat),
                   
                   establishmentMeans == 'MANAGED' & !is.na(establishmentMeans),
                   
                   year < 1950 & !is.na(year),
                   
                   is.na(year),
                   
                   coordinateUncertaintyInMeters > 1000 & 
                     !is.na(coordinateUncertaintyInMeters)
                   
                 )
                 
) %>% 
  
  as.data.frame %>%  
  
  setNames(c('NO_COORD', 'MANAGED', 
             'PRE_1950', 'NO_YEAR', 'COORD_UNCERT', 'COUNT'))


## write file to CSV
write.csv(problems, 'problems.csv', row.names = FALSE)


## clean 
GBIF.CLEAN <- GBIF.HIA.SPP.RECORDS.ALL %>% 
  
  select(one_of(
    
    c('basisOfRecord', 'lon', 'lat', 'establishmentMeans', 
      'year', 'scientificName', 'searchTaxon', 'gbifID', 
      'coordinateUncertaintyInMeters'))) %>% 
  
  filter(!is.na(lon) & !is.na(lat),
         establishmentMeans!='MANAGED' | is.na(establishmentMeans),
         year >= 1950 & !is.na(year))





#########################################################################################################################
## 3). Remove spatial outliers
#########################################################################################################################


world.temp = raster("//sci-7910/F/data/worldclim/world/0.5/bio/current/bio_01")

xy <- cellFromXY(world.temp, GBIF.CLEAN[c("lon", "lat")]) %>% 
  unique %>% 
  #length
  xyFromCell(world.temp, .)

#onland <- cells[!is.na(world.temp[cells])]  ## placeholder for what's been passed from previous argument
# vals <- gdalUtils::gdallocationinfo(
#   "//sci-7910/F/data/worldclim/world/0.5/bio/current/bio_01",
#   coords=xy,
#   valonly=TRUE, raw_output=FALSE, geoloc=T, wgs84=T)
# 


##
onland <- extract(world.temp, xy) %>% !is.na %>% xy[.,]
onland = z %>% is.na %>%  `!` # %>% xy[.,]  cells on land or not

onland.points = filter(GBIF.CLEAN, cellFromXY(world.temp, GBIF.CLEAN[c("lon", "lat")]) %in% 
                         unique(cellFromXY(world.temp, GBIF.CLEAN[c("lon", "lat")]))[onland]) 
  

##
points(onland.points[c("lon", "lat")], pch = ".")


## clip to coastline
WORLD <- readOGR("./data/base/URBAN/TM_WORLD_BORDERS-0.3.shp", layer = "TM_WORLD_BORDERS-0.3")
GBIF.LAND.POINTS <- GBIF.TRIM.POINTS[WORLD, ]


## hard to tell if the points in the ocean are on islands?
plot(WORLD)
points(GBIF.LAND.POINTS, cex = 0.05, col = "blue", pch = 19)


