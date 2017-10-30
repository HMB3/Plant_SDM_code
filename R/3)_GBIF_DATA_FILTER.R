#########################################################################################################################
####################################################  FILTER GBIF DATA ################################################## 
#########################################################################################################################


## Aim of this code is to filer the GBIF data to the most reliable recrods. Two sources of uncertainty:

## Spatial 
## Taxonomic


## Which columns overlap between GBIF and ALA/AVH? The problem here is that we don't know what each GBIF field means,

## dataProvider             (ALA), GBIF ?
## basisOfRecord            (ALA), basisOfRecord (GBIF)
## occCultivatedEscapee     (ALA), establishmentMeans (GBIF)
## year
## taxonIdentificationIssue (ALA), GBIF no equivalent
## fatal assertions         (ALA), GBIF no equivalent


## A record with the same lat/long to x decimal places, collected in the same month by the same person is a duplicate


#########################################################################################################################
## LOAD
load("./data/base/HIA_LIST/COMBO/GBIF_TRIM.RData")
dim(GBIF.TRIM)
names(GBIF.TRIM)
length(unique(GBIF.TRIM$searchTaxon))  ## has the list updated with extra species? YES!





#########################################################################################################################
## 1). MARK CULTIVATED RECORDS AND CHECK TAXONOMY
#########################################################################################################################


## Try and make this one command: can these terms be searched for across the whole data frame (i.e. any column)
GBIF.TRIM$CULTIVATED <- ifelse(grepl("garden|cultiva",   GBIF.TRIM$locality,           ignore.case = TRUE) | 
                                 grepl("garden|cultiva", GBIF.TRIM$habitat,            ignore.case = TRUE) | 
                                 grepl("garden|cultiva", GBIF.TRIM$eventRemarks,       ignore.case = TRUE) |
                                 grepl("garden|cultiva", GBIF.TRIM$cloc,               ignore.case = TRUE) |
                                 grepl("managed",        GBIF.TRIM$establishmentMeans, ignore.case = TRUE),
                               
                               "CULTIVATED", "UNKNOWN")


## How many records are knocked out by using this definition?
## This is probably a bit strict, in that for some of the fields, garden doesn't = cultivated
GBIF.CULTIVATED = subset(GBIF.TRIM, CULTIVATED == "CULTIVATED")
dim(GBIF.CULTIVATED)[1]/dim(GBIF.TRIM)[1]


#########################################################################################################################
## CHECK TAXONOMY THAT GBIF RETURNED 
## Use "Taxonstand"
GBIF.TAXO <- TPL(unique(GBIF.TRIM$scientificName), infra = TRUE,
                 corr = TRUE)
sort(names(GBIF.TAXO))


## Then join the GBIF data to the taxonomic check, using "scientificName" as the join field... 
GBIF.TRIM <- GBIF.TRIM %>%
  left_join(., GBIF.TAXO, by = c("scientificName" = "Taxon"))


## So we can filter by the agreement between "scientificName", and "New.Taxonomic.status"?
## Not sure if this completely checks the searched and returned species match, but maybe good enough...
unique(GBIF.TAXO.CHECK$New.Taxonomic.status)


## Also keep the managed records:
GBIF.UNRESOLVED <- GBIF.TAXO.CHECK %>% 
  
  ## Note that these filters are very forgiving...
  ## Unless we include the NAs, very few records are returned!
  filter(New.Taxonomic.status == 'Unresolved')


## Also keep the managed records:
# GBIF.TAXO.CHECK <- GBIF.TAXO.CHECK %>% 
#   
#   ## Note that these filters are very forgiving...
#   ## Unless we include the NAs, very few records are returned!
#   filter(New.Taxonomic.status == 'Accepted')


## Unique(GBIF.UNRESOLVED$New.Taxonomic.status)
save(GBIF.UNRESOLVED, file = paste("./data/base/HIA_LIST/GBIF/GBIF_UNRESOLVED.RData"))





#########################################################################################################################
## 2). CREATE TABLE OF PRE-CLEAN FLAGS 
#########################################################################################################################


#########################################################################################################################
## Then create a table which counts the number of records meeting each criteria:
## Note that TRUE indicates there is a problem (e.g. if a record has no lat/long, it will = TRUE)
GBIF.PROBLEMS <- with(GBIF.TRIM,
                 
                 table(
                   
                   ## Note this list could be exanded for other data types
                   ## ALA/AVH, council data, etc.
                   
                   ## No coordinates
                   is.na(lon)|is.na(lat),
                   
                   ## Also add geoclean
                   #GEOCLEAN   == 'FALSE',
                   
                   ## Also add duplicated
                   #DUPLICATED == 'TRUE',
                   
                   ## Taxon rank is genus/form?
                   #taxonRank == 'GENUS' & 'FORM',
                   
                   ## Taxonomic status
                   New.Taxonomic.status == 'Unresolved',
                   
                   ## Cultivated
                   CULTIVATED == 'CULTIVATED',
                   
                   ## Establishment means is "MANAGED", can change this:
                   establishmentMeans == 'MANAGED' & !is.na(establishmentMeans),
                   
                   ## Collected before 1950 or na.year
                   year < 1950 & !is.na(year),
                   
                   ## No year
                   is.na(year),
                   
                   ## Coordinate uncertainty is > 100 or is not NA
                   coordinateUncertaintyInMeters > 1000 & 
                     !is.na(coordinateUncertaintyInMeters)
                   
                   ## Add maybe the centre of Australia?
                   ## Lamber centre of Aus: 25.610111, 134.354806
                   
                 )
                 
) %>% 
  
  ## Create a data frame and set the names
  as.data.frame %>%  
  
  setNames(c('NO_COORD', 'TAXON_STATUS',
             'CULTIVATED', #'GEOCLEAN', 'DUPLICATED', 
             'MANAGED',  'PRE_1950', 'NO_YEAR', 
             'COORD_UNCERT', 'COUNT'))


## Print table to screen: in .Rmd file, no need to save
## Note that TRUE indicates there is a problem (e.g. no lat/long = TRUE)
kable(GBIF.PROBLEMS)


## Quickly check the total record number matches the count of problems
Total.count = sum(GBIF.PROBLEMS$COUNT)
identical(dim(GBIF.TRIM)[1], Total.count)  ## identical matches two objects


## Probably don't need this
save(GBIF.PROBLEMS,  file = paste("./data/base/HIA_LIST/GBIF/GBIF_PROBLEMS.RData"))


## Also keep the managed records:
GBIF.MANAGED <- GBIF.TRIM %>% 
  
  ## Note that these filters are very forgiving...
  ## Unless we include the NAs, very few records are returned!
  filter(establishmentMeans == 'MANAGED')


## Unique(GBIF.MANAGED$establishmentMeans)
save(GBIF.MANAGED, file = paste("./data/base/HIA_LIST/GBIF/GBIF_MANAGED.RData"))


## Add field for managed/unmanaged 
# COMBO.APNI = merge(COMBO.NICHE.CONTEXT, APNI, by = "searchTaxon", all.x = TRUE) 
# COMBO.APNI$APNI[is.na(COMBO.APNI$APNI)] <- "FALSE"
# unique(COMBO.APNI$APNI)





#########################################################################################################################
## 3). FILTER RECORDS 
#########################################################################################################################


## Filter the GBIF records using conditions which are not too restrictive
GBIF.CLEAN <- GBIF.TRIM %>% 
  
  ## Note that these filters are very forgiving...
  ## Unless we include the NAs, very few records are returned!
  filter(!is.na(lon) & !is.na(lat),
         establishmentMeans!='MANAGED' | is.na(establishmentMeans),
         year >= 1950 & !is.na(year))
         # DUPLICATED != 'TRUE',
         # GEOCLEAN   != 'FALSE')


## The table above gives the details, but worth documenting how many records are knocked out by each filter
Remaining.records = dim(GBIF.CLEAN)[1] 
Remaining.percent = dim(GBIF.CLEAN)[1]/Total.count*100
Filters.applied = "NA COORD | MANAGED/NA | < 1950/NA"
Remaining.percent ## 57% of records remain after cleaning 
gc()


## Check
dim(GBIF.CLEAN)
head(GBIF.CLEAN)


## What does this dataframe look like?
names(GBIF.CLEAN)





#########################################################################################################################
## 4). REMOVE POINTS OUTSIDE WORLDCLIM LAYERS...
#########################################################################################################################


## Can use WORLDCIM rasters to get only records where wordlclim data is. 
## At the global scale, there probably is no alterntive to using WORLDCLIM...


## First get one of the BIOCLIM variables
world.temp = raster("//sci-7910/F/data/worldclim/world/0.5/bio/current/bio_01")
plot(world.temp)


## Now get the XY centroids of the unique 1km * 1km WORLDCLIM blocks where GBIF records are found
## Get cell number(s) of WORLDCLIM raster from row and/or column numbers. Cell numbers start at 1 in the upper left corner, 
## and increase from left to right, and then from top to bottom. The last cell number equals the number of raster cells 
xy <- cellFromXY(world.temp, GBIF.CLEAN[c("lon", "lat")]) %>% 
  
  ## get the unique raster cells
  unique %>% 

  ## Get coordinates of the center of raster cells for a row, column, or cell number of WORLDCLIM raster
  xyFromCell(world.temp, .)


## Take a look at xy: NA's should be removed...
summary(xy)
str(xy)
points(xy, pch = ".", col = "red")


## For some reason, we need to convert the xy coords to a spatial points data frame, in order to avoid this error:
## 'NAs introduced by coercion to integer range'
xy <- SpatialPointsDataFrame(coords = xy, data = as.data.frame(xy),
                             proj4string = CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))


## Now extract the temperature values for the unique 1km centroids which contain GBIF data
class(xy)
z   = extract(world.temp, xy)

# Warning message:
#   In .doExtract(x, i, ..., drop = drop) :
#   some indices are invalid (NA returned)

hist(z, border = NA, col = "orange", breaks = 50, main = "", xlab = "Worldclim Annual temp")


## Then track which values of Z are on land or not
onland = z %>% is.na %>%  `!` # %>% xy[.,]  cells on land or not
summary(onland)


## Finally, filter the cleaned GBIF data to only those points on land. 
## This is achieved with the final [onland]
GBIF.LAND = filter(GBIF.CLEAN, cellFromXY(world.temp, GBIF.CLEAN[c("lon", "lat")]) %in% 
                     unique(cellFromXY(world.temp, GBIF.CLEAN[c("lon", "lat")]))[onland])


## how many records were on land?
records.ocean = dim(GBIF.CLEAN)[1] - dim(GBIF.LAND)[1]  ## 91575 records are in the ocean   
records.ocean


#########################################################################################################################
## Plot cleaned data
#plot(world.temp)
#points(GBIF.CLEAN[c("lon", "lat")], pch = ".", col = "red")


## Plot cleaned data that's in the worldclim raster
#plot(world.temp)
#points(GBIF.LAND[c("lon", "lat")], pch = ".", col = "blue")
gc()


## note that this still leaves lots of points in the Islands. So we need to decide if those are legitimate.
WORLD <- readOGR("./data/base/CONTEXTUAL/TM_WORLD_BORDERS-0.3.shp", layer = "TM_WORLD_BORDERS-0.3")
LAND  <- readOGR("./data/base/CONTEXTUAL/ne_10m_land.shp",          layer = "ne_10m_land")
plot(WORLD)
plot(LAND)



#########################################################################################################################
## save data
save(GBIF.LAND,        file = paste("./data/base/HIA_LIST/GBIF/GBIF_LAND_POINTS.RData"))
save(GBIF.TAXO.ERRORS, file = paste("./data/base/HIA_LIST/GBIF/GBIF_TAXO_ERRORS.RData", sep = ""))
gc()


## Now save .RData file for the next session
save.image("STEP_3_GBIF_CLEAN.RData")





#########################################################################################################################
## OUTSTANDING CLEANING TASKS:
#########################################################################################################################


## When should the additional filters be run in? Just after GBIF.TRIM?            - Ask John

## GBIF taxonomic errors                                                          - Use taxonstand, check with John, Dave K.

## Keep cultivated records as a separate column/file                              - Have them: check search with Rach, Linda

## Estimate native ranges: as a separate colum too?                               - Ubran polygons for AUS, USA, EU? Ask Dave Kendall

## GBIF duplicates                                                                - Check GBIF issues, but John's SDM code will get rid of more

## GBIF spatial outliers: Ocean, middle of Australia, etc. pp? Duplicated?        - Check each map, email Alex Zika RE geoclean

## Duplicates between GBIF and ALA                                                - See email from CSIRO, check in again with them





#########################################################################################################################
#####################################################  END ############################################################## 
#########################################################################################################################