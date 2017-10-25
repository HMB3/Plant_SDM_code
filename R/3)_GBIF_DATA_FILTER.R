#########################################################################################################################
####################################################  FILTER GBIF DATA ################################################## 
#########################################################################################################################


#########################################################################################################################
## 1). CREATE TABLE OF PRE-CLEAN FLAGS 
#########################################################################################################################


## Which columns overlap between GBIF and ALA/AVH? The problem here is that we don't know what each GBIF field means, and 
## the documentation is a bit dodgy.

## dataProvider             (ALA), GBIF ?
## basisOfRecord            (ALA), basisOfRecord (GBIF)
## occCultivatedEscapee     (ALA), establishmentMeans (GBIF)
## year
## taxonIdentificationIssue (ALA), GBIF no equivalent
## fatal assertions         (ALA), GBIF no equivalent


## Rach considers anything with the same lat/long to x decimal places, collected in the same month by the same person as
## a duplicate. Collection month and collector could be added, but in most cases, they will be knocked out by using 
## one record per cell.



## GBIF issues? Not that helpful...
# Abelia.geosp.t = gbif('Abelia grandiflora', args = list("hasGeospatialIssue=true"))
# Abelia.geosp.f = gbif('Abelia grandiflora', args = list("hasGeospatialIssue=false"))
# Abelia         = gbif('Abelia grandiflora', args = list("hasGeospatialIssue=false"))


#########################################################################################################################
## LOAD
load("./data/base/HIA_LIST/COMBO/GBIF_TRIM.RData")
dim(GBIF.TRIM)
names(GBIF.TRIM)
length(unique(GBIF.TRIM$searchTaxon))  ## has the list updated with extra species? YES!


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
  
  setNames(c('NO_COORD', #'GEOCLEAN', 'DUPLICATED', 
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
## 2). FILTER RECORDS 
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
## 3). REMOVE POINTS OUTSIDE WORLDCLIM LAYERS...
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
## So here we could add extra filters?
names(GBIF.LAND)





#########################################################################################################################
## save data
save(GBIF.LAND, file = paste("./data/base/HIA_LIST/GBIF/GBIF_LAND_POINTS.RData"))
gc()


## Now save .RData file for the next session
save.image("STEP_3_GBIF_CLEAN.RData")





#########################################################################################################################
## OUTSTANDING CLEANING TASKS:
#########################################################################################################################


## When should the additional filters be run in? Just after GBIF.TRIM?            -

## Keep managed records as a separate file                                        - Have them

## GBIF duplicates                                                                - Check GBIF issues, but John's SDM code will get rid of more

## GBIF species match                                                             - Species summary will take care of it. 

## GBIF spatial outliers: Ocean, middle of Australia, etc. ppp? Duplicated?       - Check each map, email Alex Zika

## GBIF taxonomic errors?                                                         - Get GBIF issues uisng GBIF ID? Or redownload with occ_search

## Duplicates between GBIF and ALA                                                - See email from CSIRO





#########################################################################################################################
#####################################################  END ############################################################## 
#########################################################################################################################