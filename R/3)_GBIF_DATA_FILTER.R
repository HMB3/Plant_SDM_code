#########################################################################################################################
####################################################  FILTER GBIF DATA ################################################## 
#########################################################################################################################


#########################################################################################################################
## 1). CREATE TABLE OF PRE-CLEAN FLAGS 
#########################################################################################################################


## Which columns overlap between GBIF and ALA/AVH? The problem here is that we don't no what each GBIF field means, and 
## the documentation is a bit dodgy.

## dataProvider             (ALA), GBIF ?
## basisOfRecord            (ALA), basisOfRecord (GBIF)
## occCultivatedEscapee     (ALA), establishmentMeans (GBIF)
## year
## taxonIdentificationIssue (ALA), GBIF no equivalent
## fatal assertions         (ALA), GBIF no equivalent

## Currently, what is missing from this workflow is checking duplicate refords, checking taxonomy and checking the 
## spatial outliers.


#########################################################################################################################
## Create a table which counts the number of records meeting all the criteria
## Note that TRUE indicates there is a problem (e.g. no lat/long = TRUE)
GBIF.PROBLEMS <- with(GBIF.TRIM,
                 
                 table(
                   
                   ## note this list could be exanded for other data types
                   ## ALA/AVH, council data, etc.
                   
                   ## no coordinates
                   is.na(lon)|is.na(lat),
                   
                   ## establishment means is "MANAGED", can change this
                   establishmentMeans == 'MANAGED' & !is.na(establishmentMeans),
                   
                   ## also add duplicated
                   
                   ## collected before 1950 or na.year
                   year < 1950 & !is.na(year),
                   
                   ## no year
                   is.na(year),
                   
                   ## coordinate uncertainty is > 100 or is not NA
                   coordinateUncertaintyInMeters > 1000 & 
                     !is.na(coordinateUncertaintyInMeters)
                   
                   ## add maybe the centre of Australia?
                   ## Lamber centre of Aus: 25.610111, 134.354806
                   
                 )
                 
) %>% 
  
  ## create a data frame and set the names
  as.data.frame %>%  
  
  setNames(c('NO_COORD', 'MANAGED', 
             'PRE_1950', 'NO_YEAR', 
             'COORD_UNCERT', 'COUNT'))


## print table to screen: in .Rmd file, no need to save
## Note that TRUE indicates there is a problem (e.g. no lat/long = TRUE)
kable(GBIF.PROBLEMS)


## quickly check the total record number matches the count of problems
Total.count = sum(GBIF.PROBLEMS$COUNT)
identical(total.records, total.count)  ## identical matches two objects


## probably don't need this
write.csv(GBIF.PROBLEMS, "./output/tables/GBIF_PROBLEMS.csv", row.names = FALSE)





#########################################################################################################################
## 2). FILTER RECORDS 
#########################################################################################################################


## Filter the GBIF records using conditions which are not too restrictive
GBIF.CLEAN <- GBIF.TRIM %>% 
  
  ## Note that these filters are very forgiving...
  ## unless we include the NAs, very few records are returned!
  filter(!is.na(lon) & !is.na(lat),
         establishmentMeans!='MANAGED' | is.na(establishmentMeans),
         year >= 1950 & !is.na(year))


## The table above gives the details, but worth documenting how many records are knocked out by each
Remaining.records = dim(GBIF.CLEAN)[1] 
Remaining.percent = dim(GBIF.CLEAN)[1]/total.records*100  
Remaining.percent ## 65% of records remain after cleaning 
gc()


## check
dim(GBIF.CLEAN)
head(GBIF.CLEAN)


## Store character string of the filters applied
Filters.applied = "NA COORD | MANAGED/NA | < 1950/NA"





#########################################################################################################################
## 3). REMOVE POINTS OUTSIDE WORLDCLIM LAYERS
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


## take a look at xy: NA's should be removed...but the numbers are too big
summary(xy)
str(xy)
points(xy, pch = ".", col = "red")


## For some reason, we need to convert the xy coords to a spatial points data frame, in order to avoid the error:
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
plot(world.temp)
points(GBIF.CLEAN[c("lon", "lat")], pch = ".", col = "red")


## Plot cleaned data that's in the worldclim raster
plot(world.temp)
points(GBIF.LAND[c("lon", "lat")], pch = ".", col = "blue")
gc()


## note that this still leaves lots of points in the Islands. So we need to decide if those are legitimate.
WORLD <- readOGR("./data/base/CONTEXTUAL/TM_WORLD_BORDERS-0.3.shp", layer = "TM_WORLD_BORDERS-0.3")
LAND  <- readOGR("./data/base/CONTEXTUAL/ne_10m_land.shp", layer = "ne_10m_land")
plot(WORLD)
plot(LAND)


## consider some more ways of spatially checking the records in the ocean...





#########################################################################################################################
## 4). CHECK EACH SPECIES FOR SPATIAL OUTLIERS
#########################################################################################################################


##



#########################################################################################################################
## save data
save(GBIF.LAND, file = paste("./data/base/HIA_LIST/GBIF/GBIF_LAND_POINTS.RData"))
gc()




#########################################################################################################################
#####################################################  END ############################################################## 
#########################################################################################################################