#########################################################################################################################
#############################################  CREATE GBIF CLEANING COLUMNS ############################################# 
#########################################################################################################################


#########################################################################################################################
## 1). CREATE PRE-CLEAN FLAGS 
#########################################################################################################################


## Which columns overlap between GBIF and ALA/AVH? The problem here is that we don't no what each GBIF field means, and 
## the documentation is a bit dodgy.

## dataProvider             (ALA), GBIF ?
## basisOfRecord            (ALA), basisOfRecord (GBIF)
## occCultivatedEscapee     (ALA), establishmentMeans (GBIF)
## year
## taxonIdentificationIssue (ALA), GBIF no equivalent
## fatal assertions         (ALA), GBIF no equivalent


## create copy and free memory
GBIF.CLEAN = GBIF.TRIM
gc()


#########################################################################################################################
## Create a table which counts the number of records meeting all the criteria
## Note that TRUE indicates there is a problem (e.g. no lat/long = TRUE)
GBIF.PROBLEMS <- with(GBIF.TRIM,
                 
                 table(
                   
                   ## no coordinates
                   is.na(lon)|is.na(lat),
                   
                   ## establishment means is "MANAGED", can change this
                   establishmentMeans == 'MANAGED' & !is.na(establishmentMeans),
                   
                   ## collected before 1950 or na.year
                   year < 1950 & !is.na(year),
                   
                   ## no year
                   is.na(year),
                   
                   ## coordinate uncertainty is > 100 or is not NA
                   coordinateUncertaintyInMeters > 1000 & 
                     !is.na(coordinateUncertaintyInMeters)
                   
                 )
                 
) %>% 
  
  ## create a data frame and set the names
  as.data.frame %>%  
  
  setNames(c('NO_COORD', 'MANAGED', 
             'PRE_1950', 'NO_YEAR', 
             'COORD_UNCERT', 'COUNT'))


## print table to screen: in .Rmd file, no need to save 
kable(GBIF.PROBLEMS)


## write file to CSV
write.csv(GBIF.PROBLEMS, "./output/tables/GBIF_PROBLEMS.csv", row.names = FALSE)





#########################################################################################################################
## 2). FILTER RECORDS 
#########################################################################################################################


## Filter the GBIF records using conditions which are not too restrictive
GBIF.CLEAN <- GBIF.CLEAN %>% 
  
  ## Note that these filters are very forgiving...
  ## unless we exclude the NAs, very few records are returned! 
  filter(!is.na(lon) & !is.na(lat),
         establishmentMeans! = 'MANAGED' | is.na(establishmentMeans),
         year >= 1950 & !is.na(year))

## The table above gives the details, but worth documenting how many records are knocked out by each
remaining.records = dim(GBIF.CLEAN)[1]/total.records*100  
remaining.records ## 65% of records remain after cleaning 
gc()





#########################################################################################################################
## 3). REMOVE SPATIAL OUTLIERS
#########################################################################################################################


## Can use WORLDCIM rasters to get only records where wordlclim data is. 
## At the global scale, there probably is no alterntive to using WORLDCLIM...


## First get one of the BIOCLIM variables
world.temp = raster("//sci-7910/F/data/worldclim/world/0.5/bio/current/bio_01")
plot(world.temp)


## Now get the XY centroids of the unique 1km * 1km WORLDCLIM blocks where GBIF records are found
xy <- cellFromXY(world.temp, GBIF.CLEAN[c("lon", "lat")]) %>% 
  unique %>% 
  #length
  cellFromXY(world.temp, .)


## take a look at xy: NA's should be removed
## but the numbers are too big
summary(xy)
str(xy)
points(xy, pch = ".", col = "red")


## For some reason, we need to convert the xy coords to a spatial points data frame
## This is to avoid:  'NAs introduced by coercion to integer range'
xy <- SpatialPointsDataFrame(coords = xy, data = df,
                             proj4string = CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))


## Now extract the temperature values for the unique 1km centroids which contain GBIF data
z   = extract(world.temp, xy)


## Then remove the NA values from the temperature values: these are on land
onland = z %>% is.na %>%  `!` # %>% xy[.,]  cells on land or not
summary(onland)


## Finally, filter the cleaned GBIF data to only those points on land 
GBIF.LAND = filter(GBIF.CLEAN, cellFromXY(world.temp, GBIF.CLEAN[c("lon", "lat")]) %in% 
                         unique(cellFromXY(world.temp, GBIF.CLEAN[c("lon", "lat")]))[onland])


## how many records were on land?
dim(GBIF.CLEAN)[1] - dim(GBIF.LAND)[1]  ## 91575 missing records   
  

#########################################################################################################################
## Plot cleaned data
plot(world.temp)
points(GBIF.CLEAN[c("lon", "lat")], pch = ".", col = "red")


## Plot cleaned data that's in the worldclim raster
plot(world.temp)
points(GBIF.LAND[c("lon", "lat")], pch = ".", col = "red")
gc()


#########################################################################################################################
## save data
save(GBIF.LAND, file = paste("./data/base/HIA_LIST/GBIF/GBIF_LAND_POINTS.RData"))





#########################################################################################################################
#####################################################  END ############################################################## 
#########################################################################################################################