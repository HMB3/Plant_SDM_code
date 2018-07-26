#########################################################################################################################
############################################# UPDATE ALA  ############################################################### 
#########################################################################################################################


#########################################################################################################################
## This code updates the ALA data


## John's ALA data :: 
load("./data/base/HIA_LIST/ALA/ALA_LAND_POINTS.RData")


## Load ALA update :: the coordinates are coming in as characters
ALA.UPDATE = read.csv("./data/base/HIA_LIST/ALA/Plants.csv", stringsAsFactors = FALSE)
length(unique(ALA.UPDATE$scientificName))  ## 30k species 




#########################################################################################################################
## 1). LOAD TAXA, ADD SEARCH TAXA, DROP COLUMNS AND COMBINE SPECIES DATAFRAMES INTO ONE
#########################################################################################################################


## What are the names?
sort(names(ALA.UPDATE))
intersect(names(ALA.UPDATE), names(ALA.LAND))
intersect(names(ALA.UPDATE), ALA.keep)


## Some coordinates are invalid? Can't look at the whole file
ALA.UPDATE$lat = as.numeric(ALA.UPDATE$decimalLatitude)
ALA.UPDATE$lon = as.numeric(ALA.UPDATE$decimalLongitude)
summary(ALA.UPDATE$lat)
summary(ALA.UPDATE$lon)


#########################################################################################################################
## Now get just the columns we want to keep. Note gc() frees up RAM
ALA.TRIM <- ALA.UPDATE %>% 
  select(one_of(ALA.keep))


## Check names
dim(ALA.TRIM)
names(ALA.TRIM)
setdiff(names(ALA.TRIM), ALA.keep)


## Just get the species we need - this should be the intersection of Ale's list and the HIA list
ALA.TRIM = ALA.TRIM[ALA.TRIM$scientificName %in% CLEAN.SPP$Binomial, ]
dim(ALA.TRIM)
length(unique(ALA.TRIM$scientificName))     ## not many records per species


#########################################################################################################################
## Test the difference between John's version and mine for a subset of species
test.new.ALA = ALA.TRIM[ALA.TRIM$scientificName %in% test.exotics, ]
test.old.ALA = ALA.LAND[ALA.LAND$scientificname %in% test.exotics, ]


## Are the counts of records per species different?
table(test.new.ALA$scientificName)
table(test.old.ALA$scientificname)



## Rename columns to match GBIF code
# names(ALA.TRIM)[names(ALA.TRIM) == 'decimalLatitude']  <- 'lat'
# names(ALA.TRIM)[names(ALA.TRIM) == 'decimalLongitude'] <- 'lon'


summary(ALA.CLEAN$lat)
summary(ALA.CLEAN$lon)
head(ALA.CLEAN$lon)



#########################################################################################################################
## 2). CREATE TABLE OF PRE-CLEAN FLAGS AND FILTER RECORDS
#########################################################################################################################


#########################################################################################################################
## Create a table which counts the number of records meeting each criteria:
## Note that TRUE indicates there is a problem (e.g. if a record has no lat/long, it will = TRUE)
ALA.TRIM.TAXO = ALA.TRIM
ALA.PROBLEMS <- with(ALA.TRIM.TAXO,
                      
                      table(
                        
                        ## Note this list is incomplete
                        
                        ## No coordinates
                        is.na(lon)|is.na(lat),
                        
                        ## Taxon rank is genus/form?
                        #taxonRank == 'GENUS' & 'FORM',
                        
                        ## Taxonomic status: consider if this is the right filter, it seems too restrictive
                        #New.Taxonomic.status == 'Unresolved',
                        
                        ## Cultivated
                        #CULTIVATED == 'CULTIVATED',
                        
                        ## Establishment means is "MANAGED" is included in the above
                        #establishmentMeans == 'MANAGED' & !is.na(establishmentMeans),
                        
                        ## Collected before 1950 or na.year
                        year < 1950 & !is.na(year),
                        
                        ## No year
                        is.na(year),
                        
                        ## Coordinate uncertainty is > 100 or is not NA
                        coordinateUncertaintyInMeters > 1000 & 
                          !is.na(coordinateUncertaintyInMeters)
                        
                        ## Other checks using coordinateCleaner
                        
                      )
                      
) %>% 
  
  ## Create a data frame and set the names
  as.data.frame %>%  
  
  setNames(c('NO_COORD',     
             #'TAXON_STATUS', 
             #'CULTIVATED', 
             'PRE_1950',     'NO_YEAR', 
             'COORD_UNCERT', 'COUNT'))


#########################################################################################################################
## Now filter the ALA records using conditions which are not too restrictive
ALA.CLEAN <- ALA.TRIM.TAXO %>% 
  
  ## Note that these filters are very forgiving...
  ## Unless we include the NAs, very few records are returned!
  filter(!is.na(lon) & !is.na(lat),
         
         ## CULTIVATED == 'CULTIVATED', 
         ## #establishmentMeans!='MANAGED' | is.na(establishmentMeans),
         ## New.Taxonomic.status == 'Accepted'
         
         year >= 1950 & !is.na(year))


## Check
names(ALA.CLEAN)





#########################################################################################################################
## 6). REMOVE POINTS OUTSIDE WORLDCLIM LAYERS...
#########################################################################################################################


## Can use WORLDCIM rasters to get only records where wordlclim data is. 
## First, get one of the BIOCLIM variables
world.temp = raster("./data/base/worldclim/world/0.5/bio/current/bio_01")
#plot(world.temp)


## Now get the XY centroids of the unique 1km * 1km WORLDCLIM blocks where ALA records are found
## Get cell number(s) of WORLDCLIM raster from row and/or column numbers. Cell numbers start at 1 in the upper left corner, 
## and increase from left to right, and then from top to bottom. The last cell number equals the number of raster cells 
xy <- cellFromXY(world.temp, ALA.CLEAN[c("lon", "lat")]) %>% 
  
  ## get the unique raster cells
  unique %>% 
  
  ## Get coordinates of the center of raster cells for a row, column, or cell number of WORLDCLIM raster
  xyFromCell(world.temp, .)


## Take a look at xy: NA's should be removed...
# summary(xy)
# str(xy)
# points(xy, pch = ".", col = "red")


## For some reason, we need to convert the xy coords to a spatial points data frame, in order to avoid this error:
## 'NAs introduced by coercion to integer range'
xy <- SpatialPointsDataFrame(coords = xy, data = as.data.frame(xy),
                             proj4string = CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))


## Now extract the temperature values for the unique 1km centroids which contain ALA data
class(xy)
z   = extract(world.temp, xy)

# Warning message:
#   In .doExtract(x, i, ..., drop = drop) :
#   some indices are invalid (NA returned)

hist(z, border = NA, col = "orange", breaks = 50, main = "", xlab = "Worldclim Annual temp")


## Then track which values of Z are on land or not
onland = z %>% is.na %>%  `!` # %>% xy[.,]  cells on land or not
#summary(onland)


## Finally, filter the cleaned ALA data to only those points on land. 
## This is achieved with the final [onland]
ALA.LAND = filter(ALA.CLEAN, cellFromXY(world.temp, ALA.CLEAN[c("lon", "lat")]) %in% 
                     unique(cellFromXY(world.temp, ALA.CLEAN[c("lon", "lat")]))[onland])


## how many records were on land?
records.ocean = dim(ALA.CLEAN)[1] - dim(ALA.LAND)[1]  ## 91575 records are in the ocean   


## Print the dataframe dimensions to screen
dim(ALA.LAND)
length(unique(ALA.LAND$searchTaxon))


## Free some memory
gc()



#########################################################################################################################
## save data
saveRDS(ALA.UPDATE, file = paste("./data/base/HIA_LIST/ALA/ALA_UPDATE_POINTS.rds"))
save.image("STEP_4_NICHES.RData")




#########################################################################################################################
## OUTSTANDING NICHE TASKS:
#########################################################################################################################


## 
#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################