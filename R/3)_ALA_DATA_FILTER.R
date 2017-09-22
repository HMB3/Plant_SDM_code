#########################################################################################################################
####################################################  FILTER ALA DATA ################################################## 
#########################################################################################################################


#########################################################################################################################
## 3). REMOVE POINTS OUTSIDE WORLDCLIM LAYERS
#########################################################################################################################


## Can use WORLDCIM rasters to get only records where wordlclim data is. 
## At the global scale, there probably is no alterntive to using WORLDCLIM...


## First get one of the BIOCLIM variables
world.temp = raster("//sci-7910/F/data/worldclim/world/0.5/bio/current/bio_01")
plot(world.temp)


## Now get the XY centroids of the unique 1km * 1km WORLDCLIM blocks where ALA records are found
## Get cell number(s) of WORLDCLIM raster from row and/or column numbers. Cell numbers start at 1 in the upper left corner, 
## and increase from left to right, and then from top to bottom. The last cell number equals the number of raster cells 
xy <- cellFromXY(world.temp, AVH.OEH.VASC.XY[c("lon", "lat")]) %>% 
  
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


## Now extract the temperature values for the unique 1km centroids which contain ALA data
class(xy)
z   = extract(world.temp, xy)


# Warning message:
#   In .doExtract(x, i, ..., drop = drop) :
#   some indices are invalid (NA returned)
hist(z, border = NA, col = "orange", breaks = 50, main = "", xlab = "Worldclim Annual temp")


## Then track which values of Z are on land or not
onland = z %>% is.na %>%  `!` # %>% xy[.,]  cells on land or not
summary(onland)


## Finally, filter the cleaned ALA data to only those points on land. 
## This is achieved with the final [onland]
ALA.LAND = filter(AVH.OEH.VASC.XY, cellFromXY(world.temp, AVH.OEH.VASC.XY[c("lon", "lat")]) %in% 
                    unique(cellFromXY(world.temp, AVH.OEH.VASC.XY[c("lon", "lat")]))[onland])


## how many records were on land?
records.ocean = dim(AVH.OEH.VASC.XY)[1] - dim(ALA.LAND)[1]  ## 91575 records are in the ocean   
records.ocean


#########################################################################################################################
## Plot cleaned data
plot(world.temp)
points(ALA.CLEAN[c("lon", "lat")], pch = ".", col = "red")


## Plot cleaned data that's in the worldclim raster
plot(world.temp)
points(ALA.LAND[c("lon", "lat")], pch = ".", col = "blue")
gc()


## note that this still leaves lots of points in the Islands. So we need to decide if those are legitimate.
WORLD <- readOGR("./data/base/CONTEXTUAL/TM_WORLD_BORDERS-0.3.shp", layer = "TM_WORLD_BORDERS-0.3")
LAND  <- readOGR("./data/base/CONTEXTUAL/ne_10m_land.shp", layer = "ne_10m_land")
plot(WORLD)
plot(LAND)


#########################################################################################################################
## save data
save(ALA.LAND, file = paste("./data/base/HIA_LIST/ALA/ALA_LAND_POINTS.RData"))
gc()



#########################################################################################################################
## OUTSTANDING CLEANING TASKS:
#########################################################################################################################


## Get GBIF duplicates 

## Get GBIF spatial outliers

## Get GBIF taxonomic errors - probably the most serious erro

## Get duplicates between GBIF and ALA - keep in touch with CSIRO


#########################################################################################################################
#####################################################  END ############################################################## 
#########################################################################################################################