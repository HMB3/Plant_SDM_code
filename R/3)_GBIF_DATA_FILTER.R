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
  
  as.data.frame %>%  
  
  setNames(c('NO_COORD', 'MANAGED', 
             'PRE_1950', 'NO_YEAR', 
             'COORD_UNCERT', 'COUNT'))


## print table to screen: in .Rmd file, no need to save 
kable(GBIF.PROBLEMS)


## write file to CSV
write.csv(GBIF.PROBLEMS, "./output/tables/GBIF_PROBLEMS.csv", row.names = FALSE)


#########################################################################################################################
## Now filter the records
GBIF.CLEAN <- GBIF.HIA.SPP.RECORDS.ALL %>% 
  
  select(one_of(
    
    c('basisOfRecord', 'lon', 'lat', 'establishmentMeans', 
      'year', 'scientificName', 'searchTaxon', 'gbifID', 
      'coordinateUncertaintyInMeters'))) %>% 
  
  filter(!is.na(lon) & !is.na(lat),
         establishmentMeans!='MANAGED' | is.na(establishmentMeans),
         year >= 1950 & !is.na(year))





#########################################################################################################################
## 3). REMOVE SPATIAL OUTLIERS
#########################################################################################################################


## first get one of the BIOCLIM variables
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


