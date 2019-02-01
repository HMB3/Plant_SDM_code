#########################################################################################################################
####################################################  FILTER GBIF DATA ################################################## 
#########################################################################################################################


#########################################################################################################################
## This code combines the GBIF data into one file, and filters to the most reliable recrods. 
## Two main sources of uncertainty:


## Taxonomic
## Spatial. 


#########################################################################################################################
## Read in all data to run the SDM code :: species lists, shapefile, rasters & tables
#source('./R/HIA_LIST_MATCHING.R')
rasterTmpFile()


#########################################################################################################################
## 1). COMBINE GBIF DATAFRAMES INTO ONE
#########################################################################################################################


## Print the species run to the screen
message('Combining GBIF occurrence data for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")


#########################################################################################################################
## Create a list of species from the downloaded files
gbif.download = list.files(GBIF_path, pattern = ".RData")
length(gbif.download)


## Now these lists are getting too long for the combine step. Restrict to just the matching strings for each run?
gbif.spp.download <- paste(GBIF.spp, "_GBIF_records.RData", sep = "")
gbif.download     = gbif.download[gbif.download %in% gbif.spp.download ] 
message('downloaded species ', length(gbif.download), ' analyzed species ', length(GBIF.spp))



#########################################################################################################################
## Combine all the taxa into a single dataframe at once
GBIF.ALL <- gbif.download %>%   
  
  ## Pipe the list into lapply
  lapply(function(x) {
    
    ## Create a character string of each .RData file
    f <- sprintf(paste0(GBIF_path, "%s"), x)
    
    ## Load each file
    d <- get(load(f))
    
    ## Now drop the columns which we don't need
    message ('Reading GBIF data for ', x)
    
    ## Check if the dataframes have data
    if (nrow(d) <= 2) {
      
      ## If the species has < 2 records, escape the loop
      print (paste ("No GBIF records for ", x, " skipping "))
      return (d)
      
    }
    
    dat <- data.frame(searchTaxon = x, d[, colnames(d) %in% gbif.keep],
                      stringsAsFactors = FALSE)
    
    if(!is.character(dat$gbifID)) {
      
      dat$gbifID <- as.character(dat$gbifID)
      
    }
    
    ## Need to print the object within the loop
    names(dat)[names(dat) == 'decimalLatitude']  <- 'lat'
    names(dat)[names(dat) == 'decimalLongitude'] <- 'lon'
    dat$searchTaxon = gsub("_GBIF_records.RData", "", dat$searchTaxon)
    return(dat)
    
  }) %>%
  
  ## Finally, bind all the rows together
  bind_rows


#########################################################################################################################
## What proportion of the dataset has no lat/lon? Need to check this so we know the latest download is working
formatC(dim(GBIF.ALL)[1], format = "e", digits = 2)
(sum(is.na(GBIF.ALL$lat))            + dim(subset(GBIF.ALL, year < 1950))[1])/dim(GBIF.ALL)[1]*100


## Almost none of the GBIF data has no scientificName. This is the right field to use for matching taxonomy
(sum(is.na(GBIF.ALL$scientificName)) + dim(subset(GBIF.ALL, scientificName == ""))[1])/dim(GBIF.ALL)[1]*100


## Now get just the columns we want to keep.
GBIF.TRIM <- GBIF.ALL %>% 
  select(one_of(gbif.keep))
names(GBIF.TRIM)
gc()


#########################################################################################################################
## Just get the newly downloaded species
GBIF.TRIM = GBIF.TRIM[GBIF.TRIM$searchTaxon %in% GBIF.spp, ]
formatC(dim(GBIF.TRIM)[1], format = "e", digits = 2)


## What are the unique species?
length(unique(GBIF.TRIM$species))
length(unique(GBIF.TRIM$searchTaxon)) 
length(unique(GBIF.TRIM$scientificName))  
 




#########################################################################################################################
## 3). FILTER RECORDS TO THOSE WITH COORDINATES, AND AFTER 1950
#########################################################################################################################


#########################################################################################################################
## Now filter the GBIF records using conditions which are not too restrictive
GBIF.CLEAN <- GBIF.TRIM %>% 
  
  ## Note that these filters are very forgiving...
  ## Unless we include the NAs, very few records are returned!
  filter(!is.na(lon) & !is.na(lat),
         
         ## CULTIVATED == 'CULTIVATED', 
         ## #establishmentMeans!='MANAGED' | is.na(establishmentMeans),
         ## New.Taxonomic.status == 'Accepted'
         
         year >= 1950 & !is.na(year))


## How many species are there?
names(GBIF.CLEAN)
length(unique(GBIF.CLEAN$searchTaxon))





#########################################################################################################################
## 4). REMOVE POINTS OUTSIDE WORLDCLIM LAYERS...
#########################################################################################################################


#########################################################################################################################
## Can use WORLDCIM rasters to get only records where wordlclim data is. 
message('Removing GBIF points in the ocean for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")


## Now get the XY centroids of the unique 1km * 1km WORLDCLIM blocks where GBIF records are found
## Get cell number(s) of WORLDCLIM raster from row and/or column numbers. Cell numbers start at 1 in the upper left corner, 
## and increase from left to right, and then from top to bottom. The last cell number equals the number of raster cells 
xy <- cellFromXY(world.grids.current, GBIF.CLEAN[c("lon", "lat")]) %>% 
  
  ## get the unique raster cells
  unique %>% 
  
  ## Get coordinates of the center of raster cells for a row, column, or cell number of WORLDCLIM raster
  xyFromCell(world.grids.current, .)


## For some reason, we need to convert the xy coords to a spatial points data frame, in order to avoid this error:
## 'NAs introduced by coercion to integer range'
xy <- SpatialPointsDataFrame(coords = xy, data = as.data.frame(xy),
                             proj4string = CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))


## Now extract the temperature values for the unique 1km centroids which contain GBIF data
message('Removing GBIF points in the ocean for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")
class(xy)
z   = raster::extract(world.grids.current[["bio_01"]], xy)
#hist(z, border = NA, col = "orange", breaks = 50, main = "", xlab = "Worldclim Annual temp")


## Then track which values of Z are on land or not
onland = z %>% is.na %>%  `!` # %>% xy[.,]  cells on land or not


## Finally, filter the cleaned GBIF data to only those points on land. 
## This is achieved with the final [onland]
GBIF.LAND = filter(GBIF.CLEAN, cellFromXY(world.grids.current, GBIF.CLEAN[c("lon", "lat")]) %in% 
                     unique(cellFromXY(world.grids.current,    GBIF.CLEAN[c("lon", "lat")]))[onland])


## how many records were on land?
records.ocean = dim(GBIF.CLEAN)[1] - dim(GBIF.LAND)[1]  ## 91575 records are in the ocean   


## Print the dataframe dimensions to screen
dim(GBIF.LAND)
length(unique(GBIF.LAND$searchTaxon))


## Add a source column
GBIF.LAND$SOURCE = 'GBIF'
unique(GBIF.LAND$SOURCE)
names(GBIF.LAND)

## Free some memory
gc()


## note that this still leaves lots of points in the Islands. So we need to decide if those are legitimate.
# WORLD <- readOGR("./data/base/CONEXTUAL/TM_WORLD_BORDERS-0.3.shp", layer = "TM_WORLD_BORDERS-0.3")
# LAND  <- readOGR("./data/base/CONTEXTUAL/ne_10m_land.shp", layer = "ne_10m_land")
# plot(WORLD)
# plot(LAND)




#########################################################################################################################
## save data
if(save_data == "TRUE") {
  
  ## save .rds file for the next session
  saveRDS(GBIF.LAND, paste0(GBIF_path, 'GBIF_LAND_', save_run, '.rds'))
 
} else {
  
  message(' skip file saving, not many species analysed')   ##
  
}

## Get rid of some memory
gc()




#########################################################################################################################
## OUTSTANDING CLEANING TASKS:
#########################################################################################################################


## Increase the taxonomic check to include all species on HIA list


#########################################################################################################################
###################################################### TBC ############################################################## 
#########################################################################################################################