#########################################################################################################################
############################################ CHECK TAXONOMY OF RETURNED DATA ############################################ 
#########################################################################################################################


#########################################################################################################################
## 1). READ IN SPECIES RANGES
#########################################################################################################################


## Load data
EURO.RANGES = read.csv("./data/base/TRAITS/EURO_RANGE_SPECIES.csv", stringsAsFactors = FALSE)
load("./data/base/HIA_LIST/GBIF/GBIF_LAND_POINTS.RData")


## Get species list which have ranges, and restrict big data frame to just those species
HIA.RANGE.SPP = intersect(HIA.SPP$Binomial, EURO.RANGES$Species)
GBIF.RANGE    = GBIF.LAND[GBIF.LAND$searchTaxon %in% HIA.RANGE.SPP, ] 


## Also when doing a subset, need to decide on the factors for choosing species:
## Mix of families
## Mix of % cultivated/non
## Spatial data available for native vs. non-native
## Mix of natives/exotic


#########################################################################################################################
## Create a list of shapefile and read them in
Fagus.sylv = readOGR("./data/base/CONTEXTUAL/RANGES/Fagus_sylvatica_EUFORGEN.shp", layer = "Fagus_sylvatica_EUFORGEN")
names(Fagus.sylv)


range.list <- list.files(path      = "./data/base/CONTEXTUAL/RANGES/",  ## include the $ to stop the XML's being included
                         pattern   = "*.shp$", full.names = TRUE,
                         recursive = TRUE,     include.dirs = FALSE)


## Create a list of all the data frames
range.shp <- lapply(range.list, function(x) {readOGR(dsn = x, 
                                                     layer = ogrListLayers(x))})

## Can access each shapefile by indexing the list
class(range.shp[[1]])
range.list = 2:5


## Plot each shapefile
lapply(range.list, function(x) {plot(range.shp[[x]], 
                                     col = "light blue", 
                                     main = range.shp[[x]]$Species)
  
  points(GBIF.RANGE[ which(GBIF.RANGE$searchTaxon == range.shp[[x]]$Species), ][, c("lon", "lat")],
         pch = ".", col = "red",
         cex = 1.3, asp = 1)
})





#########################################################################################################################
## 2). CHECK IF SPECIES ARE IN THE RANGE OR NOT
#########################################################################################################################


## Can we achieve the native calculation by just adding a column for native range? This will be 


## Convert to spatial points data frame
GBIF.RANGE.SP   = SpatialPointsDataFrame(coords = GBIF.RANGE[c("lon", "lat")], 
                                         data   = GBIF.RANGE,
                                         proj4string = CRS("+init=epsg:4326"))


## Project
Betula.pendula.range = range.shp[[2]]
CRS.new  <- CRS("+init=epsg:4326") # EPSG:3577
Betula.pendula.range  = spTransform(Betula.pendula.range, CRS.new)


## Now what is the easiest way to record the native range? By checking if points are in the polygon?
class(Betula.pendula.range)
class(Betula.pendula)
names(Betula.pendula.range)


Betula.pendula = GBIF.RANGE.SP[GBIF.RANGE.SP@data$searchTaxon == "Betula pendula", ]
Betula.over    = over(GBIF.RANGE.SP, #Betula.pendula, 
                      Betula.pendula.range)


## Not sure why this is not working...
str(Betula.over)





#########################################################################################################################
#####################################################  END ############################################################## 
#########################################################################################################################