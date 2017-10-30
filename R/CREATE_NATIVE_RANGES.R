#########################################################################################################################
############################################ CHECK TAXONOMY OF RETURNED DATA ############################################ 
#########################################################################################################################



#########################################################################################################################
## 1). CHOOSE SPECIES
#########################################################################################################################


## Also when doing a subset, need to decide on the factors for choosing species:
## Mix of families
## Mix of % cultivated/non
## Spatial data available for native vs. non-native
## Mix of natives/exotic


## The top 50
HIA.SORT = HIA.SPP[order(HIA.SPP$Number.of.growers, rev(HIA.SPP$Top_200), decreasing = TRUE), ]
View(HIA.SORT)


## From the top 50, choose
View(head(HIA.SORT, 50)) 
TEST.SPP = head(HIA.SORT$Binomial, 50) 


## What is the overlap between
intersect(renee.full$Species, TEST.SPP)
intersect(GBIF.RANGE$searchTaxon, TEST.SPP)


## Add on the species which have geographic ranges
TEST.SPP = c(TEST.SPP, "Betula pendula", "Fraxinus excelsior", "Quercus robur")
TEST.SPP


## quickly check the GBIF names:
sp.n = TEST.SPP[10]
sp.n

taxa.test = gbif(sp.n, download = TRUE)
GBIF.names = sort(names(taxa.test))


#########################################################################################################################
## Try and find terms "garden" and "cultivated" in particular columns
taxa.test$CULTIVATED <- ifelse(grepl("garden|cultiva",   taxa.test$locality,           ignore.case = TRUE) | 
                                 grepl("garden|cultiva", taxa.test$habitat,            ignore.case = TRUE) | 
                                 grepl("garden|cultiva", taxa.test$eventRemarks,       ignore.case = TRUE) |
                                 grepl("garden|cultiva", taxa.test$cloc,               ignore.case = TRUE) |
                                 grepl("managed",        taxa.test$establishmentMeans, ignore.case = TRUE),
                               
                               "CULTIVATED", "UNKNOWN")


## This is probably a bit strict, in that for some of the fields, garden doesn't = cultivated
unique(taxa.test$CULTIVATED)
test.cult = subset(taxa.test, CULTIVATED == "CULTIVATED")
dim(test.cult)[1]/dim(taxa.test)[1]
View(test.cult)



#########################################################################################################################
## GBIF issues? Not that helpful...
# Abelia.geosp.t = gbif('Abelia grandiflora', args = list("hasGeospatialIssue=true"))
# Abelia.geosp.f = gbif('Abelia grandiflora', args = list("hasGeospatialIssue=false"))
# Abelia         = gbif('Abelia grandiflora', args = list("hasGeospatialIssue=false"))


#########################################################################################################################
## 2). READ IN SPECIES RANGES
#########################################################################################################################


## Load data
EURO.RANGES = read.csv("./data/base/TRAITS/EURO_RANGE_SPECIES.csv", stringsAsFactors = FALSE)
load("./data/base/HIA_LIST/GBIF/GBIF_LAND_POINTS.RData")


## Get species list which have ranges, and restrict big data frame to just those species
HIA.RANGE.SPP = intersect(HIA.SPP$Binomial, EURO.RANGES$Species)
GBIF.RANGE    = GBIF.LAND[GBIF.LAND$searchTaxon %in% HIA.RANGE.SPP, ]
unique(GBIF.RANGE$searchTaxon)


#########################################################################################################################
## Create a list of shapefile and read them in
# Fagus.sylv = readOGR("./data/base/CONTEXTUAL/RANGES/Fagus_sylvatica_EUFORGEN.shp", layer = "Fagus_sylvatica_EUFORGEN")
# names(Fagus.sylv)
## The computer is running very slowly...


## List the .shp
range.list <- list.files(path      = "./data/base/CONTEXTUAL/RANGES/",  ## include the $ to stop the XML's being included
                         pattern   = "*.shp$", full.names = TRUE,
                         recursive = TRUE,     include.dirs = FALSE)


## Read them in
range.shp <- lapply(range.list, function(x) {readOGR(dsn = x, 
                                                     layer = ogrListLayers(x))})


## Can access each shapefile by indexing the list
class(range.shp[[1]])
shp.list = 2:5


## Plot each shapefile
lapply(range.list, function(x) {plot(range.shp[[x]], 
                                     col = "light blue", 
                                     main = range.shp[[x]]$Species)
  
  points(GBIF.RANGE[ which(GBIF.RANGE$searchTaxon == range.shp[[x]]$Species), ][, c("lon", "lat")],
         pch = ".", col = "red",
         cex = 1.3, asp = 1)
})





#########################################################################################################################
## 3). CHECK IF SPECIES ARE IN THE RANGE OR NOT
#########################################################################################################################


## Can we achieve the native calculation by just adding a column for native range? This will be 


## Convert to spatial points data frame
GBIF.RANGE.SP   = SpatialPointsDataFrame(coords = GBIF.RANGE[c("lon", "lat")], 
                                         data   = GBIF.RANGE,
                                         proj4string = CRS("+init=epsg:4326"))


## Project
Betula.pendula.range = range.shp[[2]]
Betula.pendula.range = range.shp[[3]]
Betula.pendula.range = range.shp[[4]]
Betula.pendula.range = range.shp[[5]]





CRS.new  <- CRS("+init=epsg:4326") # EPSG:3577
Betula.pendula.range  = spTransform(Betula.pendula.range, CRS.new)


## Now what is the easiest way to record the native range? By checking if points are in the polygon?
class(Betula.pendula.range)
class(Betula.pendula)
names(Betula.pendula.range)


Betula.pendula = GBIF.RANGE.SP[GBIF.RANGE.SP@data$searchTaxon == "Betula pendula", ]
Betula.over    = over(GBIF.RANGE.SP, #Betula.pendula, 
                      Betula.pendula.range)


## So over returns the species name if the point is inside that species range.
## but this would make the table too big, because you would need a column for each species
str(Betula.over)
unique(Betula.over$Species)
count(Betula.over, Species)




#########################################################################################################################
#####################################################  END ############################################################## 
#########################################################################################################################