#########################################################################################################################
############################################ CHECK TAXONOMY OF RETURNED DATA ############################################ 
#########################################################################################################################



#########################################################################################################################
## 1). CHOOSE SPECIES FOR SDM TEST
#########################################################################################################################


## When doing a subset, need to decide on the factors for choosing species:
## Mix of families
## Mix of % cultivated/non
## Spatial data available for native vs. non-native
## Mix of natives/exotic


## The top 50?
HIA.SORT = HIA.SPP[order(HIA.SPP$Number.of.growers, rev(HIA.SPP$Top_200), decreasing = TRUE), ]
HIA.SORT = head(HIA.SORT$Binomial, 50) 
intersect(renee.full$Species, HIA.SORT)


## Or, just use Renee's shortlist, and add on the the species which have geographic ranges.
## Renee and Dave have done a selection process which
test.spp = sort(unique(c(renee.full$Species, "Betula pendula", "Fraxinus excelsior", "Quercus robur", "Fagus sylvatica")))
test.spp 


## Quickly check the GBIF names:
sp.n = test.spp[1]
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
## Get just the records on the trial list
COMBO.RASTER.TRIAL = COMBO.RASTER.CONTEXT[COMBO.RASTER.CONTEXT$searchTaxon %in% test.spp, ]



#########################################################################################################################
## 2). READ IN SPECIES RANGES
#########################################################################################################################


## Load data
EURO.RANGES = read.csv("./data/base/TRAITS/EURO_RANGE_SPECIES.csv", stringsAsFactors = FALSE)
load("./data/base/HIA_LIST/GBIF/GBIF_LAND_POINTS.RData")


## Get species list which have ranges, and restrict big data frame to just those species
HIA.RANGE.SPP = intersect(HIA.SPP$Binomial, EURO.RANGES$Species)
GBIF.RANGE    = GBIF.LAND[GBIF.LAND$searchTaxon %in% test.spp, ]
names(GBIF.RANGE)
#GBIF.RANGE    = COMBO.RASTER.CONTEXT[COMBO.RASTER.CONTEXT$searchTaxon %in% test.spp, ]


## Check these species
dim(GBIF.RANGE)
unique(GBIF.RANGE$searchTaxon)  ## 900k records for 54 species


## Plot this restricted dataset, just to check that it covers the whole world, rather than just Europe
LAND  <- readOGR("./data/base/CONTEXTUAL/ne_10m_land.shp", layer = "ne_10m_land")
plot(LAND)
points(GBIF.RANGE[c("lon", "lat")], pch = ".", col = "red")




#########################################################################################################################
## Create a list of shapefile and read them in
## List the .shp's
range.list <- list.files(path      = "./data/base/CONTEXTUAL/RANGES/",  ## include the $ to stop the XML's being included
                         pattern   = "*.shp$", full.names = TRUE,
                         recursive = TRUE,     include.dirs = FALSE)


## Read them in
range.shp <- lapply(range.list, function(x) {readOGR(dsn = x, 
                                                     layer = ogrListLayers(x))})


## Can access each shapefile by indexing the list
names(range.shp)
class(range.shp[[1]])
shps = 2:5


## Plot each shapefile
lapply(shp.list, function(x) {plot(range.shp[[x]], 
                                     col = "light blue", 
                                     main = range.shp[[x]]$Species)
  
  points(GBIF.RANGE[ which(GBIF.RANGE$searchTaxon == range.shp[[x]]$Species), ][, c("lon", "lat")],
         pch = ".", col = "red",
         cex = 1.3, asp = 1)
  
})





#########################################################################################################################
## 3). CHECK IF SPECIES ARE IN THE RANGE OR NOT
#########################################################################################################################


## Can we achieve the native calculation by just adding a column for native range? Also can we do this for n species, 
## rather than one species? This is just a question of the size/time to run the over calculation. Or, can we 



## Convert to spatial points data frame
GBIF.RANGE.SP   = SpatialPointsDataFrame(coords = GBIF.RANGE[c("lon", "lat")], 
                                         data   = GBIF.RANGE,
                                         proj4string = CRS("+init=epsg:4326"))


## Individual shapefiles: should loop this
Betula.pendula.range      = range.shp[[2]]
Fagus.sylvatica.range     = range.shp[[3]]
Fraxinus.excelsior.range  = range.shp[[4]]
Quercus.robur.range       = range.shp[[5]]


## Project, and check they are the same
CRS.new  <- CRS("+init=epsg:4326") # EPSG:3577

Betula.pendula.range      = spTransform(Betula.pendula.range, CRS.new)
Fagus.sylvatica.range     = spTransform(Fagus.sylvatica.range, CRS.new)
Fraxinus.excelsior.range  = spTransform(Fraxinus.excelsior.range, CRS.new)
Quercus.robur.range       = spTransform(Quercus.robur.range, CRS.new)


## 
projection(GBIF.RANGE.SP);projection(Betula.pendula.range)


## Now what is the easiest way to record the native range? By checking if points are in the polygon?
class(Betula.pendula.range);class(Betula.pendula)
names(Betula.pendula.range);names(Betula.pendula)


## Comgine into one data frame
combined.ranges = rbind(Betula.pendula.range,
                        Fagus.sylvatica.range,
                        Fraxinus.excelsior.range,
                        Quercus.robur.range)

##
class(combined.ranges)





#########################################################################################################################
## 4). OVERLAY POLYGONS FOR MULTIPLE SPECIES
#########################################################################################################################


#########################################################################################################################
## Now what is the easiest way to record the native range? By checking if points are in the polygon?
## Use the over function to see which points are in the polygon. One at a time...
Betula.over    = over(GBIF.RANGE.SP, 
                      Betula.pendula.range)

Fagus.over     = over(GBIF.RANGE.SP, 
                      Fagus.range)

Fraxinus.over  = over(GBIF.RANGE.SP, 
                      Fraxinus.range)

Quercus.over   = over(GBIF.RANGE.SP, 
                      Quercus.range)


## Rename and reassign so the meaning is less cryptic
Betula.over = Species.over %>%
  setNames(c('Betula_pendula_range')) 
str(Species.over)

Fagus.over = Fagus.over %>%
  setNames(c('Fagus_sylvatica_range')) 
str(Fagus.over)

Fraxinus.over = Fraxinus.over %>%
  setNames(c('Fraxinus_excelsior_range')) 
str(Fraxinus.over)

Quercus.over = Quercus.over %>%
  setNames(c('Quercus_robur_range')) 
str(Quercus.over)


## What are the new columns for this dataset? If there is only one, we can't differentiate...
## Change the NAs to "outside". To reassign factor levels, we might need to actually assign the column.
test = Species.over
## test$Betula_pendula_range = levels(test$Betula_pendula_range) <- sub("Betula pendula", "INSIDE_RANGE", levels(test$Betula_pendula_range))
Species.over$Betula_pendula_range = `levels<-`(addNA(Species.over$Betula_pendula_range), c(levels(Species.over$Betula_pendula_range), "OUTSIDE_RANGE"))




unique(Species.over$Betula_pendula_range)


## So using 'over' with this data returns the species name if the point is inside that species range, and NA if the point is 
## outside the range. But would this make the table too big, because you would need a column for each species?
## Potentially, the ranges could be stored as a lookup table, not in the main table. Then a row index could be used to return
## the rows which. 


## So 58367/177219, or ~30% of the points are inside the native species range according to the polygon, and ~70% are outside...
count(Species.over, Betula_pendula_range)[2, 2]/ dim(Betula.pendula)[1] *100


## Do a visual check to see if the "species/NA" split makes sense
GBIF.RANGE = cbind.data.frame(GBIF.RANGE, Species.over)
head(GBIF.RANGE)


## Plot the Betula points inside the range
plot(LAND)
points(GBIF.RANGE[ which(GBIF.RANGE[20] == "Betula pendula"), ][, c("lon", "lat")],
       pch = ".", col = "red",
       cex = 1.3, asp = 1)


## Plot the Betula points outside the range ok, it works
plot(LAND)
points(GBIF.RANGE[ which(GBIF.RANGE[20] == "OUTSIDE_RANGE"), ][, c("lon", "lat")],
       pch = ".", col = "red",
       cex = 1.3, asp = 1)



## So, of this code works, how would I put it in the workflow? I could create a loop or function that would use "sp.n" to generate 
## the out/in column for each species. 

## A dataframe with just these species is probably easier, using the clean records and ignoring the other contextual data.
## This could then be passed to the SDM maxent and prediction code as a separate analysis. So just replace SDM.DATA with the
## native range dataframe...


## So not all analyses will be replicated all they way through. For any species subset, we can do culivated / uncultivated, for
## All scenarios. However, we can't do this for inside / outside native range, except for a small subset of species. These analyses
## Are just for exploration/publication, so probably not crucial to the industry tool.


## Now save .RData file for the next session
save(Species.over, file = paste("./data/base/HIA_LIST/GBIF/GBIF_SPECIES_RANGES.RData"))
save.image("SPECIES_RANGES.RData")





#########################################################################################################################
#####################################################  END ############################################################## 
#########################################################################################################################