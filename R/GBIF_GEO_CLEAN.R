#########################################################################################################################
################################################# FLAG SPATIAL OUTLIERS #################################################
#########################################################################################################################


## Use geoclean to Automatcially screen out the dodgy spatial records. See ::
#https://www.rdocumentation.org/packages/speciesgeocodeR/versions/1.0-4/topics/GeoClean
#https://github.com/cran/speciesgeocodeR/blob/master/inst/doc/data_cleaning_and_exploration.Rmd


## Load GBIF data
load("GEO_CLEAN_EG.RData")
# library(speciesgeocodeR)


## Subset to only the 10 test species
# GBIF.TRIM.TEST = subset(GBIF.TRIM, 
#                          searchTaxon == "Backhousia citriodora" |
#                          searchTaxon == "Calodendrum capense" |
#                          searchTaxon == "Cordyline australis" |
#                          searchTaxon == "Dianella caerulea" |
#                          searchTaxon == "Eucalyptus erythrocorys" |
#                          searchTaxon == "Ficus brachypoda" |
#                          searchTaxon == "Liriope muscari" |
#                          searchTaxon == "Lomandra longifolia" |
#                          searchTaxon == "Murraya paniculata" |
#                          searchTaxon == "Ulmus parvifolia")


## Quickly check the data
dim(GBIF.TRIM.TEST)
names(GBIF.TRIM.TEST)
unique(GBIF.TRIM.TEST$TAXON)


#########################################################################################################################
## 'Geoclean' provides several different tests to clean datasets with geographic coordinates...
GBIF.TRIM.GEO = dplyr::rename(GBIF.TRIM.TEST, identifier = TAXON,
                              XCOOR = LAT, YCOOR = LON, country = COUNTRY)


## Perform some simple test for coordinate validity. Each argument of GeoClean represents one test, 
## see ?GeoClean for details. Creates a GeoClean table with columms for each data check
tt.long <- GeoClean(GBIF.TRIM.GEO, isna = TRUE, isnumeric = TRUE, coordinatevalidity = TRUE,
                    containszero = TRUE, zerozero = TRUE, zerozerothresh = 1, latequallong = TRUE, outp = "detailed")
tt.long = cbind(GBIF.TRIM.TEST["OBS"], tt.long)


## So most of the data have invalid coordinates?
View(tt.long)


## Or, create a vectorm with FALSE if ANY of the tests fail...
tt <- GeoClean(GBIF.TRIM.GEO, isna = TRUE, isnumeric = TRUE, coordinatevalidity = TRUE,
               containszero = TRUE, zerozero = TRUE, zerozerothresh = 1, latequallong = TRUE)


## So 99.5% of the data fails the test? Seems harsh...
simple.clean <- subset(GBIF.TRIM.GEO, tt == TRUE)
dim(simple.clean)[1]/dim(GBIF.TRIM.GEO)[1]*100


#########################################################################################################################
## Map the points, are they really invalid?
GBIF.POINTS = SpatialPointsDataFrame(coords      = GBIF.TRIM.TEST[c("LON", "LAT")], 
                                     data        = GBIF.TRIM.TEST,
                                     proj4string = CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))


## Create list of species 
TAXA <- unique(as.character(GBIF.POINTS$TAXON))


#########################################################################################################################
## Then, loop over the species list and plot them
for (i in 1:length(TAXA)) {
  
  ## Need to check the OBS column matches up - or do we not need this again?
  spp.points <- GBIF.POINTS[GBIF.POINTS$TAXON == TAXA[i], ] 
  plot(aus, main = TAXA[i])
  points(spp.points$LON, spp.points$LAT, col = "red", cex = .6, pch = 19)
  
}


## EG Ficus brachypoda (no. 6 in the plots) :: most of these records are valid?
ficus.clean = subset (tt.long, identifier == "Ficus brachypoda")
View(ficus.clean)


##save.image("GEO_CLEAN_EG.RData")


#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################