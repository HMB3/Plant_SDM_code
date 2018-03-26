#########################################################################################################################
################################################# FLAG SPATIAL OUTLIERS #################################################
#########################################################################################################################


## Use geoclean to Automatcially screen out the dodgy spatial records. See ::
#https://www.rdocumentation.org/packages/speciesgeocodeR/versions/1.0-4/topics/GeoClean
#https://github.com/cran/speciesgeocodeR/blob/master/inst/doc/data_cleaning_and_exploration.Rmd


## Load GBIF data
# install.packages("speciesgeocodeR")
# install.packages("CoordinateCleaner")
# library(speciesgeocodeR)
# library(CoordinateCleaner)
load("GEO_CLEAN_EG.RData")
load("COMBO_RASTER_CONTEXT_1601_2018.RData")


## load a buffered coastline
load("./data/base/CONTEXTUAL/buffland_1deg.rda")
class(buffland)
projection(buffland)
plot(buffland)  ## Raster extract will get rid of these anyway... 


##
GBIF.TRIM.TEST  = COMBO.RASTER.CONTEXT[COMBO.RASTER.CONTEXT$searchTaxon %in% kop.spp, ]
dim(GBIF.TRIM.TEST)


## Quickly check the data
GBIF.TRIM.TEST$TAXON = as.character(GBIF.TRIM.TEST$TAXON)
dim(GBIF.TRIM.TEST)
names(GBIF.TRIM.TEST)
unique(GBIF.TRIM.TEST$TAXON)


## Do an example species
#GBIF.TRIM.TEST = subset(GBIF.TRIM.TEST, TAXON == "Ficus brachypoda")
dim(GBIF.TRIM.TEST)
unique(GBIF.TRIM.TEST$TAXON)


#########################################################################################################################
## 'Geoclean' provides several different tests to clean datasets with geographic coordinates...
GBIF.TRIM.GEO = dplyr::rename(GBIF.TRIM.TEST, identifier = TAXON,
                              XCOOR = LON, YCOOR = LAT, country = COUNTRY)


## Perform some simple test for coordinate validity. Each argument of GeoClean represents one test, 
## See ?GeoClean for details. Creates a GeoClean table with columms for each data check
tt.long <- GeoClean(GBIF.TRIM.GEO, isna = TRUE, isnumeric = TRUE, coordinatevalidity = TRUE,
                    containszero = TRUE, zerozero = TRUE, zerozerothresh = 1, 
                    latequallong = TRUE, outp = "detailed")


## AZ :: now all records pass, you seem to have some high quality data... 
mean(tt.long$summary) 


# AZ I have further devloped the cleaning functions in the CoordinateCleaner package, which runs additional tests.
# I suggest you do the following instead
GBIF.TRIM.GEO = dplyr::rename(GBIF.TRIM.TEST, 
                              species = TAXON,
                              decimallongitude = LON, 
                              decimallatitude = LAT, 
                              country = COUNTRY)


## Flag records for default problems :: ?CleanCoordinates :: ?CleanCoordinates
tt.long <- CleanCoordinates(GBIF.TRIM.GEO)#, seas.ref = buffland)


## Takes a little while for the seas test, due to the global extent. Flagging 1%, that seems reasonable
summary(tt.long)
summary(tt.long)[10]/dim(tt.long)[1]*100


## Plot this for each species to double check as a larger file, how to see the icons?
plot(tt.long)


## That looks much better. Depending on your exact study question/area, I suggest to possibly use a buffered 
## gazetter for the seas test, to avoid falsely flagging records on the coastline
## you can find one here: https://github.com/azizka/CoordinateCleaner/tree/master/extra_gazetteers
## you can also tweak the buffer area around capital and centroids using the function arguments


# Also check out my newest data cleaning tutorial for GBIF data:
# https://github.com/azizka/CoordinateCleaner/tree/master/Tutorials


## Now the outlier test (you can also run it through CleanCoordinates, but unfortunately it runs for awhile...)
## You might want to consider some tweaking of the mltpl paramters
## Why is the length of the value vector less than the 
?cc_outl
GBIF.CLEAN    <- cc_outl(GBIF.TRIM.GEO, value = "clean", mltpl = 6)
GBIF.SPAT.OUT <- cc_outl(GBIF.TRIM.GEO, value = "flags", mltpl = 6)


## Check the output ::
TEST.CLEAN = completeFun(GBIF.CLEAN, "species")
TEST.SPAT  = GBIF.SPAT.OUT[!is.na(GBIF.SPAT.OUT)]


## 
dim(TEST.CLEAN)
dim(GBIF.TRIM.GEO)[1] - dim(TEST.CLEAN)[1]

unique(TEST.SPAT)
length(TEST.SPAT)


## Create a new column for the spatial outliers
TEST.CLEAN$SPAT_OUT = "TRUE"
SPAT.OUT = TEST.CLEAN[c("OBS", "SPAT_OUT")]



#########################################################################################################################
## Now bind the GBIF data with the cleaning test columns
test.geo = cbind(GBIF.TRIM.TEST, tt.long)


## Check these columns
names(test.geo)
names(TEST.CLEAN)


## Try joining them
test.clean          = join(test.geo, SPAT.OUT)
test.clean$SPAT_OUT = ifelse(is.na(test.clean$SPAT_OUT), 
                             'FALSE', test.clean$SPAT_OUT)
unique(test.clean$SPAT_OUT)
length(test.clean$SPAT_OUT)


#########################################################################################################################
## Map the points, are they really invalid?
GBIF.POINTS = SpatialPointsDataFrame(coords      = test.clean [c("LON", "LAT")], 
                                     data        = test.clean ,
                                     proj4string = CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))


## Create list of species 
TAXA <- unique(as.character(test.geo$TAXON))


#########################################################################################################################
## Then, loop over the species list and plot them
for (i in 1:length(TAXA)) {
  
  ## Need to check the OBS column matches up - or do we not need this again?
  spp.points <- GBIF.POINTS[GBIF.POINTS$TAXON == TAXA[i], ] 
  plot(aus, main = TAXA[i])
  points(spp.points$LON, spp.points$LAT, col = "red", cex = .6, pch = 19)
  
}


## EG Ficus brachypoda (no. 6 in the plots) ::
ficus.clean = subset (GBIF.POINTS, TAXON == "Ficus brachypoda")


## Save files
save(test.geo, file = paste("./data/base/HIA_LIST/COMBO/GBIF_coordclean.RData"))
writeOGR(obj = GBIF.POINTS, dsn = "./data/base/HIA_LIST/COMBO/CLEAN_GBIF", 
         layer = "GBIF.DIANELLA.CLEAN", driver = "ESRI Shapefile")




#########################################################################################################################
#####################################################  TBC ##############################################################
#########################################################################################################################