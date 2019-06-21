#########################################################################################################################
################################################# FLAG SPATIAL OUTLIERS #################################################
#########################################################################################################################


## Use geoclean to Automatcially screen out the dodgy spatial records. See ::
#https://www.rdocumentation.org/packages/speciesgeocodeR/versions/1.0-4/topics/GeoClean
#https://github.com/cran/speciesgeocodeR/blob/master/inst/doc/data_cleaning_and_exploration.Rmd


## Can we flag herbarium records, plus the spatial outliers?


#########################################################################################################################
## 1). CHECKK DATA FOR AN EXAMPLE SPECIES
#########################################################################################################################


## Load GBIF data
load("./data/base/HIA_LIST/COMBO/COMBO_RASTER_CONTEXT_1601_2018.RData")
load("./data/base/CONTEXTUAL/urbanareas.rda")
LAND = readRDS("F:/green_cities_sdm/data/base/CONTEXTUAL/LAND_world.rds")
aus  = readRDS("F:/green_cities_sdm/data/base/CONTEXTUAL/aus_states.rds")
names(COMBO.RASTER.CONTEXT)


## Create lists
source('./R/HIA_LIST_MATCHING.R')


## HB Do an example species :: Ficus Brachypoda
GBIF.TRIM.TEST  = COMBO.RASTER.CONTEXT[COMBO.RASTER.CONTEXT$searchTaxon %in% kop.spp, ]
FICUS.TEST      = subset(GBIF.TRIM.TEST, searchTaxon == "Ficus brachypoda")

unique(GBIF.TRIM.TEST$searchTaxon)
unique(FICUS.TEST$searchTaxon)


## HB Hard to see - but we'd love to remove the point in Madagascar 
plot(LAND)
points(FICUS.TEST$lon,  FICUS.TEST$lat, cex = 0.2, col = "red", pch = 19)


#########################################################################################################################
## HB and in Australia, we'd love to remove the points on the east coast and north QLD (two in top right), 
## so it looks more like this :: https://bie.ala.org.au/species/http://id.biodiversity.org.au/node/apni/2897559
plot(aus)
points(FICUS.TEST$lon,  FICUS.TEST$lat, cex = 0.8, col = "red", pch = 19)


## HB The trouble comes from combining GBIF data with data from the Atlas of Living Australia (ALA). This is needed
## because we're modelling global niches for a mix of Australian natives and exotics. But it causes headaches! 





#########################################################################################################################
## 2). FLAG GBIF AND SPATIAL OUTLIERS FOR ONE SPECIES
#########################################################################################################################


#########################################################################################################################
## AZ I have further devloped the cleaning functions in the CoordinateCleaner package, which runs additional tests.
## I suggest you do the following instead
GBIF.TRIM.GEO = dplyr::rename(FICUS.TEST, 
                              species = searchTaxon,
                              decimallongitude = lon, 
                              decimallatitude  = lat)


## Now use ?CleanCoordinates :: need to work on a tibble, for some reason !
TIB.TEST <- as_tibble(GBIF.TRIM.GEO)
TIB.TEST


## Run a check :: I've already stripped out the records that fall outside
## the worldclim raster boundaries, so the sea test is probably not the most important for me
## Study area is the globe, but we are only projecting models onto Australia
flags  <- CleanCoordinates(TIB.TEST,
                           countrycheck     = TRUE,
                           duplicates       = TRUE,
                           outliers         = TRUE,
                           outliers.method  = "quantile",
                           outliers.mtp     = 5,
                           outliers.td      = 1000,
                           seas             = FALSE)

## Flagging < 1%, but that's just coastal records for this species - pretty dodgy
summary(flags)
summary(flags)[10]/dim(flags)[1]*100


## A plot like this for each species would be awesome, fantastic work!
plot(flags)


## AZ That looks much better. Depending on your exact study question/area, I suggest to possibly use a buffered 
## gazetter for the seas test, to avoid falsely flagging records on the coastline
## you can find one here: https://github.com/azizka/CoordinateCleaner/tree/master/extra_gazetteers
## you can also tweak the buffer area around capital and centroids using the function arguments


## HB All the records in this dataset are inside a worldclim pixel, so prob not too worried about the sea..


# Also check out my newest data cleaning tutorial for GBIF data:
# https://github.com/azizka/CoordinateCleaner/tree/master/Tutorials


## Looks good!


## Now the outlier test (you can also run it through CleanCoordinates, but unfortunately it runs for awhile...)
## You might want to consider some tweaking of the mltpl paramters
?cc_outl


## This creates a data frame of just the clean records
## Not sure where the extra rows are coming from.....???
GBIF.CLEAN    <- cc_outl(GBIF.TRIM.GEO, 
                         lon     = "decimallongitude", 
                         lat     = "decimallatitude", 
                         species = "species", 
                         method  = "quantile", mltpl = 5, tdi = 1000, 
                         value   = "clean",  verbose = T)


## This creates a vector of true/false for each record
GBIF.SPAT.OUT <- cc_outl(GBIF.TRIM.GEO, 
                         lon     = "decimallongitude", 
                         lat     = "decimallatitude", 
                         species = "species", 
                         method  = "quantile", mltpl = 5, tdi = 1000, 
                         value   = "flags",  verbose = T)


## Check the output ::
TEST.CLEAN = completeFun(GBIF.CLEAN, "species")
TEST.SPAT  = GBIF.SPAT.OUT[!is.na(GBIF.SPAT.OUT)]


## So the outlier test is not working in this form....
dim(TEST.CLEAN)
dim(GBIF.TRIM.GEO)[1] - dim(TEST.CLEAN)[1]

unique(TEST.SPAT)
length(TEST.SPAT)





#########################################################################################################################
## 3). PLOT CLEAN DATA
#########################################################################################################################


## HB I haven't removed any spatial outliers yet...
plot(LAND)
points(TEST.CLEAN$decimallongitude,  TEST.CLEAN$decimallatitude, cex = 0.2, col = "red", pch = 19)


## Still have the outliers on the eastern half of Australia
plot(aus)
points(TEST.CLEAN$decimallongitude,  TEST.CLEAN$decimallatitude, cex = 0.8, col = "red", pch = 19)





#########################################################################################################################
## 4). CLEAN WHOLE DATASET
#########################################################################################################################


#########################################################################################################################
#####################################################  TBC ##############################################################
#########################################################################################################################