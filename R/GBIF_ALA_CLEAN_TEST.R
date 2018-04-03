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


## HB Do an example species :: SPP Brachypoda
GBIF.TRIM.TEST  = COMBO.RASTER.CONTEXT[COMBO.RASTER.CONTEXT$searchTaxon %in% spp.all, ]   ## kop.spp for the test spp
SPP.TEST        = subset(GBIF.TRIM.TEST, searchTaxon == "Ficus brachypoda")

str(unique(GBIF.TRIM.TEST$searchTaxon))
unique(SPP.TEST$searchTaxon)
'Ficus brachypoda' %in% GBIF.TRIM.TEST$searchTaxon  


## HB hard to see - but we'd love to remove the point in Madagascar 
plot(LAND)
points(SPP.TEST$lon,  SPP.TEST$lat, cex = 0.2, col = "red", pch = 19)


#########################################################################################################################
## HB and in Australia, we'd love to remove the points on the east coast and north QLD (two in top right), 
## so it looks more like this :: https://bie.ala.org.au/species/http://id.biodiversity.org.au/node/apni/2897559
plot(aus)
points(SPP.TEST$lon,  SPP.TEST$lat, cex = 0.8, col = "red", pch = 19)


## HB The trouble comes from combining GBIF data with data from the Atlas of Living Australia (ALA). This is needed
## because we're modelling global niches for a mix of Australian natives and exotics. But it causes headaches! 





#########################################################################################################################
## FLAG GBIF AND SPATIAL OUTLIERS FOR ONE SPECIES
#########################################################################################################################


#########################################################################################################################
## AZ I have further devloped the cleaning functions in the CoordinateCleaner package, which runs additional tests.
## I suggest you do the following instead
GBIF.TRIM.GEO = dplyr::rename(SPP.TEST, 
                              species = searchTaxon,
                              decimallongitude = lon, 
                              decimallatitude  = lat)


## Now use ?CleanCoordinates :: need to work on a tibble, for some reason !
TIB.TEST <- as_tibble(GBIF.TRIM.GEO)
TIB.TEST


## Run a check :: I've already stripped out the records that fall outside
## the worldclim raster boundaries, so the sea test is probably not the most important for me
## Study area is the globe, but we are only projecting models onto Australia
FLAGS  <- CleanCoordinates(TIB.TEST,
                           countrycheck     = TRUE,
                           duplicates       = TRUE,
                           outliers         = TRUE,
                           outliers.method  = "quantile",
                           outliers.mtp     = 5,
                           outliers.td      = 1000,
                           seas             = FALSE)


## Flagging < 1%, but that's just coastal records for this species - pretty dodgy
summary(FLAGS)
summary(FLAGS)[10]/dim(FLAGS)[1]*100


## A plot like this for each species would be awesome, fantastic work!
plot(FLAGS)


## AZ That looks much better. Depending on your exact study question/area, I suggest to possibly use a buffered 
## gazetter for the seas test, to avoid falsely flagging records on the coastline
## you can find one here: https://github.com/azizka/CoordinateCleaner/tree/master/extra_gazetteers
## you can also tweak the buffer area around capital and centroids using the function arguments


## HB All the records in this dataset are inside a worldclim pixel, so prob not too worried about the sea..


# Also check out my newest data cleaning tutorial for GBIF data:
# https://github.com/azizka/CoordinateCleaner/tree/master/Tutorials


## Looks good! You need to convert the data frame to a tibble to make it work...


## Now the outlier test (you can also run it through CleanCoordinates, but unfortunately it runs for awhile...)
## You might want to consider some tweaking of the mltpl paramters
?cc_outl


## This creates a data frame of just the clean records
## Not sure where the extra rows are coming from.....???
GBIF.CLEAN    <- cc_outl(TIB.TEST, 
                         lon     = "decimallongitude", 
                         lat     = "decimallatitude", 
                         species = "species", 
                         method  = "quantile", mltpl = 5, tdi = 1000, 
                         value   = "clean",  verbose = T)


## This creates a vector of true/false for each record
GBIF.SPAT.OUT <- cc_outl(TIB.TEST, 
                         lon     = "decimallongitude", 
                         lat     = "decimallatitude", 
                         species = "species", 
                         method  = "quantile", mltpl = 5, tdi = 1000, 
                         value   = "FLAGS",  verbose = T)


## Check the output ::
dim(GBIF.CLEAN);length(GBIF.SPAT.OUT)


## Join data
TEST.GEO = cbind(TIB.TEST, FLAGS, GBIF.SPAT.OUT)
class(TEST.GEO)
names(TEST.GEO)
dim(TEST.GEO)


## Check the values for each field
unique(TEST.GEO$validity)
unique(TEST.GEO$equal)
unique(TEST.GEO$equal)
unique(TEST.GEO$zeros)
unique(TEST.GEO$capitals)
unique(TEST.GEO$centroids)
unique(TEST.GEO$outliers)
unique(TEST.GEO$gbif)
unique(TEST.GEO$institution)
unique(TEST.GEO$outliers)
unique(TEST.GEO$gbif)
unique(TEST.GEO$summary)
unique(TEST.GEO$GBIF.SPAT.OUT)





#########################################################################################################################
## PLOT CLEAN DATA
#########################################################################################################################


## HB I haven't removed any spatial outliers yet...
plot(LAND)
points(GBIF.CLEAN$decimallongitude,  GBIF.CLEAN$decimallatitude, cex = 0.2, col = "red", pch = 19)


## Still have the outliers on the eastern half of Australia
plot(aus)
points(GBIF.CLEAN$decimallongitude,  GBIF.CLEAN$decimallatitude, cex = 0.8, col = "red", pch = 19)





#########################################################################################################################
## 2). NOW, CLEAN WHOLE DATASET
#########################################################################################################################


#########################################################################################################################
## FLAG GBIF AND SPATIAL OUTLIERS FOR ONE SPECIES
#########################################################################################################################


#########################################################################################################################
## AZ I have further devloped the cleaning functions in the CoordinateCleaner package, which runs additional tests.
## I suggest you do the following instead
str(unique(GBIF.TRIM.TEST$searchTaxon))
GBIF.TRIM.GEO = dplyr::rename(GBIF.TRIM.TEST, 
                              species = searchTaxon,
                              decimallongitude = lon, 
                              decimallatitude  = lat)


## Now use ?CleanCoordinates :: need to work on a tibble, for some reason !
TIB.TEST <- as_tibble(GBIF.TRIM.GEO)
dim(TIB.TEST)
str(unique(TIB.TEST$species))


## Run a check :: I've already stripped out the records that fall outside
## the worldclim raster boundaries, so the sea test is probably not the most important for me
## Study area is the globe, but we are only projecting models onto Australia

#########################################################################################################################
## Don't run the outliers test here, it is slower
FLAGS  <- CleanCoordinates(TIB.TEST,
                           #countries        = "country",    ## too many flagged here...
                           countrycheck     = TRUE,
                           duplicates       = TRUE,
                           #outliers         = TRUE,
                           # outliers.method  = "quantile",
                           # outliers.mtp     = 5,
                           # outliers.td      = 1000,
                           seas             = FALSE)


## Flagging < 2.7 % seems reasonable
summary(FLAGS)
summary(FLAGS)[8]/dim(FLAGS)[1]*100
FLAGS = FLAGS[ ,!(colnames(FLAGS) == "decimallongitude" | colnames(FLAGS) =="decimallatitude")]


## A plot like this for each species would be awesome, fantastic work!
#plot(FLAGS)


#########################################################################################################################
## This creates a vector of true/false for each record
GBIF.SPAT.OUT.1 <- cc_outl(head(TIB.TEST, (dim(TIB.TEST)[1]/2)), 
                         lon     = "decimallongitude", 
                         lat     = "decimallatitude", 
                         species = "species", 
                         method  = "quantile", 
                         mltpl   = 5, 
                         tdi     = 1000, 
                         value   = "flags")#, verbose = TRUE

GBIF.SPAT.OUT.2 <- cc_outl(tail(TIB.TEST, (dim(TIB.TEST)[1]/2)), 
                           lon     = "decimallongitude", 
                           lat     = "decimallatitude", 
                           species = "species", 
                           method  = "quantile", 
                           mltpl   = 5, 
                           tdi     = 1000, 
                           value   = "flags")#, verbose = TRUE)


## Check the output ::
#dim(FLAGS);length(GBIF.SPAT.OUT)
FLAGS = bind_rows(GBIF.SPAT.OUT.1, GBIF.SPAT.OUT.2 )


#########################################################################################################################
## Join data :: exclude the decimal lat/long 
GBIF.SPAT.OUT = cbind(GBIF.TRIM.TEST, FLAGS, GBIF.SPAT.OUT)
identical(TEST.GEO$searchTaxon, TEST.GEO$species)                              ## order matches


## Check the values for each flag
summary(TEST.GEO$validity)
summary(TEST.GEO$equal)
summary(TEST.GEO$equal)
summary(TEST.GEO$zeros)
summary(TEST.GEO$capitals)
summary(TEST.GEO$centroids)
summary(TEST.GEO$gbif)
summary(TEST.GEO$institution)
summary(TEST.GEO$gbif)
summary(TEST.GEO$summary)
summary(TEST.GEO$GBIF.SPAT.OUT)





#########################################################################################################################
## PLOT CLEAN DATA
#########################################################################################################################


# plot(LAND)
# points(GBIF.CLEAN$decimallongitude,  GBIF.CLEAN$decimallatitude, cex = 0.2, col = "red", pch = 19)


## Still have the outliers on the eastern half of Australia
# plot(aus)
# points(GBIF.CLEAN$decimallongitude,  GBIF.CLEAN$decimallatitude, cex = 0.8, col = "red", pch = 19)





#########################################################################################################################
## SUBSET AND SAVE TABLES
#########################################################################################################################


## This unique ID column can be applied across the project  
TEST.GEO$OBS <- 1:nrow(TEST.GEO)
dim(TEST.GEO)[1];length(TEST.GEO$OBS)  


## So ~0.8% of the data is dodgy according to the GBIF fields or spatial outliers
## This seems ok as a median figure across the data set?
dim(subset(TEST.GEO, summary == "FALSE" | GBIF.SPAT.OUT == "FALSE"))[1]/dim(TEST.GEO)[1]*100
CLEAN.FALSE = subset(TEST.GEO, summary == "FALSE" | GBIF.SPAT.OUT == "FALSE")


## Check one species
coordyline = subset(TEST.GEO, searchTaxon == "Cordyline australis")
View(subset(coordyline, summary == "FALSE" | GBIF.SPAT.OUT == "FALSE")[, c("searchTaxon", "OBS",
                                                                           "lon",
                                                                           "lat",
                                                                           "validity", 
                                                                           "equal",
                                                                           "zeros",
                                                                           "capitals",
                                                                           "centroids",
                                                                           "duplicates",
                                                                           "gbif",
                                                                           "institution",
                                                                           "summary",
                                                                           "GBIF.SPAT.OUT")])



## Not sure why the inverse did not work :: get only the records which were not flagged as being dodgy.
## dim(subset(TEST.GEO, summary == "TRUE" | GBIF.SPAT.OUT == "TRUE"))
CLEAN.TRUE = TEST.GEO[!TEST.GEO$OBS %in% CLEAN.FALSE$OBS, ]
dim(CLEAN.TRUE)
unique(CLEAN.TRUE$summary)                                             ## works
unique(CLEAN.TRUE$GBIF.SPAT.OUT)                                       ## works


## How many species?
str(unique(CLEAN.TRUE$searchTaxon))
str(unique(TEST.GEO$searchTaxon))
'Ficus brachypoda' %in% TEST.GEO$searchTaxon  


#########################################################################################################################
##    
saveRDS(TEST.GEO,   'data/base/HIA_LIST/COMBO/CLEAN_FLAGS_HIA_SPP.rds')
saveRDS(CLEAN.TRUE, 'data/base/HIA_LIST/COMBO/CLEAN_ONLY_HIA_SPP.rds')





#########################################################################################################################
#####################################################  TBC ##############################################################
#########################################################################################################################