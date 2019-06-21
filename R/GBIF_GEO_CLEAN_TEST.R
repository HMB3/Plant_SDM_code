#########################################################################################################################
################################################# FLAG SPATIAL OUTLIERS #################################################
#########################################################################################################################


## Use geoclean to Automatcially screen out the dodgy spatial records. See ::
#https://www.rdocumentation.org/packages/speciesgeocodeR/versions/1.0-4/topics/GeoClean
#https://github.com/cran/speciesgeocodeR/blob/master/inst/doc/data_cleaning_and_exploration.Rmd


#########################################################################################################################
## 1). CHECK DATA FOR AN EXAMPLE SPECIES
#########################################################################################################################


## Load GBIF data
# install.packages("speciesgeocodeR")
# install.packages("CoordinateCleaner")
# library(speciesgeocodeR)
# library(CoordinateCleaner)
#load("GEO_CLEAN_EG.RData")
load("./data/base/HIA_LIST/COMBO/COMBO_RASTER_CONTEXT_1601_2018.RData")
load("./data/base/CONTEXTUAL/urbanareas.rda")
LAND = readRDS("F:/green_cities_sdm/data/base/CONTEXTUAL/LAND_world.rds")
aus  = readRDS("F:/green_cities_sdm/data/base/CONTEXTUAL/aus_states.rds")
names(COMBO.RASTER.CONTEXT)


## HB Do an example species :: Ficus Brachypoda
GBIF.TRIM.TEST  = COMBO.RASTER.CONTEXT#[COMBO.RASTER.CONTEXT$searchTaxon %in% kop.spp, ]
GBIF.TRIM.TEST = subset(GBIF.TRIM.TEST, searchTaxon == "Ficus brachypoda")
unique(GBIF.TRIM.TEST$searchTaxon)


## HB Hard to see - but we'd love to remove the point in Madagascar 
plot(LAND)
points(GBIF.TRIM.TEST$lon,  GBIF.TRIM.TEST$lat, cex = 0.2, col = "red", pch = 19)


## HB and in Australia, we'd love to remove the points on the east coast and north QLD (two in top right), 
## so it looks more like this :: https://bie.ala.org.au/species/http://id.biodiversity.org.au/node/apni/2897559
plot(aus)
points(GBIF.TRIM.TEST$lon,  GBIF.TRIM.TEST$lat, cex = 0.8, col = "red", pch = 19)


## HB The trouble comes from combining GBIF data with data from the Atlas of Living Australia (ALA). This is needed
## because we're modelling global niches for a mix of Australian natives and exotics. But it causes headaches! 





#########################################################################################################################
## 2). FLAG GBIF AND SPATIAL OUTLIERS
#########################################################################################################################


#########################################################################################################################
## AZ I have further devloped the cleaning functions in the CoordinateCleaner package, which runs additional tests.
## I suggest you do the following instead
GBIF.TRIM.GEO = dplyr::rename(GBIF.TRIM.TEST, 
                              species = searchTaxon,
                              decimallongitude = lon, 
                              decimallatitude  = lat)
GBIF.TRIM.GEO = GBIF.TRIM.GEO[c("species", "decimallongitude", "decimallatitude", "country")]


## HB I'm haivng trouble getting the cc_outl function working by itself 
## 427 duplicates looks too many, but 8 geographic outliers looks about right
tt.long <- CleanCoordinates(GBIF.TRIM.GEO,
                            countrycheck     = TRUE,
                            duplicates       = TRUE,
                            outliers         = TRUE,
                            outliers.method  = "quantile",
                            outliers.mtp     = 5,
                            seas             = FALSE)

## However, the outliers test causes an unequal number of rows in the output, not sure why...
# Error in data.frame(inp, validity = val, equal = equ, zeros = zer, capitals = cap,  : 
#                       arguments imply differing number of rows: 2023, 2883295

## This works, but only flags records in the sea
tt.long <- CleanCoordinates(GBIF.TRIM.GEO)


## Takes a little while for the seas test, due to the global extent. Flagging 1%
summary(tt.long)
summary(tt.long)[10]/dim(tt.long)[1]*100


## Plot this for each species to double check as a larger file, how to see the icons?
plot(tt.long)


## AZ That looks much better. Depending on your exact study question/area, I suggest to possibly use a buffered 
## gazetter for the seas test, to avoid falsely flagging records on the coastline
## you can find one here: https://github.com/azizka/CoordinateCleaner/tree/master/extra_gazetteers
## you can also tweak the buffer area around capital and centroids using the function arguments


## HB All the records in this dataset are inside a worldclim pixel, so I can probably trust the records...

# Also check out my newest data cleaning tutorial for GBIF data:
# https://github.com/azizka/CoordinateCleaner/tree/master/Tutorials


## Now the outlier test (you can also run it through CleanCoordinates, but unfortunately it runs for awhile...)
## You might want to consider some tweaking of the mltpl paramters
## Why is the length of the value vector less than the 
?cc_outl


## This creates a data frame of just the clean records
## Not sure where the extra rows are coming from....???
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
# test.clean          = join(test.geo, SPAT.OUT)
# test.clean$SPAT_OUT = ifelse(is.na(test.clean$SPAT_OUT), 
#                              'FALSE', test.clean$SPAT_OUT)
# unique(test.clean$SPAT_OUT)
# length(test.clean$SPAT_OUT)



#########################################################################################################################
## 3). MAP POINTS
#########################################################################################################################


#########################################################################################################################
## Map the points, are they really invalid?
GBIF.POINTS = SpatialPointsDataFrame(coords      = test.geo[c("LON", "LAT")], 
                                     data        = test.geo,
                                     proj4string = CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))


## Create list of species 
TAXA <- unique(as.character(test.geo$searchTaxon))


#########################################################################################################################
## Then, loop over the species list and plot them
for (i in 1:length(TAXA)) {
  
  ## Need to check the OBS column matches up - or do we not need this again?
  spp.points <- GBIF.POINTS[GBIF.POINTS$searchTaxon == TAXA[i], ] 
  plot(aus, main = TAXA[i])
  points(spp.points$lon, spp.points$lat, col = "red", cex = .6, pch = 19)
  
}


## EG Ficus brachypoda (no. 6 in the plots) ::
ficus.clean = subset (GBIF.POINTS, searchTaxon == "Ficus brachypoda")





#########################################################################################################################
## 3). SUBSET TO THE CLEAN RECORDS
#########################################################################################################################


## Just get the false records
GBIF.FALSE  = subset(test.geo, summary == "FALSE")# | GBIF.SPAT.OUT == "FALSE")


## Check one species
coordyline = subset(GBIF.FALSE, searchTaxon == "Cordyline australis")
View(subset(coordyline, summary == "FALSE")# | GBIF.SPAT.OUT == "FALSE")
     [, c("searchTaxon", 
          #"OBS",
          "decimallongitude",
          "decimallatitude",
          "validity", 
          "equal",
          "zeros",
          "capitals",
          "centroids",
          "sea",
          "gbif",
          "institution",
          "summary")])


## Not sure why the inverse did not work :: get only the records which were not flagged as being dodgy.
## dim(subset(RASTER.CLEAN, summary == "TRUE" | GBIF.SPAT.OUT == "TRUE"))
#GBIF.TRUE = test.geo[!test.geo$OBS %in% GBIF.FALSE$OBS, ]
GBIF.TRUE  = subset(test.geo, summary == "TRUE")
unique(GBIF.TRUE$summary)                                              ## works
dim(GBIF.FALSE)[1]/dim(test.geo)[1]*100                                ## 1.3% of records excluded



#########################################################################################################################
## Save the points
COMBO.RASTER.CONTEXT = GBIF.TRUE
save(COMBO.RASTER.CONTEXT, file = paste("./data/base/HIA_LIST/COMBO/COMBO_GBIF_TRUE_TEST_SPP.RData", sep = ""))


## Save files
save(test.geo, file = paste("./data/base/HIA_LIST/COMBO/GBIF_coordclean.RData"))
writeOGR(obj = GBIF.POINTS, dsn = "./data/base/HIA_LIST/COMBO/CLEAN_GBIF", 
         layer = "GBIF.DIANELLA.CLEAN", driver = "ESRI Shapefile")




#########################################################################################################################
#####################################################  TBC ##############################################################
#########################################################################################################################