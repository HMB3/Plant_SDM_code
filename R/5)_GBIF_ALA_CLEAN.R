#########################################################################################################################
################################################# FLAG SPATIAL OUTLIERS #################################################
#########################################################################################################################


## This code uses the  to Automatcially screen out the dodgy spatial records. See ::
#https://www.rdocumentation.org/packages/speciesgeocodeR/versions/1.0-4/topics/GeoClean
#https://github.com/cran/speciesgeocodeR/blob/master/inst/doc/data_cleaning_and_exploration.Rmd


## Can we flag herbarium records, plus the spatial outliers?


#########################################################################################################################
## Read in all data to run the SDM code :: species lists, shapefile, rasters & tables
#source('./R/HIA_LIST_MATCHING.R')
rasterTmpFile()


#########################################################################################################################
## 1). CHECK DATA FOR AN EXAMPLE SPECIES...
#########################################################################################################################


## Load GBIF data
# COMBO.RASTER.CONVERT = readRDS("./data/base/HIA_LIST/COMBO/COMBO_RASTER_CONTEXT_APRIL_2018.rds")
# COMBO.NICHE.CONTEXT  = readRDS("./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT_APRIL_2018_STANDARD_CLEAN.rds")
names(COMBO.RASTER.CONVERT)


## Restrict ALA data to just those species on the big list
HIA.SPP.JOIN     = CLEAN.SPP
HIA.SPP.JOIN     = dplyr::rename(HIA.SPP.JOIN, searchTaxon = Binomial)


## Set NA to blank, then sort by no. of growers to get them to the top
HIA.SPP.JOIN[is.na(HIA.SPP.JOIN)] <- 0
HIA.SPP.JOIN = HIA.SPP.JOIN[with(HIA.SPP.JOIN, rev(order(Number.of.growers))), ]
head(HIA.SPP.JOIN[, c("searchTaxon", "Number.of.growers")])
View(HIA.SPP.JOIN)


## HB do an example species :: SPP Brachypoda
GBIF.TRIM.TEST  = COMBO.RASTER.CONVERT#[COMBO.RASTER.CONVERT$searchTaxon %in% spp.all, ]      
SPP.TEST        = subset(GBIF.TRIM.TEST, searchTaxon == "Ficus hillii")                       


## Does it work?
str(unique(GBIF.TRIM.TEST$searchTaxon))
unique(SPP.TEST$searchTaxon)
'Cordyline australis' %in% GBIF.TRIM.TEST$searchTaxon  


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


# #########################################################################################################################
# ## AZ I have further devloped the cleaning functions in the CoordinateCleaner package, which runs additional tests.
# ## I suggest you do the following instead
# GBIF.TRIM.GEO = dplyr::rename(SPP.TEST,
#                               species = searchTaxon,
#                               decimallongitude = lon,
#                               decimallatitude  = lat)
# 
# 
# ## Now use ?CleanCoordinates :: need to work on a tibble, for some reason !
# TIB.TEST <- timetk::tk_tbl(GBIF.TRIM.GEO)
# TIB.TEST
# 
# 
# ## Run a check :: I've already stripped out the records that fall outside
# ## the worldclim raster boundaries, so the sea test is probably not the most important for me
# ## Study area is the globe, but we are only projecting models onto Australia
# FLAGS  <- CleanCoordinates(TIB.TEST,
#                            #countries        = "country",    ## too many flagged here...
#                            countrycheck     = TRUE,
#                            duplicates       = TRUE,
#                            capitals.rad     = 0.12,
#                            seas             = FALSE)
# 
# 
# ## Flagging < 1%, but that's just coastal records for this species - pretty dodgy
# summary(FLAGS)
# summary(FLAGS)[10]/dim(FLAGS)[1]*100
# 
# 
# ## A plot like this for each species would be awesome, fantastic work!
# plot(FLAGS)
# 
# 
# ## AZ That looks much better. Depending on your exact study question/area, I suggest to possibly use a buffered
# ## gazetter for the seas test, to avoid falsely flagging records on the coastline
# ## you can find one here: https://github.com/azizka/CoordinateCleaner/tree/master/extra_gazetteers
# ## you can also tweak the buffer area around capital and centroids using the function arguments
# 
# 
# ## HB All the records in this dataset are inside a worldclim pixel, so prob not too worried about the sea..
# 
# 
# # Also check out my newest data cleaning tutorial for GBIF data:
# # https://github.com/azizka/CoordinateCleaner/tree/master/Tutorials
# 
# 
# ## Looks good! You need to convert the data frame to a tibble to make it work...
# 
# 
# ## Now the outlier test (you can also run it through CleanCoordinates, but unfortunately it runs for awhile...)
# ## You might want to consider some tweaking of the mltpl paramters
# ?cc_outl
# 
# 
# ## This creates a data frame of just the clean records
# # ## Not sure where the extra rows are coming from.....???
# # GBIF.CLEAN    <- cc_outl(TIB.TEST,
# #                          lon     = "decimallongitude",
# #                          lat     = "decimallatitude",
# #                          species = "species",
# #                          method  = "quantile", mltpl = 5, tdi = 1000,
# #                          value   = "clean",  verbose = T)
# #
# #
# # ## This creates a vector of true/false for each record
# # GBIF.SPAT.OUT <- cc_outl(TIB.TEST,
# #                          lon     = "decimallongitude",
# #                          lat     = "decimallatitude",
# #                          species = "species",
# #                          method  = "quantile", mltpl = 5, tdi = 1000,
# #                          value   = "FLAGS",  verbose = T)
# 
# 
# ## Check the output ::
# dim(GBIF.CLEAN);length(GBIF.SPAT.OUT)
# 
# 
# ## Join data
# TEST.GEO = cbind(TIB.TEST, FLAGS) #, GBIF.SPAT.OUT)
# class(TEST.GEO)
# names(TEST.GEO)
# dim(TEST.GEO)
# 
# 
# ## Check the values for each field
# unique(TEST.GEO$validity)
# unique(TEST.GEO$equal)
# unique(TEST.GEO$equal)
# unique(TEST.GEO$zeros)
# unique(TEST.GEO$capitals)
# unique(TEST.GEO$centroids)
# unique(TEST.GEO$outliers)
# unique(TEST.GEO$gbif)
# unique(TEST.GEO$institution)
# unique(TEST.GEO$outliers)
# unique(TEST.GEO$gbif)
# unique(TEST.GEO$summary)
# unique(TEST.GEO$GBIF.SPAT.OUT)





#########################################################################################################################
## PLOT CLEAN DATA
#########################################################################################################################


# ## Just get the clean data
# TEST.CLEAN = subset(TEST.GEO, summary == "TRUE")
# unique(TEST.CLEAN$country)
# 
#   
# ## HB I haven't removed any spatial outliers yet...
# plot(LAND)
# points(TEST.CLEAN$decimallongitude,  TEST.CLEAN$decimallatitude, cex = 0.2, col = "red", pch = 19)
# 
# 
# ## Still have the outliers on the eastern half of Australia
# plot(aus)
# points(TEST.CLEAN$decimallongitude,  TEST.CLEAN$decimallatitude, cex = 0.8, col = "red", pch = 19)





#########################################################################################################################
## 2). NOW, CLEAN WHOLE DATASET
#########################################################################################################################


#########################################################################################################################
## FLAG GBIF OUTLIERS
#########################################################################################################################


#########################################################################################################################
## Rename the columns to fit the CleanCoordinates format
str(unique(COMBO.RASTER.CONVERT$searchTaxon))
GBIF.TRIM.GEO = dplyr::rename(COMBO.RASTER.CONVERT, 
                              species = searchTaxon,
                              decimallongitude = lon, 
                              decimallatitude  = lat)


## The ?CleanCoordinates function needs to work on a tibble, for some reason !
TIB.TEST <- timetk::tk_tbl(GBIF.TRIM.GEO)
dim(TIB.TEST)
str(unique(TIB.TEST$species))


## I've already stripped out the records that fall outside
## the worldclim raster boundaries, so the sea test is probably not the most important for me
## Study area is the globe, but we are only projecting models onto Australia


#########################################################################################################################
## Don't run the outliers test here, it is slower...
FLAGS  <- CleanCoordinates(TIB.TEST,
                           #countries        = "country",    ## too many flagged here...
                           capitals.rad     = 0.12,
                           countrycheck     = TRUE,
                           duplicates       = TRUE,
                           seas             = FALSE)


## save/load the flags
#saveRDS(FLAGS, 'data/base/HIA_LIST/COMBO/COMBO_GBIF_FLAGS.rds')
# FLAGS = readRDS('data/base/HIA_LIST/COMBO/COMBO_GBIF_FLAGS.rds')


## Flagging ~ 1.64%, excluding the spatial outliers. Seems reasonable?
summary(FLAGS)
FLAGS = FLAGS[ ,!(colnames(FLAGS) == "decimallongitude" | colnames(FLAGS) == "decimallatitude")]
names(FLAGS)
summary(FLAGS)[8]/dim(FLAGS)[1]*100


## A plot like this for each species would be awesome
#plot(FLAGS)





#########################################################################################################################
## 3). FLAG SPATIAL OUTLIERS
#########################################################################################################################


#########################################################################################################################
## Then split the data into n maneageable subsets to check the spatial outliers, but just for the 120 modelled species
# TIB.MILE  = TIB.TEST[TIB.TEST$species %in% spp.mile, ] 
# 
# 
# ## 8 subsets of 15 species each
# n = 8
# dim(TIB.MILE)[1]/n
# REP <- rep(1:n, each = round(dim(TIB.MILE)[1]/n, digits = 0))
# REP <- head(REP, dim(TIB.MILE)[1])
# 
# identical(dim(TIB.MILE)[1], length(REP))
# head(REP)
# tail(REP)
# 
# 
# ## If the vector is a non factoris-abubble length, make it the same
# dim(TIB.MILE)[1] - length(REP)
# REP <- c(REP, rep(n, dim(TIB.MILE)[1] - length(REP)))
# dim(TIB.MILE)[1] - length(REP)
# TIB.MILE$REP = REP
# head(TIB.MILE$REP);tail(TIB.MILE$REP)
# 
# 
# ## Could create a list of data frames :: save data to run multiple sessions
# ## 100k points shouldn't be too bad to process
# OUT <- split( TIB.MILE , f = TIB.MILE$REP )
# dim(OUT[[1]]);dim(OUT[[4]]);dim(OUT[[8]])


#save.image("STEP_5_COORD_CLEAN.RData")
#load("STEP_COORD_CLEAN.RData")


#########################################################################################################################
## Now run the spatial clean on each list element :: one at a time because it is too big
## For all elements
# GBIF.SPAT.OUT = sapply( OUT , function(x) cc_outl( x,
#                                                    lon     = "decimallongitude",
#                                                    lat     = "decimallatitude",
#                                                    species = "species",
#                                                    method  = "quantile",
#                                                    mltpl   = 5,
#                                                    tdi     = 1000,
#                                                    value   = "flags") )
# 
# 
# ## Save spatial outliers
# saveRDS(GBIF.SPAT.OUT, 'data/base/HIA_LIST/COMBO/SPAT_OUT/SPAT_OUT.rds')


#########################################################################################################################
## Join data :: exclude the decimal lat/long, check the length 
dim(GBIF.TRIM.TEST)[1];dim(FLAGS)[1]#;length(GBIF.SPAT.OUT)
names(FLAGS)[1] = c("coord_spp")
identical(COMBO.RASTER.CONVERT$searchTaxon, FLAGS$coord_spp)                                            ## order matches
identical(dim(FLAGS)[1], dim(GBIF.TRIM.GEO)[1])


TEST.GEO = cbind(COMBO.RASTER.CONVERT, FLAGS)#, GBIF.SPAT.OUT)
identical(TEST.GEO$searchTaxon, TEST.GEO$coord_spp)                                                     ## order matches


## Check the values for each flag
# summary(TEST.GEO$validity)
# summary(TEST.GEO$equal)
# summary(TEST.GEO$zeros)
# summary(TEST.GEO$capitals)
# summary(TEST.GEO$centroids)
# summary(TEST.GEO$gbif)
# summary(TEST.GEO$institution)
# summary(TEST.GEO$gbif)
# summary(TEST.GEO$summary)
#summary(TEST.GEO$GBIF.SPAT.OUT)





#########################################################################################################################
## 4). SUBSET FOR THE NICHES : JUST USE COORDCLEAN SUMMARY
#########################################################################################################################


## So ~2.6% of the data is dodgy according to the GBIF fields or spatial outliers
## This seems ok as a median figure across the data set?
#dim(subset(TEST.GEO, summary == "FALSE" | GBIF.SPAT.OUT == "FALSE"))[1]/dim(TEST.GEO)[1]*100
CLEAN.FALSE = subset(TEST.GEO, summary == "FALSE")# | GBIF.SPAT.OUT == "FALSE")


## Check one species
# coordyline = subset(TEST.GEO, searchTaxon == "Cordyline australis")
# View(subset(coordyline, summary == "FALSE") #| GBIF.SPAT.OUT == "FALSE")
#      [, c("searchTaxon", "OBS",
#           "lon",
#           "lat",
#           "validity", 
#           "equal",
#           "zeros",
#           "capitals",
#           "centroids",
#           "duplicates",
#           "gbif",
#           "institution",
#           "summary")])


## Not sure why the inverse did not work :: get only the records which were not flagged as being dodgy.
## dim(subset(TEST.GEO, summary == "TRUE" | GBIF.SPAT.OUT == "TRUE"))
#CLEAN.TRUE = TEST.GEO[!TEST.GEO$OBS %in% CLEAN.FALSE$OBS, ]
CLEAN.TRUE = subset(TEST.GEO, summary == "TRUE")
identical(dim(CLEAN.TRUE)[1], (dim(COMBO.RASTER.CONVERT)[1] - dim(subset(TEST.GEO, summary == "FALSE"))[1]))
unique(CLEAN.TRUE$summary)                                              
#unique(CLEAN.TRUE$GBIF.SPAT.OUT)                                       


## How many species?
(dim(CLEAN.TRUE)[1]/dim(TEST.GEO)[1])*100                                               ## x% of the records are retained





#########################################################################################################################
## 5). INTERSECT SPECIES RECORDS WITH LOCAL GOV AREAS AND SIGNIFICANT URBAN AREAS
#########################################################################################################################


#########################################################################################################################
## See the ABS for details :: there are 563 LGAs
## http://www.abs.gov.au/ausstats/abs@.nsf/Lookup/by%20Subject/1270.0.55.003~July%202016~Main%20Features~Local%20Government%20Areas%20(LGA)~7


## We want to know the count of species that occur in 'n' LGAs, across a range of climates. Read in LGA and SUA
projection(LGA);projection(SUA);projection(AUS)


## Convert the raster data back into a spdf
COMBO.RASTER.SP   = SpatialPointsDataFrame(coords      = CLEAN.TRUE[c("lon", "lat")], 
                                           data        = CLEAN.TRUE,
                                           proj4string = CRS.WGS.84)


## Project using a projected rather than geographic coordinate system
LGA.WGS  = spTransform(LGA, CRS.WGS.84)
SUA.WGS  = spTransform(SUA, CRS.WGS.84)
AUS.WGS  = spTransform(AUS, CRS.WGS.84)


## Remove the columns we don't need
LGA.WGS = LGA.WGS[, c("LGA_CODE16", "LGA_NAME16")] 


#########################################################################################################################
## This creates a shapefile which is just the SUA: in or out...
# IN.SUA <- SUA[ !(SUA$SUA_NAME11 %in% c("Not in any Significant Urban Area (NSW)", 
#                                        "Not in any Significant Urban Area (Vic.)",
#                                        "Not in any Significant Urban Area (Qld)",
#                                        "Not in any Significant Urban Area (SA)",
#                                        "Not in any Significant Urban Area (WA)",
#                                        "Not in any Significant Urban Area (Tas.)",
#                                        "Not in any Significant Urban Area (NT)",
#                                        "Not in any Significant Urban Area (OT)",
#                                        "Not in any Significant Urban Area (ACT)")), ]
# 
# plot(IN.SUA)
# writeOGR(obj = IN.SUA, dsn = "./data/base/CONTEXTUAL", layer = "IN.SUA", driver = "ESRI Shapefile")


## Then, we want to create a layer which is just in the urban area, or not. This would need to combine the above fields into one
# IN.SUA   = readOGR("F:/green_cities_sdm/data/base/CONTEXTUAL/INSIDE_AUS_SUA.shp", layer = "INSIDE_AUS_SUA")
# saveRDS(IN.SUA, file.path("F:/green_cities_sdm/data/base/CONTEXTUAL/", 'IN_SUA_AUS.rds'))
# 
# IN.SUA   = IN.SUA[,-(1)];
# IN.SUA   = IN.SUA[,-(2)]
# names(IN.SUA)
# IN.SUA  = spTransform(IN.SUA, CRS.WGS.84)


#########################################################################################################################
## Run join between species records and LGAs/SUAs :: Double check they are the same
projection(COMBO.RASTER.SP);projection(LGA.WGS);projection(SUA.WGS);projection(AUS.WGS)
LGA.JOIN      = over(COMBO.RASTER.SP, LGA.WGS)              ## =SUA.JOIN      = over(COMBO.RASTER.SP, SUA.WGS) 
COMBO.SUA.LGA = cbind.data.frame(COMBO.RASTER.SP, LGA.JOIN) 
#saveRDS(COMBO.SUA.LGA, file = paste("./data/base/HIA_LIST/GBIF/COMBO_SUA_LGA.rds"))
## COMBO.SUA.LGA = readRDS("./data/base/HIA_LIST/GBIF/COMBO_SUA_LGA.rds")
## str(unique(COMBO.SUA.LGA$searchTaxon))


#########################################################################################################################
## AGGREGATE THE NUMBER OF LGAs EACH SPECIES IS FOUND IN. NA LGAs ARE OUTSIDE AUS
LGA.AGG   = tapply(COMBO.SUA.LGA$LGA_NAME16, COMBO.SUA.LGA$searchTaxon, function(x) length(unique(x))) ## group LGA by species name
AUS.AGG   = aggregate(LGA_CODE16 ~ searchTaxon, data = COMBO.SUA.LGA, function(x) {sum(!is.na(x))}, na.action = NULL)
LGA.AGG   = as.data.frame(LGA.AGG)
LGA.AGG   = cbind.data.frame(AUS.AGG, LGA.AGG)
names(LGA.AGG) = c("searchTaxon", "AUS_RECORDS", "LGA_COUNT")


## Check
dim(LGA.AGG)
head(LGA.AGG)


## 
names(COMBO.SUA.LGA)
COMBO.SUA.LGA = subset(COMBO.SUA.LGA, select = -c(lon.1, lat.1))
names(COMBO.SUA.LGA)
dim(COMBO.SUA.LGA)
str(unique(COMBO.SUA.LGA$searchTaxon))





#########################################################################################################################
## 6). CREATE NICHES FOR SELECTED TAXA
#########################################################################################################################


#########################################################################################################################
## Now summarise the niches. But figure out a cleaner way of doing this
env.variables = c("Annual_mean_temp",
                  "Mean_diurnal_range",
                  "Isothermality",
                  "Temp_seasonality",
                  "Max_temp_warm_month",
                  "Min_temp_cold_month",
                  "Temp_annual_range",
                  "Mean_temp_wet_qu",
                  "Mean_temp_dry_qu",
                  "Mean_temp_warm_qu",
                  "Mean_temp_cold_qu",
                  
                  "Annual_precip",
                  "Precip_wet_month",
                  "Precip_dry_month",
                  "Precip_seasonality",
                  "Precip_wet_qu",
                  "Precip_dry_qu",
                  "Precip_warm_qu",
                  "Precip_col_qu",
                  "PET")


#########################################################################################################################
## Change the raster values here: See http://worldclim.org/formats1 for description of the interger conversion. 
## All temperature variables wer multiplied by 10, so divide by 10 to reverse it.
# COMBO.RASTER.CONVERT = as.data.table(COMBO.SUA.LGA)                           ## Check this works, also inefficient
# COMBO.RASTER.CONVERT[, (env.variables [c(1:11)]) := lapply(.SD, function(x) 
#   x / 10 ), .SDcols = env.variables [c(1:11)]]
# COMBO.RASTER.CONVERT = as.data.frame(COMBO.RASTER.CONVERT)                   ## Find another method without using data.table
# 
# 
# ## Check Looks ok?
# summary(COMBO.RASTER.CONVERT$Annual_mean_temp)  ## -23 looks too low/. Check where these are ok
# summary(COMBO.RASTER$Annual_mean_temp)
# 
# summary(COMBO.RASTER.CONVERT$Isothermality)
# summary(COMBO.RASTER$Isothermality)
# 
# 
# ## Plot a few points to see :: do those look reasonable?
# plot(LAND, col = 'grey', bg = 'sky blue')
# points(COMBO.RASTER.CONVERT[ which(COMBO.RASTER.CONVERT$Annual_mean_temp < -5), ][, c("lon", "lat")], 
#        pch = ".", col = "red", cex = 3, asp = 1)


#########################################################################################################################
## Create niche summaries for each environmental condition like this...
## Here's what the function will produce :
NICHE.DF = completeFun(COMBO.SUA.LGA, "PET")
dim(NICHE.DF)
head(niche_estimate (DF = NICHE.DF, colname = "Annual_mean_temp"))  ## including the q05 and q95


## So lets use lapply on the "SearchTaxon"
## test = run_function_concatenate(list, DF, "DF, colname = x") 
COMBO.NICHE <- env.variables %>% 
  
  ## Pipe the list into lapply
  lapply(function(x) {
    
    ## Now use the niche width function on each colname (so 8 environmental variables)
    ## Also, need to figure out how to make the aggregating column generic (species, genus, etc.)
    ## currently it only works hard-wired
    niche_estimate (DF = NICHE.DF, colname = x)
    
    ## would be good to remove the duplicate columns here
    
  }) %>% 
  
  ## finally, create one dataframe for all niches
  as.data.frame


## Remove duplicate Taxon columns and check the output :: would be great to skip these columns when running the function
names(COMBO.NICHE)
COMBO.NICHE = subset(COMBO.NICHE, select = -c(searchTaxon.1,  searchTaxon.2,  searchTaxon.3,  searchTaxon.4,
                                              searchTaxon.5,  searchTaxon.6,  searchTaxon.7,  searchTaxon.7,
                                              searchTaxon.8,  searchTaxon.9,  searchTaxon.10, searchTaxon.11,
                                              searchTaxon.12, searchTaxon.13, searchTaxon.14, searchTaxon.15,
                                              searchTaxon.16, searchTaxon.17, searchTaxon.18, searchTaxon.19))


#########################################################################################################################
## Add counts for each species, and record the total number of taxa processed
## dim(COMBO.RASTER.CONVERT);dim(CLEAN.TRUE)
COMBO.count = as.data.frame(table(COMBO.RASTER.CONVERT$searchTaxon))$Freq
identical(length(COMBO.count), dim(COMBO.NICHE )[1])

Total.taxa.processed = dim(COMBO.NICHE)[1]
COMBO.NICHE  = cbind(COMBO.count, COMBO.NICHE)
names(COMBO.NICHE)
dim(COMBO.NICHE)
#saveRDS(COMBO.NICHE, file = paste("./data/base/HIA_LIST/GBIF/COMBO_NICHE_CLEAN.rds"))


#########################################################################################################################
## Add the counts of Australian records for each species to the niche database
names(COMBO.NICHE)
names(LGA.AGG)
dim(COMBO.NICHE)
dim(LGA.AGG)


COMBO.LGA = join(COMBO.NICHE, LGA.AGG)                            ## The tapply needs to go where the niche summaries are
names(COMBO.LGA)

dim(COMBO.LGA)
head(COMBO.LGA$AUS_RECORDS)
head(COMBO.LGA$LGA_COUNT)



#########################################################################################################################
## 7). CALCULATE AREA OF OCCUPANCY RANGES 
#########################################################################################################################


## These numbers don't look that accurate: Try convert into a sp data frame and projecting into a projected system?
## Create a species list to estimate the ranges for
# spp.geo = as.character(unique(COMBO.RASTER$searchTaxon)) 
# data    = COMBO.RASTER


#########################################################################################################################
## AREA OF OCCUPANCY (AOO)
## For every species in the list: calculate the AOO
# GBIF.AOO <- spp.geo[c(1:length(spp.geo))] %>%
# 
#   ## Pipe the list into lapply
#   lapply(function(x) {
# 
#     ## Subset the the data frame
#     DF      = subset(data, searchTaxon == x)[, c("lon", "lat")]
# 
#     ## Calculate area of occupancy according the the "red" package
#     aoo (DF)
# 
#     ## Warning messages: Ask John if this is a problem
#     ## In rgdal::project(longlat, paste("+proj=utm +zone=", zone,  ... :
#     ## 3644 projected point(s) not finite
# 
#   }) %>%
# 
#   ## Finally, create one dataframe for all niches
#   as.data.frame


#########################################################################################################################
## EXTENT OF OCCURRENCE
## For every species in the list: calculate the EOO
# GBIF.EOO <- spp.geo[c(1:length(spp.geo))] %>% 
#   
#   ## Pipe the list into lapply
#   lapply(function(x) {
#     
#     ## Subset the the data frame 
#     DF      = subset(data, searchTaxon == x)[, c("searchTaxon", "lon", "lat")]
#     DF.GEO  = dplyr::rename(DF, 
#                             identifier = searchTaxon,
#                             XCOOR      = lon,
#                             YCOOR      = lat)
#     
#     ## Calculate area of occupancy according the the "red" package
#     CalcRange (DF.GEO)
#     
#     ## Warning messages: Ask John if this is a problem
#     ## In rgdal::project(longlat, paste("+proj=utm +zone=", zone,  ... :
#     ## 3644 projected point(s) not finite
#     
#   }) %>% 
#   
#   ## Finally, create one dataframe for all niches
#   as.data.frame


#########################################################################################################################
## Clean it up :: the order of species should be preserved
# GBIF.AOO = gather(GBIF.AOO)
# str(GBIF.AOO)
# str(unique(GBIF.ALA.COMBO.HIA$searchTaxon))                     ## same number of species...
# 
# 
# ## Now join on the GEOGRAPHIC RANGE
# COMBO.NICHE$AREA_OCCUPANCY = GBIF.AOO$value                     ## vectors same length so don't need to match


## AOO is calculated as the area of all known or predicted cells for the species. The resolution will be 2x2km as 
## required by IUCN. A single value in km2.

#########################################################################################################################
## Add the counts of Australian records for each species to the niche database
# names(COMBO.NICHE)
# names(LGA.AGG)
# dim(COMBO.NICHE)
# dim(LGA.AGG)
# 
# 
# COMBO.LGA = join(COMBO.NICHE, LGA.AGG)                            ## The tapply needs to go where the niche summaries are
# names(COMBO.LGA)
# 
# dim(COMBO.LGA)
# head(COMBO.LGA$AUS_RECORDS)
# head(COMBO.LGA$LGA_COUNT)





#########################################################################################################################
## 8). JOIN ON CONTEXTUAL DATA
#########################################################################################################################


#########################################################################################################################
## Now join the horticultural contextual data onto one or both tables ()
names(COMBO.RASTER.CONVERT)
names(CLEAN.SPP)
COMBO.RASTER.CONTEXT = join(COMBO.RASTER.CONVERT, HIA.SPP.JOIN)
#COMBO.RASTER.CONTEXT  = COMBO.RASTER.CONTEXT[,  c(42, 1, 65, 2:41, 43:61, 62:64, 66:78)]                         ## REDO
names(COMBO.RASTER.CONTEXT)


## Now join hort context to all the niche
COMBO.NICHE.CONTEXT = join(COMBO.LGA, HIA.SPP.JOIN)
COMBO.NICHE.CONTEXT =  COMBO.NICHE.CONTEXT[, c(2, 185, 1, 183:184, 186:198, 3:182)] 
head(COMBO.NICHE.CONTEXT$AUS_RECORDS)
head(COMBO.NICHE.CONTEXT$LGA_COUNT)


## Join on the total growers? The synonyms won't match here
names(COMBO.NICHE.CONTEXT)
COMBO.NICHE.CONTEXT$searchTaxon [! COMBO.NICHE.CONTEXT$searchTaxon  %in% TOT.GROW$searchTaxon]
COMBO.NICHE.CONTEXT$searchTaxon [! COMBO.NICHE.CONTEXT$searchTaxon  %in% HIA.list$Binomial]
COMBO.NICHE.CONTEXT$Total.growers = join(COMBO.NICHE.CONTEXT, TOT.GROW)$Total.growers


## Set NA to blank, then sort by no. of growers
COMBO.NICHE.CONTEXT$Number.of.growers[is.na(COMBO.NICHE.CONTEXT$Number.of.growers)] <- 0
COMBO.NICHE.CONTEXT = COMBO.NICHE.CONTEXT[with(COMBO.NICHE.CONTEXT, rev(order(Number.of.growers))), ]


## View the data
names(COMBO.RASTER.CONTEXT)
names(COMBO.NICHE.CONTEXT)
dim(COMBO.RASTER.CONTEXT)
dim(COMBO.NICHE.CONTEXT)


## Print the dataframe dimensions to screen
dim(CLEAN.TRUE)
dim(COMBO.NICHE.CONTEXT)
length(unique(CLEAN.TRUE$searchTaxon))
length(COMBO.NICHE.CONTEXT$searchTaxon)


#########################################################################################################################
## Save
# saveRDS(TEST.GEO,                'data/base/HIA_LIST/COMBO/CLEAN_FLAGS_HIA_SPP.rds')
# saveRDS(CLEAN.TRUE,              'data/base/HIA_LIST/COMBO/CLEAN_ONLY_HIA_SPP.rds')
# saveRDS(CLEAN.NICHE.CONTEXT,     'data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT_APRIL_2018_COORD_CLEAN.rds')
# write.csv(CLEAN.NICHE.CONTEXT,   "./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT_APRIL_2018_COORD_CLEAN.csv", row.names = FALSE)





#########################################################################################################################
#####################################################  TBC ##############################################################
#########################################################################################################################