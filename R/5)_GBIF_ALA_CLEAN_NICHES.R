#########################################################################################################################
################################################# FLAG SPATIAL OUTLIERS #################################################
#########################################################################################################################


## This code uses the CoordinateCleaner package to Automatcially screen out the dodgy spatial records. See ::
## https://github.com/azizka/CoordinateCleaner


## Can we flag herbarium records, plus the spatial outliers?


#########################################################################################################################
## Read in all data to run the SDM code :: species lists, shapefile, rasters & tables
#source('./R/HIA_LIST_MATCHING.R')
#TI.RASTER.CONVERT = readRDS("./data/base/HIA_LIST/COMBO/TI_RASTER_CONVERT.rds")
#COMBO.RASTER.CONVERT = readRDS("./data/base/HIA_LIST/COMBO/COMBO_RASTER_CONVERT_APRIL_2018.rds")
rasterTmpFile()



#########################################################################################################################
## 1). FLAG GBIF OUTLIERS
#########################################################################################################################


#########################################################################################################################
## Rename the columns to fit the CleanCoordinates format and create a tibble. A Tibble is needed for coordinate cleaner
str(unique(COMBO.RASTER.CONVERT$searchTaxon))
TIB.GBIF <- COMBO.RASTER.CONVERT %>% dplyr::rename(species          = searchTaxon,
                                                   decimallongitude = lon, 
                                                   decimallatitude  = lat) %>%
  timetk::tk_tbl()


## Add a column for unique observation so we can check the records match up
TIB.GBIF$CC.OBS <- 1:nrow(TIB.GBIF)


DF.GBIF <- COMBO.RASTER.CONVERT %>% dplyr::rename(species          = searchTaxon,
                                                  decimallongitude = lon, 
                                                  decimallatitude  = lat) 
DF.GBIF$CC.OBS <- 1:nrow(DF.GBIF)


## saveRDS(TIB.GBIF, file = paste("./data/base/HIA_LIST/GBIF/TIB_GBIF.rds"))
## DF.TEST = DF.GBIF[DF.GBIF$species %in% unique(DF.GBIF$species)[31:33], ] 
## saveRDS(DF.TEST, file = paste("./data/base/HIA_LIST/GBIF/DF_GBIF_TEST.rds"))
## I've already stripped out the records that fall outside
## the worldclim raster boundaries, so the sea test is probably not the most important
## Study area is the globe, but we are only projecting models onto Australia


#########################################################################################################################
## Don't run the outliers test here, it is slower. Also, can't run cleaning on the urban tree inventory data, because this
## removes all the records near capital cities 
FLAGS  <- CleanCoordinates(TIB.GBIF,
                           #countries        = "country",    ## too many flagged here...
                           capitals.rad     = 0.12,
                           countrycheck     = TRUE,
                           duplicates       = TRUE,
                           seas             = FALSE)


## save/load the flags
saveRDS(FLAGS, 'data/base/HIA_LIST/COMBO/ALA_GBIF_FLAGS.rds')
#FLAGS = readRDS('data/base/HIA_LIST/COMBO/ALA_GBIF_FLAGS.rds')


## Flagging ~ 1.64%, excluding the spatial outliers. Seems reasonable?
summary(FLAGS)
FLAGS = FLAGS[ ,!(colnames(FLAGS) == "decimallongitude" | colnames(FLAGS) == "decimallatitude")]
message(round(summary(FLAGS)[8]/dim(FLAGS)[1]*100, 2), " % records removed")


## A plot like this for each species would be awesome
#plot(FLAGS)





#########################################################################################################################
## 2). FLAG SPATIAL OUTLIERS
#########################################################################################################################


## This is not working, chews up all the RAM. Contact Alex Zika for more info RE the memory leak........................
## alexander.zizka@bioenv.gu.se 


## Alex Zizka ::
## The outlier function is limited in the amount of records it can process. It uses a distance matrix of all records per species, 
## which means that for instance for a species with 200,000k records will result in a 200,000x200,000 cells matrix, which will 
## probably choke your computer. The latest version of cc_outl includes a subsampling heuristic to address this problem. 
## I think this will work for you case, but it might run for a while


## Create a data frame of species name and spatial outlier
SPAT.OUT <- unique(TIB.GBIF$species)[1:2] %>%         
  
  ## pipe the list of species into lapply
  lapply(function(x) {
    
    ## Create the species df by subsetting by species
    f <- subset(TIB.GBIF, species == x) 
    
    ## Run the spatial outlier detection
    sp.flag <- cc_outl(f,
                       lon     = "decimallongitude",
                       lat     = "decimallatitude",
                       species = "species",
                       method  = "quantile",
                       mltpl   = 5,
                       tdi     = 1000,
                       value   = "flags")
    
    ## Now add attache column for species, and the flag for each record
    d = data.frame(searchTaxon = x, SPAT_OUT = sp.flag)
    d = cbind(searchTaxon = x, 
              SPAT_OUT = sp.flag, f)[c("searchTaxon", "SPAT_OUT", "CC.OBS")] 
    
    ## Remeber to explicitly return the df at the end of loop, so we can bind
    message("Is the spatial vector the same length as the DF? ", identical(dim(d)[1], dim(f)[1]))
    return(d)
    
  }) %>%
  
  ## Finally, bind all the rows together
  bind_rows



## Save data
head(SPAT.OUT)
unique(SPAT.OUT$searchTaxon)
dim(SPAT.OUT)
saveRDS(SPAT.OUT, 'data/base/HIA_LIST/COMBO/SPAT_OUT/GBIF_SPAT_OUT.rds')


#########################################################################################################################
## Join data :: exclude the decimal lat/long, check the length 
dim(TIB.GBIF)[1];dim(FLAGS)[1]#;length(GBIF.SPAT.OUT)
names(FLAGS)[1] = c("coord_spp")
identical(COMBO.RASTER.CONVERT$searchTaxon, FLAGS$coord_spp)                                            ## order matches
identical(dim(FLAGS)[1], dim(GBIF.TRIM.GEO)[1])


#########################################################################################################################
## This bind might introduce some NA's
TEST.GEO = cbind(COMBO.RASTER.CONVERT, FLAGS)#, GBIF.SPAT.OUT)
summary(TEST.GEO)
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
## 3). SUBSET FOR THE NICHES : JUST USE COORDCLEAN SUMMARY
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
length(unique(CLEAN.TRUE$searchTaxon))
#unique(CLEAN.TRUE$GBIF.SPAT.OUT)                                       


## How many species?
message(round(dim(CLEAN.TRUE)[1]/dim(TEST.GEO)[1]*100, 2), " % records retained")                                               


#########################################################################################################################
## Now bind on the urban tree inventory data. We are assuming this data is clean, after we manually fix the taxonomy
## Check the NAs
names(TI.RASTER.CONVERT)
names(CLEAN.TRUE)

CLEAN.TRUE = bind_rows(CLEAN.TRUE, TI.RASTER.CONVERT)
names(CLEAN.TRUE)
unique(CLEAN.TRUE$SOURCE) 
unique(CLEAN.TRUE$INVENTORY) 
length(unique(CLEAN.TRUE$searchTaxon))
summary(CLEAN.TRUE$Annual_mean_temp)


## Then create a unique ID column which can be used to identify outlier records 
CLEAN.TRUE$OBS <- 1:nrow(CLEAN.TRUE)
dim(CLEAN.TRUE)[1];length(CLEAN.TRUE$OBS)  
identical(dim(CLEAN.TRUE)[1], length(CLEAN.TRUE$OBS))


## How many records are added by including the tree inventories?
message("Tree inventory data increases records by ", round(dim(CLEAN.TRUE)[1]/dim(TEST.GEO)[1]*100, 2), " % ")     


## Save the clean species
## Plot new combined data to check for bias in "SAVE_BIASED_SPDF.R".....................................................
#saveRDS(TEST.GEO,    'data/base/HIA_LIST/COMBO/CLEAN_FLAGS_INV_SPP.rds')
#saveRDS(CLEAN.TRUE,  'data/base/HIA_LIST/COMBO/CLEAN_ONLY_INV_SPP.rds')


#########################################################################################################################
## Then create SPDF
CLEAN.SPDF   = SpatialPointsDataFrame(coords      = CLEAN.TRUE[c("lon", "lat")],
                                      data        = CLEAN.TRUE,
                                      proj4string = CRS.WGS.84)


## Project the SDM data into WGS
CLEAN.SPDF <- spTransform(CLEAN.SPDF, CRS("+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
projection(CLEAN.SPDF)


## Save the shapefile, to be subsampled in ArcMap
names(CLEAN.SPDF);head(CLEAN.SPDF)
#writeOGR(obj = CLEAN.SPDF, dsn = "./data/base/HIA_LIST/COMBO", layer = "CLEAN_SPDF", driver = "ESRI Shapefile")





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
## Run join between species records and LGAs/SUAs :: Double check they are the same
projection(COMBO.RASTER.SP);projection(LGA.WGS);projection(SUA.WGS);projection(AUS.WGS)
LGA.JOIN      = over(COMBO.RASTER.SP, LGA.WGS)              
COMBO.SUA.LGA = cbind.data.frame(COMBO.RASTER.SP, LGA.JOIN) 
#saveRDS(COMBO.SUA.LGA, file = paste("./data/base/HIA_LIST/GBIF/COMBO_SUA_LGA.rds"))
## COMBO.SUA.LGA = readRDS("./data/base/HIA_LIST/GBIF/COMBO_SUA_LGA.rds")
## str(unique(COMBO.SUA.LGA$searchTaxon))


## Join on the

#########################################################################################################################
## AGGREGATE THE NUMBER OF LGAs EACH SPECIES IS FOUND IN. NA LGAs ARE OUTSIDE AUS
## Could Also include the koppen zones here, within Australia :: "Koppen_aus"
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
unique(COMBO.SUA.LGA$SOURCE)





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
COMBO.count = as.data.frame(table(COMBO.SUA.LGA$searchTaxon))$Freq
identical(length(COMBO.count), dim(COMBO.NICHE)[1])

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
COMBO.RASTER.CONTEXT = join(CLEAN.TRUE, HIA.SPP.JOIN)
#COMBO.RASTER.CONTEXT  = COMBO.RASTER.CONTEXT[,  c(42, 1, 65, 2:41, 43:61, 62:64, 66:78)]                         ## REDO
names(COMBO.RASTER.CONTEXT)


## Now join hort context to all the niche
COMBO.NICHE.CONTEXT = join(COMBO.LGA, HIA.SPP.JOIN)
COMBO.NICHE.CONTEXT =  COMBO.NICHE.CONTEXT[, c(2, 185, 1, 183:184, 186:197, 3:182)] 
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
## OUTSTANDING CLEANING TASKS:
#########################################################################################################################


## Check the spatial outlier function works. Try doing it one species at a time, then stich the table back together

## Check distribution maps for each species for spatial outliers  - Keep a spreasheet of all species...

## Estimate native/naturalised ranges as a separate colum         - Rachael's package

## GBIF taxonomic errors                                          - Use TPL






#########################################################################################################################
#####################################################  TBC ##############################################################
#########################################################################################################################