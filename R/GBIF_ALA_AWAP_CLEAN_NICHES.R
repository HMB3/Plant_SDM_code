#########################################################################################################################
################################################# FLAG SPATIAL OUTLIERS #################################################
#########################################################################################################################


## This code uses the CoordinateCleaner package to Automatcially screen out the dodgy spatial records. See ::
## https://github.com/azizka/CoordinateCleaner


## It then creates a table of the species niche :

## 2). A table with one row for each species, including contextual data and species attributes (niches, traits, etc.)

## These tables are subsequently used to estimate the current global realised niche/climatic tolerance using the best 
## available data, and susequently model the niches using the maxent algorithm.


#########################################################################################################################
## Read in the three data tables
#COMBO.AWAP.CONVERT = readRDS("./data/base/HIA_LIST/COMBO/COMBO_RASTER_CONVERT_NEW_ALA_300_SPAT.rds")
rasterTmpFile()


## Check dimensions of the occurrence and inventory data tables.
length(unique(COMBO.AWAP.CONVERT$searchTaxon))
formatC(dim(COMBO.AWAP.CONVERT)[1], format = "e", digits = 2)
names(COMBO.AWAP.CONVERT)


#########################################################################################################################
## Create a unique identifier
COMBO.AWAP.CONVERT$CC.OBS <- 1:nrow(TIB.GBIF)





#########################################################################################################################
## 1). FLAG GBIF OUTLIERS
#########################################################################################################################


#########################################################################################################################
## Rename the columns to fit the CleanCoordinates format and create a tibble. 
## A Tibble is needed for running the spatial outlier cleaning
TIB.GBIF <- COMBO.AWAP.CONVERT %>% dplyr::rename(species          = searchTaxon,
                                                   decimallongitude = lon, 
                                                   decimallatitude  = lat) %>%
  timetk::tk_tbl()


## We've already stripped out the records that fall outside
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
identical(dim(FLAGS)[1], dim(TIB.GBIF)[1])
#saveRDS(FLAGS, paste0('data/base/HIA_LIST/COMBO/ALA_GBIF_FLAGS_', save_run, '.rds'))
#FLAGS = readRDS('data/base/HIA_LIST/COMBO/ALA_GBIF_FLAGS_OLD_ALA.rds')


## Flagging ~ 1.64%, excluding the spatial outliers. Seems reasonable?
summary(FLAGS)
FLAGS = FLAGS[ ,!(colnames(FLAGS) == "decimallongitude" | colnames(FLAGS) == "decimallatitude")]
message(round(summary(FLAGS)[8]/dim(FLAGS)[1]*100, 2), " % records removed")





#########################################################################################################################
## 2). FLAG SPATIAL OUTLIERS
#########################################################################################################################


## Alex Zizka ::
## The outlier function is limited in the amount of records it can process. It uses a distance matrix of all records per 
## species, which means that a species with 200k records will result in a 200,000x200,000 cells matrix, which will 
## probably choke your computer. The latest version of cc_outl includes a subsampling heuristic to address this problem. 
## I think this will work for you case, but it might run for a while


## Check the frequency table first, to see if any species are likely to hit this threshold ::
## So there are a few species with +100k records, this will be hard for the computer 
## Could split them up into species under 200k or not?
COMBO.LUT = as.data.frame(table(TIB.GBIF$species))
names(COMBO.LUT) = c("species", "FREQUENCY")
COMBO.LUT = COMBO.LUT[with(COMBO.LUT, rev(order(FREQUENCY))), ] 
head(COMBO.LUT);summary(COMBO.LUT$FREQUENCY)  ## Quercus robur, 214, and Fraxinus excelsior, 156
LUT.100K = as.character(subset(COMBO.LUT, FREQUENCY < 100000)$species)
LUT.100K = LUT.100K [order(LUT.100K)]


## Unfortunately, the cc_outl function can't handle vectors of a certain size - over 40 GB at least.
## So we are problably better off moving this code to after the data has been thinned by the SDM
## step 


## Create a data frame of species name and spatial outlier
SPAT.OUT <- unique(TIB.GBIF$species) %>%  ## LUT.100K
  
  ## pipe the list of species into lapply
  lapply(function(x) {
    
    ## Create the species df by subsetting by species
    f <- subset(TIB.GBIF, species == x)
    
    ## Run the spatial outlier detection
    message("Running spatial outlier detection for ", x)
    message(dim(f)[1], " records for ", x)
    sp.flag <- cc_outl(f,
                       lon     = "decimallongitude",
                       lat     = "decimallatitude",
                       species = "species",
                       method  = "distance",
                       #mltpl   = 5,
                       tdi     = 300,
                       value   = "flags")
    
    ## Now add attache column for species, and the flag for each record
    #d = data.frame(searchTaxon = x, SPAT_OUT = sp.flag)
    d = cbind(species = x,
              SPAT_OUT = sp.flag, f)[c("species", "SPAT_OUT", "CC.OBS")]
    
  }) %>%
  
  ## Finally, bind all the rows together
  bind_rows


## Check spatial data
identical(dim(SPAT.OUT)[1], dim(TIB.GBIF)[1])
unique(SPAT.OUT$species)
summary(SPAT.OUT$SPAT_OUT)
head(SPAT.OUT)




#########################################################################################################################
## 3). FILTER DATA TO THE CLEAN RECORDS
#########################################################################################################################


#########################################################################################################################
## Join data :: exclude the decimal lat/long, check the length 
identical(dim(TIB.GBIF)[1],dim(FLAGS)[1])#;length(GBIF.SPAT.OUT)
names(FLAGS)[1]    = c("coord_spp")
names(SPAT.OUT)[1] = c("spatial_spp")
identical(COMBO.AWAP.CONVERT$searchTaxon, FLAGS$coord_spp)                                            ## order matches


#########################################################################################################################
## Is the species column the same as the searchTaxon column?
TEST.GEO   = cbind(COMBO.AWAP.CONVERT, FLAGS)
TEST.GEO   = join(TEST.GEO, SPAT.OUT)
identical(TEST.GEO$searchTaxon, TEST.GEO$coord_spp)                                                     ## order matches
identical(TEST.GEO$searchTaxon, TEST.GEO$spatial_spp) 
identical(COMBO.AWAP.CONVERT$searchTaxon, TEST.GEO$spatial_spp)  


## Keep records which passed the GBIF and spatial test
dim(subset(TEST.GEO, summary == "TRUE" | SPAT_OUT == "TRUE"))
CLEAN.TRUE = subset(TEST.GEO, summary == "TRUE" & SPAT_OUT == "TRUE")
unique(CLEAN.TRUE$summary);unique(CLEAN.TRUE$SPAT_OUT)   

                                    
## What percentage of records are retained?
identical(dim(CLEAN.TRUE)[1], (dim(COMBO.AWAP.CONVERT)[1] - dim(subset(TEST.GEO, summary == "FALSE" | SPAT_OUT == "FALSE"))[1]))
length(unique(CLEAN.TRUE$searchTaxon))
message(round(dim(CLEAN.TRUE)[1]/dim(TEST.GEO)[1]*100, 2), " % records retained")                                               


#########################################################################################################################
## Save niches
saveRDS(CLEAN.TRUE, paste0('data/base/HIA_LIST/COMBO/CLEAN_TRUE_', save_run, '.rds'))





#########################################################################################################################
## 4). CREATE SHAPEFILES TO CHECK OUTLIERS ARCMAP
#########################################################################################################################


## Select columns
GBIF.ALA.AWAP  = select(CLEAN.TRUE,     CC.OBS, searchTaxon, scientificName, lat, lon, SOURCE, year, 
                        coordinateUncertaintyInMetres, 
                        country, locality, basisOfRecord, institutionCode, 
                        rank, Taxonomic.status, New.Taxonomic.status)


## Rename the fields so that ArcMap can handle them
GBIF.ALA.AWAP     = dplyr::rename(GBIF.ALA.AWAP, 
                                  TAXON     = searchTaxon,
                                  LAT       = lat,
                                  LON       = lon,
                                  SC_NAME   = scientificName,
                                  RANK      = rank,
                                  BASIS     = basisOfRecord,                
                                  LOCAL     = locality,                      
                                  INSTIT    = institutionCode,                
                                  COUNTRY   = country,                
                                  COORD_UN  = coordinateUncertaintyInMetres,
                                  YEAR      = year,
                                  TAX_STAT  = Taxonomic.status,
                                  NEW_STAT  = New.Taxonomic.status)
names(GBIF.ALA.AWAP)


#########################################################################################################################
## Then create SPDF
GBIF.ALA.AWAP.SPDF    = SpatialPointsDataFrame(coords      = GBIF.ALA.AWAP[c("LON", "LAT")],
                                               data        = GBIF.ALA.AWAP,
                                               proj4string = CRS.WGS.84)

## Write the shapefile out just in case
writeOGR(obj = GBIF.ALA.AWAP.SPDF, dsn = "./data/base/HIA_LIST/COMBO", layer = "CLEAN_AWAP_SPDF", driver = "ESRI Shapefile")





#########################################################################################################################
## 5). INTERSECT SPECIES RECORDS WITH LOCAL GOV AREAS AND SIGNIFICANT URBAN AREAS
#########################################################################################################################


#########################################################################################################################
## See the ABS for details :: there are 563 LGAs
## http://www.abs.gov.au/ausstats/abs@.nsf/Lookup/by%20Subject/1270.0.55.003~July%202016~Main%20Features~Local%20Government%20Areas%20(LGA)~7


## We want to know the count of species that occur in 'n' LGAs, across a range of climates. Read in LGA and SUA
length(unique(CLEAN.TRUE$searchTaxon))
projection(LGA);projection(AUS);projection(SUA.16)


## Convert the raster data back into a spdf
COMBO.AWAP.SP   = SpatialPointsDataFrame(coords      = CLEAN.TRUE[c("lon", "lat")], 
                                         data        = CLEAN.TRUE,
                                         proj4string = CRS.WGS.84)


## Project using a projected rather than geographic coordinate system
LGA.WGS  = spTransform(LGA,    CRS.WGS.84)
SUA.WGS  = spTransform(SUA.16, CRS.WGS.84)
AUS.WGS  = spTransform(AUS,    CRS.WGS.84)


## Remove the columns we don't need
LGA.WGS = LGA.WGS[, c("LGA_CODE16", "LGA_NAME16")] 
head(SUA.WGS)


#########################################################################################################################
## Run join between species records and LGAs/SUAs :: Double check they are the same
projection(COMBO.AWAP.SP);projection(LGA.WGS);projection(SUA.WGS);projection(AUS.WGS)
SUA.JOIN      = over(COMBO.AWAP.SP, SUA.WGS)              
COMBO.SUA.LGA = cbind.data.frame(COMBO.AWAP.SP, SUA.JOIN) 
saveRDS(COMBO.SUA.LGA, file = paste0('data/base/HIA_LIST/COMBO/COMBO_SUA_OVER_', save_run, '.rds'))
## COMBO.SUA.LGA = readRDS("./data/base/HIA_LIST/GBIF/COMBO_SUA_LGA.rds")
## str(unique(COMBO.SUA.LGA$searchTaxon))


#########################################################################################################################
## AGGREGATE THE NUMBER OF SUAs EACH SPECIES IS FOUND IN. NA LGAs ARE OUTSIDE AUS
SUA.AGG   = tapply(COMBO.SUA.LGA$SUA_NAME16, COMBO.SUA.LGA$searchTaxon, function(x) length(unique(x))) ## group SUA by species name
SUA.AGG   = as.data.frame(SUA.AGG)
AUS.AGG   = aggregate(SUA_NAME16 ~ searchTaxon, data = COMBO.SUA.LGA, function(x) {sum(!is.na(x))}, na.action = NULL)

SUA.AGG   = cbind.data.frame(AUS.AGG, SUA.AGG)
names(SUA.AGG) = c("searchTaxon", "AUS_RECORDS", "SUA_COUNT")


## Now create a table of all the SUA's that each species occurrs
SUA.SPP.COUNT = as.data.frame(table(COMBO.SUA.LGA[["SUA_NAME16"]], COMBO.SUA.LGA[["searchTaxon"]]))
names(SUA.SPP.COUNT) = c("SUA", "SPECIES", "SUA_COUNT")
saveRDS(SUA.SPP.COUNT, paste0('data/base/HIA_LIST/COMBO/SUA_SPP_COUNT', save_run, '.rds'))


## Check : That's ok, but we want a table of which SUA each species is actually in.
dim(SUA.AGG)
head(SUA.AGG)


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
                  "PET", 
                  
                  "Drought_freq_extr", 
                  "Drought_max_dur_extr", 
                  "Drought_max_int_extr", 
                  "Drought_max_rel_int_extr",
                  "Drought_mean_dur_extr", 
                  "Drought_mean_int_extr", 
                  "Drought_mean_rel_int_extr")


#########################################################################################################################
## Create niche summaries for each environmental condition like this...
## Here's what the function will produce :
NICHE.DF = completeFun(COMBO.SUA.LGA, "PET")
NICHE.DF = completeFun(NICHE.DF,      "Drought_freq_extr")
dim(NICHE.DF)
head(niche_estimate (DF = NICHE.DF, colname = "Drought_freq_extr"))  ## including the q05 and q95


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
## dim(COMBO.AWAP.CONVERT);dim(CLEAN.TRUE)
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
names(SUA.AGG)
dim(COMBO.NICHE)
dim(SUA.AGG)


COMBO.LGA = join(COMBO.NICHE, SUA.AGG)                            ## The tapply needs to go where the niche summaries are
names(COMBO.LGA)

dim(COMBO.LGA)
head(COMBO.LGA$AUS_RECORDS)
head(COMBO.LGA$SUA_COUNT)



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
# names(SUA.AGG)
# dim(COMBO.NICHE)
# dim(SUA.AGG)
# 
# 
# COMBO.LGA = join(COMBO.NICHE, SUA.AGG)                            ## The tapply needs to go where the niche summaries are
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
# names(COMBO.AWAP.CONVERT)
# names(CLEAN.SPP)
# COMBO.AWAP.CONTEXT = join(CLEAN.TRUE, HIA.SPP.JOIN)
# #COMBO.AWAP.CONTEXT  = COMBO.AWAP.CONTEXT[,  c(42, 1, 65, 2:41, 43:61, 62:64, 66:78)]                         ## REDO
# names(COMBO.AWAP.CONTEXT)
# 
# 
# ## Now join hort context to all the niche
# COMBO.NICHE.CONTEXT = join(COMBO.LGA, HIA.SPP.JOIN)
# #COMBO.NICHE.CONTEXT =  COMBO.NICHE.CONTEXT[, c(2, 185, 1, 183:184, 186:197, 3:182)] 
# head(COMBO.NICHE.CONTEXT$AUS_RECORDS)
# head(COMBO.NICHE.CONTEXT$LGA_COUNT)
# 
# 
# ## Join on the total growers? The synonyms won't match here
# names(COMBO.NICHE.CONTEXT)
# COMBO.NICHE.CONTEXT$searchTaxon [! COMBO.NICHE.CONTEXT$searchTaxon  %in% TOT.GROW$searchTaxon]
# COMBO.NICHE.CONTEXT$searchTaxon [! COMBO.NICHE.CONTEXT$searchTaxon  %in% HIA.list$Binomial]
# COMBO.NICHE.CONTEXT$Total.growers = join(COMBO.NICHE.CONTEXT, TOT.GROW)$Total.growers
# 
# 
# ## Set NA to blank, then sort by no. of growers
# COMBO.NICHE.CONTEXT$Number.of.growers[is.na(COMBO.NICHE.CONTEXT$Number.of.growers)] <- 0
# COMBO.NICHE.CONTEXT = COMBO.NICHE.CONTEXT[with(COMBO.NICHE.CONTEXT, rev(order(Number.of.growers))), ]
# 
# 
# ## View the data
# names(COMBO.AWAP.CONTEXT)
# names(COMBO.NICHE.CONTEXT)
# dim(COMBO.AWAP.CONTEXT)
# dim(COMBO.NICHE.CONTEXT)
# 
# 
# ## Print the dataframe dimensions to screen
# dim(CLEAN.TRUE)
# dim(COMBO.NICHE.CONTEXT)
# length(unique(CLEAN.TRUE$searchTaxon))
# length(COMBO.NICHE.CONTEXT$searchTaxon)


#########################################################################################################################
## Save
saveRDS(COMBO.LGA,             paste0('data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT_',  save_run, '.rds'))
#saveRDS(COMBO.AWAP.CONTEXT,  paste0('data/base/HIA_LIST/COMBO/COMBO_RASTER_CONTEXT_', save_run, '.rds'))




#########################################################################################################################
## OUTSTANDING CLEANING TASKS:
#########################################################################################################################


## Check the automatic cleaning results against manual cleaning results




#########################################################################################################################
#####################################################  TBC ##############################################################
#########################################################################################################################