#########################################################################################################################
################################################ CLEAN RECORDS ########################################################## 
#########################################################################################################################


## This code creates individual shapefiles for each species, so we can create a list of the spatial outliers manually
## It also makes a table of koppen zones * species occurrences


## Load the niche data and the records from step 4, after cruncing ALA and GBIF, but prior to cleaning (GBIF.ALA.COMBO.HIA)
source('./R/HIA_LIST_MATCHING.R')
load("./data/base/HIA_LIST/COMBO/GBIF_ALA_COMBO_PRE_CLEAN.RData")
load("./data/base/HIA_LIST/COMBO/COMBO_RASTER_CONTEXT_1601_2018.RData")


#########################################################################################################################
## 1). CHECK ALA/GBIF OVERLAP
#########################################################################################################################


## This column can be applied across the project  
names(GBIF.ALA.COMBO.HIA)
dim(GBIF.ALA.COMBO.HIA)                                                                           


## Check the unique IDs match
COMBO.RASTER.CONTEXT$OBS <- 1:nrow(COMBO.RASTER.CONTEXT)
dim(COMBO.RASTER.CONTEXT)[1];length(COMBO.RASTER.CONTEXT$OBS)  

head(COMBO.RASTER.CONTEXT, 10)[, c("OBS", "searchTaxon", "lat", "lon")]
head(GBIF.ALA.COMBO.HIA,   10)[, c("OBS", "searchTaxon", "lat", "lon")]


## Just get the environmetnal columns
OCC.ENV    <- select(COMBO.RASTER.CONTEXT, searchTaxon, OBS,
                     Annual_mean_temp,    Mean_diurnal_range, Isothermality,      Temp_seasonality,   Max_temp_warm_month, 
                     Mean_temp_wet_qu,    Mean_temp_dry_qu,   Temp_annual_range,  Mean_temp_warm_qu,  Mean_temp_cold_qu, 
                     Min_temp_cold_month,
                     Annual_precip,       Precip_wet_month,   Precip_dry_month,   Precip_seasonality, Precip_wet_qu,     
                     Precip_dry_qu,       Precip_warm_qu,     Precip_col_qu)


#########################################################################################################################
## Check the field names :: are ALA and GBIF ID fields mutually exclusive?
test_exotic = subset(GBIF.ALA.COMBO.HIA, searchTaxon == "Cycas revoluta")
test_native = subset(GBIF.ALA.COMBO.HIA, searchTaxon == "Syzygium floribundum")
test_spp    = subset(GBIF.ALA.COMBO.HIA, searchTaxon == "Cycas revoluta" | searchTaxon == "Syzygium floribundum")


## Are there records with both ALA and GBIF IDs? Seems not...
dim(GBIF.ALA.COMBO.HIA[!is.na(GBIF.ALA.COMBO.HIA$gbifID) & !is.na(GBIF.ALA.COMBO.HIA$id),])


## So the records don't overlap
# View(head(test_exotic, 162)[, c("searchTaxon",
#                                 "lon",
#                                 "lat",
#                                 "gbifID", 
#                                 "id")])
# 
# View(head(test_native, 521)[, c("searchTaxon",
#                                 "lon",
#                                 "lat",
#                                 "gbifID", 
#                                 "id")])


## Update
names(GBIF.ALA.COMBO.HIA)


## So, excluding the outliers is not as simple as just removing the GBIF records for native species...
## This will work for some taxa, but for others in could lead to equal bias in the model by further underestimating
## the fundamental niche.

## So the preferable option is to either use a density-based approach, and/or clean them spatially using manual inspection.
## No approach will get everything, and the scaling problems with occurrence vs spatial resolution remain. The key is to
## get sensible looking maps at the end :]


## Also, we should consider that some species (e.g. Agave attenuata) are simply so poorly sampled, that we can't model them
## well, or estimate their distribution. We should discuss how to present those, if we indeed do.


#########################################################################################################################
## 2). CREATE INDIVIDUAL SHAPEFILES WITH UNIQUE ID FOR MANUAL CLEANING
#########################################################################################################################


## Merge on the origin columns
length(unique(GBIF.ALA.COMBO.HIA$searchTaxon)) 
OCC.CONTEXT = COMBO.NICHE.CONTEXT[,c("searchTaxon", "Origin", "Plant.type", "Top_200")] 


## The whole table of records is too big, so split it up into each species...
GBIF.ALA.CLEAN               = GBIF.ALA.COMBO.HIA
GBIF.ALA.CLEAN               = join(GBIF.ALA.COMBO.HIA, OCC.CONTEXT)
dim(GBIF.ALA.CLEAN)
names(GBIF.ALA.CLEAN)


## Reorder the columns to see in ArcMap
GBIF.ALA.CLEAN  = select(GBIF.ALA.CLEAN, searchTaxon, OBS,     lat,            lon,
                         Origin,         Plant.type,  Top_200, scientificName, commonname, 
                         taxonRank,      genus,       family,  TPL_binomial,   taxo_agree, 
                         gbifID,         id,          cloc,    basisOfRecord,  locality, 
                         establishmentMeans, institutionCode, datasetName, habitat, catalogNumber, 
                         country, coordinateUncertaintyInMeters, geodeticDatum, year, month,                         
                         day, eventDate, eventID, CULTIVATED, POLY_ID, SUB_CODE_7, REG_CODE_7)


## Rename the fields so that ArcMap can handle them
GBIF.ALA.CLEAN     = dplyr::rename(GBIF.ALA.CLEAN, 
                                   TAXON     = searchTaxon,
                                   ORIGIN    = Origin,
                                   LAT       = lat,
                                   LON       = lon,
                                   TYPE      = Plant.type,
                                   T200      = Top_200,
                                   SC_NAME   = scientificName,
                                   COM_NAME  = commonname,
                                   RANK      = taxonRank,
                                   GENUS     = genus,
                                   FAMILY    = family,
                                   TPL       = TPL_binomial,
                                   TX_AGR    = taxo_agree,
                                   GBIF_ID   = gbifID,
                                   ALA_ID    = id,
                                   CLOC      = cloc,                     
                                   BASIS     = basisOfRecord,                
                                   LOCAL     = locality,                      
                                   ESTAB     = establishmentMeans,          
                                   INSTIT    = institutionCode,                
                                   DATASET   = datasetName,                      
                                   HABITAT   = habitat,                
                                   CAT_NO    = catalogNumber,                
                                   COUNTRY   = country,                
                                   COORD_UN  = coordinateUncertaintyInMeters,
                                   DATUM     = geodeticDatum,                 
                                   YEAR      = year,                          
                                   MNTH      = month,                        
                                   DAY       = day,                      
                                   DATE      = eventDate,                     
                                   EV_ID     = eventID,                      
                                   CULT      = CULTIVATED)


##
names(GBIF.ALA.CLEAN)


#########################################################################################################################
## The whole table of records is too big, so split it up into each species...
## Create spatial points dataframe
GBIF.ALA.POINTS = SpatialPointsDataFrame(coords      = GBIF.ALA.CLEAN[c("LON", "LAT")], 
                                         data        = GBIF.ALA.CLEAN,
                                         proj4string = CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))

TEST.POINTS     = SpatialPointsDataFrame(coords      = test_spp[c("lon", "lat")], 
                                         data        = test_spp,
                                         proj4string = CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))

## OK
names(GBIF.ALA.POINTS)
dim(GBIF.ALA.POINTS)


## Could also remove the species with insufficent records to run maxent?
#COMBO.RASTER.CLEAN   = COMBO.RASTER.CLEAN[!COMBO.RASTER.CLEAN$searchTaxon %in% unique(MISSING$searchTaxon), ]
#FINAL.MAXENT.LIST =  unique(COMBO.RASTER.CLEAN$searchTaxon[!COMBO.RASTER.CLEAN$searchTaxon %in% unique(MISSING$searchTaxon) ])


## Create list of species 
TAXA <- unique(GBIF.ALA.POINTS$TAXON)


#########################################################################################################################
## Then, loop over the species list and create a shapefile for each 
for (i in 1:length(TAXA)) {
  
  ## Need to check the OBS column matches up - or do we not need this again?
  tmp <- GBIF.ALA.POINTS[GBIF.ALA.POINTS$TAXON == TAXA[i], ] 
  writeOGR(tmp, dsn = "./data/base/HIA_LIST/COMBO/CLEAN_GBIF", TAXA[i], driver = "ESRI Shapefile", overwrite_layer = TRUE)
  
}


#########################################################################################################################
## Now merge the shapefiles into 5-10 big files for easier previewing





#########################################################################################################################
## 3). INTERSECT SPECIES RECORDS WITH CLIMATE ZONES
#########################################################################################################################


## Read in the Koppen shapefile :: 1975 centred data
## https://www.climond.org/Core/Authenticated/KoppenGeiger.aspx
Koppen = readOGR("F:/green_cities_sdm/data/base/CONTEXTUAL/WC05_1975H_Koppen.shp", layer = "WC05_1975H_Koppen")
head(Koppen)
str(unique(Koppen$Koppen))


## Project :: need to be the same for the spatial overlay 
CRS.new          <- CRS("+init=epsg:4326")        ## EPSG:3577
Koppen           = spTransform(Koppen, CRS.new)
GBIF.ALA.POINTS  = spTransform(GBIF.ALA.POINTS, CRS.new)
TEST.POINTS      = spTransform(TEST.POINTS, CRS.new)


## Check
projection(Koppen)
projection(GBIF.ALA.POINTS)
projection(TEST.POINTS)


## Intersect the species with the data 
KOPPEN.TAXA = over(TEST.POINTS, Koppen) 
head(KOPPEN.TAXA)
head(TEST.POINTS[1])


## Join the koppen classification onto the data
TAXA.KOPPEN.JOIN = cbind.data.frame(TEST.POINTS, KOPPEN.TAXA)
unique(TAXA.KOPPEN.JOIN$searchTaxon)


#########################################################################################################################
## How many Koppen zones each taxa falls in ... 
KOPP.TAXA   = tapply(TAXA.KOPPEN.JOIN$Koppen, TAXA.KOPPEN.JOIN$searchTaxon, function(x) length(unique(x)))
KOPP.TAXA   = as.data.frame(KOPP.TAXA)
setDT(KOPP.TAXA, keep.rownames = TRUE)[]
names(KOPP.TAXA) <-  c('searchTaxon', 'NUMB_KOPPZN')


## How many species records are in each koppen zone (total)?
KOPP.OCC        = data.frame(table(TAXA.KOPPEN.JOIN$Koppen))        ## group Koppen by occurrence
names(KOPP.OCC) <-  c('KOPPEN', 'NUMB_REC')


## First, produce a "long" format table with records (rows) * Koppen (columns)
TAXA.KOPPEN.JOIN = TAXA.KOPPEN.JOIN[with(TAXA.KOPPEN.JOIN, order(Koppen)), ]
head(TAXA.KOPPEN.JOIN, 10)[, c("Koppen", "searchTaxon", "lat", "lon", "OBS")]
TAXA.KOPPEN.JOIN = TAXA.KOPPEN.JOIN[, c("Koppen", "searchTaxon", "lat", "lon", "OBS")]


#########################################################################################################################
## Now try counting the number of records in each Koppen Zone
TAXA.KOPPEN.COUNT     = TAXA.KOPPEN.JOIN[, c("Koppen", "searchTaxon")]


## We want a table with one row per species, and a column for every koppen  
TEST.CAST             = dcast(TAXA.KOPPEN.COUNT, searchTaxon ~ Koppen)
TEST.CAST$ALL_RECORDS = rowSums(TEST.CAST[2:14])


## Then can we make another column with the % for every species
## Create a list
col.list = names(TEST.CAST)
col.list["searchTaxon"] <- NULL


## Update the
koppen.list = as.list(unique(Koppen$Koppen))


## Loop over all the Koppens and create a new column 
for(kop in koppen.list){
  
  TEST.CAST[[paste("Percent_rec_in_", kop, sep = "")]] <- TEST.CAST[[kop]]/TEST.CAST$ALL_RECORDS
  
}





## Then create a table of how many records (for all spp) are in each koppen zone
head(KOPP.OCC) 


## Then, summarise into a table of how many koppen zones each species falls in
head(KOPP.TAXA)




#########################################################################################################################
## Now try intersecting the environmental matrix with the koppen * records
test_env           = subset(OCC.ENV, searchTaxon == "Cycas revoluta" | searchTaxon == "Syzygium floribundum")
TAXA.KOP.ENV       = join(TAXA.KOPPEN.JOIN, test_env)


#########################################################################################################################
## Need a criteria for removing records :: key ecological information is in the density of records 

## Spatial outliers     - long way from where they are most dense, e.g. in a garden of central Australian homestead
## Logistical outliers  - Herbarium specimens, etc (we know the location of Australian herbaria, but not overseas)



#########################################################################################################################
## 4). MANUALLY INSPECT THE SHAPEFILES
#########################################################################################################################


#########################################################################################################################
## Need a criteria for removing records :: key ecological information is in the density of records 

## Spatial outliers     - long way from where they are most dense, e.g. in a garden of central Australian homestead
## Logistical outliers  - Herbarium specimens, etc (we know the location of Australian herbaria, but not overseas)







#########################################################################################################################
## Then Read back in a list of suspect points and remove them ................................................................





