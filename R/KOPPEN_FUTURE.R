#########################################################################################################################
################################################ CLEAN RECORDS ########################################################## 
#########################################################################################################################


#########################################################################################################################
## This code combines the table of Koppen * spp records, and uses Renee's species list to choose 10 species to trial
## The Koppen approach, outlined below :: 


## We could then group species according to their distribution in Koppen zones. For each of these major groupings, we 
## would then come up with a list of variables to use in the models, based on our understanding of the climate of those 
## zones. So, for example, we would be able to identify a group of species that are primarily Montane. We would then decide 
## on which variables to use for these. A second group of species might be primarily ari##semi-arid, and so on.
## So the point of overlaying the occurrences with the Koppen is simply to identify what zones the species falls in. We 
## would do that by calculating the % of a species' records in each zone. I have done similar things before, and I guarantee 
## it still wont be straightforward with grouping species, but it will be a start.


## Load the records from step 4
source('./R/HIA_LIST_MATCHING.R')
COMBO.RASTER.CONTEXT = readRDS('data/base/HIA_LIST/COMBO/CLEAN_ONLY_HIA_SPP.rds')
COMBO.NICHE.CONTEXT  = readRDS('data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT_APRIL_2018_COORD_CLEAN.rds')
  

## Just look at the experimental species
COMBO.POINTS = COMBO.RASTER.CONTEXT[COMBO.RASTER.CONTEXT$searchTaxon %in% spp.combo , ]
str(unique(sort(COMBO.POINTS$searchTaxon))) 
dim(COMBO.POINTS)


## Read in the Koppen shapefile :: 1975 centred data
## https://www.climond.org/Core/Authenticated/KoppenGeiger.aspx
# Koppen.1975      = readOGR("F:/green_cities_sdm/data/base/CONTEXTUAL/WC05_1975H_Koppen_Shapefile/WC05_1975H_Koppen_Kriticos_2012.shp",
#                       layer = "WC05_1975H_Koppen_Kriticos_2012")
# 
# Kopp.2030     = readOGR("H:/green_cities_sdm/data/base/CONTEXTUAL/CM10_Kop_Shp_V1.2/KOP_2030_A1B_CS.shp",
#                         layer = "KOP_2030_A1B_CS")
# 
# Kopp.2050     = readOGR("H:/green_cities_sdm/data/base/CONTEXTUAL/CM10_Kop_Shp_V1.2/KOP_2050_A1B_CS.shp",
#                         layer = "KOP_2050_A1B_CS")
# 
# Kopp.2070     = readOGR("H:/green_cities_sdm/data/base/CONTEXTUAL/CM10_Kop_Shp_V1.2/KOP_2070_A1B_CS.shp",
#                         layer = "KOP_2070_A1B_CS")
# 
# saveRDS(Kopp.1975, 'data/base/CONTEXTUAL/KOP_1975.rds')
# saveRDS(Kopp.2030, 'data/base/CONTEXTUAL/KOP_2030.rds')
# saveRDS(Kopp.2050, 'data/base/CONTEXTUAL/KOP_2050.rds')
# saveRDS(Kopp.2070, 'data/base/CONTEXTUAL/KOP_2070.rds')


## saveRDS(Kopp.1975, 'data/base/CONTEXTUAL/KOP_1975.rds')
Kopp.1975   = readRDS('data/base/CONTEXTUAL/KOP_1975.rds')
Kopp.2030   = readRDS('data/base/CONTEXTUAL/KOP_2030.rds')
Kopp.2050   = readRDS('data/base/CONTEXTUAL/KOP_2050.rds')
Kopp.2070   = readRDS('data/base/CONTEXTUAL/KOP_2070.rds')


## Check
head(Kopp.1975)
names(Kopp.2030) = c("Koppen_2030_A1B_CS")
names(Kopp.2050) = c("Koppen_2050_A1B_CS")
names(Kopp.2070) = c("Koppen_2070_A1B_CS")


## I think this will be ok :: the grouping will take care of the multipart problem
dim(Kopp.1975)
dim(Kopp.2030)
dim(Kopp.2050)
dim(Kopp.2070)


## All are the same
sort(unique(Kopp.1975$Koppen))
sort(unique(Kopp.2030$Koppen_2030_A1B_CS))
sort(unique(Kopp.2050$Koppen_2050_A1B_CS))
sort(unique(Kopp.2070$Koppen_2070_A1B_CS))





#########################################################################################################################
## 1). CREATE SPATIAL POINTS DATA FRAMES FOR OVERLAY
#########################################################################################################################


## Check the unique IDs match
head(COMBO.POINTS, 10)[, c("OBS", "searchTaxon", "lat", "lon")]


## Just get the environmental columns
# OCC.ENV    <- select(COMBO.RASTER.CONTEXT, searchTaxon, OBS,
#                      Annual_mean_temp,     Mean_diurnal_range, Isothermality,      Temp_seasonality,   Max_temp_warm_month, 
#                      Mean_temp_wet_qu,     Mean_temp_dry_qu,   Temp_annual_range,  Mean_temp_warm_qu,  Mean_temp_cold_qu, 
#                      Min_temp_cold_month,
#                      Annual_precip,        Precip_wet_month,   Precip_dry_month,   Precip_seasonality, Precip_wet_qu,     
#                      Precip_dry_qu,        Precip_warm_qu,     Precip_col_qu)


## Create spatial points dataframe
GBIF.ALA.POINTS = SpatialPointsDataFrame(coords      = COMBO.POINTS[c("lon", "lat")], 
                                         data        = COMBO.POINTS,
                                         proj4string = CRS.WGS.84)

dim(GBIF.ALA.POINTS)
length(unique(GBIF.ALA.POINTS$searchTaxon))





#########################################################################################################################
## 2). INTERSECT SPECIES RECORDS WITH CLIMATE ZONES
#########################################################################################################################


## The projections need to match for the spatial overlay 
Kopp.1975        = spTransform(Kopp.1975, CRS.WGS.84)
Kopp.2030        = spTransform(Kopp.2030, CRS.WGS.84)
Kopp.2050        = spTransform(Kopp.2050, CRS.WGS.84)
Kopp.2070        = spTransform(Kopp.2070, CRS.WGS.84)


## Check they are the same
projection(Kopp.1975)
projection(Kopp.2030)
projection(Kopp.2050)
projection(Kopp.2070)
projection(GBIF.ALA.POINTS)


#########################################################################################################################
## Now, intersect the species with the data :: these objects are huge!
KOPPEN.TAXA.1975 = over(GBIF.ALA.POINTS, Kopp.1975)
saveRDS(KOPPEN.TAXA.1975, 'data/base/CONTEXTUAL/KOPPEN_JOIN/KOPPEN_TAXA_1975.rds')
#KOPPEN.TAXA.1975 = readRDS('data/base/CONTEXTUAL/KOPPEN_JOIN/KOPPEN_TAXA_1975.rds')


KOPPEN.TAXA.2030 = over(GBIF.ALA.POINTS, Kopp.2030)
saveRDS(KOPPEN.TAXA.2030, 'data/base/CONTEXTUAL/KOPPEN_JOIN/KOPPEN_TAXA_2030.rds')
#KOPPEN.TAXA.2030 =readRDS('data/base/CONTEXTUAL/KOPPEN_JOIN/KOPPEN_TAXA_2030.rds')


KOPPEN.TAXA.2050 = over(GBIF.ALA.POINTS, Kopp.2050)
saveRDS(KOPPEN.TAXA.2050, 'data/base/CONTEXTUAL/KOPPEN_JOIN/KOPPEN_TAXA_2050.rds')
#KOPPEN.TAXA.2050 =readRDS('data/base/CONTEXTUAL/KOPPEN_JOIN/KOPPEN_TAXA_2050.rds')


KOPPEN.TAXA.2070 = over(GBIF.ALA.POINTS, Kopp.2070)
saveRDS(KOPPEN.TAXA.2070, 'data/base/CONTEXTUAL/KOPPEN_JOIN/KOPPEN_TAXA_2070.rds')
#KOPPEN.TAXA.2070 = readRDS('data/base/CONTEXTUAL/KOPPEN_JOIN/KOPPEN_TAXA_2070.rds')



## What do they look like?
dim(KOPPEN.TAXA.1975)
length(sort(unique(KOPPEN.TAXA.1975$Koppen)))

dim(KOPPEN.TAXA.2030)
length(sort(unique(KOPPEN.TAXA.2030$Koppen_2030_A1B_CS)))

dim(KOPPEN.TAXA.2050)
length(sort(unique(KOPPEN.TAXA.2050$Koppen_2050_A1B_CS)))

dim(KOPPEN.TAXA.2070)
length(sort(unique(KOPPEN.TAXA.2050$Koppen_2070_A1B_CS)))


## Then join the koppen classification onto the data
TAXA.KOPPEN.1975.JOIN = cbind.data.frame(GBIF.ALA.POINTS, KOPPEN.TAXA.1975)
str(unique(TAXA.KOPPEN.1975.JOIN$searchTaxon))

TAXA.KOPPEN.2030.JOIN = cbind.data.frame(GBIF.ALA.POINTS, KOPPEN.TAXA.2030)
str(unique(TAXA.KOPPEN.2030.JOIN$searchTaxon))

TAXA.KOPPEN.2050.JOIN = cbind.data.frame(GBIF.ALA.POINTS, KOPPEN.TAXA.2050)
str(unique(TAXA.KOPPEN.2050.JOIN$searchTaxon))

TAXA.KOPPEN.2070.JOIN = cbind.data.frame(GBIF.ALA.POINTS, KOPPEN.TAXA.2070)
str(unique(TAXA.KOPPEN.2070.JOIN$searchTaxon))


## Produce a "long" format table with records (rows) * Koppen (columns).
## Not sure the order matters, but this analysis is framed about Koppen so order by that
TAXA.KOPPEN.1975.JOIN  = TAXA.KOPPEN.1975.JOIN [with(TAXA.KOPPEN.1975.JOIN , order(Koppen)), ]
TAXA.KOPPEN.2030.JOIN  = TAXA.KOPPEN.2030.JOIN[with(TAXA.KOPPEN.2030.JOIN, order(Koppen_2030_A1B_CS)), ]
TAXA.KOPPEN.2050.JOIN  = TAXA.KOPPEN.2050.JOIN[with(TAXA.KOPPEN.2050.JOIN, order(Koppen_2050_A1B_CS)), ]
TAXA.KOPPEN.2070.JOIN  = TAXA.KOPPEN.2070.JOIN[with(TAXA.KOPPEN.2070.JOIN, order(Koppen_2070_A1B_CS)), ]


TAXA.KOPPEN.1975.COUNT = TAXA.KOPPEN.1975.JOIN[, c("Koppen", "searchTaxon")]
TAXA.KOPPEN.2030.COUNT = TAXA.KOPPEN.2030.JOIN[, c("Koppen_2030_A1B_CS", "searchTaxon")]
TAXA.KOPPEN.2050.COUNT = TAXA.KOPPEN.2050.JOIN[, c("Koppen_2050_A1B_CS", "searchTaxon")]
TAXA.KOPPEN.2070.COUNT = TAXA.KOPPEN.2070.JOIN[, c("Koppen_2070_A1B_CS", "searchTaxon")]


#########################################################################################################################
## Create a table with one row per species, and a column for every koppen  
KOPPEN.CAST.1975    = dcast(TAXA.KOPPEN.1975.COUNT, searchTaxon ~ Koppen)    ## can't use aggregation function on characters
KOPPEN.CAST.2030    = dcast(TAXA.KOPPEN.2030.COUNT, searchTaxon ~ Koppen_2030_A1B_CS)
KOPPEN.CAST.2050    = dcast(TAXA.KOPPEN.2050.COUNT, searchTaxon ~ Koppen_2050_A1B_CS)
KOPPEN.CAST.2070    = dcast(TAXA.KOPPEN.2070.COUNT, searchTaxon ~ Koppen_2070_A1B_CS)
dim(KOPPEN.CAST.1975)          ## so two koppen zones are not in the data at all....seems strange?
dim(KOPPEN.CAST.2070)


## Summ the records
koppen.list                  = as.list(names(KOPPEN.CAST.1975 [2:dim(KOPPEN.CAST.1975 )[2]]))
KOPPEN.CAST.1975$ALL_RECORDS = rowSums(KOPPEN.CAST.1975 [2:dim(KOPPEN.CAST.1975)[2]])
KOPPEN.CAST.2030$ALL_RECORDS = rowSums(KOPPEN.CAST.2030 [2:dim(KOPPEN.CAST.2030)[2]])
KOPPEN.CAST.2050$ALL_RECORDS = rowSums(KOPPEN.CAST.2050 [2:dim(KOPPEN.CAST.2050)[2]])
KOPPEN.CAST.2070$ALL_RECORDS = rowSums(KOPPEN.CAST.2070 [2:dim(KOPPEN.CAST.2070)[2]])


## Loop over all the Koppens and create a new column for each, which is the % of records in each koppen
for(kop in koppen.list){
  
  KOPPEN.CAST.1975[[paste("Percent_rec_in_", kop, sep = "")]] <- KOPPEN.CAST.1975 [[kop]]/KOPPEN.CAST.1975$ALL_RECORDS * 100
  
}


## 2030
for(kop in koppen.list){
  
  KOPPEN.CAST.2030[[paste("Percent_rec_in_", kop, sep = "")]] <- KOPPEN.CAST.2030 [[kop]]/KOPPEN.CAST.2030$ALL_RECORDS * 100
  
}


## 2050
for(kop in koppen.list){
  
  KOPPEN.CAST.2050[[paste("Percent_rec_in_", kop, sep = "")]] <- KOPPEN.CAST.2050 [[kop]]/KOPPEN.CAST.2050$ALL_RECORDS * 100
  
}


## 2070
for(kop in koppen.list){
  
  KOPPEN.CAST.2070[[paste("Percent_rec_in_", kop, sep = "")]] <- KOPPEN.CAST.2070 [[kop]]/KOPPEN.CAST.2070$ALL_RECORDS * 100
  
}


## Check
names(KOPPEN.CAST.1975)
View(KOPPEN.CAST.1975)


#########################################################################################################################
## Join to metadata
ALL.SPP  = MOD.ALL[, c("searchTaxon", "Origin", "Plant.type", "Number.of.growers", "Top_200")]


#########################################################################################################################
## Save both tables for safe keeping
write.csv(KOPPEN.MOD.SPP, "./data/base/HIA_LIST/COMBO/KOPPEN_TEST_SPP.csv", row.names = FALSE)
write.csv(KOPPEN.ALL.SPP, "./data/base/HIA_LIST/COMBO/KOPPEN_ALL_SPP.csv", row.names = FALSE)





#########################################################################################################################
## 3). CREATE A KOPPEN * SPP TABLE FOR 10 EXPERIMENTAL SPECIES 
#########################################################################################################################


# ## Check the contextual data for just the test species
# intersect(mod.compare.spp, COMBO.ALL$searchTaxon)
# MOD.COMP.NICHE = COMBO.ALL[COMBO.ALL$searchTaxon %in% intersect(mod.compare.spp, COMBO.ALL$searchTaxon), ]
# View(MOD.COMP.NICHE[, c(2:18)][with(MOD.COMP.NICHE[, c(2:18)], rev(order(COMBO.count))), ])
# 
# 
# ## So lets choose :
# ## 5 natives: 2 with many records, 3 with few
# ## 5 exotics: 2 with many records, 3 with few
# 
# 
# ## Create a list of just the species to compare model outputs
# mod.comp.spp = c('Lomandra longifolia', 'Dianella caerulea',  'Backhousia citriodora',   'Eucalyptus erythrocorys',  'Ficus brachypoda',
#                  'Cordyline australis', 'Murraya paniculata', 'Calodendrum capense',     'Liriope muscari',          'Ulmus parvifolia')
# 
# 
# ## Restrict the data frame to just the species we want
# MOD.COMP = COMBO.ALL[COMBO.ALL$searchTaxon %in% mod.comp.spp, ]
# MOD.ALL  = COMBO.ALL[COMBO.ALL$searchTaxon %in% unique(COMBO.RASTER.CONTEXT$searchTaxon), ]
# 
# MOD.SPP  = MOD.COMP[, c("searchTaxon", "Origin", "Plant.type", "Number.of.growers", "Top_200")]
# ALL.SPP  = MOD.ALL[, c("searchTaxon", "Origin", "Plant.type", "Number.of.growers", "Top_200")]
# 
# 
# ## Note that the number of records is a bit different (by one or two records). This is a versioning problem and not important, 
# ## but ideally should be the same.
# KOPPEN.MOD.SPP = join(MOD.SPP, KOPPEN.CAST)
# KOPPEN.ALL.SPP = join(ALL.SPP, KOPPEN.CAST)
# View(KOPPEN.MOD.SPP)
# View(KOPPEN.ALL.SPP)
# 
# 
# ## Create a shapefile for just the test species
# TEST.OCC = GBIF.ALA.POINTS[GBIF.ALA.POINTS$searchTaxon %in% KOPPEN.MOD.SPP$searchTaxon, ]
# dim(test)
# TEST.OCC = as.data.frame(TEST.OCC)
# 
# 
# ## Reorder the columns to see in ArcMap
# TEST.OCC  = select(TEST.OCC, searchTaxon, OBS,     lat,            lon,
#                    Origin,         Plant.type,  Top_200, scientificName, commonname, 
#                    taxonRank,      genus,       family,  TPL_binomial,   taxo_agree, 
#                    gbifID,         id,          cloc,    basisOfRecord,  locality, 
#                    establishmentMeans, institutionCode, datasetName, habitat, catalogNumber, 
#                    country, coordinateUncertaintyInMeters, geodeticDatum, year, month,                         
#                    day, eventDate, eventID, CULTIVATED, POLY_ID, SUB_CODE_7, REG_CODE_7)
# 
# 
# ## Rename the fields so that ArcMap can handle them
# TEST.OCC     = dplyr::rename(TEST.OCC, 
#                              TAXON     = searchTaxon,
#                              ORIGIN    = Origin,
#                              LAT       = lat,
#                              LON       = lon,
#                              TYPE      = Plant.type,
#                              T200      = Top_200,
#                              SC_NAME   = scientificName,
#                              COM_NAME  = commonname,
#                              RANK      = taxonRank,
#                              GENUS     = genus,
#                              FAMILY    = family,
#                              TPL       = TPL_binomial,
#                              TX_AGR    = taxo_agree,
#                              GBIF_ID   = gbifID,
#                              ALA_ID    = id,
#                              CLOC      = cloc,                     
#                              BASIS     = basisOfRecord,                
#                              LOCAL     = locality,                      
#                              ESTAB     = establishmentMeans,          
#                              INSTIT    = institutionCode,                
#                              DATASET   = datasetName,                      
#                              HABITAT   = habitat,                
#                              CAT_NO    = catalogNumber,                
#                              COUNTRY   = country,                
#                              COORD_UN  = coordinateUncertaintyInMeters,
#                              DATUM     = geodeticDatum,                 
#                              YEAR      = year,                          
#                              MNTH      = month,                        
#                              DAY       = day,                      
#                              DATE      = eventDate,                     
#                              EV_ID     = eventID,                      
#                              CULT      = CULTIVATED)
# 
# 
# ##
# TEST.OCC = SpatialPointsDataFrame(coords      = TEST.OCC[c("LON", "LAT")], 
#                                   data        = TEST.OCC,
#                                   proj4string = CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))
# 
# 
# ##
# writeOGR(obj = TEST.OCC, dsn = "./data/base/CONTEXTUAL", layer = "TEST.OCC", driver = "ESRI Shapefile")
# 
# 
# #########################################################################################################################
# ## Save both tables for safe keeping
# write.csv(KOPPEN.MOD.SPP, "./data/base/HIA_LIST/COMBO/KOPPEN_TEST_SPP.csv", row.names = FALSE)
# write.csv(KOPPEN.ALL.SPP, "./data/base/HIA_LIST/COMBO/KOPPEN_ALL_SPP.csv", row.names = FALSE)
# 
# 
# ## Now save .RData file for the next session...
# save.image("KOPPEN_JOIN.RData")





#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################