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
load("./data/base/HIA_LIST/COMBO/COMBO_RASTER_CONTEXT_1601_2018.RData")


## Read in the Koppen shapefile :: 1975 centred data
## https://www.climond.org/Core/Authenticated/KoppenGeiger.aspx
Koppen = readOGR("F:/green_cities_sdm/data/base/CONTEXTUAL/WC05_1975H_Koppen_Shapefile/WC05_1975H_Koppen_Kriticos_2012.shp", 
                 layer = "WC05_1975H_Koppen_Kriticos_2012")
head(Koppen)
str(unique(Koppen$Koppen))
#plot(Koppen)


#########################################################################################################################
## 1). CREATE SPATIAL POINTS DATA FRAMES FOR OVERLAY
#########################################################################################################################


## Check the unique IDs match
COMBO.RASTER.CONTEXT$OBS <- 1:nrow(COMBO.RASTER.CONTEXT)
dim(COMBO.RASTER.CONTEXT)[1];length(COMBO.RASTER.CONTEXT$OBS)  
head(COMBO.RASTER.CONTEXT, 10)[, c("OBS", "searchTaxon", "lat", "lon")]


## Just get the environmetnal columns
OCC.ENV    <- select(COMBO.RASTER.CONTEXT, searchTaxon, OBS,
                     Annual_mean_temp,     Mean_diurnal_range, Isothermality,      Temp_seasonality,   Max_temp_warm_month, 
                     Mean_temp_wet_qu,     Mean_temp_dry_qu,   Temp_annual_range,  Mean_temp_warm_qu,  Mean_temp_cold_qu, 
                     Min_temp_cold_month,
                     Annual_precip,        Precip_wet_month,   Precip_dry_month,   Precip_seasonality, Precip_wet_qu,     
                     Precip_dry_qu,        Precip_warm_qu,     Precip_col_qu)


## Create spatial points dataframe
GBIF.ALA.POINTS = SpatialPointsDataFrame(coords      = COMBO.RASTER.CONTEXT[c("lon", "lat")], 
                                         data        = COMBO.RASTER.CONTEXT,
                                         proj4string = CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))


##
dim(GBIF.ALA.POINTS)
length(unique(GBIF.ALA.POINTS$searchTaxon))


##
test = GBIF.ALA.POINTS[GBIF.ALA.POINTS$searchTaxon %in% mod.compare.spp, ]





#########################################################################################################################
## 2). INTERSECT SPECIES RECORDS WITH CLIMATE ZONES
#########################################################################################################################


## The projections need to match for the spatial overlay 
CRS.new          <- CRS("+init=epsg:4326")        ## EPSG:3577
Koppen           = spTransform(Koppen, CRS.new)
GBIF.ALA.POINTS  = spTransform(GBIF.ALA.POINTS, CRS.new)


## Check they are the same
projection(Koppen)
projection(GBIF.ALA.POINTS)


## Now, intersect the species with the data 
KOPPEN.TAXA = over(GBIF.ALA.POINTS, Koppen) 
dim(KOPPEN.TAXA)
head(GBIF.ALA.POINTS[1])


## Then join the koppen classification onto the data
TAXA.KOPPEN.JOIN = cbind.data.frame(GBIF.ALA.POINTS, KOPPEN.TAXA)
str(unique(TAXA.KOPPEN.JOIN$searchTaxon))


## Produce a "long" format table with records (rows) * Koppen (columns)
## Not sure the order matters, but this analysis is framed about Koppen so order on this
TAXA.KOPPEN.JOIN  = TAXA.KOPPEN.JOIN[with(TAXA.KOPPEN.JOIN, order(Koppen)), ]
TAXA.KOPPEN.COUNT = TAXA.KOPPEN.JOIN[, c("Koppen", "searchTaxon")]


#########################################################################################################################
## Create a table with one row per species, and a column for every koppen  
KOPPEN.CAST             = dcast(TAXA.KOPPEN.COUNT, searchTaxon ~ Koppen)    ## can't use aggregation function on characters
dim(KOPPEN.CAST)          ## so two koppen zones are not in the data at all....seems strange?


## Summ the records
koppen.list             = as.list(names(KOPPEN.CAST[2:dim(KOPPEN.CAST)[2]]))
KOPPEN.CAST$ALL_RECORDS = rowSums(KOPPEN.CAST[2:dim(KOPPEN.CAST)[2]])


## Loop over all the Koppens and create a new column for each, which is the % of records in each koppen
for(kop in koppen.list){
  
  KOPPEN.CAST[[paste("Percent_rec_in_", kop, sep = "")]] <- KOPPEN.CAST[[kop]]/KOPPEN.CAST$ALL_RECORDS * 100
  
}


## Check
names(KOPPEN.CAST)
View(KOPPEN.CAST)





#########################################################################################################################
## 3). CREATE A KOPPEN * SPP TABLE FOR 10 EXPERIMENTAL SPECIES 
#########################################################################################################################


## Just look at the experimental species
mod.compare.spp = trimws(sort(unique(c(renee.full$Species, MQ.glasshouse$Species))))


## Check the contextual data for just the test species
intersect(mod.compare.spp, COMBO.ALL$searchTaxon)
MOD.COMP.NICHE = COMBO.ALL[COMBO.ALL$searchTaxon %in% intersect(mod.compare.spp, COMBO.ALL$searchTaxon), ]
View(MOD.COMP.NICHE[, c(2:18)][with(MOD.COMP.NICHE[, c(2:18)], rev(order(COMBO.count))), ])


## So lets choose :
## 5 natives: 2 with many records, 3 with few
## 5 exotics: 2 with many records, 3 with few


## Create a list of just the species to compare model outputs
mod.comp.spp = c('Lomandra longifolia', 'Dianella caerulea',  'Backhousia citriodora',   'Eucalyptus erythrocorys',  'Ficus brachypoda',
                 'Cordyline australis', 'Murraya paniculata', 'Calodendrum capense',     'Liriope muscari',          'Ulmus parvifolia')


## Restrict the data frame to just the species we want
MOD.COMP = COMBO.ALL[COMBO.ALL$searchTaxon %in% mod.comp.spp, ]
MOD.ALL  = COMBO.ALL[COMBO.ALL$searchTaxon %in% unique(COMBO.RASTER.CONTEXT$searchTaxon), ]

MOD.SPP  = MOD.COMP[, c("searchTaxon", "Origin", "Plant.type", "Number.of.growers", "Top_200")]
ALL.SPP  = MOD.ALL[, c("searchTaxon", "Origin", "Plant.type", "Number.of.growers", "Top_200")]


## Note that the number of records is a bit different (by one or two records). This is a versioning problem and not important, 
## but ideally should be the same.
KOPPEN.MOD.SPP = join(MOD.SPP, KOPPEN.CAST)
KOPPEN.ALL.SPP = join(ALL.SPP, KOPPEN.CAST)
View(KOPPEN.MOD.SPP)
View(KOPPEN.ALL.SPP)


## Create a shapefile for just the test species
TEST.OCC = GBIF.ALA.POINTS[GBIF.ALA.POINTS$searchTaxon %in% KOPPEN.MOD.SPP$searchTaxon, ]
dim(test)
TEST.OCC = as.data.frame(TEST.OCC)


## Reorder the columns to see in ArcMap
TEST.OCC  = select(TEST.OCC, searchTaxon, OBS,     lat,            lon,
                   Origin,         Plant.type,  Top_200, scientificName, commonname, 
                   taxonRank,      genus,       family,  TPL_binomial,   taxo_agree, 
                   gbifID,         id,          cloc,    basisOfRecord,  locality, 
                   establishmentMeans, institutionCode, datasetName, habitat, catalogNumber, 
                   country, coordinateUncertaintyInMeters, geodeticDatum, year, month,                         
                   day, eventDate, eventID, CULTIVATED, POLY_ID, SUB_CODE_7, REG_CODE_7)


## Rename the fields so that ArcMap can handle them
TEST.OCC     = dplyr::rename(TEST.OCC, 
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
TEST.OCC = SpatialPointsDataFrame(coords      = TEST.OCC[c("LON", "LAT")], 
                                  data        = TEST.OCC,
                                  proj4string = CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))


##
writeOGR(obj = TEST.OCC, dsn = "./data/base/CONTEXTUAL", layer = "TEST.OCC", driver = "ESRI Shapefile")


#########################################################################################################################
## Save both tables for safe keeping
write.csv(KOPPEN.MOD.SPP, "./data/base/HIA_LIST/COMBO/KOPPEN_TEST_SPP.csv", row.names = FALSE)
write.csv(KOPPEN.ALL.SPP, "./data/base/HIA_LIST/COMBO/KOPPEN_ALL_SPP.csv", row.names = FALSE)


## Now save .RData file for the next session...
save.image("KOPPEN_JOIN.RData")





#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################