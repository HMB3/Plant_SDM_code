#########################################################################################################################
################################################ CLEAN RECORDS ########################################################## 
#########################################################################################################################


#########################################################################################################################
## This code combines the table of Koppen * SPP records, and uses Renee's species list to choose 10 species to trial
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
Koppen = readOGR("F:/green_cities_sdm/data/base/CONTEXTUAL/WC05_1975H_Koppen.shp", layer = "WC05_1975H_Koppen")
head(Koppen)
str(unique(Koppen$Koppen))
#plot(Koppen)


## Make a table of the Koppen descriptions
n.samples = 30
table.columns = c("Koppen", "Description")


## Create big table to store the data, with n rows = n samples
table.length         = length(table.columns)
KOP.LUT              = matrix(0, n.samples, table.length)
colnames(KOP.LUT)    = table.columns 
KOP.LUT              = as.data.frame(KOP.LUT)
KOP.LUT$Koppen       = as.character(sort(unique(Koppen$Koppen)))
str(KOP.LUT)


## Loop over all the Koppens and create a new column
## descr.list = as.character(unique(KOP.LUT$Description))
# for(i in 1:30){
#   
#   for(descr in descr.list){
#     
#     KOP.LUT[i, "Description"] = descr
#     
#   }
#   
# }


## Fill the table
KOP.LUT[1, "Koppen"]
KOP.LUT[1, "Description"] = "Tropical rain forest" 

KOP.LUT[2, "Koppen"]
KOP.LUT[2, "Description"] = "Tropical monsoon" 

KOP.LUT[3, "Koppen"]
KOP.LUT[3, "Description"] = "Tropical wet and dry savanna" 

KOP.LUT[4, "Koppen"]
KOP.LUT[4, "Description"] = "Hot semi-arid " 

KOP.LUT[5, "Koppen"]
KOP.LUT[5, "Description"] = "Cold semi-arid"

KOP.LUT[6, "Koppen"]
KOP.LUT[6, "Description"] = "Hot desert climate"

KOP.LUT[7, "Koppen"]
KOP.LUT[7, "Description"] = "Cold desert climate"

KOP.LUT[8, "Koppen"]
KOP.LUT[8, "Description"] = "Humid subtropical"

KOP.LUT[9, "Koppen"]
KOP.LUT[9, "Description"] = "Temperate oceanic"

KOP.LUT[10, "Koppen"]
KOP.LUT[10, "Description"] = "Subpolar oceanic"


##
KOP.LUT[11, "Koppen"]
KOP.LUT[11, "Description"] = "Hot-summer Mediterranean" 

KOP.LUT[12, "Koppen"]
KOP.LUT[12, "Description"] = "Warm-summer Mediterranean" 

KOP.LUT[13, "Koppen"]
KOP.LUT[13, "Description"] = "Cold-summer Mediterranean" 

KOP.LUT[14, "Koppen"]
KOP.LUT[14, "Description"] = "Monsoon-influenced humid subtropical" 

KOP.LUT[15, "Koppen"]
KOP.LUT[15, "Description"] = "Monsoon-influenced temperate oceanic"

KOP.LUT[16, "Koppen"]
KOP.LUT[16, "Description"] = "Monsoon-influenced subpolar oceanic"

KOP.LUT[17, "Koppen"]
KOP.LUT[17, "Description"] = "Hot-summer humid continental"

KOP.LUT[18, "Koppen"]
KOP.LUT[18, "Description"] = "Warm-summer humid continental"

KOP.LUT[19, "Koppen"]
KOP.LUT[19, "Description"] = "Subarctic climate"

KOP.LUT[20, "Koppen"]
KOP.LUT[20, "Description"] = "Cold subarctic climate"


##
KOP.LUT[21, "Koppen"]
KOP.LUT[21, "Description"] = "Mediterranean-influenced hot-summer humid continental" 

KOP.LUT[22, "Koppen"]
KOP.LUT[22, "Description"] = "Mediterranean-influenced warm-summer humid continental" 

KOP.LUT[23, "Koppen"]
KOP.LUT[23, "Description"] = "Mediterranean-influenced subarctic" 

KOP.LUT[24, "Koppen"]
KOP.LUT[24, "Description"] = "Mediterranean-influenced cold" 

KOP.LUT[25, "Koppen"]
KOP.LUT[25, "Description"] = "Monsoon-influenced hot-summer humid continental"

KOP.LUT[26, "Koppen"]
KOP.LUT[26, "Description"] = "Monsoon-influenced warm-summer humid continental"

KOP.LUT[27, "Koppen"]
KOP.LUT[27, "Description"] = "Monsoon-influenced subarctic climate"

KOP.LUT[28, "Koppen"]
KOP.LUT[28, "Description"] = "Monsoon-influenced cold subarctic"

KOP.LUT[29, "Koppen"]
KOP.LUT[29, "Description"] = "Ice cap"

KOP.LUT[30, "Koppen"]
KOP.LUT[30, "Description"] = "Tundra"


## check
head(KOP.LUT, 30)


##
write.csv(KOP.LUT, "./data/base/HIA_LIST/COMBO/KOPPEN_LUT.csv", row.names = FALSE)





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


## Check the field names :: are ALA and GBIF ID fields mutually exclusive?
test_spp    = subset(COMBO.RASTER.CONTEXT, searchTaxon == "Cycas revoluta" | searchTaxon == "Syzygium floribundum")


## Create spatial points dataframe
GBIF.ALA.POINTS = SpatialPointsDataFrame(coords      = COMBO.RASTER.CONTEXT[c("lon", "lat")], 
                                         data        = COMBO.RASTER.CONTEXT,
                                         proj4string = CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))

TEST.POINTS     = SpatialPointsDataFrame(coords      = test_spp[c("lon", "lat")], 
                                         data        = test_spp,
                                         proj4string = CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))





#########################################################################################################################
## 2). INTERSECT SPECIES RECORDS WITH CLIMATE ZONES
#########################################################################################################################


## Projection need to be the same for the spatial overlay 
CRS.new          <- CRS("+init=epsg:4326")        ## EPSG:3577
Koppen           = spTransform(Koppen, CRS.new)
GBIF.ALA.POINTS  = spTransform(GBIF.ALA.POINTS, CRS.new)
TEST.POINTS      = spTransform(TEST.POINTS, CRS.new)


## Check they are the same
projection(Koppen)
projection(GBIF.ALA.POINTS)
projection(TEST.POINTS)


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


## Loop over all the Koppens and create a new column
for(kop in koppen.list){
  
  KOPPEN.CAST[[paste("Percent_rec_in_", kop, sep = "")]] <- KOPPEN.CAST[[kop]]/KOPPEN.CAST$ALL_RECORDS * 100
  
}


## Check
names(KOPPEN.CAST)
View(KOPPEN.CAST)


#########################################################################################################################
## Now try intersecting the environmental matrix with the koppen * records
# test_env           = subset(OCC.ENV, searchTaxon == "Cycas revoluta" | searchTaxon == "Syzygium floribundum")
# TAXA.KOP.ENV       = join(TAXA.KOPPEN.JOIN, test_env)





#########################################################################################################################
## 3). CREATE A KOPPEN * SPP TABLE FOR 10 EXPERIMENTAL SPECIES 
#########################################################################################################################


## Just look at the experimental species
mod.compare.spp = trimws(sort(unique(c(renee.full$Species, MQ.glasshouse$Species))))


## Check the contextual data for just the test species
intersect(mod.compare.spp, COMBO.ALL$searchTaxon)
MOD.COMP.NICHE = COMBO.ALL[COMBO.ALL$searchTaxon %in% intersect(mod.compare.spp, COMBO.ALL$searchTaxon), ]
View(MOD.COMP.NICHE[, c(2:18)][with(MOD.COMP.NICHE[, c(2:18)], rev(order(COMBO.count))), ])
#View(MOD.COMP.NICHE[, c(2:18)][with(MOD.COMP.NICHE[, c(2:18)], rev(order(Number.of.growers))), ])


## So lets choose :
## 5 natives: 2 with many records, 3 with few
## 5 exotics: 2 with many records, 3 with few


## Create a list of just the species to compare model outputs
mod.comp.spp = c('Lomandra longifolia', 'Dianella caerulea',  'Backhousia citriodora',   'Eucalyptus erythrocorys',  'Ficus brachypoda',
                 'Cordyline australis', 'Murraya paniculata', 'Calodendrum capense',     'Liriope muscari',          'Ulmus parvifolia')

##
MOD.COMP = COMBO.ALL[COMBO.ALL$searchTaxon %in% mod.comp.spp, ]
View(MOD.COMP[, c(2:18)][with(MOD.COMP[, c(2:18)], rev(order(COMBO.count))), ])
MOD.SPP = MOD.COMP[, c("searchTaxon", "Origin", "Plant.type", "Number.of.growers", "Top_200")]


## Note that the number of records is a bit different (by one or two records). This is a versioning problem and not important, 
## but ideally should be the same
KOPPEN.MOD.SPP = join(MOD.SPP, KOPPEN.CAST)
View(KOPPEN.MOD.SPP)


#########################################################################################################################
## Save the table for safe keeping
write.csv(KOPPEN.MOD.SPP, "./data/base/HIA_LIST/COMBO/KOPPEN_TEST_SPP.csv", row.names = FALSE)





#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################