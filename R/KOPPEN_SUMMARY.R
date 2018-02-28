#########################################################################################################################
################################################ CLEAN RECORDS ########################################################## 
#########################################################################################################################


## This code creates individual shapefiles for each species, so we can create a list of the spatial outliers manually.
## It also makes a table of koppen zones * species occurrences

#  which Koppen zones a given species is primarily found in (by overlaying the occurrences with the Koppen data). We could then group species according to their distribution in Koppen zones. For each of these major groupings, we would then come up with a list of variables to use in the models, based on our understanding of the climate of those zones. So, for example, we would be able to identify a group of species that are primarily Montane. We would then decide on which variables to use for these. A second group of species might be primarily arid/semi-arid, and so on.
# 
# So the point of overlaying the occurrences with the Koppen is simply to identify what zones the species falls in. We would do that by calculating the % of a species' records in each zone. I have done similar things before, and I guarantee it still wont be straightforward with grouping species, but it will be a start.


## Load the niche data and the records from step 4, after cruncing ALA and GBIF, but prior to cleaning (GBIF.ALA.COMBO.HIA)
## Load Koppen zones 
source('./R/HIA_LIST_MATCHING.R')
#load("./data/base/HIA_LIST/COMBO/GBIF_ALA_COMBO_PRE_CLEAN.RData")
load("./data/base/HIA_LIST/COMBO/COMBO_RASTER_CONTEXT_1601_2018.RData")


## Read in the Koppen shapefile :: 1975 centred data
## https://www.climond.org/Core/Authenticated/KoppenGeiger.aspx
Koppen = readOGR("F:/green_cities_sdm/data/base/CONTEXTUAL/WC05_1975H_Koppen.shp", layer = "WC05_1975H_Koppen")
head(Koppen)
str(unique(Koppen$Koppen))
plot(Koppen)


#########################################################################################################################
## 1). CHECK ALA/GBIF OVERLAP
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


## Read in the Koppen shapefile :: 1975 centred data
## https://www.climond.org/Core/Authenticated/KoppenGeiger.aspx
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
KOPPEN.TAXA = over(GBIF.ALA.POINTS, Koppen) 
dim(KOPPEN.TAXA)
head(GBIF.ALA.POINTS[1])


## Join the koppen classification onto the data
TAXA.KOPPEN.JOIN = cbind.data.frame(GBIF.ALA.POINTS, KOPPEN.TAXA)
str(unique(TAXA.KOPPEN.JOIN$searchTaxon))


## First, produce a "long" format table with records (rows) * Koppen (columns)
## Not sure the order matters, but this analysis is framed about Koppen
TAXA.KOPPEN.JOIN  = TAXA.KOPPEN.JOIN[with(TAXA.KOPPEN.JOIN, order(Koppen)), ]
TAXA.KOPPEN.COUNT = TAXA.KOPPEN.JOIN[, c("Koppen", "searchTaxon")]


#########################################################################################################################
## We want a table with one row per species, and a column for every koppen  
KOPPEN.CAST             = dcast(TAXA.KOPPEN.COUNT, searchTaxon ~ Koppen)    ## can't use aggregation function on characters
dim(KOPPEN.CAST)          ## so two koppen zones are not in the data at all....seems strange?


## Loop over all the Koppens and create a new column :: 
koppen.list             = as.list(names(KOPPEN.CAST[2:dim(KOPPEN.CAST)[2]]))
KOPPEN.CAST$ALL_RECORDS = rowSums(KOPPEN.CAST[2:dim(KOPPEN.CAST)[2]])


## Loop over all the Koppens and create a new column
for(kop in koppen.list){
  
  KOPPEN.CAST[[paste("Percent_rec_in_", kop, sep = "")]] <- KOPPEN.CAST[[kop]]/KOPPEN.CAST$ALL_RECORDS
  
}


##
names(KOPPEN.CAST)View(KOPPEN.CAST)



#########################################################################################################################
## Now try intersecting the environmental matrix with the koppen * records
test_env           = subset(OCC.ENV, searchTaxon == "Cycas revoluta" | searchTaxon == "Syzygium floribundum")
TAXA.KOP.ENV       = join(TAXA.KOPPEN.JOIN, test_env)





#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################



