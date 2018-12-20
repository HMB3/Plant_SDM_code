#########################################################################################################################
########################################  CREATE SPECIES LISTS ########################################################## 
#########################################################################################################################


## This code takes the raw list of plants with 25 or more growers supplied by Matt Plumber and Anthony Manea, and
## then cleans the list as best as possible in R to use the species binomial as the unit for downloading and analysis.

## All the cleaining methods will throw up some anomalies, which need to be tracked, and checked with the team for how
## each case is treated (see outstanding tasks at the bottom).


#########################################################################################################################
## Setup for project 


## Load only the packages needed for the analysis
p <- c('ff',    'things',    'raster',        'dismo',        'sp',           'latticeExtra', 'data.table', 
       'rgdal', 'rgeos',     'gdalUtils',     'rmaxent',      'readr',        'plyr',         'dplyr',        
       'tidyr', 'readr',     'rnaturalearth', 'rasterVis',    'RColorBrewer', 'latticeExtra', 'parallel',     
       'ALA4R', 'stringr',   'Taxonstand',    'CoordinateCleaner', 'gsubfn',  'PerformanceAnalytics',
       'rvest', 'magrittr',  'devtools',      'ggplot2',      'reshape2',     'rmarkdown', 'flexdashboard', 'shiny', 'rgbif',
       'ENMeval', 'tibble',  'ncdf4',         'Cairo', 'velox', 'taxonlookup', 'kgc', 'maptools', 'DataCombine', 'mgcv', 'rsq')


## Require packages
sapply(p, require, character.only = TRUE)
source_gist('26e8091f082f2b3dd279',             filename = 'polygonizer.R')
source_gist('c6a1cb61b8b6616143538950e6ec34aa', filename = 'hatch.R')


## Source functions, and set temporary directory (for both raster files and just generally)
source('./R/GREEN_CITIES_FUNCTIONS.R')
source('./R/MAXENT_FUNCTIONS.R')
source('./R/MAPPING_FUNCTIONS.R')
source('./R/HIA_CLEAN_MATCHING.R')
rasterOptions(tmpdir = file.path('./RTEMP')) 





#########################################################################################################################
## 1). READ IN SPATIAL DATA
#########################################################################################################################


#########################################################################################################################
## Set coordinate system definitions :: best to minimise the number of projection used in this project
## Also get rid of the '+init=ESRI", which are not compatible with some systems 
CRS.MOL      <- CRS('+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs')
CRS.MOL.SDM  <- CRS('+init=ESRI:54009 +proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs84=0,0,0')
CRS.WGS.84   <- CRS("+init=epsg:4326")
CRS.AUS.ALB  <- CRS("+init=EPSG:3577")
ALB.CONICAL  <- CRS('+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')


#########################################################################################################################
## Read in spatial data once, rather than in each script
aus           = readRDS("./data/base/CONTEXTUAL/aus_states.rds")
LAND          = readRDS("./data/base/CONTEXTUAL/LAND_world.rds")
SUA_2016      = readRDS("./data/base/CONTEXTUAL/SUA/SUA_2016_AUST.rds")
SUA_2016      = SUA_2016[order(SUA_2016$SUA_NAME16),]
SUA_2016      = SUA_2016 %>%
  spTransform(ALB.CONICAL)


# SUA_2011      = readOGR('data/base/CONTEXTUAL/SUA/SUA_2011_AUST.shp')
Koppen        = readRDS('data/base/CONTEXTUAL/Koppen_1975.rds')
Koppen_aus    = readRDS('data/base/CONTEXTUAL/KOPPEN_AUS.rds')
Koppen_zones  = unique(readOGR('data/base/CONTEXTUAL/WC05_1975H_Koppen_Shapefile/WC05_1975H_Koppen_Kriticos_2012.shp')@data[, 1:2])
Koppen_1975   = raster('data/Koppen_1000m_Mollweide54009.tif')
AUS_RAIN      = readRDS('data/base/CONTEXTUAL/BOM/BOM_RAIN_AGG.rds')

IN.SUA        = readRDS("./data/base/CONTEXTUAL/IN_SUA_AUS.rds")
SUA.16        = readRDS("./data/base/CONTEXTUAL/SUA/SUA_2016_AUST.rds")
LGA           = readRDS("./data/base/CONTEXTUAL/LGA.rds")
AUS           = readRDS("./data/base/CONTEXTUAL/aus_states.rds")
ALL.SUA.POP   = read.csv("./data/base/CONTEXTUAL/ABS_SUA_POP.csv", stringsAsFactors = FALSE)
URB.POP       = read.csv("./data/base/CONTEXTUAL/ABS_URBAN_CENTRE_POP.csv", stringsAsFactors = FALSE)


## Create centroids for SUAs
# class(SUA.16)
# plot(SUA.16)
# writeSpatialShape(SUA.16, "SUA.16")
# cents <- coordinates(SUA.16)
# cents <- SpatialPointsDataFrame(coords = cents, data = SUA.16@data, 
#                                 proj4string = CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"))
# points(cents, col = "Blue")
# writeSpatialShape(cents, "cents")
# 
# SUA_centroids <- coordinates(SUA.16)
# SUA_location  = data.frame(SUA_centroids)
# SUA_location  = cbind(SUA.16$SUA_NAME16, SUA_location)
# names(SUA_location)   = c("SUA","rndCoord.lon", "rndCoord.lat")
# SUA_location$rndCoord.lon = RoundCoordinates(SUA_location$rndCoord.lon)
# SUA_location$rndCoord.lat = RoundCoordinates(SUA_location$rndCoord.lat)
# 
# 
# points(centroids, pch = 3, col = "Red")
# Kop.loc <- data.frame(SUA_location, ClimateZ = LookupCZ(SUA_location))
# Kop.loc = Kop.loc[c("SUA", "ClimateZ")]


## Load template rasters
template.raster = raster("./data/template_hasData.tif")
template.cells  = readRDS("./data/hasData_cells.rds")
load("./data/base/CONTEXTUAL/urbanareas.rda")


#########################################################################################################################
## Read in SUA shapefile and convert columns to numeric and character
# SUA_2016  = readOGR("./data/base/CONTEXTUAL/SUA/SUA_2016_AUST.shp",
#                    layer = "SUA_2016_AUST", stringsAsFactors = FALSE)
# SUA_2016$SUA_CODE16 <- as.numeric(as.character(SUA_2016$SUA_CODE16))
# SUA_2016$SUA_NAME16 <- as.character(SUA_2016$SUA_NAME16)
# class(SUA_2016$SUA_CODE16)
# class(SUA_2016$SUA_NAME16)
# head(SUA_2016)
# 
# 
# unique(SUA_2016$SUA_NAME16)
# saveRDS(SUA_2016 , file = paste("./data/base/CONTEXTUAL/SUA/SUA_2016_AUST.rds"))


###################################################################################################################
## Rasterize the SUA shapefile
# areal_unit = readRDS("./data/base/CONTEXTUAL/SUA/SUA_2016_AUST.rds") %>%
#   spTransform(ALB.CONICAL)
# 
# ##  Rasterize shapefile
# message('rasterizing SUA shapefile')
# areal_unit = areal_unit[order(areal_unit$SUA_NAME16),]
# f <- tempfile()
# 
# writeOGR(areal_unit, tempdir(), basename(f), 'ESRI Shapefile')
# template <- raster('./output/maxent/SUA_TREES_ANALYSIS/Acacia_retinodes/full/Acacia_retinodes_ac85bi30.tif')
# 
# areal_unit_rast <- gdal_rasterize(
#   normalizePath(paste0(f, '.shp')), 
#   
#   'H:/green_cities_sdm/data/base/CONTEXTUAL/SUA/SUA_2016_AUST.tif', tr=res(template),
#   te = c(bbox(template)), a = 'SUA_CODE16', a_nodata = 0, init = 0, ot = 'UInt16', output_Raster = TRUE)
# 
# areal_unit_vec <- c(areal_unit_rast[])
# summary(areal_unit_vec)


## Save RDS for areal_unit_rast and areal_unit_vec
# saveRDS(areal_unit_rast, file = paste("./data/base/CONTEXTUAL/SUA/SUA_2016_RAST.rds"))
# saveRDS(areal_unit_vec, file = paste("./data/base/CONTEXTUAL/SUA/SUA_2016_VEC.rds"))

SUA_2016_rast = readRDS("./data/base/CONTEXTUAL/SUA/SUA_2016_RAST.rds")
SUA_2016_vec  = readRDS("./data/base/CONTEXTUAL/SUA/SUA_2016_VEC.rds")
summary(SUA_2016_vec)




#########################################################################################################################
## 2). READ IN RASTER COMPONENTS
#########################################################################################################################


#########################################################################################################################
## Now create the variables needed to access current environmental conditions + their names in the functions
sdm.predictors <- c("Annual_mean_temp",    "Mean_diurnal_range",  "Isothermality",      "Temp_seasonality",  
                    "Max_temp_warm_month", "Min_temp_cold_month", "Temp_annual_range",  
                    "Mean_temp_wet_qu",    "Mean_temp_dry_qu",    
                    "Mean_temp_warm_qu",   "Mean_temp_cold_qu",   
                    
                    "Annual_precip",       "Precip_wet_month",   "Precip_dry_month",    "Precip_seasonality",  
                    "Precip_wet_qu",       "Precip_dry_qu",      
                    "Precip_warm_qu",      "Precip_col_qu")

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

## A-priori worldclim predictors
sdm.select     <- c("Annual_mean_temp", "Temp_seasonality",    "Max_temp_warm_month", "Min_temp_cold_month",
                    "Annual_precip",    "Precip_seasonality",  
                    "Precip_wet_month", "Precip_dry_month")


grid.names = c('Annual_mean_temp',    'Mean_diurnal_range',  'Isothermality',      'Temp_seasonality', 
               'Max_temp_warm_month', 'Min_temp_cold_month', 'Temp_annual_range',  'Mean_temp_wet_qu',
               'Mean_temp_dry_qu',    'Mean_temp_warm_qu',   'Mean_temp_cold_qu',  'Annual_precip',
               'Precip_wet_month',    'Precip_dry_month',    'Precip_seasonality', 'Precip_wet_qu',
               'Precip_dry_qu',       'Precip_warm_qu',      'Precip_col_qu')


#########################################################################################################################
## Create a raster stack of current global environmental conditions
world.grids.current = stack(
  file.path('./data/base/worldclim/world/0.5/bio/current',
            sprintf('bio_%02d', 1:19)))


## Create a raster stack of current Australian environmental conditions
inventory.grids.current = stack(
  file.path('./data/base/worldclim/aus/1km/bio/current/WGS/', 
            sprintf('bio_%02d.tif', 1:19))) 


## Also get the PET raster
PET               = raster("./data/base/worldclim/world/1km/pet_he_yr1.tif")


## Create a velox raster of the world raster stack
# vx.world <- velox(world.grids.current)
# writeRaster(vx.world, './data/base/worldclim/world/0.5/bio/current/vx_world.tif')
# vx.world  = raster('./data/base/worldclim/world/0.5/bio/current/vx_world.tif')
# i  <- match(sdm.predictors, sdm.predictors)
# ff <- file.path('./data/base/worldclim/world/0.5/bio/current',
#                 sprintf('bio_%02d.tif', i))
# env.grids.current = stack(sub('0.5', '1km', ff))


## Name the grids :: these should be indentical
# names(world.grids.current) <- sdm.predictors
# identical(names(world.grids.current),sdm.predictors)


h <- read_html('http://www.worldclim.org/cmip5_30s') 
gcms <- h %>% 
  html_node('table') %>% 
  html_table(header = TRUE) %>% 
  filter(rcp85 != '')

id.50 <- h %>% 
  html_nodes(xpath = "//a[text()='bi']/@href") %>% 
  as.character %>% 
  grep('85bi50', ., value = TRUE) %>% 
  basename %>% 
  sub('\\.zip.', '', .)

id.70 <- h %>% 
  html_nodes(xpath = "//a[text()='bi']/@href") %>% 
  as.character %>% 
  grep('85bi70', ., value = TRUE) %>% 
  basename %>% 
  sub('\\.zip.', '', .)


## Work around because the 2030 data comes from CC in Aus, not worldclim
id.30 = gsub("50", "30", id.50)


## Create the IDs
gcms.30 <- cbind(gcms, id.30)
gcms.30$GCM = sub(" \\(#\\)", "", gcms$GCM)

gcms.50 <- cbind(gcms, id.50)
gcms.50$GCM = sub(" \\(#\\)", "", gcms$GCM)  ## sub replaces first instance in a string, gsub = global

gcms.70 <- cbind(gcms, id.70)
gcms.70$GCM = sub(" \\(#\\)", "", gcms$GCM) 


## Now create the scenario lists across which to loop
gcms.50 ; gcms.70 ; gcms.30


## Just get the 6 models picked by CSIRO for Australia, for 2030, 2050 and 2070
scen_2030 = c("mc85bi30", "no85bi30", "ac85bi30", "cc85bi30", "gf85bi30", "hg85bi30")
scen_2050 = c("mc85bi50", "no85bi50", "ac85bi50", "cc85bi50", "gf85bi50", "hg85bi50")
scen_2070 = c("mc85bi70", "no85bi70", "ac85bi70", "cc85bi70", "gf85bi70", "hg85bi70")


## Then create a stack of current environmental conditions outside the function, and an Australia shapefile for the mapping later...
# aus <- ne_states(country = 'Australia') %>% 
#   subset(!grepl('Island', name))
# 
# shapefile = aus


#########################################################################################################################
## Create a raster stack of current Australian environmental conditions, and divide the current environmental grids by 10
aus.grids.current <- stack(
  file.path('./data/base/worldclim/aus/1km/bio/current',   ## ./data/base/worldclim/aus/1km/bio
            sprintf('bio_%02d.tif', 1:19)))

for(i in 1:11) {
  
  ## simple loop
  message(i)
  aus.grids.current[[i]] <- aus.grids.current[[i]]/10
  
}   





#########################################################################################################################
## 3). READ IN SPECIES LISTS
#########################################################################################################################


#########################################################################################################################
## Read in the niche data
COMBO.NICHE.CONTEXT = readRDS("./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT_APRIL_2018_STANDARD_CLEAN.rds")
CLEAN.NICHE.CONTEXT = readRDS("./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT_APRIL_2018_COORD_CLEAN.rds")
MAXENT.CHECK        = read.csv("./output/maxent/CHECK_SPP_MAPS_BIAS_0310_2018.csv", stringsAsFactors = FALSE)
OVERALL.LOSS        = read.csv("./output/tables/OVERALL_LOSS.csv",                  stringsAsFactors = FALSE)
APPENDIX            = read.csv("./data/base/HIA_LIST/COMBO/Appendix_table.csv",     stringsAsFactors = FALSE)
AUS.FROGS           = read.csv("./data/base/Frog/AUS_FROGS.csv",                    stringsAsFactors = FALSE)
AUS.FROGS           = AUS.FROGS$Scientific.Name

HOLLOW.SPP          = read.csv("./data/base/HIA_LIST/Tree_hollows/Hollow_species.csv", stringsAsFactors = FALSE)
HOLLOW.SPP          = HOLLOW.SPP$Hollow_species

str(unique(COMBO.NICHE.CONTEXT$searchTaxon))  ## long enough
background          = readRDS("./data/base/HIA_LIST/COMBO/SDM_DATA_CLEAN_052018.rds")


# writeOGR(obj    = background, 
#          dsn    = "./data/base/HIA_LIST/COMBO", 
#          layer  = paste0('background_SUA'),
#          driver = "ESRI Shapefile", overwrite_layer = TRUE)


#########################################################################################################################
## Horticultural lists ::
## This evergreen list (HIA.list) derives from all species and varieties sold anywhere in Australia in the last 5 years. Anthony Manea cleaned 
## up the data and cross-linked to growth form and exotic/native status and derived a list of ~1000 species that are the 
## Most commonly sold, covering the right ratio of growth forms, regional representation and native/exotic
HIA.list            = read.csv("./data/base/HIA_LIST/HIA/GREEN_CITIES_DRAFT_LIST_2709_2017.csv", stringsAsFactors = FALSE)
HIA.RAW             = read.csv("./data/base/HIA_LIST/HIA/HIA_ORIGINAL_RAW.csv",                  stringsAsFactors = FALSE)
CLEAN.list          = read.csv("./data/base/HIA_LIST/HIA/HIA.CLEAN.csv",                         stringsAsFactors = FALSE)
GROWING             = read.csv("./data/base/HIA_LIST/HIA/database_aus_sp_growing.csv",           stringsAsFactors = FALSE)
MOD_2               = read.csv("./data/base/HIA_LIST/HIA/MOD2_LIST.csv",                         stringsAsFactors = FALSE)
Manuel              = read.csv("./data/base/HIA_LIST/MANUEL/species_module2.csv",                stringsAsFactors = FALSE)
MOD.2.3             = read.csv("./data/base/HIA_LIST/HIA/MODULE_2_3.csv",                        stringsAsFactors = FALSE)
MOD.3               = read.csv("./data/base/HIA_LIST/HIA/MODULE_2_3.csv",                        stringsAsFactors = FALSE)
APNI                = readRDS("./data/base/HIA_LIST/ALA/APNI_LIST.rds")

RISK.LIST           = read.csv("./data/base/HIA_LIST/HIA/RISK_LIST.csv",                         stringsAsFactors = FALSE)
RISK.BINOMIAL.CLEAN = read.csv("./data/base/HIA_LIST/HIA/RISK_BINOMIAL_DF.csv",                  stringsAsFactors = FALSE)
MAXENT.RATING       = read.csv("./output/maxent/MAXENT_RATING_26_2018.csv",                      stringsAsFactors = FALSE)
#MXT.CHECK           = read.csv("./output/maxent/CHECK_SPP_MAPS_BIAS_0310_2018.csv",              stringsAsFactors = FALSE)
MXT.CHECK           = read.csv("./output/maxent/MAXNET_ORIGIN_RESULTS.csv",                      stringsAsFactors = FALSE)
native.good.models  = subset(MXT.CHECK, Origin == "Native" & Check.map <=2)$searchTaxon
native.sua          = subset(MXT.CHECK, Origin == "Native")$searchTaxon

INV.CHECK           = read.csv("./output/maxent/TREES_INVENTORY_RESULTS.csv",                    stringsAsFactors = FALSE)
INV.SPP             = INV.CHECK$searchTaxon 


URBAN.FOREST        = read.csv("./data/base/HIA_LIST/URBAN/URBAN_FOREST.csv",                    stringsAsFactors = FALSE)
URBAN.FOR.SPP       = URBAN.FOREST$Species 
TI.LIST             = read.csv("./data/base/HIA_LIST/COMBO/ALE_TREE_SPP_LIST.csv",               stringsAsFactors = FALSE)
TI.XY               = read.csv("./data/base/HIA_LIST/COMBO/ALE_TREE_SPP_XY.csv",                 stringsAsFactors = FALSE)
SPP.BIAS            = read.csv("./output/maxent/SPP_BOUNDARY_BIAS.csv", stringsAsFactors = FALSE)


## Check Ale's data
TI.LIST = na.omit(TI.LIST)
TI.LIST = TI.LIST[with(TI.LIST, rev(order(Plantings))), ]
head(TI.LIST$searchTaxon, 10)


## Rename tree inventory data
TI.XY = TI.XY[c("searchTaxon", "POINT_X", "POINT_Y", "FILENAME")]
names(TI.XY)[names(TI.XY) == 'FILENAME'] <- 'INVENTORY'
names(TI.XY)[names(TI.XY) == 'POINT_X']  <- 'lon'
names(TI.XY)[names(TI.XY) == 'POINT_Y']  <- 'lat'


## Remove the gunk
TI.XY$INVENTORY = gsub("_trees.shp",      "", TI.XY$INVENTORY)
TI.XY$INVENTORY = gsub("_trees_.shp",     "", TI.XY$INVENTORY)
TI.XY$INVENTORY = gsub("_trees_sign.shp", "", TI.XY$INVENTORY)
unique(TI.XY$INVENTORY)
TI.XY$SOURCE = 'INVENTORY'


## Filter the XY tree inventory data to just the cleaned species
TI.XY  = TI.XY[TI.XY$searchTaxon %in% unique(TI.LIST$searchTaxon), ]
length(unique(TI.XY$searchTaxon))
TI.XY = na.omit(TI.XY)


## The list of species with checked maxent maps
SPP.BIAS      = merge(SPP.BIAS, COMBO.NICHE.CONTEXT[c("searchTaxon", "Origin")])
SPP.BIAS      = subset(SPP.BIAS, AUS_BOUND_BIAS == "TRUE" & Origin == "Native")
SPP.BIAS      = subset(SPP.BIAS)$searchTaxon



#########################################################################################################################
## European garden nursery lists
EURO.SPP            = read.csv("./data/base/HIA_LIST/URBAN/Euro_garden_flora_spp.csv",            stringsAsFactors = FALSE)
EURO.NURSE          = read.csv("./data/base/HIA_LIST/URBAN/Euro_garden_flora_spp_nurse.csv",      stringsAsFactors = FALSE)
EURO.NURSE.LOC      = read.csv("./data/base/HIA_LIST/URBAN/Euro_garden_flora_nurse_location.csv", stringsAsFactors = FALSE)


#########################################################################################################################
## Experimental trait lists ::
#HEAT.RISK  = read.csv("./data/base/HIA_LIST/RENEE/MOD3_HEAT_RANKS_072018.csv",                   stringsAsFactors = FALSE)
#TRAIT.SPP  = read.csv("./data/base/HIA_LIST/RENEE/RankingTraits_Control_latest.csv",              stringsAsFactors = FALSE)
TRAIT.SPP  = read.csv("./data/base/HIA_LIST/RENEE/RankingTraits_ALL.csv",                          stringsAsFactors = FALSE)
MOD.3.SPP  = read.csv("./data/base/HIA_LIST/RENEE/MOD3_OCT18.csv",                            stringsAsFactors = FALSE)
MOD.3.CHK  = join(MOD.3.SPP, MXT.CHECK[c("searchTaxon", "Check.map", "Origin", "genus", "order", "group",
                                         "Plant.type",  "Plantings", "COMBO.count", "AUS_RECORDS",
                                         "Total.growers",  "Number.of.States")], type = "full") 
MOD.3.CHK  = MOD.3.CHK [with(MOD.3.CHK , order(Check.map)), ]
#write.csv(MOD.3.CHK, "./data/base/HIA_LIST/RENEE/MOD3_OCT18_MAP.csv", row.names = FALSE)
## Now find the match between the trait species and the trait species... 
#colnames(HEAT.RISK)[colnames(HEAT.RISK)=="Species"] <- "searchTaxon"
colnames(TRAIT.SPP)[colnames(TRAIT.SPP)=="Species"] <- "searchTaxon"
#TRAIT.SPP = merge(TRAIT.SPP, HEAT.RISK, by = "searchTaxon")
trait.spp = unique(TRAIT.SPP$searchTaxon)
trait_spp = gsub(" ", "_", trait.spp)


#########################################################################################################################
## Experimental lists :: the next round of experiments will be... 
top.200              = read.csv("./data/base/HIA_LIST/HIA/HIA_TOP_200_1309_2017.csv",       stringsAsFactors = FALSE)
renee.full           = read.csv("./data/base/HIA_LIST/HIA/RENEE_FULL_LIST.csv",             stringsAsFactors = FALSE) 
renee.taxa           = read.csv("./data/base/HIA_LIST/HIA/RENEE_TAXA.csv",                  stringsAsFactors = FALSE)
renee.50             = read.csv("./data/base/HIA_LIST/HIA/RENEE_TOP_50.csv",                stringsAsFactors = FALSE)
MQ.glasshouse        = read.csv("./data/base/HIA_LIST/HIA/MQ_glasshouse.csv",               stringsAsFactors = FALSE)
Manuel.experimental  = read.csv("./data/base/HIA_LIST/HIA/Manuel_experimental_species.csv", stringsAsFactors = FALSE)
Manuel.group         = read.csv("./data/base/HIA_LIST/MANUEL/SUA_by_SPP.csv",               stringsAsFactors = FALSE)
TREE.NETWORK         = read.csv("./data/base/HIA_LIST/MANUEL/COMBO_APNI.csv",               stringsAsFactors = FALSE)
NURSE.MATCH          = read.csv("./data/base/HIA_LIST/MANUEL/nurseries.csv",                stringsAsFactors = FALSE)
campbelltown         = read.csv("./data/base/HIA_LIST/HIA/campbelltown_species.csv",        stringsAsFactors = FALSE)


## The species for the SUA analysis
SUA.spp = read.csv("./output/maxent/MAXENT_SUA_SPP.csv", stringsAsFactors = FALSE)
SUA.spp = SUA.spp[order(SUA.spp$searchTaxon),] 
SUA.spp = unique(SUA.spp$searchTaxon)
SUA_spp = gsub(" ", "_", SUA.spp)


##
intersect(campbelltown$Species, CLEAN.NICHE.CONTEXT$searchTaxon)
intersect(NURSE.MATCH$species, CLEAN.NICHE.CONTEXT$searchTaxon)
new.spp  = trimws(unique((sort(c(NURSE.MATCH$species, campbelltown$Species)))))
camp.spp = trimws(campbelltown$Species)
camp_spp = gsub(" ", "_", camp.spp) 




#########################################################################################################################
## 4). CLEAN THE HIA LIST
#########################################################################################################################


#########################################################################################################################
## Create a list of the raw HIA list, but removing the weird characters...
## Just use "TRIM(CELL)" in excel
## trim <- function (x) gsub("^\\s+|\\s+$", "", x) ##  HIA.list$Species <- trim(HIA.list$Species)
RAW.HIA.SPP = gsub("  ",     " ", HIA.list$Species)
RAW.HIA.SPP = gsub(" $",     "",  HIA.list$Species, perl = TRUE)
RAW.HIA.SPP = gsub("    $",  "",  HIA.list$Species, perl = TRUE)
length(RAW.HIA.SPP)


## also, add the "Top 200" species in here
spp.200          = top.200[c("Species")]
spp.200$Species  <- sub('(^\\S+ \\S+).*', '\\1', spp.200$Species) # \\s = white space; \\S = not white space

spp.200$Species  = gsub("  ",     " ", spp.200$Species)
spp.200$Species  = gsub(" $",     "",  spp.200$Species, perl = TRUE)
spp.200$Species  = gsub("    $",  "",  spp.200$Species, perl = TRUE)
spp.200$Top_200  = "TRUE"
spp.200          = dplyr::rename(spp.200, Binomial = Species)


## Just get the species renee selected that are not on the top 1000 or 200
renee.list       = renee.taxa[c("Species", "Growth_Form")]


#########################################################################################################################
## Merge the ~1000 with the top 200
## This merge won't get the ones that match to binomial
HIA.list$Binomial <- sub('(^\\S+ \\S+).*', '\\1',   HIA.list$Species) # \\s = white space; \\S = not white space
CLEAN.list$Binomial <- sub('(^\\S+ \\S+).*', '\\1', CLEAN.list$Species) # \\s = white space; \\S = not white space


#########################################################################################################################
## Now sum the number of growers across multiple varieties
## Just get the count for each unique binomial
n <- tapply(HIA.list$Number.of.growers, HIA.list$Binomial, sum, na.rm = TRUE)
HIA.list$Number.of.growers.total <- n[HIA.list$Binomial]
TOT.GROW = HIA.list[c("Binomial",
                      "Plant.type",
                      "Origin",
                      "Number.of.growers.total",
                      "Number.of.States")]
names(TOT.GROW) = c("searchTaxon", "Plant.type", "Origin", "Total.growers", "Number.of.States")
TOT.GROW        = TOT.GROW[!duplicated(TOT.GROW[,c('searchTaxon')]),] 


#########################################################################################################################
## Now sum the number of growers across multiple varieties
## Just get the count for each unique binomial
n.spp <- tapply(CLEAN.list$Number.of.growers, CLEAN.list$Binomial, sum, na.rm = TRUE)
CLEAN.list$Number.of.growers.total <- n.spp[CLEAN.list$Binomial]
CLEAN.GROW = CLEAN.list[c("Binomial",
                         "Number.of.growers.total")]
names(CLEAN.GROW ) = c("searchTaxon", "Total.growers")
CLEAN.GROW         = CLEAN.GROW [!duplicated(CLEAN.GROW [,c('searchTaxon')]),] 


## Join the total growers to the NICHE data
COMBO.NICHE.CONTEXT = join(COMBO.NICHE.CONTEXT, TOT.GROW)
COMBO.NICHE.CONTEXT =  COMBO.NICHE.CONTEXT[, c(1:14, 199, 15:198)] 
names(COMBO.NICHE.CONTEXT[1:15])


CLEAN.NICHE.CONTEXT = join(CLEAN.NICHE.CONTEXT, TOT.GROW)
names(CLEAN.NICHE.CONTEXT[1:15])


## Check this reduces the number of Top 200 missing from this list
HIA.list = merge(HIA.list, spp.200, by = "Binomial", all.x = TRUE) 
HIA.list$Top_200[is.na(HIA.list$Top_200)] <- "FALSE"
HIA.list$Origin <- gsub(" ",  "", HIA.list$Origin)


## check
str(HIA.list)
head(HIA.list)
unique(HIA.list$Top_200)
unique(HIA.list$Origin)
length(unique(HIA.list$Binomial)) ## 654 unique binomials
names(HIA.list)





#########################################################################################################################
## 5). CREATE LIST OF BINOMIALS : might not need this now
#########################################################################################################################


## First, taxa with "spp". Return to these later...
DRAFT.HIA.TAXA         = HIA.list
DRAFT.HIA.TAXA         = DRAFT.HIA.TAXA[DRAFT.HIA.TAXA$Species %in% DRAFT.HIA.TAXA$Species[!grepl("spp.", DRAFT.HIA.TAXA$Species)], ]
dim(DRAFT.HIA.TAXA)    ## 948 species after we cut out the "spp."


## Remove weird characters...
DRAFT.HIA.TAXA$Species = gsub(" x",     "",  DRAFT.HIA.TAXA$Species)
DRAFT.HIA.TAXA$Species = gsub("NA",     "",  DRAFT.HIA.TAXA$Species)
DRAFT.HIA.TAXA$Species = gsub("  ",     " ", DRAFT.HIA.TAXA$Species)
DRAFT.HIA.TAXA$Species = gsub(" $",     "",  DRAFT.HIA.TAXA$Species, perl = TRUE)
DRAFT.HIA.TAXA$Species = gsub("    $",  "",  DRAFT.HIA.TAXA$Species, perl = TRUE)


#########################################################################################################################
## Now create a table of how many varieties each species has
# length(unique(HIA.list$Binomial)) 
# length(unique(sub('(^\\S+ \\S+).*', '\\1', DRAFT.HIA.TAXA$Species)))


## Create another binomial column
DRAFT.HIA.TAXA$Binomial <- sub('(^\\S+ \\S+).*', '\\1', DRAFT.HIA.TAXA$Species) # \\s = white space; \\S = not white space


## And count how many varieties each taxa has?
## Before, this needs to change so that every variety that matches a binomial (e.g. magnolia grandiflora) is added to the
## Number of varieties. Also, can we take the variety with the highest number of growers? There are 8 different varieties
## of magnolia, currently I'm not getting the most popular variety.


## Some spp/varieties are not working..............................................................................................
HIA.VARIETY <- 
  DRAFT.HIA.TAXA$Binomial[DRAFT.HIA.TAXA$Binomial != DRAFT.HIA.TAXA$Species] %>% 
  table %>% 
  as.data.frame %>% 
  setNames(c('Binomial', 'No.of.Varieties')) %>% 
  full_join(DRAFT.HIA.TAXA) %>% 
  dplyr::select(Species, Binomial, No.of.Varieties, Plant.type:WA, Number.of.growers.total, Number.of.States, Origin, Top_200)

HIA.VARIETY %>% 
  filter(Binomial==Species)


#######################################################################################################################
## Which taxa have the second word capitalised, and are not captured by eliminating the "spp"?
grep('^\\S+ [A-Z]', HIA.VARIETY$Species, val = TRUE)


## Now for GBIF, just get the unique species...
HIA.SPP = HIA.VARIETY[!duplicated(HIA.VARIETY["Binomial"]),]
HIA.SPP = dplyr::rename(HIA.SPP, HIA.Taxa = Species)


## Reorder by species
HIA.SPP = HIA.SPP[with(HIA.SPP, order(Binomial)), ] 
#View(HIA.SPP)


#######################################################################################################################
## Now create list of HIA species. 610 unique species, mius the corrections, etc. 
spp            = unique(as.character(HIA.SPP$Binomial))
length(spp)   


########################################################################################################################
## Try using taxonlookup to check the taxonomy
HIA.SPP.LOOKUP = lookup_table(HIA.SPP[["Binomial"]], by_species = TRUE) ## convert rows to column and merge
HIA.SPP.LOOKUP = setDT(HIA.SPP.LOOKUP , keep.rownames = TRUE)[]
HIA.SPP.LOOKUP = dplyr::rename(HIA.SPP.LOOKUP, Binomial = rn)
head(HIA.SPP.LOOKUP) ## Can merge on the bilogical data here...
## write.csv(HIA.SPP.LOOKUP, "./data/base/HIA_LIST/HIA/HIA_BINOMIAL_LOOKUP.csv",   row.names = FALSE)
## write.csv(CLEAN.SPP,      "./data/base/HIA_LIST/HIA/CLEAN_BINOMIAL_LOOKUP.csv", row.names = FALSE)
## write.csv(CLEAN.SPP,      "./data/base/HIA_LIST/HIA/CLEAN_BINOMIAL_LOOKUP.csv", row.names = FALSE)

length(intersect(TI.LIST$searchTaxon, CLEAN.SPP$Binomial))
ALL.INV.EV = intersect(TI.LIST$searchTaxon, CLEAN.SPP$Binomial)
# write.csv(ALL.INV.EV, "./data/base/HIA_LIST/HIA/CLEAN_BINOMIAL_LOOKUP.csv", row.names = FALSE)





#########################################################################################################################
## 6). MANUSCRIPT LIST
#########################################################################################################################


#########################################################################################################################
## Now just get the species with good models and enough data
MODEL.CHECK = merge(CLEAN.NICHE.CONTEXT, MAXENT.CHECK)
MODEL.CHECK = MODEL.CHECK[, c("searchTaxon",
                              "COMBO.count",
                              "AUS_RECORDS",
                              "Total.growers",
                              "Origin",
                              "Plant.type",
                              "Top_200",                                                                        
                              "check.map")]


## Now how many trees do we have? 68 modelled, another 86 with ok data.
## This will increase with Alessandro's data. Re-create niche file, and re-run the criteria.
## Let's get another 40 exotic trees with good models
checked.trees = subset(MODEL.CHECK, Total.growers >= 25 & check.map <= 2 & Plant.type == "Tree")$searchTaxon ## 67
subset(MODEL.CHECK, Total.growers >= 25 & check.map <= 2 & Plant.type == "Tree" & Origin == "Exotic")$searchTaxon


#########################################################################################################################
## What is the intersection betwen the evergreen list and the Tree inventories?
TREE.HIA.SPP = intersect(subset(TI.LIST, Plantings > 50)$searchTaxon, CLEAN.SPP$Binomial) ##
TREE.200.SPP = intersect(subset(TI.LIST, Plantings > 500)$searchTaxon, CLEAN.SPP$Binomial) 
##intersect(head(TI.LIST, 600)$searchTaxon, CLEAN.SPP$Binomial)
summary(TI.LIST[TI.LIST$searchTaxon %in% TREE.HIA.SPP, ])


## Create a list of synonyms using the TPL search
synonyms     = c("Fraxinus angustifolia", "Ficus microcarpa", "Callistemon viminalis", "Schinus molle", 
                 "Alnus acuminata",       "Schinus areira",   "Ficus microcarpa",      "Alnus jorullensis", 
                 "Quercus coccinea",      "Libidibia ferrea", "Fraxinus oxycarpa",     "Prunus blireiana",
                 "Cupressus leylandii",   "Elaeocarpus angustifolius",                 "Eucalyptus globulus",
                 "Corymbia ficifolia",    "Eucalyptus bicolor",                        "Ulmus minor")


## Get the intersection of the synonyms
TREE.HIA.SPP = sort(unique(c(TREE.HIA.SPP, synonyms)))
TREE_HIA_SPP = gsub(" ", "_", TREE.HIA.SPP)
TREE_HIA_SPP = gsub("_trees_sign_.shp", "", TREE.HIA.SPP)


## Replace the weird shapefile name for hobart
TI.XY$INVENTORY = gsub("_trees_sign_.shp", "", TI.XY$INVENTORY)
head(TI.XY)
unique(TI.XY$INVENTORY)


## Only process biased species from the target list, and remove those which are not really biased...........
SPP.BIAS    = intersect(SPP.BIAS, TREE.HIA.SPP)
SPP.BIAS    = setdiff(SPP.BIAS, c("Banksia integrifolia", "Brachychiton acerifolius", "Callistemon viminalis",
                                  "Corymbia maculata", "Elaeocarpus grandis", "Eucalyptus camaldulensis",
                                  "Flindersia schottiana", "Leptospermum petersonii", "Lophostemon suaveolens",
                                  "Stenocarpus sinuatus", "Syzygium luehmannii", "Acacia pendula"))
SPP_BIAS      = gsub(" ", "_", SPP.BIAS)


sessionInfo()## Reverse the model list
TREE.HIA.REV = sort(TREE.HIA.SPP, decreasing = TRUE) 
TREE.HIA.REV = sort(TREE.HIA.SPP, decreasing = TRUE) 


## Create a table of the species modelled - join this on in step 8
EVERGREEN = CLEAN.SPP[c("Binomial", "Number.of.growers", "Number.of.States", "Origin")]
names(EVERGREEN )[names(EVERGREEN) == 'Binomial'] <- 'searchTaxon'
TREE.EVERGREEN = merge(EVERGREEN, TI.LIST)
TREE.EVERGREEN = TREE.EVERGREEN[TREE.EVERGREEN$searchTaxon %in% TREE.HIA.SPP, ]
TREE.EVERGREEN = TREE.EVERGREEN [with(TREE.EVERGREEN, rev(order(Plantings))), ]
TREE.EVERGREEN = TREE.EVERGREEN[c("searchTaxon", "Origin", "Plantings", "Number.of.States", "Number.of.growers")]
TREE.EVERGREEN = TREE.EVERGREEN[!duplicated(TREE.EVERGREEN[,c('searchTaxon')]),]


## What does the dataset look like?
TREE.EVERGREEN$Origin = trimws(TREE.EVERGREEN$Origin)
dim(TREE.EVERGREEN)
head(TREE.EVERGREEN)
round(with(TREE.EVERGREEN, table(Origin)/sum(table(Origin))*100), 1)


## How many trees have already been modelled well
length(intersect(TREE.EVERGREEN$searchTaxon, checked.trees))
write.csv(TREE.EVERGREEN, "./data/base/HIA_LIST/COMBO/TREE_EVERGREEN.csv", row.names = FALSE)



## What is the difference between the SUA list and the greater list?
length(setdiff(HIA.SPP$Binomial, TREE.EVERGREEN$searchTaxon))
length(setdiff(CLEAN.SPP$Binomial, TREE.EVERGREEN$searchTaxon))
length(intersect(CLEAN.SPP$Binomial, HIA.SPP$Binomial))




#########################################################################################################################
## 7). CHECK TAXONOMY AND ORIGIN OF THE MANUSCRIPT SPECIES
#########################################################################################################################


## The problems is the mis-match between what we searched, and what GBIF returned. 

## 1). Create the inital list by combing the planted trees with evergreen list
##     TREE.HIA.SPP = intersect(subset(TI.LIST, Plantings > 50)$searchTaxon, CLEAN.SPP$Binomial) 
##     This is about 400 species
##     Origin
##     NA     Exotic Native 
##     20.5   25.8   53.8 

## 2). Clean this list using the GBIF backbone taxonomy :: use the "species" column in from the GBIF "species lookup" tool
##     https://www.gbif.org/tools/species-lookup

## 3). Run the GBIF "species" list through the TPL taxonomy. Take "New" Species and Genus as the "searchTaxon"



# ## Run taxonstand for the manuscript tree species
# TREE.TAXO <- TPL(unique(TREE.EVERGREEN$searchTaxon), infra = TRUE,
#                  corr = TRUE, repeats = 100)  ## to stop it timing out...
# 
# 
# ## Most of the species are resolved 
# sort(names(TREE.TAXO))
# table(TREE.TAXO$New.Taxonomic.status)


## Unique(TREE.EVERGREEN$New.Taxonomic.status)
#saveRDS(TREE.EVERGREEN, file = paste("./data/base/HIA_LIST/COMBO/TREE_LIST.rds"))
#SUA.SPP             = read.csv("./data/base/HIA_LIST/COMBO/SUA_SPP_TPL.csv",                     stringsAsFactors = FALSE)
SUA.SPP             = read.csv("./data/base/HIA_LIST/COMBO/GBIF_SPP.csv",                     stringsAsFactors = FALSE)
SUA.SPP             = sort(unique(SUA.SPP$GBIF_SPP))
SUA_SPP             = gsub(" ", "_", SUA.SPP)
setdiff(TREE.HIA.SPP, SUA.SPP)


HIA.GBIF             = read.csv("./data/base/HIA_LIST/COMBO/HIA_SPP_GBIF_MATCH.csv",      stringsAsFactors = FALSE)
CLEAN.GBIF           = read.csv("./data/base/HIA_LIST/COMBO/CLEAN_SPP_GBIF_MATCH.csv",    stringsAsFactors = FALSE)
HIA.GBIF             = sort(unique(HIA.GBIF$HIA_GBIF_SPP_MATCH))
CLEAN.GBIF           = sort(unique(CLEAN.GBIF$CLEAN_GBIF_SPP_MATCH))
length(HIA.GBIF);length(CLEAN.GBIF)


#########################################################################################################################
## only run TPL checks once, they are slow
# TPL.SUA   <- TPL(unique(SUA.SPP), infra = TRUE, corr = TRUE, repeats = 100)
# TPL.HIA   <- TPL(unique(HIA.GBIF), infra = TRUE, corr = TRUE, repeats = 100)
# TPL.CLEAN <- TPL(unique(CLEAN.GBIF), infra = TRUE, corr = TRUE, repeats = 100)
# TPL.INV   <- TPL(unique(ALL.INV.EV), infra = TRUE, corr = TRUE, repeats = 100)
# View(TPL.SUA[c("Taxon", "Taxonomic.status", "New.Taxonomic.status", "New.Genus", "New.Species")])
# saveRDS(TPL.SUA, file = paste("./data/base/HIA_LIST/GBIF/TPL_SUA.rds"))
# saveRDS(TPL.HIA, file = paste("./data/base/HIA_LIST/HIA/TPL_HIA.rds"))
# saveRDS(TPL.CLEAN, file = paste("./data/base/HIA_LIST/HIA/TPL_CLEAN.rds"))
# saveRDS(TPL.INV, file = paste("./data/base/HIA_LIST/HIA/TPL_INV.rds"))
#write.csv(TPL.SUA, "./data/base/HIA_LIST/GBIF/TPL_SUA.csv", row.names = FALSE)


#########################################################################################################################
##  Create species list from the TPL taxonomy
## Note that TPL is not always right for Australian species :: Eucalyptus largiflorens is the accepted Australian name
## bicolor is not recognised by the ALA. Here we are 
TPL.SUA   = readRDS("./data/base/HIA_LIST/GBIF/TPL_SUA.rds")
TPL.HIA   = readRDS("./data/base/HIA_LIST/HIA/TPL_HIA.rds")
TPL.CLEAN = readRDS("./data/base/HIA_LIST/HIA/TPL_CLEAN.rds")
TPL.INV   = readRDS("./data/base/HIA_LIST/HIA/TPL_CLEAN.rds")

TPL.SPP   = paste(TPL.SUA$New.Genus, TPL.SUA$New.Species, sep = " ")  ## using the latest binomial names from TPL
TPL_SPP   = gsub(" ", "_", TPL.SPP)
setdiff(TREE.HIA.SPP, SUA.SPP)

TPL.HIA   = paste(TPL.HIA$New.Genus, TPL.HIA$New.Species, sep = " ")
TPL_HIA   = gsub(" ", "_", TPL.HIA)

TPL.CLEAN = paste(TPL.CLEAN$New.Genus, TPL.CLEAN$New.Species, sep = " ")
TPL_HIA   = gsub(" ", "_", TPL.CLEAN)


length(TPL.SUA);length(TPL.SPP);length(TPL.CLEAN)


#########################################################################################################################
## 8). CREATE LIST OF INVENTORY SPECIES
#########################################################################################################################


## Losers
losers  = c("Cryptocarya laevigata", "Syzygium wilsonii", "Dysoxylum fraserianum")
winners = c("Alectryon oleifolius", "Flindersia australis", "Flindersia maculata", "Grevillea baileyana")
gard.aus = c("Lophostemon confertus", "Eucalyptus leucoxylon", "Cupaniopsis anacardioides", 
             "Tristaniopsis laurina", "Jacaranda mimosifolia", "Pongamia pinnata", "Melia azedarach")


#########################################################################################################################
## Create a count of the horticultural species
TI.LUT = as.data.frame(table(TI.XY$searchTaxon))
names(TI.LUT) = c("searchTaxon", "FREQUENCY")
TI.LUT = TI.LUT[with(TI.LUT, rev(order(FREQUENCY))), ] 
head(TI.LUT);dim(TI.LUT)

length(intersect(CLEAN.SPP$Binomial, TI.LUT $searchTaxon))
length(intersect(TI.LUT$searchTaxon, unique(APNI$searchTaxon)))
length(intersect(CLEAN.SPP$Binomial, unique(APNI$searchTaxon)))



## Intersection between Manuel's list and all the HIA species
WPW.spp = unique(c(TPL.HIA, TPL.CLEAN, ALL.INV.EV))
WPW.MANUEL = setdiff(Manuel$Species, WPW.spp)


#########################################################################################################################
## Join on the native data and the APNI
CLEAN.ORIGIN    = dplyr::rename(CLEAN.SPP, 
                                searchTaxon = Binomial)[c("searchTaxon", "Origin", "Plant.type", "Number.of.States")]


## Pipe the dataset into here
TI.COUNT.NAT  <- TI.LUT %>% 
  
  join(., APNI) %>% 
  
  join(., CLEAN.ORIGIN) %>% 
  
  join(., TOT.GROW) %>% 
  
  join(., MAXENT.CHECK[c("searchTaxon", "check.map")]) 
  
head(TI.COUNT.NAT)
#write.csv(TI.COUNT.NAT, './data/base/HIA_LIST/COMBO/TREE_INVENTORY_COUNT.csv', row.names = FALSE)
#length(intersect(MAXENT.CHECK$))
length(intersect(MAXENT.CHECK$searchTaxon, TI.LUT$searchTaxon))





#########################################################################################################################
## OUTSTANDING LIST TASKS:
#########################################################################################################################


## Remove gunk from this file............................................................................................

## Increase the taxonomic check to include all species on HIA list


#########################################################################################################################
############################################  END OF HIA LIST CODE ###################################################### 
#########################################################################################################################