#########################################################################################################################
######################################## READ SPATIAL DATA FOR HIA MODELLING ############################################ 
#########################################################################################################################


## This code reads in the vector and raster data that is used to create species distribution models for Module 1 of the 
## WPW project. This data is also used to create results for the Stoten publication ::


#########################################################################################################################
## SETUP
#########################################################################################################################


## save.image("UPDATE_DATA.RData") 
## load("UPDATE_DATA.RData") 


#########################################################################################################################
## Load only the packages needed for the analysis
## lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)
p <- c('ff',      'things',    'raster',        'dismo',             'sp',           'latticeExtra',          'data.table',
       'rgdal',   'rgeos',     'gdalUtils',     'rmaxent',           'readr',        'plyr',                  'dplyr',
       'tidyr',   'readr',     'rnaturalearth', 'rasterVis',         'RColorBrewer', 'latticeExtra',          'parallel',
       'ALA4R',   'stringr',   'Taxonstand',    'CoordinateCleaner', 'gsubfn',       'PerformanceAnalytics',  'utf8',
       'rvest',   'magrittr',  'devtools',      'ggplot2',           'reshape2',     'rmarkdown',             'flexdashboard', 
       'shiny',   'rgbif',     'ENMeval',       'tibble',            'ncdf4',        'Cairo',                 'taxonlookup',  
       'kgc',     'betareg',   'gridExtra',     'grid',              'lattice',      'ConR',                  'writexl',
       'ggmap',   'DataCombine', 'exactextractr', 'sf')


## Require packages
sapply(p, require, character.only = TRUE)
source_gist('26e8091f082f2b3dd279',             filename = 'polygonizer.R')
source_gist('c6a1cb61b8b6616143538950e6ec34aa', filename = 'hatch.R')
devtools::source_gist('306e4b7e69c87b1826db',   filename = 'diverge0.R')


## Source functions, and set temporary directory (for both raster files and generally)
source('./R/WPW_GENERAL_FUNCTIONS.R')
source('./R/WPW_MAXENT_FUNCTIONS.R')
source('./R/WPW_MAPPING_FUNCTIONS.R')
rasterOptions(tmpdir = './RTEMP')




#########################################################################################################################
## 1). READ IN SHAPEFILES
#########################################################################################################################


#########################################################################################################################
## Set coordinate system definitions :: best to minimise the number of projection used in this project
## Also get rid of the '+init=ESRI", which are not compatible with some systems (e.g. libraries on linux)
CRS.MOL      <- CRS('+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs')
CRS.MOL.SDM  <- CRS('+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs84=0,0,0')
CRS.WGS.84   <- CRS("+init=epsg:4326")
CRS.AUS.ALB  <- CRS("+init=EPSG:3577")
ALB.CONICAL  <- CRS('+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')
ALB.CON      <- '+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'
sp_epsg54009 <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs84=0,0,0"


#########################################################################################################################
## Read in spatial data once, rather than in each script
## The SUA shapefiles are created by the code below
LAND          = readRDS("./data/base/CONTEXTUAL/LAND_world.rds")
SUA_2016      = readRDS("./data/base/CONTEXTUAL/SUA/SUA_2016_AUST.rds")
SUA_2016      = SUA_2016[order(SUA_2016$SUA_NAME16),]
SUA_2016      = SUA_2016 %>%
  spTransform(ALB.CONICAL)


## Read in the Local government areas
LGA           = readRDS("./data/base/CONTEXTUAL/LGA.rds")
POA_2016      = readOGR('data/base/CONTEXTUAL/SUA/POA_2016_AUST.shp')
AUS           = readRDS("./data/base/CONTEXTUAL/aus_states.rds")
ALL.SUA.POP   = read.csv("./data/base/CONTEXTUAL/ABS_SUA_POP.csv", stringsAsFactors = FALSE)
URB.POP       = read.csv("./data/base/CONTEXTUAL/ABS_URBAN_CENTRE_POP.csv", stringsAsFactors = FALSE)


#########################################################################################################################
## Check how the Koppen zones were calculated
Koppen_zones     = unique(readOGR('data/base/CONTEXTUAL/WC05_1975H_Koppen_Shapefile/WC05_1975H_Koppen_Kriticos_2012.shp')@data[, 1:2])
Koppen_shp       = readOGR('data/base/CONTEXTUAL/WC05_1975H_Koppen_Shapefile/WC05_1975H_Koppen_Kriticos_2012.shp')
Koppen_1975_1km  = raster('data/world_koppen/Koppen_1000m_Mollweide54009.tif')
# Koppen_1975_2km  = raster('data/world_koppen/Koppen_2000m_Mollweide54009.tif')
# Koppen_1975_3km  = raster('data/world_koppen/Koppen_3000m_Mollweide54009.tif')
# Koppen_1975_5km  = raster('data/world_koppen/Koppen_5000m_Mollweide54009.tif')
# Koppen_1975_10km = raster('data/world_koppen/Koppen_10000m_Mollweide54009.tif')





#########################################################################################################################
## READ IN SDM POINTS
#########################################################################################################################


#########################################################################################################################
## Read in the background data :: this was generated by running the set of 3.7k tree inventory species through the whole code
## 


## Re-create the background points using Alessandro's latest data........................................................
## change BG points to be the SDM_SPAT_OCC_BG_ALL data set.
# hollow.points     = readRDS("./data/ANALYSIS/SDM_SPAT_OCC_BG_HOLLOW_SPP.rds")%>%
#   spTransform(sp_epsg54009)
# TI.points = readRDS("./data/ANALYSIS/SDM_SPAT_OCC_BG_TREE_INVENTORY_SPP.rds") %>%
#   spTransform(sp_epsg54009)


## Bind the rows of the hollow and inventory SDM tables together 
# background.points = rbind(hollow.points, TI.points) 


## 
# length(unique(background.points$searchTaxon))
# dim(background.points )
# unique(background.points$SOURCE)
#background.points.df = as.data.frame(background.points)



#########################################################################################################################
## Read in the spatial points for all HIA species 
# SDM.SPAT.OCC.BG.ALL = readRDS(paste0(DATA_path, 'SDM_SPAT_OCC_BG_ALL.rds'))
# length(unique(SDM.SPAT.OCC.BG.ALL$searchTaxon))
# length(intersect(unique(SDM.SPAT.OCC.BG.ALL$searchTaxon), WPW.spp)) 





#########################################################################################################################
## 2). CREATE RASTER OF SIGNIFCANT URBAN AREAS
#########################################################################################################################


## Urban areas are used to intersect with the SDM maps.
## This code has already been run, by reading shapefiles and converting to RDS............................................


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



###################################################################################################################
## Read in the city data
# CITY.POINTS   = SpatialPointsDataFrame(coords      = TI.CITY[c("lon", "lat")],
#                                        data        = TI.CITY,
#                                        proj4string = CRS(projection(world.grids.current)))
# 
# writeOGR(obj    = CITY.POINTS,
#          dsn    = "./data/base/Global_tree_inventories",
#          layer  = paste0('INVENTORY_CITY_POINTS'),
#          driver = "ESRI Shapefile", overwrite_layer = TRUE)
# saveRDS(CITY.POINTS, file = paste("./data/base/Global_tree_inventories/INVENTORY_CITY_POINTS.rds"))

# CITY.CLIMATE = read.csv("./data/base/Global_tree_inventories/Tree_inventory_city_climate_koppen_table.csv", stringsAsFactors = FALSE)


#########################################################################################################################
## Create Koppen climate zone for centroid of each SUA
# cents <- coordinates(CITY.POINTS)
# cents <- SpatialPointsDataFrame(coords = cents, data = CITY.POINTS@data,
#                                 proj4string = CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"))
# 
# plot(CITY.POINTS)
# points(cents, col = "Blue")
# writeSpatialShape(cents, "cents")

# city_centroids <- coordinates(CITY.POINTS)
# city_location  = data.frame(city_centroids)
# city_location  = cbind(CITY.POINTS$CityCountr, city_location)
# names(city_location)   = c("CityCountr","rndCoord.lon", "rndCoord.lat")
# city_location$rndCoord.lon = RoundCoordinates(city_location$rndCoord.lon)
# city_location$rndCoord.lat = RoundCoordinates(city_location$rndCoord.lat)
# 
# Kop.city.loc <- data.frame(city_location, ClimateZ = LookupCZ(city_location))
# Kop.city.loc = Kop.city.loc[c("CityCountr", "ClimateZ")]
# 
# 
# ## Join the city worldclim data to the 
# CITY.CLIMATE = join(Kop.city.loc, CITY.CLIMATE)
# head(CITY.CLIMATE)
# write.csv(CITY.CLIMATE, "./data/base/Global_tree_inventories/Tree_inventory_city_climate_koppen_table.csv", row.names = FALSE)


#########################################################################################################################
## Read in the files created by the code above for Section 2).
SUA_2016_rast = readRDS("./data/base/CONTEXTUAL/SUA/SUA_2016_RAST.rds")
SUA_2016_vec  = readRDS("./data/base/CONTEXTUAL/SUA/SUA_2016_VEC.rds")
summary(SUA_2016_vec)





#########################################################################################################################
## 3). READ IN RASTER COMPONENTS NEEDED FOR ANALYSIS
#########################################################################################################################


#########################################################################################################################
## Create the variables needed to access current environmental conditions + their names in the functions
## Names of all the worldclim variables used to extract the raster data
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


## Names of the 15 worldclim predictors that have reliable data
bs.predictors <- c("Annual_mean_temp",    "Mean_diurnal_range",  "Isothermality",      "Temp_seasonality",  
                    "Max_temp_warm_month", "Min_temp_cold_month", "Temp_annual_range",  
                    "Mean_temp_warm_qu",   "Mean_temp_cold_qu",   
                    
                    "Annual_precip",       "Precip_wet_month",   "Precip_dry_month",    "Precip_seasonality",  
                    "Precip_wet_qu",       "Precip_dry_qu")


## Names of the AWAP variables :: used to create niches for extreme weather events
awap.variables = c("Annual_mean_temp",
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
                   "Drought_mean_rel_int_extr",
                   "HWF", "HWA", "HWM", "HWD", "HWN",
                   "HW_CUM_ALL", "HW_CUM_AV",  "HW_CUM_HOT")


## Names of all the worldclim predictors :: uneeded. This is the same as 
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


#########################################################################################################################
## Create a velox raster to improve performance
# world.grids.velox <- velox(world.grids.current)
# writeRaster(world.grids.velox,  '/data/base/worldclim/world/0.5/bio/current/world_velox_1km.tif', datatype = 'INT2S', overwrite = TRUE)

## Create a raster stack of current Australian environmental conditions
## This is used to calculate the niches for the urban tree inventories
inventory.grids.current = stack(
  file.path('./data/base/worldclim/aus/1km/bio/current/WGS/', 
            sprintf('bio_%02d.tif', 1:19))) 


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
## Read in average of all 6 scenarios Australian environmental conditions
# aus.temp.2030 = raster('./data/base/worldclim/aus/1km/bio/BIO1_2030_MN.tif')/10
# aus.temp.2050 = raster('./data/base/worldclim/aus/1km/bio/BIO1_2050_MN.tif')/10
# aus.temp.2070 = raster('./data/base/worldclim/aus/1km/bio/BIO1_2070_MN.tif')/10
# 
# aus.rain.2030 = raster('./data/base/worldclim/aus/1km/bio/BIO12_2030_MN.tif')
# aus.rain.2050 = raster('./data/base/worldclim/aus/1km/bio/BIO12_2050_MN.tif')
# aus.rain.2070 = raster('./data/base/worldclim/aus/1km/bio/BIO12_2070_MN.tif')


#########################################################################################################################
## Also get the the soil rasters
PET  = raster("./data/base/worldclim/world/1km/pet_he_yr1.tif")
twi  = raster("./data/base/ACLEP/TWI.tif")
tpi  = raster("./data/base/ACLEP/TPI.tif")
soil = stack(file.path('./data/base/ACLEP', sprintf('PC%d.tif', 1:3)))


#########################################################################################################################
## Try resampling the soil rasters :: consider which re-sampling method is best 
# z <- file.path('./data/base/ACLEP', sprintf('PC%d.tif', 1:3))
# system.time(z2 <- gdalUtils::gdalwarp(z[1], f <- tempfile(fileext = '.tif'), 
#                                       te=c(bbox(aus.grids.current)), tr=c(1000, 1000),
#                                       t_srs=proj4string(aus.grids.current),
#                                       r='near', multi=TRUE, output_Raster=TRUE))

# system.time(z3 <- gdalUtils::gdalwarp(z[1], f <- tempfile(fileext = '.tif'), 
#                                       te=c(bbox(aus.grids.current)), tr=c(1000, 1000),
#                                       t_srs=proj4string(aus.grids.current),
#                                       r='cubic', multi=TRUE, output_Raster=TRUE))

# system.time(gdalUtils::gdalwarp(z[1], file.path('./data/base/ACLEP', sprintf('PC%d_aea.tif', 1:3)), 
#                                 te=c(bbox(aus.grids.current)), tr=c(1000, 1000),
#                                 t_srs=proj4string(aus.grids.current),
#                                 r='bilinear', multi=TRUE))
# 
# cellStats(raster(file.path('./data/base/ACLEP', sprintf('PC%d.tif', 1))), 'range', na.rm=TRUE)




#########################################################################################################################
## 4). DOWNLOAD RASTER DATA FROM NCIS
#########################################################################################################################


## Download raster data from NCIS
# data_path <- 'my_data/pan'
# s <-lapply(1:12, function(m) {
#   f <- sprintf('http://dapds00.nci.org.au/thredds/fileServer/rr9/ANUClimate/ANUClimate_v1-0_pan-evaporation_monthly-mean_0-01deg_1976-2005/00000000/ANUClimate_v1-0_pan-evaporation_monthly-mean_0-01deg_1976-2005_00000000_%02d.nc', m)
#   f_out <- file.path(data_path, basename(f))
#   download.file(f, f_out, mode='wb')
#   raster(f_out)
# }) %>% stack





#########################################################################################################################
## 5). PREPARE TEMPLATE RASTER GRIDS FOR SDM ANALYSIS
#########################################################################################################################


## This code has already been run


#########################################################################################################################
## Use GDAL to create a raster which = 1 where bio_01 has data (i.e. land), and NA where there is no data
## Also note that gdalwarp is much faster, and the trs ='+init=esri:54009' argument does not work here


# ## 1km
# template.raster.1km <- gdalwarp("data/base/worldclim/world/0.5/bio/current/bio_01",
#                                 tempfile(fileext = '.tif'),
#                                 t_srs = '+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs',
#                                 output_Raster = TRUE,
#                                 tr = c(1000, 1000),
#                                 r = "near", dstnodata = '-9999')
# xres(template.raster.1km)


## 2km
# template.raster.2km <- gdalwarp("data/base/worldclim/world/0.5/bio/current/bio_01",
#                                 tempfile(fileext = '.tif'),
#                                 t_srs = '+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs',
#                                 output_Raster = TRUE,
#                                 tr = c(2000, 2000),
#                                 r = "near", dstnodata = '-9999')


# ## Assign all the cells to be NA
# template.raster.1km  <- !is.na(template.raster.1km)
# writeRaster(template.raster.1km,  'data/template_has_data_1km.tif',   datatype = 'INT2S', overwrite = TRUE)

# template.cells.1km  <- Which(template.raster.1km  == 1, cells = TRUE, na.rm=TRUE)    ## Is this is incomplete, might be causing the error?
# saveRDS(template.cells.1km,  'data/has_data_cells_1km.rds')


#########################################################################################################################
## Load template rasters
template.raster.1km    = raster("./data/world_koppen/template_hasData.tif")
template.raster.1km.84 = raster("./data/world_koppen/template_1km_WGS84.tif")
# template.raster.2km  = raster("./data/template_has_data_2km.tif")
# template.raster.3km  = raster("./data/template_has_data_3km.tif")
# template.raster.5km  = raster("./data/template_has_data_5km.tif")
# template.raster.10km = raster("./data/template_has_data_10km.tif")
# xres(template.raster.1km);xres(template.raster.5km);xres(template.raster.10km)





#########################################################################################################################
## CREATE NAMES FOR THE GLOBAL CIRCULATION MODELS
#########################################################################################################################


#########################################################################################################################
## Create the names for the GCMs by scraping the worldclim website
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


#########################################################################################################################
## Create the IDs : A work around because the 2030 data comes from the climate change in Australia website, not worldclim
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
## See the publication for why we choose this
scen_2030 = c("mc85bi30", "no85bi30", "ac85bi30", "cc85bi30", "gf85bi30", "hg85bi30")
scen_2050 = c("mc85bi50", "no85bi50", "ac85bi50", "cc85bi50", "gf85bi50", "hg85bi50")
scen_2070 = c("mc85bi70", "no85bi70", "ac85bi70", "cc85bi70", "gf85bi70", "hg85bi70")





#########################################################################################################################
## 6). COMBINE THE NICHE RUNS TOGETHER
#########################################################################################################################


# ## Create a list of all dataframes with the extension from this run
# COMBO.NICHE.list  = list.files(DATA_path, pattern = 'COMBO_NICHE_CONTEXT_EVERGREEN',  full.names = TRUE, recursive = TRUE)
# SDM.TABLE.list = list.files(DATA_path, pattern = 'SDM_SPAT_OCC_BG_', full.names = TRUE)
# INV.TABLE.list = list.files(DATA_path, pattern = 'CLEAN_INV_EVERGREEN', full.names = TRUE)
# INV.TABLE.list = list.files(DATA_path, pattern = 'CLEAN_INV_TREE_INVENTORY', full.names = TRUE)
# RAS.TABLE.list = list.files(DATA_path, pattern = 'COMBO_RASTER_CONVERT_TREE_INVENTORY', full.names = TRUE)


# #########################################################################################################################
# ## Now combine the niche tables for each species into one table
# COMBO.NICHE.ALL <- COMBO.NICHE.list[1:8] %>%
# 
#   ## pipe the list into lapply
#   lapply(function(x) {
# 
#     ## create the character string
#     f <- paste0(x)
# 
#     ## load each .csv file
#     d <- readRDS(f)
#     d
# 
#   }) %>%
# 
#   ## finally, bind all the rows together
#   rbind()
# 
# 
# ## Update this
# str(COMBO.NICHE.ALL)
# dim(COMBO.NICHE.ALL)
# 
# 
# ## Make sure the Species are unique
# COMBO.NICHE.ALL = COMBO.NICHE.ALL[!duplicated(COMBO.NICHE.ALL[,c('searchTaxon')]),]
# dim(COMBO.NICHE.ALL)
# length(unique(COMBO.NICHE.ALL$searchTaxon))
# 
# 
#########################################################################################################################
## Now combine the raster tables for each species into one table
# SDM.TABLE.ALL <- INV.TABLE.list %>%
# 
#   ## Pipe the list into lapply
#   lapply(function(x) {
# 
#     ## Create the character string
#     f <- paste0(x)
# 
#     ## Load each file
#     message('Reading file for ', x)
#     d <- readRDS(f)
#     message(nrow(d), ' Records for ', x)
#     return(d)
# 
#   }) %>%
# 
#   ## Finally, bind all the rows together
#   bind_rows()
# 
# ## This is a summary of maxent output for current conditions
# dim(SDM.TABLE.ALL)
# names(SDM.TABLE.ALL)[1:10]


##
# length(unique(SDM.TABLE.ALL$searchTaxon))
# 
# 
# # #########################################################################################################################
# # ## Save the niche and raster data
# # saveRDS(COMBO.NICHE.ALL,  paste0(DATA_path, 'COMBO_NICHE_ALL.rds'))
# saveRDS(SDM.TABLE.ALL, paste0(DATA_path,    'CLEAN_INV_ALL_EVREGREEN_3004_2018.rds'))





#########################################################################################################################
## OUTSTANDING SPATIAL TASKS:
#########################################################################################################################


## 1). Re-create the background points using Alessandro's latest data....................................................

## 2). Clean up the variable names so that they don't repeat. EG one name for Aus, one for world, etc.



## Suggestions: 
##  Add more rasters and vectors in future :: EG soil layers, etc.






#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################