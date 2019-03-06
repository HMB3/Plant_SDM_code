#########################################################################################################################
######################################## READ SPATIAL DATA FOR HIA MODELLING ############################################ 
#########################################################################################################################


## This code reads in the vector and raster data that is used to create species distribution models for Module 1 of the 
## WPW project. This data is also used to create results for the Stoten publication ::


#########################################################################################################################
## SETUP
#########################################################################################################################


#########################################################################################################################
## Load only the packages needed for the analysis
## lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)
p <- c('ff',    'things',    'raster',        'dismo',        'sp',           'latticeExtra', 'data.table',
       'rgdal', 'rgeos',     'gdalUtils',     'rmaxent',      'readr',        'plyr',         'dplyr',
       'tidyr', 'readr',     'rnaturalearth', 'rasterVis',    'RColorBrewer', 'latticeExtra', 'parallel',
       'ALA4R', 'stringr',   'Taxonstand',    'CoordinateCleaner', 'gsubfn',  'PerformanceAnalytics',
       'rvest', 'magrittr',  'devtools',      'ggplot2',      'reshape2',     'rmarkdown', 'flexdashboard', 'shiny', 'rgbif',
       'ENMeval', 'tibble',  'ncdf4',         'Cairo', 'taxonlookup', 'kgc', 'maptools', 'DataCombine', 'mgcv', 'rsq', 'utf8',
       'betareg', 'hydroTSM', 'bomrang', 'gridExtra', 'grid', 'lattice', 'ConR')


## Require packages
sapply(p, require, character.only = TRUE)
source_gist('26e8091f082f2b3dd279',             filename = 'polygonizer.R')
source_gist('c6a1cb61b8b6616143538950e6ec34aa', filename = 'hatch.R')
devtools::source_gist('306e4b7e69c87b1826db',   filename = 'diverge0.R')


## Source functions, and set temporary directory (for both raster files and just generally)
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
AUS           = readRDS("./data/base/CONTEXTUAL/aus_states.rds")
ALL.SUA.POP   = read.csv("./data/base/CONTEXTUAL/ABS_SUA_POP.csv", stringsAsFactors = FALSE)
URB.POP       = read.csv("./data/base/CONTEXTUAL/ABS_URBAN_CENTRE_POP.csv", stringsAsFactors = FALSE)


#########################################################################################################################
## Check how the Koppen zones were calculated
Koppen_zones     = unique(readOGR('data/base/CONTEXTUAL/WC05_1975H_Koppen_Shapefile/WC05_1975H_Koppen_Kriticos_2012.shp')@data[, 1:2])
Koppen_1975_1km  = raster('data/Koppen_1000m_Mollweide54009.tif')
# Koppen_1975_2km  = raster('data/Koppen_2000m_Mollweide54009.tif')
# Koppen_1975_3km  = raster('data/Koppen_3000m_Mollweide54009.tif')
# Koppen_1975_5km  = raster('data/Koppen_5000m_Mollweide54009.tif')
# Koppen_1975_10km = raster('data/Koppen_10000m_Mollweide54009.tif')




#########################################################################################################################
## Read in the background data :: this was generated by running the set of 3.8k HIA species through the whole code

## Re-create the background points using Alessandro's latest data........................................................
## Also, need to include as many inventory points as possible
background.points = readRDS("./data/ANALYSIS/background_points.rds") %>%
  spTransform(sp_epsg54009)
length(unique(background.points $searchTaxon))
dim(background.points )
unique(background.points $SOURCE)
background.points.df = as.data.frame(background.points)


## There aint much Inventory data compared to the other sources!
## this is ok, except for some Australian species which have 
## Loads of inventory data, so it becomes hard to sample proportional to the source of the occurrence dataset
round(with(background.points.df, table(SOURCE)/sum(table(SOURCE))), 3)





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


## Names of all 19 worldclim predictors that could be used to run SDMs
sdm.predictors <- c("Annual_mean_temp",    "Mean_diurnal_range",  "Isothermality",      "Temp_seasonality",  
                    "Max_temp_warm_month", "Min_temp_cold_month", "Temp_annual_range",  
                    "Mean_temp_wet_qu",    "Mean_temp_dry_qu",    
                    "Mean_temp_warm_qu",   "Mean_temp_cold_qu",   
                    
                    "Annual_precip",       "Precip_wet_month",   "Precip_dry_month",    "Precip_seasonality",  
                    "Precip_wet_qu",       "Precip_dry_qu",      
                    "Precip_warm_qu",      "Precip_col_qu")


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


## Names of A-priori worldclim predictors :: these have been chosen to capture means, variability and extremes
sdm.select     <- c("Annual_mean_temp", "Temp_seasonality",    "Max_temp_warm_month", "Min_temp_cold_month",
                    "Annual_precip",    "Precip_seasonality",  
                    "Precip_wet_month", "Precip_dry_month")


## Names of all the worldclim predictors :: uneeded
grid.names = c('Annual_mean_temp',    'Mean_diurnal_range',  'Isothermality',      'Temp_seasonality', 
               'Max_temp_warm_month', 'Min_temp_cold_month', 'Temp_annual_range',  'Mean_temp_wet_qu',
               'Mean_temp_dry_qu',    'Mean_temp_warm_qu',   'Mean_temp_cold_qu',  'Annual_precip',
               'Precip_wet_month',    'Precip_dry_month',    'Precip_seasonality', 'Precip_wet_qu',
               'Precip_dry_qu',       'Precip_warm_qu',      'Precip_col_qu')


#########################################################################################################################
## Create a raster stack of current global environmental conditions
world.grids.current.1km = stack(
  file.path('./data/base/worldclim/world/0.5/bio/current',
            sprintf('bio_%02d', 1:19)))


# ## 5km
# world.grids.current.5km = stack(
#   file.path('./data/base/worldclim/world/0.5/bio/current/5km',
#             sprintf('bio_%02d_10km.tif', 1:19)))
# 
# ## 10km
# world.grids.current.10km = stack(
#   file.path('./data/base/worldclim/world/0.5/bio/current/10km',
#             sprintf('bio_%02d_10km.tif', 1:19)))


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
## Also get the PET raster
PET               = raster("./data/base/worldclim/world/1km/pet_he_yr1.tif")





#########################################################################################################################
## 4). PREPARE TEMPLATE RASTER GRIDS FOR SDM ANALYSIS
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

## 5km
# template.raster.3km <- gdalwarp("data/base/worldclim/world/0.5/bio/current/bio_01",
#                                 tempfile(fileext = '.tif'),
#                                 t_srs = '+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs',
#                                 output_Raster = TRUE,
#                                 tr = c(3000, 3000),
#                                 r = "near", dstnodata = '-9999')


# ## 5km
# template.raster.5km <- gdalwarp("data/base/worldclim/world/0.5/bio/current/bio_01",
#                                 tempfile(fileext = '.tif'),
#                                 t_srs = '+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs',
#                                 output_Raster = TRUE,
#                                 tr = c(5000, 5000),
#                                 r = "near", dstnodata = '-9999')
# xres(template.raster.5km)
# 
# 
# ## 10km
# template.raster.10km <- gdalwarp("data/base/worldclim/world/0.5/bio/current/bio_01",
#                                  tempfile(fileext = '.tif'),
#                                  t_srs = '+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs',
#                                  output_Raster = TRUE,
#                                  tr = c(10000, 10000),
#                                  r = "near", dstnodata = '-9999')
# xres(template.raster.10km)


# ## Assign all the cells to be NA
# template.raster.1km  <- !is.na(template.raster.1km)
# template.raster.2km  <- !is.na(template.raster.2km)
# template.raster.3km  <- !is.na(template.raster.3km)
# template.raster.5km  <- !is.na(template.raster.5km)
# template.raster.10km <- !is.na(template.raster.10km)


# writeRaster(template.raster.1km,  'data/template_has_data_1km.tif',   datatype = 'INT2S', overwrite = TRUE)
# writeRaster(template.raster.2km,  'data/template_has_data_2km.tif',   datatype = 'INT2S', overwrite = TRUE)
# writeRaster(template.raster.3km,  'data/template_has_data_3km.tif',   datatype = 'INT2S', overwrite = TRUE)
# writeRaster(template.raster.5km,  'data/template_has_data_5km.tif',   datatype = 'INT2S', overwrite = TRUE)
# writeRaster(template.raster.10km, 'data/template_has_data_10km.tif',  datatype = 'INT2S', overwrite = TRUE)


# template.cells.1km  <- Which(template.raster.1km  == 1, cells = TRUE, na.rm=TRUE)    ## Is this is incomplete, might be causing the error?
# template.cells.5km  <- Which(template.raster.5km  == 1, cells = TRUE, na.rm=TRUE)    ## ..................................................
# template.cells.10km <- Which(template.raster.10km == 1, cells = TRUE, na.rm=TRUE)


## Save the cells to file
# saveRDS(template.cells.1km,  'data/has_data_cells_1km.rds')
# saveRDS(template.cells.5km,  'data/has_data_cells_5km.rds')
# saveRDS(template.cells.10km, 'data/has_data_cells_10km.rds')


#########################################################################################################################
## Load template rasters
template.raster.1km  = raster("./data/template_hasData.tif")
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
## OUTSTANDING SPATIAL TASKS:
#########################################################################################################################


## 1). Re-create the background points using Alessandro's latest data....................................................

## 2). Clean up the variable names so that they don't repeat. EG one name for Aus, one for world, etc.



## Suggestions: 
##  Add more rasters and vectors in future :: EG soil layers, etc.






#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################