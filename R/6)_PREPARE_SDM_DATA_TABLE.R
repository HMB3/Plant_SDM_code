#########################################################################################################################
#############################################  PREPARE DATA FOR MAXENT ################################################## 
#########################################################################################################################


#########################################################################################################################
## This code takes a table of all species occurrences (rows) and environmental values (columns), and prepares them for
## SDM analysis


#########################################################################################################################
## There are a few key facors we would like to vary across the model runs:

## Environmnetal variables : run with the same core set, and also with variable selection
## RECORDS                 : ALL, CULTIVATED & UNCULT
## RANGES                  : ALL, NATIVE RANGE & NON-NATIVE
## GCMs/RCPs               : the lowest prioritym but probably do them all in the end
## UHI 
## WUE


#########################################################################################################################
## Load packages
p <- c('ff',    'things',         'raster',    'dismo',        'sp',           'latticeExtra', 'data.table', 
       'rgdal', 'rgeos',          'gdalUtils', 'rmaxent',      'readr',        'dplyr',        'tidyr',
       'readr', 'rnaturalearth',  'rasterVis', 'RColorBrewer', 'latticeExtra', 'parallel')


## Require packages
sapply(p, require, character.only = TRUE)
load("./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT_1601_2018.RData")
source('./R/GREEN_CITIES_FUNCTIONS.R')
source('./R/MAPPING_FUNCTIONS.R')
source('./R/HIA_LIST_MATCHING.R')





#########################################################################################################################
## 1). SELECT WORLDCLIM VARIABLES FOR STANDARD MODELS
#########################################################################################################################


#########################################################################################################################
## WHICH WORLDCLIM VARIABLES TO USE?


## Ignoring edaphic and topographic variables for now

# BIO1  = Annual Mean Temperature                                     ## 
# BIO2  = Mean Diurnal Range (Mean of monthly (max temp - min temp))  
# BIO3  = Isothermality (BIO2/BIO7) (* 100)
# BIO4  = Temperature Seasonality (standard deviation *100)           ##
# BIO5  = Max Temperature of Warmest Month                            ## 
# BIO6  = Min Temperature of Coldest Month                            ##
# BIO7  = Temperature Annual Range (BIO5-BIO6)
# BIO8  = Mean Temperature of Wettest Quarter
# BIO9  = Mean Temperature of Driest Quarter
# BIO10 = Mean Temperature of Warmest Quarter
# BIO11 = Mean Temperature of Coldest Quarter
# BIO12 = Annual Precipitation                                        ##
# BIO13 = Precipitation of Wettest Month                              ##
# BIO14 = Precipitation of Driest Month                               ##
# BIO15 = Precipitation Seasonality (Coefficient of Variation)        ##
# BIO16 = Precipitation of Wettest Quarter
# BIO17 = Precipitation of Driest Quarter
# BIO18 = Precipitation of Warmest Quarter
# BIO19 = Precipitation of Coldest Quarter


#########################################################################################################################
## Here we want an a-priori framework for choosing the variables. Following Guisan and zimmerman (2000), we want
## to choose variables that directly inlfuence plant growth and reprodcution. In our case, that means proxies of:

## Nutrients, Soil air and water, heat sum and PAR

## Additionally, the data-driven approach can be useful (create this as markdown/html) E.G.

# Bradie, J. and B. Leung (2017). 
# "A quantitative synthesis of the importance of variables used in MaxEnt species distribution models." 
# Journal of Biogeography 44(6): 1344-1361.

# "Over all MaxEnt models published, the ability to discriminate occurrence from reference sites was high 
# (average AUC = 0.92). Much of this discriminatory ability was due to temperature and precipitation variables. 
# Further, variability (temperature) and extremes (minimum precipitation) were the most predictive. More generally, 
# the most commonly tested variables were not always the most predictive, with, for instance, ‘distance to water’ 
# infrequently tested, but found to be very important when it was. Thus, the results from this study summarize the 
# MaxEnt SDM literature, and can aid in variable selection by identifying underutilized, but potentially important 
# variables, which could be incorporated in future modelling efforts."


## So for the simple Worldclim variables, we will adopt this approach. For both temperature and precipitation, we can
## choose measures of:

## Average across the time period
## variability across the time period
## Extremes (e.g. most extreme months)


#########################################################################################################################
## Read in Raster data
load("./data/base/HIA_LIST/COMBO/COMBO_RASTER_CONTEXT_1601_2018.RData")
dim(COMBO.RASTER.CONTEXT)    ##
names(COMBO.RASTER.CONTEXT)


# #########################################################################################################################
# ## Now quantify the correlations between these variables are related within this set seems wise.
# ## Try the correlation for everything
# combo.subset <- COMBO.RASTER.CONTEXT %>%
#   
#   ## Just get the sdm.predictors
#   select(one_of(sdm.predictors.all)) %>%
#   as.data.frame()
# 
# dim(combo.subset)
# head(combo.subset)
# 
# 
# #########################################################################################################################
# ## Create a pearson correlation matrix between all a-priori analysis variables
# correlations <- cor(combo.subset, use = "pairwise.complete.obs") 
# 
#   
# ## Turn into a upper triangle
# upperTriangle <- upper.tri(correlations, diag = F)
# 
# 
# ## Take a copy of the original cor-mat
# correlations.upperTriangle <- correlations
# 
# 
# ## Set everything not in upper triangle to NA
# correlations.upperTriangle[!upperTriangle]<-NA                   
# 
# 
# ## Use melt to reshape the matrix into triplets, and na.omit to get rid of the NA rows
# correlations.table <- na.omit(melt(correlations.upperTriangle, value.name = "correlationCoef")) 
# 
# 
# ## Rename columns
# colnames(correlations.table) <- c("Var1", "Var2", "Correlation")  
# 
# 
# ## Reorder by absolute correlation
# correlations.table = correlations.table[order(-abs(correlations.table["Correlation"])),]  
# 
# 
# #########################################################################################################################
# ## Have a look at the matrix, and also the ordered list of variable combinations:
# print(kable(correlations.upperTriangle))
# print(kable(correlations.table, row.names = FALSE))
# 
# 
# #########################################################################################################################
# ## Try a 'chart correlation', showing the histograms:
# chart.Correlation(sp.subset, 
#                   histogram = TRUE, pch = 19, main = "Worldclim variables")
# 
# 
# #########################################################################################################################
# ## Do a pairs plot of all the variables
# ## scatterplots for all variables
# pairs(Fagus.vars,
#       lower.panel = panel.cor,
#       col = "blue", pch = 19, cex = 0.7, main = "Worldclim variables")





#########################################################################################################################
## 2). PREPARE DATA TABLE FOR SDM ANALYSIS
#########################################################################################################################


## Watch the CRS.......................................................................................................
## 
## Create an empty raster with the desired properties, using raster(raster(x))
template.raster <- raster(raster("./data/base/worldclim/world/0.5/bio/current/bio_01")) %>% 
  projectRaster(res = 1000, crs = CRS('+init=ESRI:54009'))


#########################################################################################################################
## Just get the columns needed for modelling: this would include cultivated/non, and inside/outside
## Get an ID column here too, and use it to join back on the other columns to the unique cells data
COMBO.RASTER.CONTEXT$OBS <- 1:nrow(COMBO.RASTER.CONTEXT)
dim(COMBO.RASTER.CONTEXT)[1];length(COMBO.RASTER.CONTEXT$OBS)


COMBO.RASTER.ALL = COMBO.RASTER.CONTEXT[,c("searchTaxon", 
                                           "OBS", 
                                           "lon",
                                           "lat",
                                           "Annual_mean_temp", 
                                           "Mean_diurnal_range", 
                                           "Isothermality", 
                                           "Temp_seasonality", 
                                           "Max_temp_warm_month", 
                                           "Min_temp_cold_month", 
                                           "Temp_annual_range", 
                                           "Mean_temp_warm_qu", 
                                           "Mean_temp_cold_qu", 
                                           "Annual_precip", 
                                           "Precip_wet_month", 
                                           "Precip_dry_month", 
                                           "Precip_seasonality", 
                                           "Precip_wet_qu", 
                                           "Precip_dry_qu")]

# COMBO.RASTER.ALL    <- select(COMBO.RASTER.CONTEXT, searchTaxon, lon, lat, 
#                               Annual_mean_temp, Mean_diurnal_range, Isothermality, Temp_seasonality, Max_temp_warm_month, Min_temp_cold_month, 
#                               Temp_annual_range, Mean_temp_warm_qu, Mean_temp_cold_qu, Annual_precip, 
#                               Precip_wet_month, Precip_dry_month, Precip_seasonality, Precip_wet_qu, Precip_dry_qu)


#########################################################################################################################
## Create a spatial points object, and change to a projected system to calculate distance more accurately 
coordinates(COMBO.RASTER.ALL)    <- ~lon+lat
proj4string(COMBO.RASTER.ALL)    <- '+init=epsg:4326'

##
COMBO.RASTER.ALL    <- spTransform(
  COMBO.RASTER.ALL, CRS('+init=ESRI:54009'))


## Now split using the data using the species column, and get the unique occurrence cells
COMBO.RASTER.SPLIT.ALL <- split(COMBO.RASTER.ALL, COMBO.RASTER.ALL$searchTaxon)
occurrence_cells_all   <- lapply(COMBO.RASTER.SPLIT.ALL, function(x) cellFromXY(template.raster, x))
str(occurrence_cells_all)   ## this is a list of dataframes, where the number of rows for each being the species table


#########################################################################################################################
## Now get just one record within each 10*10km cell. This step should eliminate most, if not all, duplicate records
## A simple alternative to the extract problem could be to run this process before the extract?
SDM.DATA.ALL <- mapply(function(x, cells) {
  x[!duplicated(cells), ]
}, COMBO.RASTER.SPLIT.ALL, occurrence_cells_all, SIMPLIFY = FALSE) %>% do.call(rbind, .)


## Check to see we have 19 variables + the species for the standard predictors, and 19 for all predictors
str(SDM.DATA.ALL)


## Now save/load .RData file for the next session
#save.image("STEP_6_SDM_ALL_DATA.RData")
#load("STEP_6_SDM.RData")


## Check experimental taxa again
'Swainsona formosa'  %in% SDM.DATA.ALL$searchTaxon
'Templetonia retusa' %in% SDM.DATA.ALL$searchTaxon 
'Dodonaea baueri'    %in% SDM.DATA.ALL$searchTaxon 
'Platanus hispanica' %in% SDM.DATA.ALL$searchTaxon 
'Kennedia beckxiana' %in% SDM.DATA.ALL$searchTaxon  


#########################################################################################################################
## Now join back on the contextual columns for data cleaning
SDM.DATA.ALL.CHECK = merge(SDM.DATA.ALL, COMBO.RASTER.CONTEXT, by = "OBS", all.y = FALSE)   ##
dim(SDM.DATA.ALL.CHECK);dim(SDM.DATA.ALL)


## Probably create a separate file for cleaning the records?
COMBO_check_records(taxa.list = test.spp[84],                   ## c(names(COMBO.RASTER.CONTEXT))
                    columns   = c("searchTaxon",
                                  "lat",
                                  "lon",
                                  "locality",
                                  "country",
                                  "taxo_agree",
                                  "CULTIVATED"),  
                    DF        = SDM.DATA.ALL.CHECK)




## Save big tables to keep memory spare
save(template.raster,  file = paste("./data/base/HIA_LIST/COMBO/SDM_TEMPLATE_RASTER.RData"))
save(SDM.DATA.ALL,     file = paste("./data/base/HIA_LIST/COMBO/HIA_SDM_DATA_ALL_VAR.RData"))


## Remove the other data
rm(COMBO.RASTER.ALL)
rm(COMBO.RASTER.SPLIT.ALL)
save.image("STEP_6_PREPARE_SDM.RData")





#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################