#########################################################################################################################
#############################################  FIT MAXENT TO GBIF DATA ################################################## 
#########################################################################################################################


## This code takes a table of species occurrences (rows) and environmental values "columns", and runs a maxent model 
## for each species. 


## The brief again:

## The first module will focus on fifty plant species identified in the project’s Target Species List, and will develop maps 
## that demonstrate each species’ suitability to both current and future climates across Australia.
 
## These maps will be used to demonstrate how well or poorly a particular species will be able to tolerate future conditions 
## in urban centres across Australia as the climate changes, based on our current understanding of species’ climatic 
## requirements.

## There are several points to consider here:

## Which environmnetal variables are used
## Which GCMs and RCPs
## Can the maxent settings be the same for all species?


#########################################################################################################################
## Load packages
p <- c('ff',    'things',         'raster',    'dismo',        'sp',           'latticeExtra', 'data.table', 
       'rgdal', 'rgeos',          'gdalUtils', 'rmaxent',      'readr',        'dplyr',        'tidyr',
       'readr', 'rnaturalearth',  'rasterVis', 'RColorBrewer', 'latticeExtra', 'parallel')


## Require packages
sapply(p, require, character.only = TRUE)
source('./R/GREEN_CITIES_FUNCTIONS.R')
source('./R/SDM_FUNCTIONS.R')


#########################################################################################################################
## 1). SELECT WORLDCLIM VARIABLES 
#########################################################################################################################


#########################################################################################################################
## WHICH WORLDCLIM VARIABLES TO USE?


## Ignoring edaphic variables for now

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
## Extremes (e.g. warmest/coldest month)


## These will change... 
sdm.predictors <- c("Annual_mean_temp", "Temp_seasonality",    "Max_temp_warm_month", "Min_temp_cold_month",
                    "Annual_precip",    "Precip_seasonality",  "Precip_wet_qu",       "Precip_dry_qu",       
                    "Precip_wet_month", "Precip_dry_month")



#########################################################################################################################
## Read in Raster data
load("./data/base/HIA_LIST/COMBO/COMBO_RASTER_CONTEXT.RData")
dim(COMBO.RASTER.CONTEXT)    ## 19 million rows, the latest dataset 
names(COMBO.RASTER.CONTEXT)


#########################################################################################################################
## Now quantify the correlations between these variables are related within this set seems wise.
## Take a subset of the data for one well recorded species.
sp.n = "Fagus sylvatica"
sp.subset <- subset(COMBO.RASTER.CONTEXT, searchTaxon == sp.n) %>%
  
  ## just get the sdm.predictors
  select(one_of(sdm.predictors)) %>%
  as.data.frame()

dim(sp.subset)
head(sp.subset)


## Try the correlation for everything
combo.subset <- COMBO.RASTER.CONTEXT %>%
  
  ## just get the sdm.predictors
  select(one_of(sdm.predictors)) %>%
  as.data.frame()

dim(combo.subset)
head(combo.subset)


#########################################################################################################################
## Create a pearson correlation matrix between all a-priori analysis variables
correlations <- cor(combo.subset, use = "pairwise.complete.obs") 

  
## Turn into a upper triangle
upperTriangle <- upper.tri(correlations, diag = F)


## Take a copy of the original cor-mat
correlations.upperTriangle <- correlations


## Set everything not in upper triangle o NA
correlations.upperTriangle[!upperTriangle]<-NA                   


## Use melt to reshape the matrix into triplets, na.omit to get rid of the NA rows
correlations.table <- na.omit(melt(correlations.upperTriangle, value.name = "correlationCoef")) 


## Rename columns
colnames(correlations.table) <- c("Var1", "Var2", "Correlation")  


## Reorder by absolute correlation
correlations.table = correlations.table[order(-abs(correlations.table["Correlation"])),]  


#########################################################################################################################
# FAGUS.COR <- sp.subset
# FAGUS.COR %>%
#   
#   ## First get the pearson correlations...
#   cor(sp.subset) %>%
#   
#   ## Not all these steps are necessary
#   FAGUS.COR[lower.tri(FAGUS.COR, diag = TRUE)] <- NA %>%
#   
#   ## Create data table
#   as.data.frame(as.table(FAGUS.COR)) %>%
# 
#   ## Remove NAs
#   na.omit(FAGUS.COR) %>%
#   
#   ## Order the table by the absolute correlation
#   #FAGUS.COR[order(-abs(FAGUS.COR$Freq)),] %>%
#   FAGUS.COR[order(-abs(FAGUS.COR["Freq"])),] %>%
#   
#   ## Rename the columns
#   dplyr::rename(FAGUS.COR, 
#                 Variable      = Var1,
#                 Variable_2    = Var2,
#                 Pearson_R2    = Freq)  %>%
#   
#   ## finally spit out a table
#   as.data.frame()


#########################################################################################################################
## Have a look at the matrix, and also the ordered list of variable combinations:
print(kable(correlations.upperTriangle))
print(kable(correlations.table, row.names = FALSE))


#########################################################################################################################
## Try a 'chart correlation', showing the histograms:
chart.Correlation(sp.subset, 
                  histogram = TRUE, pch = 19, main = "Worldclim variables")


#########################################################################################################################
## Do a pairs plot of all the variables
## scatterplots for all variables
pairs(Fagus.vars,
      lower.panel = panel.cor,
      col = "blue", pch = 19, cex = 0.7, main = "Worldclim variables")


## Save correlations to file
save(correlations.upperTriangle, file = paste("./output/tables/variable_selection/Worldclim_select_cormatrix.RData", sep = ""))
save(correlations.table,         file = paste("./output/tables/variable_selection/Worldclim_select_cortable.RData",  sep = ""))




#########################################################################################################################
## 2). PREPARE DATA TABLE FOR SDM ANALYSIS
#########################################################################################################################


## Chose the minimum possible predictors
sdm.predictors <- c("Annual_mean_temp", "Temp_seasonality",    "Max_temp_warm_month", "Min_temp_cold_month",
                    "Annual_precip",    "Precip_seasonality",  "Precip_wet_month", "Precip_dry_month")


## Create an empty raster with the desired properties, using raster(raster(x))
template.raster <- raster(raster("//sci-7910/F/data/worldclim/world/0.5/bio/current/bio_01")) %>% 
  projectRaster(res = 1000, crs = CRS('+init=ESRI:54009'))


## Just get the columns needed for modelling...This changes based on the ecological framework
COMBO.RASTER.SELECT <- select(COMBO.RASTER.CONTEXT, 
                              searchTaxon, 
                              lon, lat, 
                              Annual_mean_temp, Temp_seasonality,   Max_temp_warm_month, Min_temp_cold_month,
                              Annual_precip,    Precip_seasonality, Precip_wet_month,    Precip_dry_month)


## Create a spatial points object, and change to a projected system to calculate distance more accurately 
coordinates(COMBO.RASTER.SELECT) <- ~lon+lat
proj4string(COMBO.RASTER.SELECT) <- '+init=epsg:4326'
COMBO.RASTER.SELECT <- spTransform(
  COMBO.RASTER.SELECT, CRS('+init=ESRI:54009'))


## Now split using the data using the species column, and get the unique occurrence cells
COMBO.RASTER.SPLIT <- split(COMBO.RASTER.SELECT, COMBO.RASTER.SELECT$searchTaxon)
occurrence_cells   <- lapply(COMBO.RASTER.SPLIT, function(x) cellFromXY(template.raster, x))
str(occurrence_cells)  ## this is a list of dataframes, where the number of rows for each being the species table


## Now get just one record within each 10*10km cell. This step should eliminate most, if not all, duplicate records
## A simple alternative to the extract problem could be to run this process before the extract?
SDM.DATA <- mapply(function(x, cells) {
  x[!duplicated(cells), ]
}, COMBO.RASTER.SPLIT, occurrence_cells, SIMPLIFY = FALSE) %>% do.call(rbind, .)


## Use the analysis data
str(SDM.DATA)





#########################################################################################################################
## 2). RUN MODELS FOR MULTIPLE SPECIES USING TABLE OF RECORDS (ROWS) AND ENVIRONMENT (COLUMNS)
#########################################################################################################################


#########################################################################################################################
## Create species subsets for analysis
spp.all  <- unique(COMBO.RASTER.CONTEXT$searchTaxon)
str(spp.all)                 ## 6782


## Just the top 200 species
top.200.spp     = subset(COMBO.RASTER.CONTEXT, Top_200 == TRUE)["searchTaxon"] 
top.200.spp     = as.character(unique(top.200.spp$searchTaxon))
str(top.200.spp)

## Just the species with > 25 growers
spp.25     = subset(COMBO.RASTER.CONTEXT, Number.of.growers >= 25)["searchTaxon"]
spp.25     = as.character(unique(spp.25$searchTaxon))
str(spp.25)


########################################################################################################################
## We can run Maxent from a cluster of cores on the local computer. Here we send (i.e. export) all the necessary ingredients 
## to the cluster. So that's the:


## template.raster raster, 
## data frame and 
## maxent function

## Use these functions to debug
## debugonce(FIT_MAXENT)
## undebug(FIT_MAXENT)

## Use this function to find all functions relating to a search term/topic (e.g. 'debug')
## apropos('debug')
## apropos('read')


## 100 species takes about 4 hours...
cl <- makeCluster(6)
clusterExport(cl, c('template.raster', 'SDM.DATA', 'FIT_MAXENT'))
clusterEvalQ(cl, {
  
  require(ff)
  require(rgeos)
  require(sp)
  require(raster)
  require(rJava)
  require(dismo)
  require(things)
  
})


## Now use 'lapply' to run maxent for multiple species
## Note that running the code in parallel causes problems
lapply(spp.25[1:length(spp.25)], function(x) { # for serial, parLapply(cl, species[1:8], function(x) { # for parallel 
  
  ## Print the taxa being processed to screen
  message('Doing ', x)
  
  ## Subset the records to only the taxa being processed
  occurrence <- subset(SDM.DATA, searchTaxon == x)
  
  ## Now get the background points. These can come from anywhere in the whole dataset,
  ## other than the species used.
  background <- subset(SDM.DATA, searchTaxon != x)
  
  ## The create a vector of the sdm.predictors used. 
  ## This should be based on an ecological framework! 
  sdm.predictors <- sdm.predictors # vector of used sdm.predictors 

  ## Finally fit the models using FIT_MAXENT
  ## There is no switch in the function to skip outputs that exist.
  ## Given all the changes likely to be made to the models, this could be wise...
  FIT_MAXENT(occ                     = occurrence, 
             bg                      = background, 
             sdm.predictors          = sdm.predictors, 
             name                    = x, 
             outdir                  = 'output/maxent/baseline', 
             template.raster,
             min_n                   = 20,   ## This should be higher...
             max_bg_size             = 100000,
             background_buffer_width = 200000,
             shapefiles              = TRUE,
             features                = 'lpq',
             replicates              = 5,
             responsecurves          = TRUE)
  
})



## Potential errors...

## this is because you don't have maxent correctly installed yet 
# Error in .local(x, p, ...) : args not understood:
#   replicates = 5, responsecurves = TRUE, threshold = FALSE, hinge = FALSE

## Ignore this error, doesn't mean anything
# In addition: Warning messages:
#   1: In writeOGR(swd_occ, outdir_sp, "occ_swd", "ESRI Shapefile", overwrite_layer = TRUE) :
#   Field names abbreviated for ESRI Shapefile driver
# 2: In writeOGR(swd_bg, outdir_sp, "bg_swd", "ESRI Shapefile", overwrite_layer = TRUE) :
#   Field names abbreviated for ESRI Shapefile driver

stopCluster(cl)


## Look at the output...





#########################################################################################################################
## 3). PROJECT MODELS
#########################################################################################################################


#########################################################################################################################
## Here we needs to choose RCPs and emission senarios. This data should be stored locally, so we can access from



# models         <- list.files('F:/output', '^maxent_fitted\\.rds$', recursive = TRUE, full = TRUE)
# models_by_type <- split(models, sub('plants_|chordata_|nonchordata_', '', 
#                                     basename(dirname(dirname(models)))))
# 
# # Project to current climate, all of Australia, for each of three predictor 
# # sets: clim+soil, clim+soil+weathering, clim+soil+weathering+topography
# lapply(names(models_by_type), function(x) {
#   
#   clim <- sprintf(
#     
#     'f:/data/narclim_from_c_drive/albers/1km/BIOCLIM/epoch_1990_2009/aus/bioclim_1990-2009_aus_%s_albers.tif', 
#     c('p02', 'p04', 'p05', 'p06', 'p13', 'p14', 'p15'))
#   
#   soil <- c('c:/data/csiro/soil_spectra/1km/aus/PC1.tif',
#             'c:/data/csiro/soil_spectra/1km/aus/PC2.tif',
#             'c:/data/csiro/soil_spectra/1km/aus/PC3.tif')
#   
#   wii <- 'c:/data/weathering/wii_1km_albers.tif'
#   tpi <- 'f:/data/narclim/CSIRO_topo_AA/CSIRO_TPI_eMast_albers.tif'
#   twi <- 'f:/data/narclim/CSIRO_topo_AA/TWI_eMAST_albers.tif'
#   
#   preds <- switch(
#     
#     x,
#     clim_soil = c(clim, soil),
#     clim_soil_topo_wii = c(clim, soil, wii, tpi, twi),
#     clim_soil_wii = c(clim, soil, wii))
#   
#   s <- stack(preds)
#   
#   scen <- 'current'
#   names(s) <- sub('.*(p\\d{2})_?.*', '\\1', names(s))
#   names(s) <- gsub("_eMast_albers|_1km_albers|_eMAST_albers|CSIRO_","",names(s))
#   
#   # Identify which cells have data for all sdm.predictors.
#   # locs contains cell numbers and corresponding coordinates of non-NA cells.
#   locs  <- which(!is.na(sum(s)[]))
#   locs  <- cbind(cell = locs, xyFromCell(s, locs))
#   cells <- locs[, 'cell']
#   r     <- raster(s)
#   
#   # Using ff_matrix objects (this will make calculating the mean and sd a lot easier).
#   s_ff <- ff(vmode = "double", dim = c(length(cells), nlayers(s)),
#              filename = ff_swd <- tempfile(fileext = '.ff'))
#   # fill ff_matrix with data
#   for(i in 1:nlayers(s)) {
#     
#     s_ff[, i] <- s[[i]][][cells]
#     
#   }
#   
#   colnames(s_ff) <- names(s)
#   rm(s); gc()
#   lapply(models_by_type[[x]], function(model) {
#     
#     species <- basename(dirname(model))
#     outfile <- sprintf('%s/%s_%s_1000m_prediction.ff', 
#                        dirname(model),
#                        gsub(' +', '_', species), scen)
#     
#     if(!file.exists(extension(outfile, '.tif'))) {
#       
#       r_pred <- r
#       m <- readRDS(model)$me_full
#       message('Doing species ', species)
#       preds_ff <- ff(vmode = "double", dim = c(length(cells), 1),
#                      filename = outfile)
#       finalizer(preds_ff) <- 'close'
#       
#       preds_ff[, 1] <- 
#         round(rmaxent::project(m, s_ff[, seq_len(ncol(s_ff))])$prediction_logistic * 1000)
#       
#       r_pred[cells] <- preds_ff[, 1]
#       writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag = -9999)
#       #saveRDS(preds_ff, extension(outfile, 'rds')) 
#       close(preds_ff)
#       delete(preds_ff)
#       rm(preds_ff, r_pred)
#       
#     }
#     
#   })
#   
#   delete(s_ff)
#   rm(s_ff)
#   gc()
#   
# })





#########################################################################################################################
## OUTSTANDING SDM TASKS:
#########################################################################################################################


## Create a switch to skip files that exist

## Further mapping and cleaning of GBIF data needed for the important species

## Can maxent setting be the same for all species?  

## Which GCMs and RCPs? Need layers for 2030, 2050, etc.

## Create maps using prediction functions?

## 



#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################