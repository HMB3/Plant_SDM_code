#########################################################################################################################
#############################################  FIT MAXENT TO GBIF DATA ################################################## 
#########################################################################################################################


#########################################################################################################################
## Load packages
p <- c('ff',    'things',         'raster',    'dismo',        'sp',           'latticeExtra', 'data.table', 
       'rgdal', 'rgeos',          'gdalUtils', 'rmaxent',      'readr',        'dplyr',        'tidyr',
       'readr', 'rnaturalearth',  'rasterVis', 'RColorBrewer', 'latticeExtra', 'parallel')


## Require packages
sapply(p, require, character.only = TRUE)


#########################################################################################################################
## 1). PREPARE TABLE FOR ANALYSIS
#########################################################################################################################


#########################################################################################################################
## Create list of species from the GBIF data...
load("./data/base/HIA_LIST/COMBO/COMBO_RASTER_CONTEXT.RData")
names(COMBO.RASTER.CONTEXT)
species  <- unique(COMBO.RASTER.CONTEXT$searchTaxon)


## Create an empty raster with the desired properties, using raster(raster(x))
## Also change the coordinate system to a projected system, so that distance can be calculated... 
template.raster <- raster(raster("//sci-7910/F/data/worldclim/world/0.5/bio/current/bio_01")) %>% 
  projectRaster(res = 1000, crs = CRS('+init=ESRI:54009'))


## Just get the columns needed for modelling...
COMBO.RASTER.CONTEXT <- select(COMBO.RASTER.CONTEXT, 
                               searchTaxon, 
                               lon, lat, Mean_diurnal_range:Temp_seasonality)


## Create a spatial points object
coordinates(COMBO.RASTER.CONTEXT) <- ~lon+lat
proj4string(COMBO.RASTER.CONTEXT) <- '+init=epsg:4326'
COMBO.RASTER.CONTEXT <- spTransform(
  COMBO.RASTER.CONTEXT, CRS('+init=ESRI:54009'))


## Now split using the data using the species column, and get the unique occurrence cells
COMBO.RASTER.SPLIT <- split(COMBO.RASTER.CONTEXT, COMBO.RASTER.CONTEXT$searchTaxon)
occurrence_cells     <- lapply(COMBO.RASTER.SPLIT, function(x) cellFromXY(template.raster, x))
str(occurrence_cells)  ## this is a list of dataframes, where the number of rows for each being the species table


## Now get just one record within each 10*10km cell. This step should eliminate most, if not all, duplicate records
## A simple alternative to the extract problem could be to use this instead, before the extract?
SDM.DATA <- mapply(function(x, cells) {
  x[!duplicated(cells), ]
}, COMBO.RASTER.SPLIT, occurrence_cells, SIMPLIFY = FALSE) %>% do.call(rbind, .)


## Use the analysis data
str(SDM.DATA)





#########################################################################################################################
## 2). RUN MODELS FOR GENERIC SPECIES USING TABLE OF RECORDS (ROWS) AND ENVIRONMENT (COLUMNS)
#########################################################################################################################



#########################################################################################################################
## We can run Maxent from a cluster. Here we send all the necessary ingredients to the cluster. So that's the  
cl <- makeCluster(5)
clusterExport(cl, c('template', 'SDM.DATA', 'FIT_MAXENT'))
clusterEvalQ(cl, {
  
  require(ff)
  require(rgeos)
  require(sp)
  require(raster)
  require(rJava)
  require(dismo)
  require(things)
  
})


## Use 'lapply' to run maxent for multiple species   
lapply(species[1], function(x) { # for serial, parLapply(cl, species[1:8], function(x) { # for parallel 
  
  ## Print the taxa being processed to screen
  message('Doing ', x)
  
  ## Split the records to the taxa being processed
  occurrence <- subset(SDM.DATA, searchTaxon == x)
  
  ## The background points can come from anywhere in the whole dataset
  background <- subset(SDM.DATA, searchTaxon != x)
  
  ## Create a vector of the predictors used. This should be based on an ecological framework! 
  predictors <- c("Mean_diurnal_range", "Isothermality", "Temp_seasonality") # vector of used predictors 

  ## Finally fit the models using  
  FIT_MAXENT(occ                     = occurrence, 
             bg                      = background, 
             predictors              = predictors, 
             name                    = x, 
             outdir                  = 'output/maxent', 
             template                = template,
             min_n                   = 20,
             max_bg_size             = 100000,
             background_buffer_width = 200000,
             shapefiles              = TRUE,
             features                = 'lpq',
             replicates              = 5,
             responsecurves          = TRUE)
  
})



## Errors... 
# Error in .local(x, p, ...) : args not understood:
#   replicates = 5, responsecurves = TRUE, threshold = FALSE, hinge = FALSE

# In addition: Warning messages:
#   1: In writeOGR(swd_occ, outdir_sp, "occ_swd", "ESRI Shapefile", overwrite_layer = TRUE) :
#   Field names abbreviated for ESRI Shapefile driver
# 2: In writeOGR(swd_bg, outdir_sp, "bg_swd", "ESRI Shapefile", overwrite_layer = TRUE) :
#   Field names abbreviated for ESRI Shapefile driver

stopCluster(cl)


## Look at the output...





#########################################################################################################################
## 2). PROJECT MODELS
#########################################################################################################################


## 
models         <- list.files('F:/output', '^maxent_fitted\\.rds$', recursive = TRUE, full = TRUE)
models_by_type <- split(models, sub('plants_|chordata_|nonchordata_', '', 
                                    basename(dirname(dirname(models)))))

# Project to current climate, all of Australia, for each of three predictor 
# sets: clim+soil, clim+soil+weathering, clim+soil+weathering+topography
lapply(names(models_by_type), function(x) {
  
  clim <- sprintf(
    
    'f:/data/narclim_from_c_drive/albers/1km/BIOCLIM/epoch_1990_2009/aus/bioclim_1990-2009_aus_%s_albers.tif', 
    c('p02', 'p04', 'p05', 'p06', 'p13', 'p14', 'p15'))
  
  soil <- c('c:/data/csiro/soil_spectra/1km/aus/PC1.tif',
            'c:/data/csiro/soil_spectra/1km/aus/PC2.tif',
            'c:/data/csiro/soil_spectra/1km/aus/PC3.tif')
  
  wii <- 'c:/data/weathering/wii_1km_albers.tif'
  tpi <- 'f:/data/narclim/CSIRO_topo_AA/CSIRO_TPI_eMast_albers.tif'
  twi <- 'f:/data/narclim/CSIRO_topo_AA/TWI_eMAST_albers.tif'
  
  preds <- switch(
    
    x,
    clim_soil = c(clim, soil),
    clim_soil_topo_wii = c(clim, soil, wii, tpi, twi),
    clim_soil_wii = c(clim, soil, wii))
  
  s <- stack(preds)
  
  scen <- 'current'
  names(s) <- sub('.*(p\\d{2})_?.*', '\\1', names(s))
  names(s) <- gsub("_eMast_albers|_1km_albers|_eMAST_albers|CSIRO_","",names(s))
  
  # Identify which cells have data for all predictors.
  # locs contains cell numbers and corresponding coordinates of non-NA cells.
  locs  <- which(!is.na(sum(s)[]))
  locs  <- cbind(cell = locs, xyFromCell(s, locs))
  cells <- locs[, 'cell']
  r     <- raster(s)
  
  # Using ff_matrix objects (this will make calculating the mean and sd a lot easier).
  s_ff <- ff(vmode = "double", dim = c(length(cells), nlayers(s)),
             filename = ff_swd <- tempfile(fileext = '.ff'))
  # fill ff_matrix with data
  for(i in 1:nlayers(s)) {
    
    s_ff[, i] <- s[[i]][][cells]
    
  }
  
  colnames(s_ff) <- names(s)
  rm(s); gc()
  lapply(models_by_type[[x]], function(model) {
    
    species <- basename(dirname(model))
    outfile <- sprintf('%s/%s_%s_1000m_prediction.ff', 
                       dirname(model),
                       gsub(' +', '_', species), scen)
    
    if(!file.exists(extension(outfile, '.tif'))) {
      
      r_pred <- r
      m <- readRDS(model)$me_full
      message('Doing species ', species)
      preds_ff <- ff(vmode = "double", dim = c(length(cells), 1),
                     filename = outfile)
      finalizer(preds_ff) <- 'close'
      
      preds_ff[, 1] <- 
        round(rmaxent::project(m, s_ff[, seq_len(ncol(s_ff))])$prediction_logistic * 1000)
      
      r_pred[cells] <- preds_ff[, 1]
      writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag = -9999)
      #saveRDS(preds_ff, extension(outfile, 'rds')) 
      close(preds_ff)
      delete(preds_ff)
      rm(preds_ff, r_pred)
      
    }
    
  })
  
  delete(s_ff)
  rm(s_ff)
  gc()
  
})





#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################