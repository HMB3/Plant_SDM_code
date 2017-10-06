#########################################################################################################################
#############################################  FIT MAXENT TO GBIF DATA ################################################## 
#########################################################################################################################


#########################################################################################################################
## Load packages
p <- c('ff', 'things', 'raster', 'dismo', 'sp', 'latticeExtra', 
       'rgdal', 'rgeos', 'gdalUtils', 'rmaxent', 'readr', 'dplyr', 'tidyr',
       'readr', 'rnaturalearth', 'rasterVis', 'RColorBrewer', 'latticeExtra')
sapply(p, require, character.only = TRUE)





#########################################################################################################################
## 1). READ IN GBIF DATA
#########################################################################################################################

## Group: 'plants', 'chordata', or 'nonchordata'
group <- 'chordata'

load("./data/base/HIA_LIST/COMBO/COMBO_RASTER_CONTEXT.RData") 

# occ <- readRDS(c(plants='', 
#                   chordata='data/occurrence/chordata/shp_cleaned_by_linda_and_bionet_etc_added/occ_minimal_albers_chordata_phase1.rds',
#                   nonchordata='')[group])


# overwrite CRS to EPSG:3577 so its expression is consistent with bg data
#proj4string(occ) <- '+init=epsg:3577'
occ_by_sp <- split(occ, occ$searchTaxon)





#########################################################################################################################
## 2). READ IN RASTER DATA
#########################################################################################################################


#Load template raster to determine cell numbers for points
r1000 <- raster('//SCI-7910/f/data/narclim/albers/1km/BIOCLIM/epoch_1990-2009/aus/bioclim_1990-2009_aus_p02_albers.tif')

# gdalUtils::gdalwarp('C:/data/weathering/wii_1km/wii_1km/wii_oz2_w1k3',
#                     'C:/data/weathering/wii_1km_albers.tif',
#                     te=c(bbox(r1000)),
#                     tr=c(1000, 1000),
#                     r='bilinear',
#                     t_srs='epsg:3577')
wii <- raster('C:/data/weathering/wii_1km_albers.tif')

# gdalUtils::gdalwarp('F:/data/narclim/CSIRO_topo_AA',
#                     'F:/data/narclim/CSIRO_topo_AA/CSIRO_TPI_eMast_albers.tif',
#                     te=c(bbox(r1000)),
#                     tr=c(1000, 1000),
#                     r='bilinear',
#                     t_srs='epsg:3577')
TPI <- raster('//sci-7910/f/data/narclim/CSIRO_topo_AA/CSIRO_TPI_eMast_albers.tif')

# gdalUtils::gdalwarp('F:/data/narclim/CSIRO_topo_AA',
#                     'F:/data/narclim/CSIRO_topo_AA/TWI_eMAST_albers.tif',
#                     te=c(bbox(r1000)),
#                     tr=c(1000, 1000),
#                     r='bilinear',
#                     t_srs='epsg:3577')
TWI <- raster('//sci-7910/f/data/narclim/CSIRO_topo_AA/TWI_eMAST_albers.tif')


#########################################################################################################################
## Create vector of predictor paths
predictors_current <- c(
  '//sci-7910/f/data/narclim_from_c_drive/albers/1km/BIOCLIM/epoch_1990_2009/aus/bioclim_1990-2009_aus_p02_albers.tif',
  '//sci-7910/f/data/narclim_from_c_drive/albers/1km/BIOCLIM/epoch_1990_2009/aus/bioclim_1990-2009_aus_p04_albers.tif', 
  '//sci-7910/f/data/narclim_from_c_drive/albers/1km/BIOCLIM/epoch_1990_2009/aus/bioclim_1990-2009_aus_p05_albers.tif', 
  '//sci-7910/f/data/narclim_from_c_drive/albers/1km/BIOCLIM/epoch_1990_2009/aus/bioclim_1990-2009_aus_p06_albers.tif', 
  '//sci-7910/f/data/narclim_from_c_drive/albers/1km/BIOCLIM/epoch_1990_2009/aus/bioclim_1990-2009_aus_p13_albers.tif', 
  '//sci-7910/f/data/narclim_from_c_drive/albers/1km/BIOCLIM/epoch_1990_2009/aus/bioclim_1990-2009_aus_p14_albers.tif', 
  '//sci-7910/f/data/narclim_from_c_drive/albers/1km/BIOCLIM/epoch_1990_2009/aus/bioclim_1990-2009_aus_p15_albers.tif',
  'c:/data/csiro/soil_spectra/1km/aus/PC1.tif',
  'c:/data/csiro/soil_spectra/1km/aus/PC2.tif',
  'c:/data/csiro/soil_spectra/1km/aus/PC3.tif',
  'c:/data/weathering/wii_1km_albers.tif', 
  '//sci-7910/f/data/narclim/CSIRO_topo_AA/CSIRO_TPI_eMast_albers.tif', 
  '//sci-7910/f/data/narclim/CSIRO_topo_AA/TWI_eMAST_albers.tif')

#predictor_set <- c(rep('clim', 7), rep('soil', 3), 'wii', 'topo', 'topo')


# Create ff_matrix objects for current, all-Australia grids
# It's pretty quick for 1000 m environmental data, but don't rerun this if
# already run. Instead, load the data with
# readRDS('C:/Users/MQ20156218/Documents/R_projects/sdms_1km_2016/data/current_all_Australia_ff_1000m.rds')
# Have commented out the next few # lines so we don't accidentally run them
# again.

# s <- stack(predictors_current)
# names(s) <- sub('.*_(p\\d{2})_.*', '\\1', names(s))
# names(s)[names(s)=='wii_1km_albers'] <- 'wii'
# names(s)[names(s)=='CSIRO_TPI_eMast_albers'] <- 'TPI'
# names(s)[names(s)=='TWI_eMAST_albers'] <- 'TWI'
# s_ff <- ff(vmode="double", dim=c(ncell(s), nlayers(s)),
#            filename='data/current_all_Aus_1000m.ff')
#  
# for(i in 1:nlayers(s)) {
#  s_ff[,i] <- s[[i]][]
# }
# colnames(s_ff) <- names(s)
# saveRDS(s_ff, file='data/current_all_Australia_ff_1000m.rds')
# rm(s); gc()
# s_ff <- readRDS('data/current_all_Australia_ff_1000m.rds')





#########################################################################################################################
#source('fit_maxent_background_buffer.R') # load in the function to fit Maxent models
# coerce bg to a SPDF so we can write it out as a shapefile
bg <- SpatialPointsDataFrame(bg, data.frame(id = seq_len(length(bg))))

lapply(c('clim_soil', 'clim_soil_wii', 'clim_soil_topo_wii'), function(x) {
  
  outdir <- file.path('//sci-7910/f/sos_output/phase1', paste(group, x, sep='_'))
  
  if(!dir.exists(outdir)) dir.create(outdir)
  
  preds <- switch(x,
         clim_soil     = c(sprintf('p%02d', c(2, 4, 5, 6, 13, 14, 15)), paste0('PC', 1:3)),
         clim_soil_wii = c(sprintf('p%02d', c(2, 4, 5, 6, 13, 14, 15)), paste0('PC', 1:3), 'wii'),
         clim_soil_wii = c(sprintf('p%02d', c(2, 4, 5, 6, 13, 14, 15)), paste0('PC', 1:3), 'wii', 'TPI', 'TWI'))
  
  mapply(function(occ, name) {
    
    message('Doing ', name)
    
    name <- gsub('\\*', '', name) # remove asterisks from names
    
    fit_maxent2(occ   =occ, bg=bg, predictors=s_ff[, preds], name=name, 
                outdir=outdir, 
                template=r1000, shapefiles=TRUE, features='lpq',
                replicates=5)
    
  }, occ_by_sp, names(occ_by_sp)) 
  
})





#########################################################################################################################
## 3). PROJECT MODELS
#########################################################################################################################

library(ff)
library(raster)
library(data.table)
library(dismo)
library(things)

models <- list.files('F:/output', '^maxent_fitted\\.rds$', recursive=TRUE, full=TRUE)

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
    clim_soil=c(clim, soil),
    clim_soil_topo_wii=c(clim, soil, wii, tpi, twi),
    clim_soil_wii=c(clim, soil, wii))
  
  s <- stack(preds)
  
  scen <- 'current'
  names(s) <- sub('.*(p\\d{2})_?.*', '\\1', names(s))
  names(s) <- gsub("_eMast_albers|_1km_albers|_eMAST_albers|CSIRO_","",names(s))
  
  # Identify which cells have data for all predictors.
  # locs contains cell numbers and corresponding coordinates of non-NA cells.
  locs <- which(!is.na(sum(s)[]))
  locs <- cbind(cell=locs, xyFromCell(s, locs))
  cells <- locs[, 'cell']
  r <- raster(s)
  
  # Using ff_matrix objects (this will make calculating the mean and sd a lot easier).
  s_ff <- ff(vmode="double", dim=c(length(cells), nlayers(s)),
             filename=ff_swd <- tempfile(fileext='.ff'))
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
      preds_ff <- ff(vmode="double", dim=c(length(cells), 1),
                     filename=outfile)
      finalizer(preds_ff) <- 'close'
      
      preds_ff[, 1] <- 
        round(rmaxent::project(m, s_ff[, seq_len(ncol(s_ff))])$prediction_logistic * 1000)
      
      r_pred[cells] <- preds_ff[, 1]
      writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999)
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
## 1). FUTURE PREDICTIONS
#########################################################################################################################


#########################################################################################################################
## Create vector of predictor data directories
dirs <- grep('\\d{4}$',
             list.dirs('f:/data/narclim_from_c_drive/albers/1km/BIOCLIM/Interpolated complete'),
             value = TRUE)

## Loop over climate scenarios, holding each swd in memory while looping over SDMs
lapply(dirs, function(d) {
  
  require(ff)
  
  clim <- list.files(d, '\\.tif$', full = T)
  soil <- c('c:/data/csiro/soil_spectra/1km/narclim/PC1.tif',
            'c:/data/csiro/soil_spectra/1km/narclim/PC2.tif',
            'c:/data/csiro/soil_spectra/1km/narclim/PC3.tif')
  wii  <- 'c:/data/weathering/wii_1km_albers_narclim.tif'
  tpi  <- 'f:/data/narclim/CSIRO_topo_AA/CSIRO_TPI_eMast_albers_narclim.tif'
  twi  <- 'f:/data/narclim/CSIRO_topo_AA/TWI_eMAST_albers_narclim.tif'
  
  preds <- switch(
    x,
    clim_soil = c(clim, soil),
    clim_soil_topo_wii = c(clim, soil, wii, tpi, twi),
    clim_soil_wii = c(clim, soil, wii))
  
  s <- stack(preds)
  gcm <- basename(dirname(dirname(d)))
  rcm <- basename(dirname(d))
  yr <- basename(d)
  
  scen <- paste(gcm, rcm, yr, sep = '_')
  names(s) <- sub('.*(p\\d{2})_?.*', '\\1', names(s))
  names(s) <- gsub("_eMast_albers|_1km_albers|_eMAST_albers|CSIRO_|_narclim","",names(s))
  
  #########################################################################################################################
  ## create empty ff_matrix
  
  s_ff <- ff(vmode = "double", dim = c(length(cells), nlayers(s)),
             filename = ff_swd <- tempfile(fileext = '.ff'))
  
  ## fill ff_matrix with data
  for(i in 1:nlayers(s)) {
    
    s_ff[, i] <- s[[i]][][cells]
    
  }
  
  colnames(s_ff) <- names(s)
  rm(s); gc()
  lapply(models, function(model) {
    
    species <- basename(dirname(model))
    outfile <- sprintf('%s/%s_%s_1000m_prediction.ff', 
                       dirname(model),
                       gsub(' +', '_', species), scen)
    
    if(!file.exists(extension(outfile, '.tif'))) {
      
      r_pred <- r1000
      m <- readRDS(model)$me_full
      message('Doing species ', species)
      preds_ff <- ff(vmode = "double", dim = c(length(cells), 1),
                     filename = outfile)
      
      finalizer(preds_ff) <- 'close'
      
      preds_ff[, 1] <- 
        
        round(rmaxent::project(m, s_ff[, seq_len(ncol(s_ff))])$prediction_logistic * 1000)
      
      r_pred[cells] <- preds_ff[, 1]
      writeRaster(r_pred, extension(outfile, '.tif'), datatype = 'INT2S', NAflag =-9999)
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
## 
ff <- list.files('f:/output', '\\.tif$', full = TRUE, recursive = TRUE)

#occ <- file.path(dirname(ff), 'occ.rds')
#spp_list <- read_csv('data/refugia2.csv')
#aus <- ne_countries(country='Australia')
#aus_albers <- spTransform(ne_countries(country='Australia'), CRS('+init=epsg:3577'))


#########################################################################################################################
## Plot 
sapply(ff, function(f) {
  
  nm <- gsub('_', ' ', basename(dirname(f)))
  
  f_png <- sprintf('%s/%s_current.png', dirname(f), basename(dirname(f)))
  
  if (!file.exists(f_png)) {
    
    message('Making map for species: ', nm)  
    png(f_png, 7, 4, units = 'in', res = 300)
    par(mfrow = c(1, 2), mar = c(0, 0, 0, 0), oma = c(1, 1, 3, 1))  
    xy <- list(NULL, readRDS(file.path(dirname(f), 'occ.rds')))
    
    r <- raster(f)/1000
    
    p <- levelplot(stack(r, r), col.regions = colorRampPalette(brewer.pal(11, 'Spectral')), 
                   at = seq(0, 1, len = 100), margin = FALSE, scales = list(draw = FALSE), 
                   colorkey = list(height = 0.6), main = list(nm, font = 4), names.attr = c('', '')) + 
      
      layer(sp.points(xy[[panel.number()]], pch = 20, cex = 0.7, col = 1), data = list(xy = xy))
    
    print(p)
    dev.off() 
    
  }
  
})



#########################################################################################################################
## 
models_by_variable_set <- split(models, dirname(dirname(models)))

results <- setNames(lapply(models_by_variable_set, function(mm) {
  
  ff <- sub('maxent_fitted.rds', 'xval/maxentResults.csv', mm)
  
  results <- sapply(ff, function(f) {
    
    d <- read.csv(f)
    d <- d[nrow(d), ]
    list(test_auc = d[, 'Test.AUC'],
         contribution = unlist(d[, grep('contribution', names(d))]),
         permutation = unlist(d[, grep('permutation', names(d))]))
    
  })
  
  nm <- basename(dirname(dirname(colnames(results))))
  
  list(test_auc = setNames(do.call(c, results[1, ]), nm), 
       contribution = setNames(do.call(cbind.data.frame, results[2, ]), nm),
       permutation = setNames(do.call(cbind.data.frame, results[3, ]), nm))
  
}), basename(names(models_by_variable_set)))



#########################################################################################################################
## Getting number of records per species from each set of predictors

models_by_variable_set <- split(models, dirname(dirname(models)))

results2 <- setNames(lapply(models_by_variable_set, function(mm) {
  
  ff <- sub('maxent_fitted.rds', 'xval/maxentResults.csv', mm)
  
  results <- sapply(ff, function(f) {
    
    d <- read.csv(f)
    d <- d[nrow(d), ]
    d2 <- read.csv(sub('xval', 'full', f))
    
    list(test_auc=d[, 'Test.AUC'],
         contribution=unlist(d[, grep('contribution', names(d))]),
         permutation=unlist(d[, grep('permutation', names(d))]),
         #n_records_xval=unlist(d[, grep('X.Training.samples', names(d))]),
         n_records_full=unlist(d2[, grep('X.Training.samples', names(d2))]))
    
  })
  
  nm <- basename(dirname(dirname(colnames(results))))
  
  list(test_auc=setNames(do.call(c, results[1, ]), nm), 
       contribution=setNames(do.call(cbind.data.frame, results[2, ]), nm),
       permutation=setNames(do.call(cbind.data.frame, results[3, ]), nm),
       #n_training_records_mean_xval=setNames(do.call(c, results['n_records_xval', ]), nm),
       n_training_records_full=setNames(do.call(c, results['n_records_full', ]), nm))
  
}), basename(names(models_by_variable_set)))



n <- lapply(split(results2, sub('_clim.*', '', names(results2))), function(x) {
  out <- do.call(cbind.data.frame, lapply(x, '[[', 'n_training_records_full'))
  names(out) <- sub('^[^_]+_', '', names(out))
  data.frame(species=gsub('_', ' ', row.names(out)), out)
})

mapply(function(x, nm) 
  write.csv(x, sprintf('n_training_records_fullmodel_%s.csv', nm), row.names=F), 
  n, names(n))
  
