## This code can be executed in step 6 OR 7
## Create a list of all worldclim variables
# pred_names <- c(
#   'Annual_mean_temp',   ## To select a column with the cursor, hold ctrl+alt and use up or down arrow
#   'Mean_diurnal_range',
#   'Isothermality',
#   'Temp_seasonality', 
#   'Max_temp_warm_month', 
#   'Min_temp_cold_month', 
#   'Temp_annual_range',
#   'Mean_temp_wet_qu', 
#   'Mean_temp_dry_qu', 
#   'Mean_temp_warm_qu',
#   'Mean_temp_cold_qu',
#   'Annual_precip', 
#   'Precip_wet_month', 
#   'Precip_dry_month', 
#   'Precip_seasonality', 
#   'Precip_wet_qu', 
#   'Precip_dry_qu', 
#   'Precip_warm_qu', 
#   'Precip_col_qu') 
# 
# 
# ## Create new files in the 1km directory 
# i  <- match(sdm.select, pred_names)
# ff <- file.path('./data/base/worldclim/world/0.5/bio/current',
#                 sprintf('bio_%02d', i))
# dir.create(dirname(sub('0.5', '1km', ff)[1]), recursive = TRUE)
# 
# 
# ## Run a loop to project just the variables we will use
# lapply(ff, function(f) {
#   message(f)
#   gdalwarp(f, sub('0.5', '1km', f), tr = c(1000, 1000),
#            t_srs = '+init=esri:54009', r = 'bilinear', 
#            multi = TRUE)
# })
# 
# 
# ## Create a raster stack of the projected grids
# env.grids.current = stack(sub('0.5', '1km', ff))
# names(env.grids.current) <- pred_names[i]

#########################################################################################################################
## GET RANDOM BACKGROUND POINTS AND THEN FIT MAXENT WITH STANDARD VARIABLES
#########################################################################################################################


## Select background points randomly, then fit maxent
FIT_MAXENT_RAND_BG <- function(occ,
                               sdm.predictors,
                               # sdm.predictors is a vector of enviro conditions that you want to include
                               name,
                               outdir,
                               template.raster,
                               # template.raster is an empty raster with extent, res and projection
                               # of final output rasters. It is used to reduce
                               # occurrences to a single point per cell.
                               min_n,
                               # min_n is the minimum number of records (unique cells)
                               # required for a model to be fit
                               max_bg_size,
                               background_buffer_width,
                               shapefiles,
                               features,
                               replicates, # number of cross-validation replicates
                               responsecurves,
                               rep_args,
                               full_args) {
  
  ########################################################################
  ## First, stop if the outdir file exists,
  if(!file.exists(outdir)) stop('outdir does not exist :(', call. = FALSE)
  outdir_sp <- file.path(outdir, gsub(' ', '_', name))
  
  ## Also Skip species which have already been modelled: again, can this be made into an argument?
  if(file.exists(outdir_sp)) {
    
    print (paste ('Skip', name, ', Model results exist for this species'))
    
  } else {
    
    ## If the file doesn't exist, split out the features
    if(!file.exists(outdir_sp)) dir.create(outdir_sp)
    features <- unlist(strsplit(features, ''))
    
    ## Make sure user features are allowed: don't run the model if the
    ## features have been incorrectly specified in the main argument
    ## l: linear
    ## p: product
    ## q: quadratic
    ## h: hinge
    ## t:
    if(length(setdiff(features, c('l', 'p', 'q', 'h', 't'))) > 1)
      stop("features must be a vector of one or more of ',
           'l', 'p', 'q', 'h', and 't'.")
    
    ## Aggregate
    b <- aggregate(gBuffer(occ, width = background_buffer_width, byid = TRUE))
    
    #####################################################################
    ## Get unique cell numbers for species occurrences
    ## Can we make the template raster 10km?
    cells <- cellFromXY(template.raster, occ)
    
    ## Clean out duplicate cells and NAs (including points outside extent of predictor data)
    ## Note this will get rid of a lot of duplicate records not filtered out by GBIF columns, etc.
    not_dupes <- which(!duplicated(cells) & !is.na(cells))
    occ       <- occ[not_dupes, ]
    cells     <- cells[not_dupes]
    message(nrow(occ), ' occurrence records (unique cells).')
    
    ## Skip species that have less than a minimum number of records: eg 20 species
    if(nrow(occ) < min_n) {
      
      print (paste ('Fewer occurrence records than the number of cross-validation ',
                    'replicates for species ', name,
                    ' Model not fit for this species'))
      
    } else {
      
      #####################################################################
      ## Now subset the background records to the raster cells with data
      f <- tempfile()
      writeOGR(SpatialPolygonsDataFrame(b, data.frame(id = seq_along(b))), 
               tempdir(), basename(f), 'ESRI Shapefile')
      
      ## Create a raster of cells with data
      bg_cells <- gdal_rasterize(extension(f, 'shp'), 
                                 tempfile(fileext = '.tif'),
                                 te = c(bbox(template.raster)),
                                 tr = res(template.raster),
                                 burn = 1, init = 0, ot = 'Byte', 
                                 output_Raster = TRUE) %>% 
                                   
        .[[1]] %>% # gdal_rasterize returns a RasterBrick but Which needs a RasterLayer, so we grab the first (only) layer
        Which(. == 1, cells = TRUE) %>% 
        intersect(template.cells)
      
      bg <- xyFromCell(template.raster, bg_cells)
      
      ## Reduce background sample if it's larger than max_bg_size
      if (nrow(bg) > max_bg_size) {
        
        message(nrow(bg), ' target species background records, reduced to random ',
                max_bg_size, '.')
        
        bg <- bg[sample(nrow(bg), max_bg_size), ]  ##  
        
      } else {
        
        message(nrow(bg), ' target species background records.')
        
      }
      
      ## Extract environmental values at the background points
      bg <- extract(env.grids.current, bg, sp = TRUE)
      bg <- SpatialPointsDataFrame(SpatialPoints(bg), data.frame(cell=bg_cells))
      
      #####################################################################
      ## Save shapefiles for future reference
      if(shapefiles) {
        
        suppressWarnings({
          
          writeOGR(SpatialPolygonsDataFrame(b, data.frame(ID = seq_len(length(b)))),
                   outdir_sp, 'bg_buffer', 'ESRI Shapefile', overwrite_layer = TRUE)
          writeOGR(bg,  outdir_sp, 'bg',   'ESRI Shapefile', overwrite_layer = TRUE)
          writeOGR(occ, outdir_sp, 'occ',  'ESRI Shapefile', overwrite_layer = TRUE)
          
        })
        
      }
      
      ## Save the background and occurrence points as .rds objects
      saveRDS(bg,  file.path(outdir_sp, 'bg.rds'))
      saveRDS(occ, file.path(outdir_sp, 'occ.rds'))
      
      #####################################################################
      ## sdm.predictors: the s_ff object containing sdm.predictors to use in the model
      ## Sample sdm.predictors at occurrence and background points
      #####################################################################
      swd_occ <- occ[, sdm.predictors]
      saveRDS(swd_occ, file.path(outdir_sp, 'occ_swd.rds'))
      
      swd_bg <- bg[, sdm.predictors]
      saveRDS(swd_bg, file.path(outdir_sp, 'bg_swd.rds'))
      
      ## Save shapefiles of the occurrence and background points
      if(shapefiles) {
        
        writeOGR(swd_occ, outdir_sp,  'occ_swd', 'ESRI Shapefile', overwrite_layer = TRUE)
        writeOGR(swd_bg,  outdir_sp,  'bg_swd',  'ESRI Shapefile', overwrite_layer = TRUE)
        
      }
      
      #####################################################################
      ## Combine occurrence and background data
      swd <- as.data.frame(rbind(swd_occ@data, swd_bg@data))
      saveRDS(swd, file.path(outdir_sp, 'swd.rds'))
      pa <- rep(1:0, c(nrow(swd_occ), nrow(swd_bg)))
      
      ## Now check the features arguments are correct
      off <- setdiff(c('l', 'p', 'q', 't', 'h'), features)
      
      ## 
      if(length(off) > 0) {
        
        off <- c(l = 'linear=false',    p = 'product=false', q = 'quadratic=false',
                 t = 'threshold=false', h = 'hinge=false')[off]
        
      }
      
      off <- unname(off)
      
      if(replicates > 1) {
        
        if(missing(rep_args)) rep_args <- NULL
        
        ## Run MAXENT for X cross validation data splits of swd : so 5 replicaes, 0-4
        ## EG xval = cross validation : "OUT_DIR\Acacia_boormanii\xval\maxent_0.html"
        me_xval <- maxent(swd, pa, path = file.path(outdir_sp, 'xval'),
                          args = c(paste0('replicates=', replicates),
                                   #paste0('biasfile=', dens),
                                   'responsecurves=true',
                                   'outputformat=logistic', # "biasfile = dens.ras"
                                   off, paste(names(rep_args), rep_args, sep = '=')))
        
      }
      
      ## Runs the full maxent model - presumably using all the data in swd
      if(missing(full_args)) full_args <- NULL
      me_full <- maxent(swd, pa, path = file.path(outdir_sp, 'full'),
                        args = c(off, paste(names(full_args), full_args, sep = '='),
                                 'responsecurves=true',
                                 'outputformat=logistic'))
      
      ## Save the full model?
      saveRDS(list(me_xval = me_xval, me_full = me_full, swd = swd, pa = pa), 
              file.path(outdir_sp, 'full', 'maxent_fitted.rds'))
      
      #####################################################################
      ## Save the chart corrleation file too for the variable set
      save_name = gsub(' ', '_', name)
      png(sprintf('%s/%s/full/%s_%s.png', outdir,
                  save_name, save_name, "predictor_correlation"),
          3236, 2000, units = 'px', res = 300)
      
      ## set margins
      par(mar   = c(3, 3, 5, 3),  ## b, l, t, r
          #mgp   = c(9.8, 2.5, 0),
          oma   = c(1.5, 1.5, 1.5, 1.5))
      
      ## Add detail to the response plot
      chart.Correlation(swd_occ@data,
                        histogram = TRUE, pch = 19) 
      #cex.lab = 2, cex.axis = 1.5,
      #main = paste0("Predictor corrleation matrix for ", spp))
      
      ## Finish the device
      dev.off()
      
    }
    
  }
  
}