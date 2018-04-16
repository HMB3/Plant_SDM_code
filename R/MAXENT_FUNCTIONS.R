#########################################################################################################################
############################################# FUNCTIONS FOR RUNNING SDMs AND MAPPING #################################### 
#########################################################################################################################


## flag issues with ..........................................................................


#########################################################################################################################
## MAXENT FUNCTIONS
#########################################################################################################################


## Arguments to run maxent line by line
# occ                     = occurrence
# bg                      = background
# name                    = spp
# outdir                  = 'output/maxent/SET_VAR_COORDCLEAN'
# sdm.predictors          = sdm.select
# env.grids               = env.grids.current
# template.raster         = template.raster
# template.cells          = template.cells
# min_n                   = 20   ## This should be higher...
# max_bg_size             = 100000 ## need a min bg size?
# background_buffer_width = 200000
# shapefiles              = TRUE
# features                = 'lpq'
# replicates              = 5
# responsecurves          = TRUE


## Selection line by line 
# occ             = swd_occ
# bg              = swd_bg
# path            = outdir
# species_column  = "species"
# replicates      = replicates
# response_curves = TRUE
# logistic_format = TRUE
# cor_thr         = 0.85
# pct_thr         = 5
# k_thr           = 5
# features        ='lpq'  # change these as necessary (or cor_thr = cor_thr, etc from FIT_MAXENT_SIMP)
# quiet           = FALSE
# type            = "PI"


#########################################################################################################################
## GET RANDOM BACKGROUND POINTS AND THEN FIT MAXENT WITH STANDARD VARIABLES
#########################################################################################################################


## Select background points randomly, then fit maxent
FIT_MAXENT_RAND_BG <- function(occ,
                               name,
                               outdir,
                               sdm.predictors,
                               env.grids,
                               template.raster,
                               template.cells,
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
  
  ## There is no minimum number of background records.................................................................
  
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
    cells <- cellFromXY(template.raster, occ)
    
    ## Clean out duplicate cells and NAs (including points outside extent of predictor data)
    ## Note this will get rid of a lot of duplicate records not filtered out by GBIF columns, etc.
    not_dupes <- which(!duplicated(cells) & !is.na(cells))
    occ       <- occ[not_dupes, ]
    cells     <- cells[not_dupes]
    
    ##
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
      
      ## There is no minimum number of background records.................................................................
      
      ## What if the species has very few records surrounding its record?
      if (nrow(bg) > max_bg_size) {
        
        message(nrow(bg), ' target species background records, reduced to random ',
                max_bg_size, '.')
        
        bg <- bg[sample(nrow(bg), max_bg_size), ]  ##  
        
      } else {
        
        message(nrow(bg), ' target species background records.')
        
      }
      
      ## Extract environmental values at the background points :: this was not working for me...
      ## projection(env.grids.current)
      bg <- SpatialPointsDataFrame(coords      = as.data.frame(bg), 
                                   data        = as.data.frame(bg),
                                   proj4string = CRS('+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'))
      
      ## Is the cell = bg_cells necessary?
      bg <- extract(env.grids, bg) ## sp = TRUE)
      ## bg <- SpatialPointsDataFrame(SpatialPoints(bg), data.frame(cell = bg_cells)) 
      
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
      ## sdm.predictors: the s_ff object containing sdm.predictors to use 
      ## in the model: sample sdm.predictors at occurrence and background points
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
      
      ## Get the maxent model
      # m <- readRDS(sprintf('%s/%s/full/maxent_fitted.rds', outdir, name)) 
      # m <- m$me_full  ## class(m);View(m)
      
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
      
      ########################################################################################################################
      ## Another .png for the global records: str(LAND$long)
      LAND       = readRDS("F:/green_cities_sdm/data/base/CONTEXTUAL/LAND_world.rds")
      
      png(sprintf('%s/%s/full/%s_%s.png', outdir,
                  name, name, "global_records"),
          16180, 10000, units = 'px', res = 600)
      
      ## How do we locate bad records in the dataset after spotting them?
      plot(LAND, 
           lwd = 0.5, asp = 1, axes = TRUE, cex.axis = 3.5,
           col = 'darkolivegreen3', bg = 'lightblue', cex.lab = 3)
      
      points(occ, pch = ".", cex = 3.5, col = "red", cex.lab = 3, cex.main = 4, cex.axis = 2, 
             main = paste0("Global occurrences for ", name), 
             xlab = "", ylab = "", asp = 1)
      
      ## Title 
      title(paste0("Global points for ", name),
            cex.main = 4,   font.main = 4, col.main = "blue")
      
      ## Finsh the device
      dev.off()
      
      ########################################################################################################################
      ## Another PNG for the background points....
      png(sprintf('%s/%s/full/%s_%s.png', outdir,
                  name, name, "background_records"),
          16180, 10000, units = 'px', res = 600)
      
      ## How do we locate bad records in the dataset after spotting them?
      plot(LAND,  
           lwd = 0.5, asp = 1, axes = TRUE, cex.axis = 3.5,
           col = 'darkolivegreen3', bg = 'lightblue', cex.lab = 3)
      
      points(bg, pch = ".", cex = 1.6, col = "blue", cex.lab = 3, cex.main = 4, cex.axis = 2, 
             main = paste0("Global occurrences for ", name), 
             xlab = "", ylab = "", asp = 1)
      
      ## Title 
      title(paste0("Bacground points for ", name),
            cex.main = 4,   font.main = 4, col.main = "blue")
      
      ## Finish the device
      dev.off() 
      
      ########################################################################################################################
      ## Plot the models: can two plots be combined into one?
      ## Make these unique names, and they can be searched in windows. Otherwise, we can just click into each subfolder. 
      ## To sort, names would need to be: spp + unique_extension
      png(sprintf('%s/%s/full/%s_current.png', outdir,
                  name, name, "variable_contribution"),
          3236, 2000, units = 'px', res = 300)
      
      ## Set the margins
      # par(mgp      = c(10, 4, 0), 
      #     oma      = c(1.5, 1.5, 1.5, 1.5),
      #     font.lab = 2)
      
      plot(me_full, col = "blue", pch = 19, cex.lab = 2, cex.axis = 5, cex.main = 2, 
           main   = paste0("Variables for ", name), 
           xlab   = "Maxent contribution (%)")
      
      ## Finish the device
      dev.off()
      
      ## Plot the response curves too
      png(sprintf('%s/%s/full/%s_%s.png', outdir,
                  name, name, "response_curves"),
          3236, 2000, units = 'px', res = 300)
      
      ## Add detail to the response plot
      response(me_full, pch = 19, cex.lab = 2, cex.axis = 1.5, lwd = 2) 
      
      ## Finish the device
      dev.off()
      
    }
    
  }
  
}





#########################################################################################################################
## GET TARGETTED BACKGROUND POINTS AND THEN FIT MAXENT WITH STANDARD VARIABLES
#########################################################################################################################


##
FIT_MAXENT_TARG_BG <- function(occ,
                               bg, # A Spatial points data frame (SPDF) of candidate background points
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
    
    
    #####################################################################
    ## Skip species that have less than a minimum number of records: eg 20 species
    if(nrow(occ) < min_n) {
      
      print (paste ('Fewer occurrence records than the number of cross-validation ',
                    'replicates for species ', name,
                    ' Model not fit for this species'))
      
    } else {
      
      #####################################################################
      ## Now subset the background records to the buffered polygon
      system.time(o <- over(bg, b))
      bg <- bg[which(!is.na(o)), ]
      bg_cells <- cellFromXY(template.raster, bg)
      
      ## Clean out duplicates and NAs (including points outside extent of predictor data)
      bg_not_dupes <- which(!duplicated(bg_cells) & !is.na(bg_cells))
      bg <- bg[bg_not_dupes, ]
      bg_cells <- bg_cells[bg_not_dupes]
      
      ## To incorporate bias in the background records, we can sample from the bias file
      # bg_bias <- xyFromCell(dens$kde, sample(which(!is.na(values(dens$kde))), 10000,
      #                                        prob = values(dens$kde)[!is.na(values(dens$kde))],
      #                                                                             replace = TRUE))
      
      ## Reduce background sample if it's larger than max_bg_size
      if (nrow(bg) > max_bg_size) {
        
        message(nrow(bg), ' target species background records, reduced to random ',
                max_bg_size, '.')
        
        bg <- bg[sample(nrow(bg), max_bg_size), ]  ## Change this to use 
        
      } else {
        
        message(nrow(bg), ' target species background records.')
        
      }
      
      #####################################################################
      ## Save objects for future reference
      if(shapefiles) {
        
        suppressWarnings({
          
          writeOGR(SpatialPolygonsDataFrame(b, data.frame(ID = seq_len(length(b)))),
                   outdir_sp, 'bg_buffer', 'ESRI Shapefile', overwrite_layer = TRUE)
          writeOGR(bg,  outdir_sp, 'bg',   'ESRI Shapefile', overwrite_layer = TRUE)
          writeOGR(occ, outdir_sp, 'occ',  'ESRI Shapefile', overwrite_layer = TRUE)
          
        })
        
      }
      
      ## Save the background and occurrence points as objects
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
        
        ## Run MAXENT for x cross validation data splits of swd : so 5 replicaes, 0-4
        ## EG xval = cross validation : "OUT_DIR\Acacia_boormanii\xval\maxent_0.html"
        me_xval <- maxent(swd, pa, path = file.path(outdir_sp, 'xval'),
                          args = c(paste0('replicates=', replicates),
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
      
      ########################################################################################################################
      ## Another .png for the global records: str(LAND$long)
      LAND       = readRDS("F:/green_cities_sdm/data/base/CONTEXTUAL/LAND_world.rds")
      occ_land   = occ %>% 
        spTransform(CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
      bg_land   = bg %>% 
        spTransform(CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
      
      png(sprintf('%s/%s/full/%s_%s.png', maxent_path,
                  save_name, save_name, "global_records"),
          16180, 10000, units = 'px', res = 600)
      
      ## How do we locate bad records in the dataset after spotting them?
      plot(LAND, 
           lwd = 0.5, asp = 1, axes = TRUE, cex.axis = 3.5,
           col = 'darkolivegreen3', bg = 'lightblue', cex.lab = 3)
      
      points(occ_land, pch = ".", cex = 3.5, col = "red", cex.lab = 3, cex.main = 4, cex.axis = 2, 
             main = paste0("Global occurrences for ", name), 
             xlab = "", ylab = "", asp = 1)
      
      ## Title 
      title(paste0("Global points for ", name),
            cex.main = 4,   font.main = 4, col.main = "blue")
      
      ## Finsh the device
      dev.off()
      
      ########################################################################################################################
      ## Another PNG for the background points....
      #if(!file.exists(sprintf('%s/%s/full/%s_%s.png', maxent_path, name, name, "background_records"))) {
      png(sprintf('%s/%s/full/%s_%s.png', outdir,
                  name, name, "background_records"),
          16180, 10000, units = 'px', res = 600)
      
      ## How do we locate bad records in the dataset after spotting them?
      plot(LAND,  
           lwd = 0.5, asp = 1, axes = TRUE, cex.axis = 3.5,
           col = 'darkolivegreen3', bg = 'lightblue', cex.lab = 3)
      
      points(bg, pch = ".", cex = 1.6, col = "blue", cex.lab = 3, cex.main = 4, cex.axis = 2, 
             main = paste0("Bacground points for ", name), 
             xlab = "", ylab = "", asp = 1)
      
      ## Title 
      title(paste0("Bacground points for ", name),
            cex.main = 4,   font.main = 4, col.main = "blue")
      
      ## Finish the device
      dev.off() 
      
      # } else {
      #   
      #   message("Background records maps exists for ", name)
      #   
      # }
      
      ########################################################################################################################
      ## Plot the models: can two plots be combined into one?
      ## Make these unique names, and they can be searched in windows. Otherwise, we can just click into each subfolder. 
      ## To sort, names would need to be: spp + unique_extension
      m <- readRDS(sprintf('%s/%s/full/maxent_fitted.rds', outdir, save_name)) 
      m <- m$me_full
      
      png(sprintf('%s/%s/full/%s_current.png', outdir,
                  name, name, "variable_contribution"),
          3236, 2000, units = 'px', res = 300)
      
      ## Set the margins
      # par(mgp      = c(10, 4, 0), 
      #     oma      = c(1.5, 1.5, 1.5, 1.5),
      #     font.lab = 2)
      
      plot(m, col = "blue", pch = 19, cex.lab = 2, cex.axis = 5, cex.main = 2, 
           main   = paste0("Variables for ", name), 
           xlab   = "Maxent contribution (%)")
      
      ## Finish the device
      dev.off()
      
      ## Plot the response curves too
      png(sprintf('%s/%s/full/%s_%s.png', outdir,
                  name, name, "response_curves"),
          3236, 2000, units = 'px', res = 300)
      
      ## Add detail to the response plot
      response(m, pch = 19, cex.lab = 2, cex.axis = 1.5, lwd = 2) 
      
      ## Finish the device
      dev.off()
      
      #####################################################################
      ## Save fitted model object, and the model-fitting data.
      #       if(replicates > 1) {
      # 
      #         saveRDS(list(me_xval = me_xval, me_full = me_full, swd = swd, pa = pa),
      #                 file.path(outdir_sp, 'maxent_fitted.rds'))
      # 
      #       } else {
      # 
      #         saveRDS(list(me_xval = NA, me_full = me_full, swd = swd, pa = pa),
      #                 file.path(outdir_sp, 'maxent_fitted.rds'))
      # 
      # }
      
    }
    
  }
  
}








#########################################################################################################################
## GET BACKGROUND POINTS AND THEN FIT MAXENT WITH BACKWARDS SELECTION 
#########################################################################################################################


FIT_MAXENT_SELECTION <- function(occ, 
                                 bg, # A Spatial points data frame (SPDF) of candidate background points
                                 sdm.predictors, 
                                 # sdm.predictors is a vector of enviro conditions that you want to include
                                 name  = spp, 
                                 outdir, 
                                 template.raster, 
                                 # template.raster is an empty raster with extent, res and projection
                                 # of final output rasters. It is used to reduce
                                 # occurrences to a single point per cell.
                                 min_n,
                                 # min_n is the minimum number of records (unique cells)
                                 # required for a model to be fit
                                 max_bg_size, # 
                                 background_buffer_width,
                                 shapefiles, 
                                 features, 
                                 replicates, # number of cross-validation replicates
                                 cor_thr, 
                                 pct_thr, 
                                 k_thr, 
                                 responsecurves) {
  
  ########################################################################
  ## First, stop if the outdir file exists, 
  #browser()
  if(!file.exists(outdir)) stop('outdir does not exist :(', call. = FALSE)
  outdir_sp <- file.path(outdir, gsub(' ', '_', name))
  
  ## Also Skip species which have already been modelled: again, can this be made into an argument?
  if(file.exists(outdir_sp)) {
    
    print (paste ('Skip', name, ', Model results exist for this species'))
    
  } else {
    
    ## If the file doesn't exist, split out the features 
    if(!file.exists(outdir_sp)) dir.create(outdir_sp)
    
    ## Select background records from within 200km of the target species occurrence records
    b <- aggregate(gBuffer(occ, width = background_buffer_width, byid = TRUE))
    
    
    #####################################################################
    ## Get unique cell numbers for species occurrences
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
      ## Now subset bg to the buffer polygon
      
      ## Describe the background procedure here....................................................................
      
      ## Within the 200km buffer around the occurrence records for each species, select up to 100k records from any
      ## species in the total dataset :: trying to create the same bias in both the occurrence and background data
      system.time(o <- over(bg, b))
      bg <- bg[which(!is.na(o)), ]
      bg_cells <- cellFromXY(template.raster, bg)
      
      ## Clean out duplicates and NAs (including points outside extent of predictor data)
      ## So we are using unique cells to select background records
      bg_not_dupes <- which(!duplicated(bg_cells) & !is.na(bg_cells)) 
      bg <- bg[bg_not_dupes, ]
      bg_cells <- bg_cells[bg_not_dupes]
      
      ## Reduce background sample if it's larger than max_bg_size
      if (nrow(bg) > max_bg_size) { 
        
        message(nrow(bg), ' target species background records, reduced to random ', 
                max_bg_size, '.')
        
        bg <- bg[sample(nrow(bg), max_bg_size), ]
        
      } else {
        
        message(nrow(bg), ' target species background records.')
        
      }
      
      #####################################################################
      ## Save objects for future reference
      if(shapefiles) {
        
        suppressWarnings({
          
          writeOGR(SpatialPolygonsDataFrame(b, data.frame(ID = seq_len(length(b)))), 
                   outdir_sp, 'bg_buffer', 'ESRI Shapefile', overwrite_layer = TRUE)  
          writeOGR(bg,  outdir_sp, 'bg',   'ESRI Shapefile', overwrite_layer = TRUE)  
          writeOGR(occ, outdir_sp, 'occ',  'ESRI Shapefile', overwrite_layer = TRUE) 
          
        })
        
      }
      
      ## Save the background and occurrence points as objects
      saveRDS(bg,  file.path(outdir_sp, 'bg.rds'))
      saveRDS(occ, file.path(outdir_sp, 'occ.rds'))
      
      #####################################################################
      ## Now get the predictors at the occurrence and background points 
      swd_occ <- occ[, sdm.predictors]
      swd_occ$species <- name           
      saveRDS(swd_occ, file.path(outdir_sp, 'occ_swd.rds'))
      
      swd_bg <- bg[, sdm.predictors]
      swd_bg$species <- name
      saveRDS(swd_bg, file.path(outdir_sp, 'bg_swd.rds'))
      
      ## Save shapefiles of the occurrence and background points
      if(shapefiles) {
        writeOGR(swd_occ, outdir_sp,  'occ_swd', 'ESRI Shapefile', overwrite_layer = TRUE)
        writeOGR(swd_bg,  outdir_sp,  'bg_swd',  'ESRI Shapefile', overwrite_layer = TRUE)
      }
      
      #####################################################################
      ## coerce to regular data.frames
      swd_occ <- as.data.frame(swd_occ)
      swd_occ$lon <- NULL
      swd_occ$lat <- NULL
      swd_bg <- as.data.frame(swd_bg)
      swd_bg$lon <- NULL
      swd_bg$lat <- NULL
      
      #####################################################################
      ## Run simplify rmaxent::simplify
      m <- rmaxent::simplify(
        
        swd_occ, 
        swd_bg,
        path            = outdir, 
        species_column  = "species",
        replicates      = replicates,
        response_curves = TRUE, 
        logistic_format = TRUE, 
        cor_thr         = cor_thr, 
        pct_thr         = pct_thr, 
        k_thr           = k_thr, 
        features        = features,  ## change these as necessary (or cor_thr = cor_thr, etc from FIT_MAXENT_SIMP)
        quiet           = FALSE)
      
    }
    
  }
  
}





#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################