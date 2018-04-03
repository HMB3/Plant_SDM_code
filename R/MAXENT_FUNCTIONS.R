#########################################################################################################################
############################################# FUNCTIONS FOR RUNNING SDMs AND MAPPING #################################### 
#########################################################################################################################


## flag issues with ..........................................................................


#########################################################################################################################
## MAXENT FUNCTIONS
#########################################################################################################################


## Arguments to run maxent line by line
# occ                     = occurrence
# bg                      = background,
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
      
      ## What if the species has vert few records surrounding its record?
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
      #if(!file.exists(sprintf('%s/%s/full/%s_%s.png', outdir_sp, name, name, "global_records"))) {
      #LAND       = readRDS("F:/green_cities_sdm/data/base/CONTEXTUAL/LAND_world.rds")
      
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
      
      # } else {
      #   
      #   message("Global records maps exists for ", name)
      #   
      # }
      
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
             main = paste0("Global occurrences for ", name), 
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
      #if(!file.exists(sprintf('%s/%s/full/%s_current.png', maxent_path, name, name, "variable_contribution"))) {
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
      
      # } else {
      #   
      #   message("Variable contribution plot exists for ", name)
      #   
      # }
      
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
      if(!file.exists(sprintf('%s/%s/full/%s_%s.png', maxent_path, species, species, "global_records"))) {
        
        png(sprintf('%s/%s/full/%s_%s.png', maxent_path,
                    species, species, "global_records"),
            16180, 10000, units = 'px', res = 600)
        
        ## How do we locate bad records in the dataset after spotting them?
        plot(LAND, 
             lwd = 0.5, asp = 1, axes = TRUE, cex.axis = 3.5,
             col = 'darkolivegreen3', bg = 'lightblue', cex.lab = 3)
        
        points(occ, pch = ".", cex = 3.5, col = "red", cex.lab = 3, cex.main = 4, cex.axis = 2, 
               main = paste0("Global occurrences for ", species), 
               xlab = "", ylab = "", asp = 1)
        
        ## Title 
        title(paste0("Global points for ", species),
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
               main = paste0("Global occurrences for ", name), 
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
        m <- readRDS(sprintf('%s/%s/full/maxent_fitted.rds', outdir, name)) 
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
## PROJECT MAXENT MODELS
#########################################################################################################################


#########################################################################################################################
## Create maxent maps for a given time period 
## x = scen_2030[1]
## species = kop_spp[3]
## time_slice = 30
## maxent_path  = "./output/maxent/SET_VAR_CLEAN"
## climate_path = "./data/base/worldclim/aus/0.5/bio"
## grid_names    = grid.names
## current_grids = env.grids.current
project_maxent_grids = function(scen_list, species_list, maxent_path, climate_path, grid_names, time_slice, current_grids) {
  
  ## First, run a loop over each scenario:    
  lapply(scen_list, function(x) {
    
    ## Create a raster stack for each of the 6 GCMs, not for each species
    s <- stack(sprintf('%s/20%s/%s/%s%s.tif', climate_path, time_slice, x, x, 1:19))
    
    ## Rename both the current and future environmental stack...
    names(s) <- names(current_grids) <- grid_names 
    
    ########################################################################################################################
    ## Divide the 11 temperature rasters by 10: NA values are the ocean
    ## s[[1:11]] <- s[[1:11]]/10 ## That code doesn't work
    message('20', time_slice, ' rasters / 10 ', x)
    for(i in 1:11) {
      
      ## Simple loop
      message(i)
      s[[i]] <- s[[ i]]/10
      
    }
    
    ## Then apply each GCM to each species
    lapply(species_list, function(species) {
      
      ## First, check if the maxent model exists
      if(file.exists(sprintf('%s/%s/full/maxent_fitted.rds', maxent_path, species))) {
        message('Doing ', species)
        
        ## Then, check if the species projection has already been run...
        if(!file.exists(sprintf('%s/%s/full/%s_%s.tif', maxent_path, species, species, x))) {
          
          ## Assign the scenario name (to use later in the plot)
          scen_name = eval(parse(text = sprintf('gcms.%s$GCM[gcms.%s$id == x]', time_slice, time_slice)))           
          
          ########################################################################################################################
          ## Now read in the SDM model calibrated on current conditions  ## maxent_fitted.rds
          #m <- readRDS(sprintf('%s/%s/full/model.rds', maxent_path, species)) 
          m <- readRDS(sprintf('%s/%s/full/maxent_fitted.rds', maxent_path, species)) 
          m <- m$me_full  ## class(m);View(m)
          
          ## Read in the occurrence points used to create the SDM :: need the transform to plot later
          occ <- readRDS(sprintf('%s/%s/occ.rds', maxent_path, species)) %>%
            spTransform(CRS('+init=epsg:4326'))
          
          ## And if the current raster doesn't exist, create it
          f_current <- sprintf('%s/%s/full/%s_current.tif', maxent_path, species, species)
          
          if(!file.exists(f_current)) {
            
            ## Report which prediction is in progress :: m$me_full, m$me_full@presence
            message('Running current prediction for ', species) 
            
            pred.current <- rmaxent::project(
              m, current_grids[[colnames(m@presence)]])$prediction_logistic
            writeRaster(pred.current, f_current, overwrite = TRUE)
            
          } else {
            
            pred.current = raster(sprintf('%s/%s/full/%s_current.tif', 
                                          maxent_path, species, species))
          }
          
          ########################################################################################################################
          ## If the future raster doesn't exist, create it 
          f_future <- sprintf('%s/%s/full/%s_%s.tif', 
                              maxent_path, species, species, x)
          
          if(!file.exists(f_future)) {
            
            ## Report which prediction is in progress
            message('Running future prediction for ', species, ' ', x) 
            
            ## Error in sum(sapply(lambdas, nrow)) : invalid 'type' (list) of argument
            pred.future <- rmaxent::project(
              m, s[[colnames(m@presence)]])$prediction_logistic
            writeRaster(pred.future, f_future, overwrite = TRUE)
            
            ## Now create the empty panel just before plotting
            empty <- init(pred.current, function(x) NA) 
            
            ## Extents don't match between current and future data .................................................................
            ## workaround for difference in raster extents :: annoying, current rasters are smaller than future
            ex          = extent(pred.current)
            pred.f      = crop(pred.future, ex)
            
            ########################################################################################################################
            ## Use the levelplot function to make a multipanel output: occurrence points, current raster and future raster
            png(sprintf('%s/%s/full/%s_%s.png', maxent_path, species, species, x),      
                11, 4, units = 'in', res = 300)
            
            ## Need an empty frame
            print(levelplot(stack(empty,
                                  pred.current, 
                                  pred.f, quick = TRUE), margin = FALSE,
                            
                            ## Create a colour scheme using colbrewer: 100 is to make it continuos
                            ## Also, make it a one-directional colour scheme
                            scales      = list(draw = FALSE), 
                            at = seq(0, 1, length = 100),
                            col.regions = colorRampPalette(rev(brewer.pal(11, 'Spectral'))),
                            
                            ## Give each plot a name: the third panel is the GCM
                            names.attr = c('Occurrence', 'Current', sprintf('%s, 20%s, RCP8.5', scen_name, time_slice)),
                            colorkey   = list(height = 0.5, width = 3), xlab = '', ylab = '',
                            main       = list(gsub('_', ' ', species), font = 4, cex = 2)) +
                    
                    ## Plot the Aus shapefile with the occurrence points for reference
                    ## Can the points be made more legible for both poorly and well recorded species?
                    layer(sp.polygons(shapefile)) +
                    layer(sp.points(occ, pch = 20, cex = 0.4, 
                                    col = c('red', 'transparent', 'transparent')[panel.number()]), data = list(occ = occ)))
            dev.off()
            
          }
          
        } else {
          
          message(species, ' ', x, ' skipped - prediction already run') ## Skip species which have already been projected
          
        }
        
      } else {
        
        message(species, ' ', x, ' skipped - SDM not yet run')          ## Skip species with no existing maxent model
        
      }
      
    })
    
  })
  
}





#########################################################################################################################
## COMBINE GCM FUNCTIONS
#########################################################################################################################


#########################################################################################################################
## Loop over directories, species and one threshold for each, also taking a time_slice argument. Next, make the lists generic too
## pervious version in R/old/model_combine.R
combine_gcm_threshold = function(DIR_list, species_list, maxent_path, thresholds, percentiles, time_slice, area_occ) {
  
  ## How can the shapefiles be read in once, not for each species?...................................................
  
  ###################################################################################################################
  ## Read in shapefiles :: this should be done outside the loop
  aus        = readRDS("F:/green_cities_sdm/data/base/CONTEXTUAL/aus_states.rds")
  LAND       = readRDS("F:/green_cities_sdm/data/base/CONTEXTUAL/LAND_world.rds")
  areal_unit = readRDS("F:/green_cities_sdm/data/base/CONTEXTUAL/SUA.rds")
  areal_unit = areal_unit[order(areal_unit$SUA_NAME11),]
  
  ## Loop over each directory
  lapply(DIR_list, function(DIR) { 
    
    ## And each species - although we don't want all possible combinations. use mapply in the function call
    lapply(species_list, function(species) {
      
      # for (slice in time_slice) {
      
      ###################################################################################################################
      ## Create a list of the rasters in each directory, then take the mean. Tidy this up with %, etc
      message('Calcualting mean of GCMs for ', species)
      
      ## Check if the mean GCM raster exists
      f_mean = sprintf('%s/%s/full/%s_20%s_suitability_mean.tif', maxent_path, species, species, time_slice)
      
      ## The mean of the GCMs doesn't exist, create it
      if(!file.exists(f_mean)) {
        
        raster.list       = list.files(as.character(DIR), pattern = sprintf('bi%s.tif', time_slice), full.names = TRUE)  
        suit              = stack(raster.list)
        suit.list         = unstack(suit)
        combo_suit_mean   = mean(suit)                            ## plot(mean.suit)
        
        writeRaster(combo_suit_mean , sprintf('%s/%s/full/%s_20%s_suitability_mean.tif', 
                                              maxent_path, species, species, time_slice), overwrite = TRUE)
        
      } else {
        
        ## Create another level without the mean calculation
        raster.list = list.files(as.character(DIR), pattern = sprintf('bi%s.tif', time_slice), full.names = TRUE)  
        suit        = stack(raster.list)
        suit.list   = unstack(suit)
        
      }
      
      #########################################################################################################################
      ## Then, create rasters that meet habitat suitability criteria thresholds, determined by the rmaxent function
      for (thresh in thresholds) {
        
        for (percent in percentiles) { 
          
          ## Can we set 0 to NA in the rasters before running the calculations?
          ## Check if the combined suitability raster exists
          f_max_train_suit <- sprintf('%s/%s/full/%s_20%s%s%s.tif', maxent_path,
                                      species, species, time_slice, "_Max_train_sensit_above_", thresh)
          
          ## If it doesn't exist, create the suitability raster
          if(!file.exists(f_max_train_suit)) {
            
            ## Print the species being analysed
            message('doing ', species, ' | Max train sensit > ', thresh, ' for 20', time_slice)
            
            ## Read in the current suitability raster
            f_current <- raster(sprintf('%s/%s/full/%s_current.tif', 
                                        maxent_path, species, species))
            
            ## First, create a simple function to threshold each of the rasters in raster.list
            thresh_greater  = function (x) {x > thresh}
            percent_greater = function (x) {x > percent}
            
            ## Then apply this to just the current suitability raster. These functions use the : 
            ## Maximum training sensitivity plus specificity Logistic threshold
            ## 10th percentile training presence training omission
            current_suit_thresh  = thresh_greater(f_current)
            current_suit_percent = percent_greater(f_current) 
            
            ## Now, apply these functions to the list of rasters (6 GCMs) for each species
            ## Also, the 'init' function initializes a raster object with values
            ## suit_ras_greater    = reduce(suit.list, thresh_above_fun, .init = suit.list[[1]] > thresh)
            
            ## Check the logic of removing the zeros.....................................................
            ## If we are just using the combined rasters without calculating the difference, don't worry about the zeros   
            
            #########################################################################################################################
            ## First, calculate the cells which are greater that the: 
            ## Maximum training sensitivity plus specificity Logistic threshold
            message('Running thresholds for ', species, ' | 20', time_slice, ' combined suitability > ', thresh)
            
            suit_ras1_thresh   = thresh_greater(suit.list[[1]])   ## Abbreviate these...
            suit_ras2_thresh   = thresh_greater(suit.list[[2]])
            suit_ras3_thresh   = thresh_greater(suit.list[[3]])
            suit_ras4_thresh   = thresh_greater(suit.list[[4]])   
            suit_ras5_thresh   = thresh_greater(suit.list[[5]])
            suit_ras6_thresh   = thresh_greater(suit.list[[6]])
            
            ## Then calculate the cells which are greater than the 10th percentile training presence training omission
            suit_ras1_percent  = percent_greater(suit.list[[1]])
            suit_ras2_percent  = percent_greater(suit.list[[2]])
            suit_ras3_percent  = percent_greater(suit.list[[3]])
            suit_ras4_percent  = percent_greater(suit.list[[4]])
            suit_ras5_percent  = percent_greater(suit.list[[5]])
            suit_ras6_percent  = percent_greater(suit.list[[6]])
            
            #########################################################################################################################
            ## Then sum them up: All the threshholds
            combo_suit_thresh   =  Reduce("+", list(suit_ras1_thresh, suit_ras2_thresh, suit_ras3_thresh,
                                                    suit_ras4_thresh, suit_ras5_thresh, suit_ras6_thresh))
            
            ## All the percentiles
            combo_suit_percent  =  Reduce("+", list(suit_ras1_percent, suit_ras2_percent, suit_ras3_percent,
                                                    suit_ras4_percent, suit_ras5_percent, suit_ras6_percent))
            
            #########################################################################################################################
            ## For each species, create a binary raster with cells > 4 GCMs above the maxent threshold = 1, and cells with < 4 GCMs = 0. 
            message('Calculating change for ', species, ' | 20', time_slice, ' combined suitability > ', thresh)
            
            ## Functions for thresholding rasters
            #band_ras     <- function(x) {ifelse(x >=  4, 4, ifelse(x > 0 & x < 4, 3, x)) }
            band_4           <- function(x) {ifelse(x >=  4, 1, 0) }
            combo_suit_4GCM  <- calc(combo_suit_thresh, fun = band_4)
            
            ## Plot to check
            plot(current_suit_thresh, main = gsub('_', ' ', (sprintf('%s current max_train_sensit > %s', species, thresh))))
            plot(combo_suit_percent,  main = gsub('_', ' ', (sprintf('%s future 10th percentile > %s',   species, percent))))
            plot(combo_suit_thresh,   main = gsub('_', ' ', (sprintf('%s future max_train_sensit > %s',  species, thresh))))
            plot(combo_suit_4GCM,     main = gsub('_', ' ', (sprintf('%s 4+ GCMs > %s',  species, thresh))))
            
            #########################################################################################################################
            ## For each species, calculate the projected rainfall and temperature increase and decreas for each GCM? Could plot this
            ## as part of the final figure.
            
            #########################################################################################################################
            ## Then using this GCM consensus, calculate whether the species is likely to be present in each SUA.
            ## Decide on a threshold of % area (10?) of the SUA that needs to be occupied, for each species to be considered present. 
            message('Running zonal stats for ', species, ' | 20', time_slice, ' combined suitability > ', thresh)
            
            ## Check the order of lists match, species, SUAs, areas need to match up ................................................
            
            ## Should we also create a simple count of the no. of cells per SUA with where > 4 GCMs met the suitability threshold?
            ## Need to fix this so that it has the same order as the shapefile
            # sp.count <- spatialEco::zonal.stats(x = areal_unit, y = combo_suit_4GCM, stat = sum, trace = TRUE, plot = TRUE) 
            # z.med   <- spatialEco::zonal.stats(x = areal_unit, y = combo_suit_4GCM, stat = median, trace = TRUE, plot = TRUE)
            
            #########################################################################################################################
            ## Then, extract the values of the presence raster for each areal unit: generates a list: 32 seconds...
            ext  <- extract(combo_suit_4GCM, areal_unit, method = 'simple')
            
            ## A function to tabulate the raster values by aerial unit, returning a data frame
            tabFunc <- function(indx, extracted, region, regname) {
              
              dat<-as.data.frame(table(extracted[[indx]]))
              dat$name<-region[[regname]][[indx]]
              return(dat)
              
            }
            
            ## Run through each areal unit and calculate a table of the count of raster cells
            tabs <- lapply(seq(ext), tabFunc, ext, areal_unit, "SUA_NAME11")
            tabs <- do.call("rbind", tabs )
            
            ## Can we get the count here, so we don't need to run zonal stats as well
            # tabs.count      = tabs
            # tabs.count$Freq = ifelse(tabs.count$Var1 == 0, 0, tabs.count$Freq) 
            # sp.count <- aggregate(tabs.count$Freq, by = list(Category = tabs$name), FUN = sum)
            # names(sp.count) = c('SUA_NAME11', 'COUNT')
            # head(sp.count)
            
            #########################################################################################################################
            ## Now mutate the table
            PERECENT.AREA <- tabs %>%
              
              group_by(name) %>%                                          ## Group by region
              mutate(totcells = sum(Freq),                                ## How many cells overall?
                     percent.area = round(100 * Freq / totcells, 2)) %>%  ## Cells /total cells
              
              dplyr::select(-c(Freq, totcells)) %>%                       ## There is a select func in raster so need to specify
              spread(key = Var1, value = percent.area, fill = 0)     %>%  ## Make wide format
              as.data.frame()
            
            ## Rename and create a column for whether or not the species occupies that area 
            names(PERECENT.AREA) =  c('SUA_NAME11', 'Absent', 'Present') 
            PERECENT.AREA$species_present = ifelse(PERECENT.AREA$Present >= area_occ, 1, 0)
            head(PERECENT.AREA)
            
            ## Create a table with the columns: REGION, AREA SPECIES, STATS (for each time slice)
            ## Change this so that we can summarise across all species as suggested by Linda 
            GCM.AREA.SUMMARY <- data.frame(SUA         = areal_unit$SUA_NAME11, 
                                           AREA_SQKM   = areal_unit$AREA_SQKM,
                                           SPECIES     = species,
                                           PERIOD      = time_slice,
                                           AREA_THRESH = area_occ,
                                           MAX_TRAIN   = thresh,
                                           #CELL_COUNT = sp.count$COUNT,
                                           PERCENT_AREA = PERECENT.AREA$Present,
                                           PRESENT      = PERECENT.AREA$species_present)
            
            ## Rename columns using sprintf, so we can include the suitability threshold and the time slice
            #' names(GCM.AREA.SUMMARY) <-  c('SUA',
            #'                               'AREA_SQKM',
            #'                               'SPECIES',
            #'                               #'CELL_COUNT',
            #'                               sprintf("Percent_area_where_4GCMs > thresh in 20%s",   thresh, time_slice),
            #'                               sprintf('Species_present 4GCMs > %s in 20%s',   thresh, time_slice))
            View(GCM.AREA.SUMMARY) ## unique(GCM.AREA.SUMMARY$SPECIES)
            
            #########################################################################################################################
            ## Then save the table of SUA results for all species to a datafile...
            ## This would be the file to loop over to create a summary of species per SUA
            write.csv(GCM.AREA.SUMMARY, sprintf('%s/%s/full/%s_20%s_%s%s.csv', maxent_path,
                                                species, species, time_slice, area_occ, "pc_area_SUA_summary"), row.names = FALSE)
            
            #########################################################################################################################
            #########################################################################################################################
            ## If they don't exist, write the rasters for each species/threshold
            if(!file.exists(f_max_train_suit)) {
              
              ## Write the current suitability raster
              message('Writing ', species, ' current', ' max train > ', thresh) 
              writeRaster(current_suit_thresh, sprintf('%s/%s/full/%s_%s%s.tif', maxent_path,
                                                       species, species, "current_suit_above_", thresh), overwrite = TRUE) 
              
              ## Write the combined suitability raster, thresholded using the maximum training value
              message('Writing ', species, ' | 20', time_slice, ' max train > ', thresh) 
              writeRaster(combo_suit_thresh, sprintf('%s/%s/full/%s_20%s%s%s.tif', maxent_path,
                                                     species, species, time_slice, "_Max_train_sensit_above_", thresh), overwrite = TRUE)
              
              ## Write the combined suitability raster, thresholded using the percentile value
              message('Writing ', species, ' | 20', time_slice, ' 10th percentile > ', percent) 
              writeRaster(combo_suit_percent, sprintf('%s/%s/full/%s_20%s%s%s.tif', maxent_path,
                                                      species, species, time_slice, "_10_percentile_omiss_above_", percent), overwrite = TRUE)
              
              ## Write the combined future raster with > 4 GCMs above the maximum training value
              message('Writing ', species, ' | 20', time_slice, ' 10th percentile > ', percent) 
              writeRaster(combo_suit_4GCM, sprintf('%s/%s/full/%s_20%s%s%s.tif', maxent_path,
                                                   species, species, time_slice, "_4GCMs_above_", thresh), overwrite = TRUE)
              
            } else {
              
              message(species, ' 20', time_slice, ' combined suitability > ', thresh, ' skipped - already exists')   ## 
              
            }
            
            ########################################################################################################################
            ## Now create the empty panel just before plotting, and read in the occurrence and background points and original model
            
            
            ## And change the projection for the operators...........................................................................
            ## '+init=esri:54009'
            empty <- init(combo_suit_thresh, function(x) NA)
            
            occ <- readRDS(sprintf('%s/%s/occ_swd.rds', maxent_path, species)) %>%
              spTransform(CRS('+init=epsg:4326'))  
            
            bg <- readRDS(sprintf('%s/%s/bg_swd.rds', maxent_path, species)) %>%
              spTransform(CRS('+init=epsg:4326'))
            
            #m <- readRDS(sprintf('%s/%s/full/model.rds', maxent_path, species)) 
            m <- readRDS(sprintf('%s/%s/full/maxent_fitted.rds', maxent_path, species)) 
            m <- m$me_full  ## class(m);View(m)
            
            ## Use the 'levelplot' function to make a multipanel output: occurences, percentiles and thresholds
            message('Writing figure for ', species, ' | 20', time_slice, ' > ', thresh) 
            
            png(sprintf('%s/%s/full/%s_20%s_suitability_above_%s.png', maxent_path,
                        species, species, time_slice, thresh),
                11, 4, units = 'in', res = 300)
            
            ## Need an empty frame
            print(levelplot(stack(empty,                ## needs to have a different colour scale,
                                  combo_suit_percent,
                                  combo_suit_thresh), margin = FALSE,
                            
                            ## Create a colour scheme using colbrewer: 100 is to make it continuos
                            ## Also, make it a one-directional colour scheme
                            scales      = list(draw = FALSE), 
                            at = seq(0, 6, length = 100),
                            col.regions = colorRampPalette(rev(brewer.pal(11, 'Spectral'))),
                            
                            ## Give each plot a name
                            names.attr = c('Aus occurrences', 
                                           sprintf('20%s GCM 10thp train omission > %s', time_slice, percent), 
                                           sprintf('20%s GCM Max train logis > %s',      time_slice, thresh)),
                            
                            colorkey   = list(height = 0.5, width = 3), xlab = '', ylab = '',
                            main       = list(gsub('_', ' ', species), font = 4, cex = 2)) +
                    
                    ## Plot the Aus shapefile with the occurrence points for reference
                    ## Can we assign different shapefiles to different panels, rather than to them all?
                    
                    layer(sp.polygons(aus)) +
                    layer(sp.points(occ, pch = 20, cex = 0.4, 
                                    col = c('red', 'transparent', 'transparent')[panel.number()]), data = list(occ = occ)))
            
            ## Finish the device
            dev.off()
            
            ########################################################################################################################
            ## Another .png for the global records: str(LAND$long)
            if(!file.exists(sprintf('%s/%s/full/%s_%s.png', maxent_path, species, species, "global_records"))) {
              
              png(sprintf('%s/%s/full/%s_%s.png', maxent_path,
                          species, species, "global_records"),
                  16180, 10000, units = 'px', res = 600)
              
              ## How do we locate bad records in the dataset after spotting them?
              plot(LAND, 
                   lwd = 0.5, asp = 1, axes = TRUE, cex.axis = 3.5,
                   col = 'darkolivegreen3', bg = 'lightblue', cex.lab = 3)
              
              points(occ, pch = ".", cex = 3.5, col = "red", cex.lab = 3, cex.main = 4, cex.axis = 2, 
                     main = paste0("Global occurrences for ", species), 
                     xlab = "", ylab = "", asp = 1)
              
              ## Title 
              title(paste0("Global points for ", species),
                    cex.main = 4,   font.main = 4, col.main = "blue")
              
              ## Finsh the device
              dev.off()
              
            } else {
              
              message("Global records maps exists for ", species)
              
            }
            
            ########################################################################################################################
            ## Another PNG for the background points....
            if(!file.exists(sprintf('%s/%s/full/%s_%s.png', maxent_path, species, species, "background_records"))) {
              
              png(sprintf('%s/%s/full/%s_%s.png', maxent_path,
                          species, species, "background_records"),
                  16180, 10000, units = 'px', res = 600)
              
              ## How do we locate bad records in the dataset after spotting them?
              plot(LAND,  
                   lwd = 0.5, asp = 1, axes = TRUE, cex.axis = 3.5,
                   col = 'darkolivegreen3', bg = 'lightblue', cex.lab = 3)
              
              points(bg, pch = ".", cex = 1.6, col = "blue", cex.lab = 3, cex.main = 4, cex.axis = 2, 
                     main = paste0("Global occurrences for ", species), 
                     xlab = "", ylab = "", asp = 1)
              
              ## Title 
              title(paste0("Bacground points for ", species),
                    cex.main = 4,   font.main = 4, col.main = "blue")
              
              ## Finish the device
              dev.off() 
              
            } else {
              
              message("Background records maps exists for ", species)
              
            }
            
            ########################################################################################################################
            ## Plot the models: can two plots be combined into one?
            ## Make these unique names, and they can be searched in windows. Otherwise, we can just click into each subfolder. 
            ## To sort, names would need to be: spp + unique_extension
            if(!file.exists(sprintf('%s/%s/full/%s_current.png', maxent_path, species, species, "variable_contribution"))) {
              
              png(sprintf('%s/%s/full/%s_current.png', maxent_path,
                          species, species, "variable_contribution"),
                  3236, 2000, units = 'px', res = 300)
              
              ## Set the margins
              # par(mgp      = c(10, 4, 0), 
              #     oma      = c(1.5, 1.5, 1.5, 1.5),
              #     font.lab = 2)
              
              plot(m, col = "blue", pch = 19, cex.lab = 2, cex.axis = 5, cex.main = 2, 
                   main   = paste0("Variables for ", species), 
                   xlab   = "Maxent contribution (%)")
              
              ## Finish the device
              dev.off()
              
              ## Plot the response curves too
              png(sprintf('%s/%s/full/%s_%s.png', maxent_path,
                          species, species, "response_curves"),
                  3236, 2000, units = 'px', res = 300)
              
              ## Add detail to the response plot
              response(m, pch = 19, cex.lab = 2, cex.axis = 1.5, lwd = 2) 
              #ylab   = "",
              #main   = paste0(species, " responses"))
              
              ## Finish the device
              dev.off()
              
            } else {
              
              message("Variable contribution plot exists for ", species)
              
            }
            
            ## current.list = as.data.frame(env.grids.current)
            ## test = similarity(occ, env.grids.current, full = FALSE)
            
          } else {
            
            message(species, ' 20', time_slice, ' combined suitability > ', thresh, ' skipped - already exists')   ## 
            
          }
          
        }
        
      }
      
    })
    
  })
  
}





#########################################################################################################################
## ANOMALY FUNCTIONS
#########################################################################################################################


#########################################################################################################################
## For each scenario, calcualte the mean annual temperature and annual rainfall anomaly
calculate.anomaly = function(scen_list, time_slice, climate_path) {
  
  #########################################################################################################################
  ## Create current rasters :: can the raster creation happen just once?
  ## And read in the current data. Can the specific "current" be removed withouth makin a seconf path argument?
  message('current rasters / 10')
  
  env.grids.current       = stack(file.path(sprintf('%s/current/bio_%02d.tif', climate_path, 1:19)))
  env.grids.current[[1]]  = env.grids.current[[1]]/10
  current.bio1            = env.grids.current[[1]]
  current.bio12           = env.grids.current[[12]]
  
  ## First, run a loop over each scenario :: 
  lapply(scen_list, function(x) {
    
    ## Assign the scenario name (to use later in the plot)
    scen_name = eval(parse(text = sprintf('gcms.%s$GCM[gcms.%s$id == x]', time_slice, time_slice)))
    
    #########################################################################################################################
    ## Create a raster stack for each 2050 GCM - also an empty raster for the final plot
    s <- stack(sprintf('%s/20%s/%s/%s%s.tif', climate_path, time_slice, x, x, 1:19))
    
    ########################################################################################################################
    ## Create future rasters
    message('20', time_slice, ' rasters / 10')
    
    s[[1]]       = s[[1]]/10
    future.bio1  = s[[1]]
    future.bio12 = s[[12]]
    
    ## Now subtract the current from the future...
    temp.anomaly = future.bio1  - current.bio1
    rain.anomaly = future.bio12 - current.bio12
    
    plot(temp.anomaly, main = paste0("BIO1  anomaly ", x))
    plot(rain.anomaly, main = paste0("BIO12 anomaly ", x))
    
    ########################################################################################################################
    ## Write the rasters for each species/threshold
    message('Writing BIO1/12 anomaly rasters for 20', time_slice, ' ', x)
    
    if(!file.exists(sprintf('%s/anomalies/%s_%s.tif', climate_path, "BIO1_anomaly",  x))) {
      
      ## Temp anomalies
      writeRaster(temp.anomaly, sprintf('%s/anomalies/%s_%s.tif', climate_path, "BIO1_anomaly",  x), overwrite = TRUE)
      
    } else {
      
      message('Skipped ', "BIO1 anomaly for ", x, ', already exists')
      
    }
    
    if(!file.exists(sprintf('%s/anomalies/%s_%s.tif', climate_path, "BIO12_anomaly",  x))) {
      
      ## Rain anomalies
      writeRaster(rain.anomaly, sprintf('%s/anomalies/%s_%s.tif', climate_path, "BIO12_anomaly", x), overwrite = TRUE)
      
    } else {
      
      message('Skipped ', "BIO12 anomaly for ", x, ', already exists')
      
    }
    
  })
  
} 





#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################