#########################################################################################################################
#############################################  FUNCTIONS FOR HORT AUS SDMS ############################################## 
#########################################################################################################################


#########################################################################################################################
## MAXENT FUNCTIONS
#########################################################################################################################

## Need an example of actually runnin the code


#########################################################################################################################
## GET BACKGROUND POINTS AND FIT MAXENT 
#########################################################################################################################


FIT_MAXENT <- function(occ, 
                       bg, # A Spatial points data frame (SPDF) of candidate background points
                       sdm.predictors, 
                       # sdm.predictors is a vector of column names for sdm.predictors
                       # that you want to include
                       name, 
                       outdir, 
                       template.raster, 
                       # template.raster is an empty raster with extent, res and projection
                       # of final output rasters. It is used to reduce
                       # occurrences to a single point per cell.
                       min_n                   = 20,
                       # min_n is the minimum number of records (unique cells)
                       # required for a model to be fit
                       max_bg_size             = 100000, 
                       background_buffer_width = 200000,
                       shapefiles              = TRUE, 
                       features, 
                       replicates, # number of cross-validation replicates
                       responsecurves          = TRUE, 
                       rep_args, 
                       full_args) {
  
  ########################################################################
  ## sdm.predictors: the s_ff object containing sdm.predictors to use in the model
  
  ## First, stop if the outdir file exists, 
  if(!file.exists(outdir)) stop('outdir does not exist :(', call. = FALSE)
  outdir_sp <- file.path(outdir, gsub(' ', '_', name))
  
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
    
    warning('Fewer occurrence records than the number of cross-validation ',
            'replicates for species ', name, 
            '. Model not fit for this species.')
    
    return(NULL)
    
  }
  
  #####################################################################
  ## Now subset bg to the buffer polygon
  system.time(o <- over(bg, b))
  bg <- bg[which(!is.na(o)), ]
  bg_cells <- cellFromXY(template.raster, bg)
  
  ## Clean out duplicates and NAs (including points outside extent of predictor data)
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
  
  bg_cells <- cellFromXY(template.raster, bg) # can probably move this into if {}
  
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
  ## Sample sdm.predictors at occurrence and background points
  swd_occ <- occ[, sdm.predictors]
  saveRDS(swd_occ, file.path(outdir_sp, 'occ_swd.rds'))
  
  swd_bg <- bg[, sdm.predictors]
  saveRDS(swd_bg, file.path(outdir_sp, 'bg_swd.rds'))
  
  ## Save shapefiles of the...
  if(shapefiles) {
    
    writeOGR(swd_occ, outdir_sp,  'occ_swd', 'ESRI Shapefile', overwrite_layer = TRUE)
    writeOGR(swd_bg,  outdir_sp,   'bg_swd', 'ESRI Shapefile', overwrite_layer = TRUE)
    
  }
  
  #####################################################################
  ## Combine occurrence and bg SWD data
  swd <- as.data.frame(rbind(swd_occ@data, swd_bg@data))
  saveRDS(swd, file.path(outdir_sp, 'swd.rds'))
  pa <- rep(1:0, c(nrow(swd_occ), nrow(swd_bg)))
  
  ## Now fit model?
  off <- setdiff(c('l', 'p', 'q', 't', 'h'), features)
  
  ## 
  if(length(off) > 0) {
    
    off <- c(l = 'linear=false',    p = 'product=false', q = 'quadratic=false',
             t = 'threshold=false', h = 'hinge=false')[off]
    
  }
  
  off <- unname(off)
  
  if(replicates > 1) {
    
    if(missing(rep_args)) rep_args <- NULL
    
    ## This runs the MAXENT. This is where the argument errors are coming in
    # Error in .local(x, p, ...) : args not understood:
    #   replicates = 5, responsecurves = TRUE, threshold = FALSE, hinge = FALSE
    me_xval <- maxent(swd, pa, path = file.path(outdir_sp, 'xval'), 
                      args = c(paste0('replicates=', replicates),
                               'responsecurves=true', 
                               'outputformat=logistic',
                               off, paste(names(rep_args), rep_args, sep = '=')))
    
  }
  
  ## And this is the same, but with a different argument 
  if(missing(full_args)) full_args <- NULL
  me_full <- maxent(swd, pa, path = file.path(outdir_sp, 'full'), 
                    args = c(off, paste(names(full_args), full_args, sep = '='),
                             'responsecurves=true',
                             'outputformat=logistic'))
  
  #####################################################################
  ## Save fitted model object, and the model-fitting data.
  # if (file.exists (file)) {
  #   
  #   print (paste ("file exists for genera", gen.n, "skipping"))
  #   next
  
  
  if(replicates > 1) {
    
    saveRDS(list(me_xval = me_xval, me_full = me_full, swd = swd, pa = pa), 
            file.path(outdir_sp, 'maxent_fitted.rds'))
    
  } else {
    
    saveRDS(list(me_xval = NA, me_full = me_full, swd = swd, pa = pa), 
            file.path(outdir_sp, 'maxent_fitted.rds'))
    
  }
  
}





#########################################################################################################################
## GET BACKGROUND POINTS AND FIT MAXENT FOR ALL VARIABLES, USE MODEL SELECTION
#########################################################################################################################


FIT_MAXENT_SELECT <- function(occ, 
                              bg, # A Spatial points data frame (SPDF) of candidate background points
                              sdm.predictors.all, 
                              # sdm.predictors.all is a vector of column names for sdm.predictors.all
                              # that you want to include
                              name, 
                              outdir, 
                              template.raster, 
                              # template.raster is an empty raster with extent, res and projection
                              # of final output rasters. It is used to reduce
                              # occurrences to a single point per cell.
                              min_n                   = 20,
                              # min_n is the minimum number of records (unique cells)
                              # required for a model to be fit
                              max_bg_size             = 100000, 
                              background_buffer_width = 200000,
                              shapefiles              = TRUE, 
                              features, 
                              replicates, # number of cross-validation replicates
                              responsecurves          = TRUE, 
                              rep_args, 
                              full_args) {
  
  
  ## Where does the model selection code come in here?
  
  
  ########################################################################
  ## sdm.predictors.all: the s_ff object containing sdm.predictors.all to use in the model
  
  ## First, stop if the outdir file exists, 
  if(!file.exists(outdir)) stop('outdir does not exist :(', call. = FALSE)
  outdir_sp <- file.path(outdir, gsub(' ', '_', name))
  
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
    
    warning('Fewer occurrence records than the number of cross-validation ',
            'replicates for species ', name, 
            '. Model not fit for this species.')
    
    return(NULL)
    
  }
  
  #####################################################################
  ## Now subset bg to the buffer polygon
  system.time(o <- over(bg, b))
  bg <- bg[which(!is.na(o)), ]
  bg_cells <- cellFromXY(template.raster, bg)
  
  ## Clean out duplicates and NAs (including points outside extent of predictor data)
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
  
  bg_cells <- cellFromXY(template.raster, bg) # can probably move this into if {}
  
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
  ## Sample sdm.predictors.all at occurrence and background points
  swd_occ <- occ[, sdm.predictors.all]
  swd_bg <- bg[, sdm.predictors.all]

  #####################################################################
  ## Here is where we need to reduce the predictors to a candidate set
  sdm.predictors.all = rmaxent::simplify(swd_occ, swd_bg, path, ## Don't think we want this one?
                                         species_column  = "searchTaxon", 
                                         response_curves = FALSE,
                                         logistic_format = TRUE, 
                                         type            = "PI", 
                                         cor_thr         = 0.7, 
                                         pct_thr         = 5, 
                                         k_thr           = 4,
                                         quiet           = FALSE)
  
  
  #####################################################################
  ## Recreate occ and bg with new predictors
  swd_occ <- occ[, sdm.predictors.all]
  swd_bg <- bg[, sdm.predictors.all]
  
  ## Then save them...
  saveRDS(swd_occ, file.path(outdir_sp, 'occ_swd.rds'))
  saveRDS(swd_bg, file.path(outdir_sp, 'bg_swd.rds'))
  
  
  ## Save shapefiles of the...
  if(shapefiles) {
    
    writeOGR(swd_occ, outdir_sp,  'occ_swd', 'ESRI Shapefile', overwrite_layer = TRUE)
    writeOGR(swd_bg,  outdir_sp,   'bg_swd', 'ESRI Shapefile', overwrite_layer = TRUE)
    
  }
  
  #####################################################################
  ## Combine occurrence and bg SWD data
  swd <- as.data.frame(rbind(swd_occ@data, swd_bg@data))
  saveRDS(swd, file.path(outdir_sp, 'swd.rds'))
  pa <- rep(1:0, c(nrow(swd_occ), nrow(swd_bg)))
  
  ## Now fit model?
  off <- setdiff(c('l', 'p', 'q', 't', 'h'), features)
  
  ## 
  if(length(off) > 0) {
    
    off <- c(l = 'linear=false',    p = 'product=false', q = 'quadratic=false',
             t = 'threshold=false', h = 'hinge=false')[off]
    
  }
  
  off <- unname(off)
  
  if(replicates > 1) {
    
    if(missing(rep_args)) rep_args <- NULL
    
    ## This runs the MAXENT. This is where the argument errors are coming in
    # Error in .local(x, p, ...) : args not understood:
    #   replicates = 5, responsecurves = TRUE, threshold = FALSE, hinge = FALSE
    me_xval <- maxent(swd, pa, path = file.path(outdir_sp, 'xval'), 
                      args = c(paste0('replicates=', replicates),
                               'responsecurves=true', 
                               'outputformat=logistic',
                               off, paste(names(rep_args), rep_args, sep = '=')))
    
  }
  
  ## And this is the same, but with a different argument 
  if(missing(full_args)) full_args <- NULL
  me_full <- maxent(swd, pa, path = file.path(outdir_sp, 'full'), 
                    args = c(off, paste(names(full_args), full_args, sep = '='),
                             'responsecurves=true',
                             'outputformat=logistic'))
  
  #####################################################################
  ## Save fitted model object, and the model-fitting data.
  # if (file.exists (file)) {
  #   
  #   print (paste ("file exists for genera", gen.n, "skipping"))
  #   next
  
  
  if(replicates > 1) {
    
    saveRDS(list(me_xval = me_xval, me_full = me_full, swd = swd, pa = pa), 
            file.path(outdir_sp, 'maxent_fitted.rds'))
    
  } else {
    
    saveRDS(list(me_xval = NA, me_full = me_full, swd = swd, pa = pa), 
            file.path(outdir_sp, 'maxent_fitted.rds'))
    
  }
  
}





#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################