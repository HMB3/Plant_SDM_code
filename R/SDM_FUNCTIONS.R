#########################################################################################################################
#############################################  FUNCTIONS FOR HORT AUS SDMS ############################################## 
#########################################################################################################################


#########################################################################################################################
## MAXENT FUNCTIONS
#########################################################################################################################


#########################################################################################################################
## Big function for getting background points...and fitting maxent? 
fit_maxent2 <- function(occ, bg, predictors, name, outdir, template, 
                        max_bg_size = 100000, 
                        shapefiles = TRUE, 
                        features, replicates, 
                        responsecurves = TRUE, 
                        rep_args, full_args) {
  
  ########################################################################
  ## predictors: the s_ff object containing predictors to use in the model
  require(ff)
  require(rgeos)
  require(sp)
  require(raster)
  require(rJava)
  require(dismo)
  require(things)
  
  ## First, stop if the outdir file exists, 
  if(!file.exists(outdir)) stop('outdir does not exist :(', call. = FALSE)
  outdir_sp <- file.path(outdir, gsub(' ', '_', name))
  
  ## What? 
  if(!file.exists(outdir_sp)) dir.create(outdir_sp)
  features <- unlist(strsplit(features, ''))
  
  ## What?
  if(length(setdiff(features, c('l', 'p', 'q', 'h', 't'))) > 1)
    stop("features must be a vector of one or more of ',
         'l', 'p', 'q', 'h', and 't'.")
  
  ## aggregate
  b <- aggregate(gBuffer(occ, width = 200000, byid = TRUE))
  
  #####################################################################
  ## Get unique cell numbers for species occurrences
  cells <- cellFromXY(template, occ)
  
  ## Clean out duplicates and NAs (including points outside extent of predictor data)
  not_dupes <- which(!duplicated(cells) & !is.na(cells)) 
  occ <- occ[not_dupes, ]
  cells <- cells[not_dupes]
  message(nrow(occ), ' occurrence records (unique cells).')
  
  ## skip species with < 20 records
  if(length(occ) < min.spp) {
    
    warning('Fewer occurrence records than the number of cross-validation ',
            'replicates for species ', name, 
            '. Model not fit for this species.')
    
    return(NULL)
    
  }
  
  #####################################################################
  ## Now subset bg to the buffer polygon
  system.time(o <- over(bg, b))
  bg <- bg[which(!is.na(o)), ]
  bg_cells <- cellFromXY(template, bg)
  
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
  
  bg_cells <- cellFromXY(template, bg) # can probably move this into if {}
  
  #####################################################################
  ## Save objects for future reference
  if(shapefiles) {
    
    suppressWarnings({
      
      writeOGR(SpatialPolygonsDataFrame(b, data.frame(ID=seq_len(length(b)))), 
               outdir_sp, 'bg_buffer', 'ESRI Shapefile', overwrite_layer = TRUE)  
      writeOGR(bg,  outdir_sp, 'bg',   'ESRI Shapefile', overwrite_layer = TRUE)  
      writeOGR(occ, outdir_sp, 'occ',  'ESRI Shapefile', overwrite_layer = TRUE) 
      
    })
    
  }
  
  saveRDS(bg, file.path(outdir_sp, 'bg.rds'))
  saveRDS(occ, file.path(outdir_sp, 'occ.rds'))
  
  #####################################################################
  ## Sample predictors at occurrence and background points
  swd_occ <- predictors[cells, ]
  swd_occ <- SpatialPointsDataFrame(
    coordinates(occ), as.data.frame(swd_occ), proj4string=crs(occ))
  saveRDS(swd_occ, file.path(outdir_sp, 'occ_swd.rds'))
  
  swd_bg <- predictors[bg_cells, ]
  swd_bg <- SpatialPointsDataFrame(coordinates(bg), as.data.frame(swd_bg), 
                                   proj4string=crs(bg))
  saveRDS(swd_bg, file.path(outdir_sp, 'bg_swd.rds'))
  
  ## save shapefiles?
  if(shapefiles) {
    
    writeOGR(swd_occ, outdir_sp, 'occ_swd', 'ESRI Shapefile', overwrite_layer = TRUE)
    writeOGR(swd_bg, outdir_sp,   'bg_swd', 'ESRI Shapefile', overwrite_layer = TRUE)
    
  }
  
  #####################################################################
  ## Combine occ and bg SWD data
  swd <- as.data.frame(rbind(swd_occ@data, swd_bg@data))
  saveRDS(swd, file.path(outdir_sp, 'swd.rds'))
  pa <- rep(1:0, c(nrow(swd_occ), nrow(swd_bg)))
  
  ## Now fit model?
  off <- setdiff(c('l', 'p', 'q', 't', 'h'), features)
  
  ## 
  if(length(off) > 0) {
    
    off <- c(l = 'linear = FALSE',    p = 'product = FALSE', q = 'quadratic = FALSE',
             t = 'threshold = FALSE', h = 'hinge = FALSE')[off]
  }
  
  off <- unname(off)
  
  if(replicates > 1) {
    
    if(missing(rep_args)) rep_args <- NULL
    
    ## What?
    me_xval <- maxent(swd, pa, path = file.path(outdir_sp, 'xval'), 
                      args=c(paste0('replicates=', replicates),
                             'responsecurves=TRUE',
                             off, paste(names(rep_args), rep_args, sep='=')))
    
  }
  
  ## 
  if(missing(full_args)) full_args <- NULL
  me_full <- maxent(swd, pa, path = file.path(outdir_sp, 'full'), 
                    args = c(off, paste(names(full_args), full_args, sep= '='),
                           'responsecurves = TRUE'))
  
  #####################################################################
  ## Save fitted model object, and the model-fitting data.
  saveRDS(list(me_xval = me_xval, me_full = me_full, swd = swd, pa = pa), 
          file.path(outdir_sp, 'maxent_fitted.rds'))
  
  ## Don't return anything?
  return(invisible(NULL))
  
}




#########################################################################################################################
## 



