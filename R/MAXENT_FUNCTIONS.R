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
# outdir                  = 'output/maxent/SET_VAR_KOPPEN'
# sdm.predictors          = sdm.select
# env.grids               = env.grids.current
# template.raster         = template.raster
# template.cells          = template.cells
# min_n                   = 20               ## This should be higher...
# max_bg_size             = 100000           ## need a min bg size?
# background_buffer_width = 200000
# shapefiles              = TRUE
# features                = 'lpq'
# replicates              = 5
# responsecurves          = TRUE
# Koppen                  = Koppen_1975


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
# FIT_MAXENT_RAND_BG <- function(occ,
#                                name,
#                                outdir,
#                                sdm.predictors,
#                                env.grids,
#                                template.raster,
#                                template.cells,
#                                # template.raster is an empty raster with extent, res and projection
#                                # of final output rasters. It is used to reduce
#                                # occurrences to a single point per cell.
#                                min_n,
#                                # min_n is the minimum number of records (unique cells)
#                                # required for a model to be fit
#                                max_bg_size,
#                                background_buffer_width,
#                                shapefiles,
#                                features,
#                                replicates, # number of cross-validation replicates
#                                responsecurves,
#                                rep_args,
#                                full_args) {
#   
#   ## There is no minimum number of background records.................................................................
#   
#   ########################################################################
#   ## First, stop if the outdir file exists,
#   if(!file.exists(outdir)) stop('outdir does not exist :(', call. = FALSE)
#   outdir_sp <- file.path(outdir, gsub(' ', '_', name))
#   
#   ## Also Skip species which have already been modelled: again, can this be made into an argument?
#   if(file.exists(outdir_sp)) {
#     
#     print (paste ('Skip', name, ', Model results exist for this species'))
#     
#   } else {
#     
#     ## If the file doesn't exist, split out the features
#     if(!file.exists(outdir_sp)) dir.create(outdir_sp)
#     features <- unlist(strsplit(features, ''))
#     
#     ## Make sure user features are allowed: don't run the model if the
#     ## features have been incorrectly specified in the main argument
#     ## l: linear
#     ## p: product
#     ## q: quadratic
#     ## h: hinge
#     ## t:
#     if(length(setdiff(features, c('l', 'p', 'q', 'h', 't'))) > 1)
#       stop("features must be a vector of one or more of ',
#            'l', 'p', 'q', 'h', and 't'.")
#     
#     ## Aggregate
#     b <- aggregate(gBuffer(occ, width = background_buffer_width, byid = TRUE))
#     
#     #####################################################################
#     ## Get unique cell numbers for species occurrences
#     cells <- cellFromXY(template.raster, occ)
#     
#     ## Clean out duplicate cells and NAs (including points outside extent of predictor data)
#     ## Note this will get rid of a lot of duplicate records not filtered out by GBIF columns, etc.
#     not_dupes <- which(!duplicated(cells) & !is.na(cells))
#     occ       <- occ[not_dupes, ]
#     cells     <- cells[not_dupes]
#     
#     ##
#     message(nrow(occ), ' occurrence records (unique cells).')
#     
#     ## Skip species that have less than a minimum number of records: eg 20 species
#     if(nrow(occ) < min_n) {
#       
#       print (paste ('Fewer occurrence records than the number of cross-validation ',
#                     'replicates for species ', name,
#                     ' Model not fit for this species'))
#       
#     } else {
#       
#       #####################################################################
#       ## Now subset the background records to the raster cells with data
#       f <- tempfile()
#       writeOGR(SpatialPolygonsDataFrame(b, data.frame(id = seq_along(b))), 
#                tempdir(), basename(f), 'ESRI Shapefile')
#       
#       ## Create a raster of cells with data
#       message('creating background cells for ', name)
#       bg_cells <- gdal_rasterize(extension(f, 'shp'), 
#                                  tempfile(fileext = '.tif'),
#                                  te = c(bbox(template.raster)),
#                                  tr = res(template.raster),
#                                  burn = 1, init = 0, ot = 'Byte', 
#                                  output_Raster = TRUE) %>% 
#         
#         .[[1]] %>% # gdal_rasterize returns a RasterBrick but Which needs a RasterLayer, so we grab the first (only) layer
#         Which(. == 1, cells = TRUE) %>% 
#         intersect(template.cells)
#       
#       bg <- xyFromCell(template.raster, bg_cells)
#       
#       ## Reduce background sample if it's larger than max_bg_size
#       ## There is no minimum number of background records.................................................................
#       
#       ## What if the species has very few records surrounding its record?
#       if (nrow(bg) > max_bg_size) {
#         
#         message(nrow(bg), ' target species background records, reduced to random ',
#                 max_bg_size, '.')
#         
#         bg <- bg[sample(nrow(bg), max_bg_size), ]  ##  
#         
#       } else {
#         
#         message(nrow(bg), ' target species background records.')
#         
#       }
#       
#       ## Extract environmental values at the background points :: this was not working for me...
#       ## Is this due to a projection problem? points are mollweide, 
#       ## projection(env.grids)
#       bg <- SpatialPointsDataFrame(coords      = as.data.frame(bg), 
#                                    data        = as.data.frame(bg),
#                                    proj4string = CRS('+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'))
#       
#       ## Is the cell = bg_cells necessary? This bit takes too long....
#       ## For 160 spp this would take about 4 days
#       message('extracting cells for ', name)
#       bg <- extract(env.grids, bg) ## sp = TRUE)   ## In .doExtract(x, i, drop = drop) : some indices are invalid (NA returned)
#       ## bg <- SpatialPointsDataFrame(SpatialPoints(bg), data.frame(cell = bg_cells)) 
#       
#       #####################################################################
#       ## Save shapefiles for future reference
#       if(shapefiles) {
#         
#         suppressWarnings({
#           
#           writeOGR(SpatialPolygonsDataFrame(b, data.frame(ID = seq_len(length(b)))),
#                    outdir_sp, 'bg_buffer', 'ESRI Shapefile', overwrite_layer = TRUE)
#           writeOGR(bg,  outdir_sp, 'bg',   'ESRI Shapefile', overwrite_layer = TRUE)
#           writeOGR(occ, outdir_sp, 'occ',  'ESRI Shapefile', overwrite_layer = TRUE)
#           
#         })
#         
#       }
#       
#       ## Save the background and occurrence points as .rds objects
#       saveRDS(bg,  file.path(outdir_sp, 'bg.rds'))
#       saveRDS(occ, file.path(outdir_sp, 'occ.rds'))
#       
#       #####################################################################
#       ## sdm.predictors: the s_ff object containing sdm.predictors to use 
#       ## in the model: sample sdm.predictors at occurrence and background points
#       swd_occ <- occ[, sdm.predictors]
#       saveRDS(swd_occ, file.path(outdir_sp, 'occ_swd.rds'))
#       
#       swd_bg <- bg[, sdm.predictors]
#       saveRDS(swd_bg, file.path(outdir_sp, 'bg_swd.rds'))
#       
#       ## Save shapefiles of the occurrence and background points
#       if(shapefiles) {
#         
#         writeOGR(swd_occ, outdir_sp,  'occ_swd', 'ESRI Shapefile', overwrite_layer = TRUE)
#         writeOGR(swd_bg,  outdir_sp,  'bg_swd',  'ESRI Shapefile', overwrite_layer = TRUE)
#         
#       }
#       
#       #####################################################################
#       ## Combine occurrence and background data
#       swd <- as.data.frame(rbind(swd_occ@data, swd_bg@data))
#       saveRDS(swd, file.path(outdir_sp, 'swd.rds'))
#       pa <- rep(1:0, c(nrow(swd_occ), nrow(swd_bg)))
#       
#       ## Now check the features arguments are correct
#       off <- setdiff(c('l', 'p', 'q', 't', 'h'), features)
#       
#       ## 
#       if(length(off) > 0) {
#         
#         off <- c(l = 'linear=false',    p = 'product=false', q = 'quadratic=false',
#                  t = 'threshold=false', h = 'hinge=false')[off]
#         
#       }
#       
#       off <- unname(off)
#       
#       if(replicates > 1) {
#         
#         if(missing(rep_args)) rep_args <- NULL
#         
#         ## Run MAXENT for X cross validation data splits of swd : so 5 replicaes, 0-4
#         ## EG xval = cross validation : "OUT_DIR\Acacia_boormanii\xval\maxent_0.html"
#         me_xval <- maxent(swd, pa, path = file.path(outdir_sp, 'xval'),
#                           args = c(paste0('replicates=', replicates),
#                                    #paste0('biasfile=', dens),
#                                    'responsecurves=true',
#                                    'outputformat=logistic', # "biasfile = dens.ras"
#                                    off, paste(names(rep_args), rep_args, sep = '=')))
#         
#       }
#       
#       ## Runs the full maxent model - presumably using all the data in swd
#       if(missing(full_args)) full_args <- NULL
#       me_full <- maxent(swd, pa, path = file.path(outdir_sp, 'full'),
#                         args = c(off, paste(names(full_args), full_args, sep = '='),
#                                  'responsecurves=true',
#                                  'outputformat=logistic'))
#       
#       ## Save the full model?
#       saveRDS(list(me_xval = me_xval, me_full = me_full, swd = swd, pa = pa), 
#               file.path(outdir_sp, 'full', 'maxent_fitted.rds'))
#       
#       #####################################################################
#       ## Save the chart corrleation file too for the variable set
#       save_name = gsub(' ', '_', name)
#       png(sprintf('%s/%s/full/%s_%s.png', outdir,
#                   save_name, save_name, "predictor_correlation"),
#           3236, 2000, units = 'px', res = 300)
#       
#       ## set margins
#       par(mar   = c(3, 3, 5, 3),  ## b, l, t, r
#           #mgp   = c(9.8, 2.5, 0),
#           oma   = c(1.5, 1.5, 1.5, 1.5))
#       
#       ## Add detail to the response plot
#       chart.Correlation(swd_occ@data,
#                         histogram = TRUE, pch = 19) 
#       #cex.lab = 2, cex.axis = 1.5,
#       #main = paste0("Predictor corrleation matrix for ", spp))
#       
#       ## Finish the device
#       dev.off()
#       
#       ########################################################################################################################
#       ## Another .png for the global records: str(LAND$long)
#       LAND       = readRDS("F:/green_cities_sdm/data/base/CONTEXTUAL/LAND_world.rds")
#       occ_land   = occ %>% 
#         spTransform(CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
#       bg_land    = bg %>% 
#         spTransform(CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
#       
#       png(sprintf('%s/%s/full/%s_%s.png', maxent_path,
#                   save_name, save_name, "global_records"),
#           16180, 10000, units = 'px', res = 600)
#       
#       ## How do we locate bad records in the dataset after spotting them?
#       plot(LAND, 
#            lwd = 0.5, asp = 1, axes = TRUE, cex.axis = 3.5,
#            col = 'darkolivegreen3', bg = 'lightblue', cex.lab = 3)
#       
#       points(occ_land, pch = ".", cex = 3.5, col = "red", cex.lab = 3, cex.main = 4, cex.axis = 2, 
#              main = paste0("Global occurrences for ", name), 
#              xlab = "", ylab = "", asp = 1)
#       
#       ## Title 
#       title(paste0("Global points for ", name),
#             cex.main = 4,   font.main = 4, col.main = "blue")
#       
#       ## Finsh the device
#       dev.off()
#       
#       ########################################################################################################################
#       ## Another PNG for the background points....
#       #if(!file.exists(sprintf('%s/%s/full/%s_%s.png', maxent_path, name, name, "background_records"))) {
#       png(sprintf('%s/%s/full/%s_%s.png', outdir,
#                   name, name, "background_records"),
#           16180, 10000, units = 'px', res = 600)
#       
#       ## How do we locate bad records in the dataset after spotting them?
#       plot(LAND,  
#            lwd = 0.5, asp = 1, axes = TRUE, cex.axis = 3.5,
#            col = 'darkolivegreen3', bg = 'lightblue', cex.lab = 3)
#       
#       points(bg, pch = ".", cex = 1.6, col = "blue", cex.lab = 3, cex.main = 4, cex.axis = 2, 
#              main = paste0("Bacground points for ", name), 
#              xlab = "", ylab = "", asp = 1)
#       
#       ## Title 
#       title(paste0("Bacground points for ", name),
#             cex.main = 4,   font.main = 4, col.main = "blue")
#       
#       ## Finish the device
#       dev.off() 
#       
#       # } else {
#       #   
#       #   message("Background records maps exists for ", name)
#       #   
#       # }
#       
#       ########################################################################################################################
#       ## Plot the models: can two plots be combined into one?
#       ## Make these unique names, and they can be searched in windows. Otherwise, we can just click into each subfolder. 
#       ## To sort, names would need to be: spp + unique_extension
#       m <- readRDS(sprintf('%s/%s/full/maxent_fitted.rds', outdir, save_name)) 
#       m <- m$me_full
#       
#       png(sprintf('%s/%s/full/%s_current.png', outdir,
#                   name, name, "variable_contribution"),
#           3236, 2000, units = 'px', res = 300)
#       
#       ## Set the margins
#       # par(mgp      = c(10, 4, 0), 
#       #     oma      = c(1.5, 1.5, 1.5, 1.5),
#       #     font.lab = 2)
#       
#       plot(m, col = "blue", pch = 19, cex.lab = 2, cex.axis = 5, cex.main = 2, 
#            main   = paste0("Variables for ", name), 
#            xlab   = "Maxent contribution (%)")
#       
#       ## Finish the device
#       dev.off()
#       
#       ## Plot the response curves too
#       png(sprintf('%s/%s/full/%s_%s.png', outdir,
#                   name, name, "response_curves"),
#           3236, 2000, units = 'px', res = 300)
#       
#       ## Add detail to the response plot
#       response(m, pch = 19, cex.lab = 2, cex.axis = 1.5, lwd = 2) 
#       
#       ## Finish the device
#       dev.off()
#       
#       #####################################################################
#       ## Save fitted model object, and the model-fitting data.
#       #       if(replicates > 1) {
#       # 
#       #         saveRDS(list(me_xval = me_xval, me_full = me_full, swd = swd, pa = pa),
#       #                 file.path(outdir_sp, 'maxent_fitted.rds'))
#       # 
#       #       } else {
#       # 
#       #         saveRDS(list(me_xval = NA, me_full = me_full, swd = swd, pa = pa),
#       #                 file.path(outdir_sp, 'maxent_fitted.rds'))
#       # 
#       # }
#       
#     }
#     
#   }
#   
# }





#########################################################################################################################
## GET TARGETTED BACKGROUND POINTS AND THEN FIT MAXENT WITH STANDARD VARIABLES
#########################################################################################################################


##
fit_maxent_targ_bg_kopp <- function(occ,
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
                                    background_buffer_width, # ignored if background_method='random'
                                    #background_method, # 'random' or 'targetgroup'
                                    Koppen,
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
  
  if(!missing('Koppen')) {
    if(!is(Koppen, 'RasterLayer'))
      stop('Koppen must be a RasterLayer, and should be in the same coordinate system as template.raster')  
  }
  
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
    
    ## Subset the background records to the 200km buffered polygon
    message(name, ' creating background cells')
    system.time(o <- over(bg, b))
    bg <- bg[which(!is.na(o)), ]
    bg_cells <- cellFromXY(template.raster, bg)
    
    ## Clean out duplicates and NAs (including points outside extent of predictor data)
    bg_not_dupes <- which(!duplicated(bg_cells) & !is.na(bg_cells))
    bg <- bg[bg_not_dupes, ]
    bg_cells <- bg_cells[bg_not_dupes]
    
    ## Find which of these cells fall within the Koppen-Geiger zones that the species occupies
    if(!missing('Koppen')) {
      
      ## Crop the Kopppen raster to the extent of the occurrences, and snap it
      message(name, ' intersecting background cells with Koppen zones')
      Koppen_crop <- crop(Koppen, occ, snap = 'out')
      
      ## Only extract and match those cells that overlap between koppen_cropp, occ and bg 
      zones               <- raster::extract(Koppen_crop, occ)
      cells_in_zones_crop <- Which(Koppen_crop %in% zones, cells = TRUE)
      cells_in_zones      <- cellFromXY(Koppen, xyFromCell(Koppen_crop, cells_in_zones_crop))
      bg_cells            <- intersect(bg_cells, cells_in_zones)
      i                   <- cellFromXY(template.raster, bg)
      bg                  <- bg[which(i %in% bg_cells), ]
      
    }
    
    ## Reduce background sample if it's larger than max_bg_size
    if (nrow(bg) > max_bg_size) {
      
      message(nrow(bg), ' target species background records, reduced to random ',
              max_bg_size, '.')
      
      bg <- bg[sample(nrow(bg), max_bg_size), ]  ## Change this to use 
      
    } else {
      
      message(nrow(bg), ' target species background records.')
      
    }
    
    #####################################################################
    ## Save occ and bg shapefiles objects for future reference
    save_name = gsub(' ', '_', name)
    if(shapefiles) {
      
      suppressWarnings({
        
        message(name, ' writing occ and bg shapefiles')
        writeOGR(SpatialPolygonsDataFrame(b, data.frame(ID = seq_len(length(b)))),
                 outdir_sp, paste0(save_name, '_bg_buffer'), 'ESRI Shapefile', overwrite_layer = TRUE)
        writeOGR(bg,  outdir_sp, paste0(save_name, '_bg'),   'ESRI Shapefile', overwrite_layer = TRUE)
        writeOGR(occ, outdir_sp, paste0(save_name, '_occ'),  'ESRI Shapefile', overwrite_layer = TRUE)
        
      })
      
    }
    
    ## Save the background and occurrence points as objects
    saveRDS(bg,  file.path(outdir_sp, paste0(save_name, '_bg.rds')))
    saveRDS(occ, file.path(outdir_sp, paste0(save_name, '_occ.rds')))
    
    #####################################################################
    ## sdm.predictors: the s_ff object containing sdm.predictors to use in the model
    ## Sample sdm.predictors at occurrence and background points
    #####################################################################
    swd_occ <- occ[, sdm.predictors]
    saveRDS(swd_occ, file.path(outdir_sp, paste0(save_name,'_occ_swd.rds')))
    
    swd_bg <- bg[, sdm.predictors]
    saveRDS(swd_bg, file.path(outdir_sp, paste0(save_name, '_bg_swd.rds')))
    
    ## Save shapefiles of the occurrence and background points
    if(shapefiles) {
      
      writeOGR(swd_occ, outdir_sp,  paste0(save_name, '_occ_swd'), 'ESRI Shapefile', overwrite_layer = TRUE)
      writeOGR(swd_bg,  outdir_sp,  paste0(save_name, '_bg_swd'),  'ESRI Shapefile', overwrite_layer = TRUE)
      
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
      message(name, ' running xval maxent')
      me_xval <- maxent(swd, pa, path = file.path(outdir_sp, 'xval'),
                        args = c(paste0('replicates=', replicates),
                                 'responsecurves=true',
                                 'outputformat=logistic',
                                 off, paste(names(rep_args), rep_args, sep = '=')))
      
    }
    
    ## Runs the full maxent model - using all the data in swd
    ## This uses DISMO to output standard files, but the names can't be altered
    if(missing(full_args)) full_args <- NULL
    message(name, ' running full maxent')
    me_full <- maxent(swd, pa, path = file.path(outdir_sp, 'full'),
                      args = c(off, paste(names(full_args), full_args, sep = '='),
                               'responsecurves=true',
                               'outputformat=logistic'))
    
    ## Save the full model
    saveRDS(list(me_xval = me_xval, me_full = me_full, swd = swd, pa = pa, 
                 koppen_gridcode=as.character(Koppen_zones$Koppen[match(unique(zones), Koppen_zones$GRIDCODE)])), 
            file.path(outdir_sp, 'full', 'maxent_fitted.rds'))
    
    #####################################################################
    ## Save the chart corrleation file too for the variable set
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







#########################################################################################################################
## GET BACKGROUND POINTS AND THEN FIT MAXENT WITH BACKWARDS SELECTION 
#########################################################################################################################


fit_maxent_targ_bs <- function(sdm.predictors, 
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
  if(!file.exists(outdir_sp)) dir.create(outdir_sp)
  features <- unlist(strsplit(features, ''))
  
  #####################################################################
  ## Read in the occ and bg points from the targetted SDM step
  message('Reading previously created occurrence and background data from targetted SDM for ', spp)
  save_spp = gsub(' ', '_', spp)
  occ     <- readRDS(sprintf('%s%s/%s_occ.rds', maxent_path, save_spp, save_spp))
  bg      <- readRDS(sprintf('%s%s/%s_bg.rds',  maxent_path, save_spp, save_spp))
  
  swd_occ <- readRDS(sprintf('%s%s/%s_occ_swd.rds', maxent_path, save_spp, save_spp))
  swd_bg  <- readRDS(sprintf('%s%s/%s_bg_swd.rds',  maxent_path, save_spp, save_spp))
  
  
  #####################################################################
  ## Get unique cell numbers for species occurrences
  # cells <- cellFromXY(template.raster, occ)
  # 
  # ## Clean out duplicate cells and NAs (including points outside extent of predictor data)
  # ## Note this will get rid of a lot of duplicate records not filtered out by GBIF columns, etc.
  # not_dupes <- which(!duplicated(cells) & !is.na(cells)) 
  # occ       <- occ[not_dupes, ]
  # cells     <- cells[not_dupes]
  # message(nrow(occ), ' occurrence records (unique cells).')
  
  ## Skip species that have less than a minimum number of records: eg 20 species
  if(nrow(occ) < min_n) {
    
    print (paste ('Fewer occurrence records than the number of cross-validation ',
                  'replicates for species ', name, 
                  ' Model not fit for this species'))
    
  } else {
    
    #####################################################################
    ## Now subset bg to the buffer polygon
    ## Within the 200km buffer around the occurrence records for each species, select up to 100k records from any
    ## species in the total dataset :: trying to create the same bias in both the occurrence and background data
    # system.time(o <- over(bg, b))
    # bg <- bg[which(!is.na(o)), ]
    # bg_cells <- cellFromXY(template.raster, bg)
    # 
    # ## Clean out duplicates and NAs (including points outside extent of predictor data)
    # ## So we are using unique cells to select background records
    # bg_not_dupes <- which(!duplicated(bg_cells) & !is.na(bg_cells)) 
    # bg <- bg[bg_not_dupes, ]
    # bg_cells <- bg_cells[bg_not_dupes]
    # 
    # ## Reduce background sample if it's larger than max_bg_size
    # if (nrow(bg) > max_bg_size) { 
    #   
    #   message(nrow(bg), ' target species background records, reduced to random ', 
    #           max_bg_size, '.')
    #   
    #   bg <- bg[sample(nrow(bg), max_bg_size), ]
    #   
    # } else {
    #   
    #   message(nrow(bg), ' target species background records.')
    #   
    # }
    
    #####################################################################
    ## Save data for future use
    saveRDS(bg,  file.path(outdir_sp, paste0('bg.rds')))
    saveRDS(occ, file.path(outdir_sp, paste0(save_spp, '_occ.rds')))
    
    saveRDS(swd_bg,  file.path(outdir_sp, paste0('bg.rds')))
    saveRDS(swd_occ, file.path(outdir_sp, paste0('swd.rds')))

    
    #####################################################################
    ## Coerce the "species with data" (SWD) files to regular data.frames
    ## This is needed to use the simplify function 
    swd_occ     <- as.data.frame(swd_occ)
    swd_occ$lon <- NULL
    swd_occ$lat <- NULL
    swd_bg      <- as.data.frame(swd_bg)
    swd_bg$lon  <- NULL
    swd_bg$lat  <- NULL
    
    swd_occ$searchTaxon <- spp
    swd_bg$searchTaxon  <- spp
    
    #####################################################################
    ## Run simplify rmaxent::simplify
    m <- local_simplify(
      swd_occ, 
      swd_bg,
      path            = outdir, 
      species_column  = "searchTaxon",
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




#########################################################################################################################
## RUN BACKWARDS SELECTION
#########################################################################################################################


## Simplify the 
local_simplify = function (occ, bg, path, species_column = "species", response_curves = TRUE, 
                           logistic_format = TRUE, type = "PI", cor_thr, pct_thr, k_thr, 
                           features = "lpq", replicates = 1, quiet = TRUE) 
{
  if (!species_column %in% colnames(occ)) 
    stop(species_column, " is not a column of `occ`", call. = FALSE)
  if (!species_column %in% colnames(bg)) 
    stop(species_column, " is not a column of `bg`", call. = FALSE)
  if (missing(path)) {
    save <- FALSE
    path <- tempdir()
  }
  else save <- TRUE
  features <- unlist(strsplit(gsub("\\s", "", features), ""))
  if (length(setdiff(features, c("l", "p", "q", "h", "t"))) > 
      1) 
    stop("features must be a vector of one or more of ',\n         'l', 'p', 'q', 'h', and 't'.")
  off <- setdiff(c("l", "p", "q", "t", "h"), features)
  if (length(off) > 0) {
    off <- c(l = "linear=FALSE", p = "product=FALSE", q = "quadratic=FALSE", 
             t = "threshold=FALSE", h = "hinge=FALSE")[off]
  }
  off <- unname(off)
  occ_by_species <- split(occ, occ[[species_column]])
  bg_by_species <- split(bg, bg[[species_column]])
  if (!identical(sort(names(occ_by_species)), sort(names(bg_by_species)))) {
    stop("The same set of species names must exist in occ and bg")
  }
  type <- switch(type, PI = "permutation.importance", PC = "contribution", 
                 stop("type must be either \"PI\" or \"PC\".", call. = FALSE))
  args <- off
  if (replicates > 1) 
    args <- c(args, paste0("replicates=", replicates))
  if (isTRUE(response_curves)) 
    args <- c(args, "responsecurves=TRUE")
  if (isTRUE(logistic_format)) 
    args <- c(args, "outputformat=logistic")
  f <- function(name) {
    if (!quiet) 
      message("\n\nDoing ", name)
    name_ <- gsub(" ", "_", name)
    swd <- rbind(occ_by_species[[name]], bg_by_species[[name]])
    swd <- swd[, -match(species_column, names(swd))]
    if (ncol(swd) < k_thr) 
      stop("Initial number of variables < k_thr", call. = FALSE)
    pa <- rep(1:0, c(nrow(occ_by_species[[name]]), nrow(bg_by_species[[name]])))
    vc <- usdm::vifcor(swd, maxobservations = nrow(swd), 
                       th = cor_thr)
    vif <- methods::slot(vc, "results")
    k <- nrow(vif)
    exclude <- methods::slot(vc, "excluded")
    if (!isTRUE(quiet) & length(exclude) > 0) {
      message("Dropped due to collinearity: ", paste0(exclude, 
                                                      collapse = ", "))
    }
    if (k < k_thr) 
      stop(sprintf("Number of uncorrelated variables (%s) < k_thr (%s). %s", 
                   k, k_thr, "Reduce k_thr, increase cor_thr, or find alternative predictors."), 
           call. = FALSE)
    swd_uncor <- swd[, as.character(vif$Variables)]
    d <- file.path(path, name_, if (replicates > 1) 
      "xval"
      else "full")
    m <- dismo::maxent(swd_uncor, pa, args = args, path = d)
    if (isTRUE(save)) 
      saveRDS(m, file.path(d, "maxent_fitted.rds"))
    pct <- m@results[grep(type, rownames(m@results)), , 
                     drop = FALSE]
    pct <- pct[, ncol(pct)]
    pct <- sort(pct)
    names(pct) <- sub(paste0("\\.", type), "", names(pct))
    if (min(pct) >= pct_thr || length(pct) == k_thr) {
      if (replicates > 1) {
        d <- file.path(path, name_, "full")
        m <- dismo::maxent(swd_uncor, pa, args = grep("replicates", 
                                                      args, value = TRUE, invert = TRUE), path = d)
      }
      if (isTRUE(save)) {
        saveRDS(m, file.path(d, "maxent_fitted.rds"))
      }
      return(m)
    }
    while (min(pct) < pct_thr && length(pct) > k_thr) {
      candidates <- vif[vif$Variables %in% names(pct)[pct == 
                                                        pct[1]], ]
      drop <- as.character(candidates$Variables[which.max(candidates$VIF)])
      message("Dropping ", drop)
      swd_uncor <- swd_uncor[, -match(drop, colnames(swd_uncor))]
      if (!quiet) 
        message(sprintf("%s variables: %s", ncol(swd_uncor), 
                        paste0(colnames(swd_uncor), collapse = ", ")))
      m <- dismo::maxent(swd_uncor, pa, args = args, path = d)
      pct <- m@results[grep(type, rownames(m@results)), 
                       , drop = FALSE]
      pct <- pct[, ncol(pct)]
      pct <- sort(pct)
      names(pct) <- sub(paste0("\\.", type), "", names(pct))
    }
    if (replicates > 1) {
      d <- file.path(path, name_, "full")
      m <- dismo::maxent(swd_uncor, pa, args = grep("replicates", 
                                                    args, value = TRUE, invert = TRUE), path = d)
    }
    if (isTRUE(save)) {
      saveRDS(m, file.path(d, "maxent_fitted.rds"))
    }
    return(m)
  }
  lapply(names(occ_by_species), f)
}





#########################################################################################################################
## FILTER RECORDS 
#########################################################################################################################


## 
filterByProximity <- function(xy, dist, mapUnits = F) {
  
  # xy can be either a SpatialPoints or SPDF object, or a matrix
  # dist is in km if mapUnits=F, in mapUnits otherwise
  
  if (!mapUnits) {
    
    d <- spDists(xy, longlat = T)
    
  }
  
  if (mapUnits) {
    
    d <- spDists(xy,longlat = F)
    
  }
  
  diag(d) <- NA
  close   <- (d <= dist)
  diag(close) <- NA
  closePts <- which(close,arr.ind=T)
  discard  <- matrix(nrow = 2,ncol = 2)
  
  if (nrow(closePts) > 0) {
    
    while (nrow(closePts) > 0) {
      
      if ((!paste(closePts[1,1],closePts[1,2],sep='_') %in% paste(discard[,1],discard[,2],sep='_')) & 
          (!paste(closePts[1,2],closePts[1,1],sep='_') %in% paste(discard[,1],discard[,2],sep='_'))) {
        discard <- rbind(discard, closePts[1,])
        closePts <- closePts[-union(which(closePts[,1] == closePts[1,1]), which(closePts[,2] == closePts[1,1])),]
        
      }
      
    }
    
    discard <- discard[complete.cases(discard),]
    return(xy[-discard[,1],])
    
  }
  
  if (nrow(closePts) == 0) {
    return(xy)
    
  }
  
}





#########################################################################################################################
## CALCULATE TSS 
#########################################################################################################################


## TSS
maxtss2 <- function(x) {
  
  # x: a directory containing the cross-validated Maxent output
  # x = 'output/maxent/BIAS_TSS'
  ## Notes:
  ## tss  = sens + spec - 1
  ## sens = 1 - omission 
  ## spec = 1 - FPA
  ## tss  = 1 - omission + 1 - fpa - 1
  ##      = 1 - (omission + fpa)
  ##     
  
  ff <- list.files(x, 'omission\\.csv$', full.names = TRUE)
  
  max_tss <- sapply(ff, function(f) {
    
    d <- read.csv(f)
    i <- which.min(d$Test.omission + d$Fractional.area)
    
    c(max_tss = 1 - min(d$Test.omission + d$Fractional.area),
      thr     = d$Corresponding.logistic.value[i])
    
  })
  
  out <- t(max_tss)
  rownames(out) <- basename(rownames(out))
  
  list(max_tss      = out, 
       max_tss_mean = mean(out[, 'max_tss']), 
       max_tss_sd   = sd(out[, 'max_tss']))
  
}




#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################