#########################################################################################################################
#############################################  FUNCTIONS FOR HORT AUS SDMS ############################################## 
#########################################################################################################################


#########################################################################################################################
## MAXENT FUNCTIONS
#########################################################################################################################

## Need an example of actually runnin the code


#########################################################################################################################
## GET BACKGROUND POINTS AND THEN FIT MAXENT USING DISMO FUNCTIONS 
#########################################################################################################################


FIT_MAXENT <- function(occ, 
                       bg, # A Spatial points data frame (SPDF) of candidate background points
                       sdm.predictors, 
                       # sdm.predictors is a vector of enviro conditions that you want to include
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
      
      #bg_cells <- cellFromXY(template.raster, bg) ## returns just the unique cells which have background points 
      
      #####################################################################
      ## Here is where we need to reduce the predictors to a candidate set
      sdm.predictors = COR_VARIABLES(occ, bg,
                                     path,
                                     species_column  = "species",
                                     type            = "PI",
                                     cor_thr         = 0.7,
                                     quiet           = FALSE)
      
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
      
      ## Have a look at the correlation structure
      chart.Correlation(swd_occ@data, 
                        histogram = TRUE, pch = 19, 
                        main = paste0("Predictor subset <0.6 correlation for ", x))
      
      ## Save shapefiles of the occurrence and background points
      if(shapefiles) {
        
        writeOGR(swd_occ, outdir_sp,  'occ_swd', 'ESRI Shapefile', overwrite_layer = TRUE)
        writeOGR(swd_bg,  outdir_sp,  'bg_swd',  'ESRI Shapefile', overwrite_layer = TRUE)
        
      }
      
      #####################################################################
      ## Combine occurrence and bg SWD data
      swd <- as.data.frame(rbind(swd_occ@data, swd_bg@data))
      saveRDS(swd, file.path(outdir_sp, 'swd.rds'))
      pa <- rep(1:0, c(nrow(swd_occ), nrow(swd_bg)))
      
      ## Now check the features and cross validation arguments are correct
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
      
      bg_cells <- cellFromXY(template.raster, bg) ## tail(bg$searchTaxon)
      
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
      ## Also, you need the 'species' column here...
      occ <- occ[, c('species', 'searchTaxon',  sdm.predictors.all)]
      bg  <- bg[,  c('species', 'searchTaxon',  sdm.predictors.all)]
      
      swd_occ <- occ[, c('species', 'searchTaxon',  sdm.predictors.all)]
      swd_bg  <- bg[,  c('species', 'searchTaxon',  sdm.predictors.all)]
      
      ## names(swd_occ);names(swd_bg)
      
      #####################################################################
      ## Here is where we need to reduce the predictors to a candidate set
      ## The same set of species names must exist in occ and bg...really? 
      ## Not how the rest of code has been set up...
      
      ## background <- subset(SDM.DATA.ALL, searchTaxon != x) conflicts with:
      test = HIA_SIMPLIFY(occ, bg,
                          path            = outdir,
                          species_column  = "species",
                          response_curves = FALSE,
                          logistic_format = TRUE,
                          type            = "PI",
                          cor_thr         = 0.7,
                          pct_thr         = 5,
                          k_thr           = 4,
                          quiet           = FALSE)
      
      
      #####################################################################
      ## Re-create occ and bg with new predictors: is 
      # swd_occ <- swd_occ[, sdm.predictors.all]
      # swd_bg  <- swd_bg[, sdm.predictors.all]
      
      ## Then save them...
      saveRDS(swd_occ, file.path(outdir_sp, 'occ_swd.rds'))
      saveRDS(swd_bg,  file.path(outdir_sp, 'bg_swd.rds'))
      
      ## Save shapefiles of the...
      if(shapefiles) {
        
        writeOGR(swd_occ, outdir_sp,  'occ_swd', 'ESRI Shapefile', overwrite_layer = TRUE)
        writeOGR(swd_bg,  outdir_sp,  'bg_swd',  'ESRI Shapefile', overwrite_layer = TRUE)
        
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
    
  }
  
}





#########################################################################################################################
## GET UNCORRELATED VARIABLES
#########################################################################################################################


## This should replace John's code
COR_VARIABLES = function (occ, bg, path, species_column = "species", 
                          type = "PI", cor_thr,
                          quiet = TRUE) {
  
  if (missing(path)) {
    
    save <- FALSE
    path <- tempdir()
    
  }
  
  else save <- TRUE
  occ_by_species <- split(occ, occ[[species_column]])
  bg_by_species  <- split(bg, bg[[species_column]])
  
  ## head(occ)
  ## head(occ_by_species)
  ## str(bg_by_species)
  
  ## This code breaks on my data, because the background points are taken from points that are not the species 
  ## background <- subset(SDM.DATA.ALL, searchTaxon != x) ## what species is 
  
  if (!identical(sort(names(occ_by_species)), sort(names(bg_by_species)))) {    ##
    print(paste0("The same set of species names must exist in occ and bg"))
  }
  
  ## 
  type <- switch(type, PI = "permutation.importance", PC = "contribution", 
                 stop("type must be either \"PI\" or \"PC\".", call. = FALSE))
  
  ## Explain
  lapply(names(occ_by_species), function(name) {
    if (!quiet) 
      message("\n\nDoing ", name)
    
    ## The problem is in here: conflict between these lines and the main background argument:
    ## background <- subset(SDM.DATA.ALL, searchTaxon != x)
    ## Not sure why this dataframe needs to be matched?
    name_ <- gsub(" ", "_", name)
    swd <- rbind(occ_by_species[[name]], bg_by_species[[name]])
    swd <- swd[, -match(species_column, names(swd))]
    
    ##
    # if (ncol(swd) < k_thr) 
    #   stop("Initial number of variables < k_thr")
    
    ## what does round do?
    round(cor(as.data.frame(swd)[2:20], use = "pairwise"), 2)
    swd.cor = as.data.frame(swd)[2:20]
    pa <- rep(1:0, c(nrow(occ_by_species[[name]]), nrow(bg_by_species[[name]])))  ## problem here
    ok <- as.character(usdm::vifcor(swd.cor, maxobservations = nrow(swd), 
                                    th = cor_thr)@results$Variables)              ## ok
    
    ## return the set of predictors which are less correlated than the threshold (e.g. 0.7)
    return(ok)
    
  })
  
}





#########################################################################################################################
## RUN BACKWARDS SELECTION
#########################################################################################################################


## This should replace John's code
HIA_SIMPLIFY = function (occ, bg, path, species_column = "species", response_curves = FALSE, 
                         logistic_format = TRUE, type = "PI", cor_thr, pct_thr, k_thr, 
                         quiet = TRUE) 
  
{
  
  if (missing(path)) {
    
    save <- FALSE
    path <- tempdir()
    
  }
  
  else save <- TRUE
  occ_by_species <- split(occ, occ[[species_column]])
  bg_by_species  <- split(bg, bg[[species_column]])
  
  names(occ);names(bg)
  names(occ_by_species);names(bg_by_species)
  
  ## This code breaks on my data, because the background points are taken from points that are not the species 
  ## background <- subset(SDM.DATA.ALL, searchTaxon != x) ## what species is 
  if (!identical(sort(names(occ_by_species)), sort(names(bg_by_species)))) {    ##
    print(paste0("The same set of species names must exist in occ and bg"))
  }
  
  ## 
  type <- switch(type, PI = "permutation.importance", PC = "contribution", 
                 stop("type must be either \"PI\" or \"PC\".", call. = FALSE))
  
  ## Maxent args?
  args <- c("threshold=false", "hinge=false")
  if (isTRUE(response_curves)) 
    args <- c(args, "responsecurves=TRUE")
  
  if (isTRUE(logistic_format)) 
    args <- c(args, "outputformat=logistic")
  
  
  ## Explain
  lapply(names(occ_by_species), function(name) {
    if (!quiet) 
      message("\n\nDoing ", name)
    
    ## The problem is in here: conflict between these lines and the main background argument:
    ## background <- subset(SDM.DATA.ALL, searchTaxon != x)
    ## Not sure why this dataframe needs to be matched?
    name_ <- gsub(" ", "_", name)
    swd <- rbind(occ_by_species[[name]], bg_by_species[[name]])
    swd <- swd[, -match(species_column, names(swd))]
    ## str(swd)
    
    ##
    if (ncol(swd) < k_thr) 
      stop("Initial number of variables < k_thr")
    
    ## What does round do? This would not work by itself
    round(cor(as.data.frame(swd)[2:20], use = "pairwise"), 2)
    swd.cor = as.data.frame(swd)[2:20]
    pa <- rep(1:0, c(nrow(occ_by_species[[name]]), nrow(bg_by_species[[name]])))  ## problem here
    ok <- as.character(usdm::vifcor(swd.cor, maxobservations = nrow(swd), 
                                    th = cor_thr)@results$Variables)              ## ok
    
    ## why did I need to change from a spatial points to DF to get it working
    # Error in (function (classes, fdef, mtable)  : 
    #             unable to find an inherited method for function ‘maxent’ for signature ‘"SpatialPointsDataFrame", "integer"’
    
    swd_uncor     <- as.data.frame(swd[, ok])                     ##
    swd_uncor$lon <- NULL
    swd_uncor$lat <- NULL
    str(swd_uncor)
    
    d <- file.path(path, name_, "full")
    m <- dismo::maxent(swd_uncor, pa, args = args, path = d)  ## Seems to work, but why is lat/long being used?
    
    ##
    if (isTRUE(save)) 
      saveRDS(m, file.path(d, "model.rds"))
    
    pct <- m@results[grep(type, rownames(m@results)), ]
    pct <- sort(pct[pct > 0])
    
    names(pct) <- sub(paste0("\\.", type), "", names(pct))
    
    if (min(pct) >= pct_thr || length(pct) <= k_thr) {
      if (isTRUE(save)) {
        d_out <- file.path(path, name_, "final")
        dir.create(d_out)
        file.copy(list.files(d, full.names = TRUE), 
                  d_out, recursive = TRUE)
        
      }
      
      return(m)
      
    }
    
    ## Explain
    while (min(pct) < pct_thr && length(pct) > k_thr) {
      
      message("Dropping ", names(pct)[1])
      swd_uncor <- swd_uncor[, -match(names(pct)[1], names(swd_uncor))]
      tmp <- tempfile()
      
      if (!quiet) 
        message(sprintf("%s variables: %s", ncol(swd_uncor), 
                        paste0(colnames(swd_uncor), collapse = ", ")))
      m <- dismo::maxent(swd_uncor, pa, args = args, path = tmp)
      pct <- m@results[grep(type, rownames(m@results)), 
                       ]
      pct <- sort(pct)
      names(pct) <- sub(paste0("\\.", type), "", names(pct))
      
    }
    
    ## Explain
    if (isTRUE(save)) {
      
      d_out <- file.path(path, name_, "final")
      file.copy(tmp, file.path(path, name_), recursive = TRUE)
      file.rename(file.path(path, name_, basename(tmp)), 
                  d_out)
      saveRDS(m, file.path(path, name_, "final/model.rds"))
      
    }
    
    ## return the set of predictors which are less correlated
    return(m)
    
  })
  
}



#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################