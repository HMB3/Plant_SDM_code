#########################################################################################################################
#############################################  FUNCTIONS FOR HORT AUS SDMS ############################################## 
#########################################################################################################################


#########################################################################################################################
## THRESHOLD FUNCTIONS
########################################################################################################################


#########################################################################################################################
## MAXENT FUNCTIONS
#########################################################################################################################

# ## Arguments to run maxent line by line
# occ                     = occurrence
# bg                      = background
# sdm.predictors          = sdm.predictors
# name                    = spp
# outdir                  = 'output/maxent/STD_VAR_ALL'
# template.raster         = template.raster
# min_n                   = 20   ## This should be higher...
# max_bg_size             = 100000
# background_buffer_width = 200000
# shapefiles              = TRUE
# features                = 'lpq'
# replicates              = 5
# cor_thr                 = 0.7
# pct_thr                 = 5
# k_thr                   = 5
# responsecurves          = TRUE
# 
# 
# ## selection line by line 
# occ             = swd_occ
# bg              = swd_bg
# path            = outdir
# species_column  = "species"
# replicates      = replicates
# response_curves = TRUE
# logistic_format = TRUE
# cor_thr         = 0.7
# pct_thr         = 5
# k_thr           = 5
# features        ='lpq'  # change these as necessary (or cor_thr = cor_thr, etc from FIT_MAXENT_SIMP)
# quiet           = FALSE
# type            = "PI"


#########################################################################################################################
## GET BACKGROUND POINTS AND THEN FIT MAXENT WITH BACKWARDS SELECTION 
#########################################################################################################################

##
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
                                 max_bg_size, 
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
      
      #####################################################################
      ## Here is where we need to reduce the predictors to a candidate set
      
      
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
      ## sdm.predictors: the s_ff object containing sdm.predictors to use in 
      ## the model. Sample sdm.predictors at occurrence and background points
      #####################################################################
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
      ## Run simplify
      #debugonce(RMAXENT_SIMPLIFY)
      ##debugonce(f)
      m <- RMAXENT_SIMPLIFY(
        
        swd_occ, 
        bg,
        path            = outdir, 
        species_column  = "species",
        #name            = spp,
        replicates      = replicates,
        response_curves = TRUE, 
        logistic_format = TRUE, 
        cor_thr         = 0.7, 
        pct_thr         = 5, 
        k_thr           = 5, 
        features        ='lpq',  ## change these as necessary (or cor_thr = cor_thr, etc from FIT_MAXENT_SIMP)
        quiet           = FALSE)
      
    }
    
  }
  
}



#########################################################################################################################
## TWEAK OF RMAXENT::SIMPLIFY
#########################################################################################################################


## Not sure why this doesn't work out of the box
RMAXENT_SIMPLIFY = function (occ, bg, path, 
                             species_column  = "species", 
                             #name            = spp,
                             response_curves = TRUE, 
                             logistic_format = TRUE, 
                             type = "PI", 
                             cor_thr, 
                             pct_thr, 
                             k_thr, 
                             features = "lpq", 
                             replicates = 1, 
                             quiet = TRUE) 
  
{
  
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
    
    ## debug at ./R/SDM_FUNCTIONS.R#297: swd <- rbind(occ_by_species[[name]], bg_by_species[[name]])
    name_ <- gsub(" ", "_", name)
    swd   <- rbind(occ_by_species[[name]], bg_by_species[[name]]) 
    
    # Error in rbind2(..1, r) : 
    #   no method for coercing this S4 class to a vector
    
    swd   <- swd[, -match(species_column, names(swd))]
    
    # Error in rbind2(..1, r) : 
    #   no method for coercing this S4 class to a vector
    if (ncol(swd) < k_thr) 
      stop("Initial number of variables < k_thr", call. = FALSE)
    
    pa <- rep(1:0, c(nrow(occ_by_species[[name]]), nrow(bg_by_species[[name]])))
    vc <- usdm::vifcor(swd, maxobservations = nrow(swd), 
                       th = cor_thr)
    
    vif <- slot(vc, "results")
    k <- nrow(vif)
    
    exclude <- slot(vc, "excluded")
    if (!isTRUE(quiet) & length(exclude) > 0) {
      
      message("Dropped due to collinearity: ", paste0(exclude, 
                                                      collapse = ", "))
    }
    
    if (k < k_thr) 
      stop(sprintf("Number of uncorrelated variables (%s) < k_thr (%s). %s", 
                   k, k_thr, "Reduce k_thr and/or cor_thr, or find alternative predictors."), 
           call. = FALSE)
    
    swd_uncor <- swd[, vif$Variables]
    d <- file.path(path, name_, if (replicates > 1) ## Is this ok?
      "xval"
      else "full")
    
    m <- dismo::maxent(swd_uncor, pa, args = args, path = d) ## takes ages because backwards selection is happening...
    if (isTRUE(save)) 
      saveRDS(m, file.path(d, "model.rds"))
    
    pct <- m@results[grep(type, rownames(m@results)), , drop = FALSE]
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
        
        saveRDS(m, file.path(d, "model.rds"))
        
      }
      
      return(m)
    }
    
    ## error occurs between here
    ## Error in FUN(X[[i]], ...) : 
    ## cannot coerce type 'closure' to vector of type 'character'
    while (min(pct) < pct_thr && length(pct) > k_thr) {
      
      if (sum(pct == pct[1]) > 1) {
        
        candidates <- subset(vif, Variables %in% names(pct)[pct == 
                                                              pct[1]])
        drop <- as.character(candidates$Variables[which.max(candidates$VIF)])
        
      }
      
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
      
      ## And here. But only in the loop. Cam be run line by line...
      
    }
    
    if (replicates > 1) {
      
      ## Run maxent for each replicate, doing the backwards selection each time?
      d <- file.path(path, name_, "full")
      m <- dismo::maxent(swd_uncor, pa, args = grep("replicates", 
                                                    args, value = TRUE, invert = TRUE), path = d)
      
    }
    
    if (isTRUE(save)) {
      
      saveRDS(m, file.path(d, "model.rds"))
      
    }
    
    return(m)
    
  }  ## function ends here
  lapply(names(occ_by_species), f) 
  ## Iterate over colnames of a df (i.e. environmental variables, names(occ_by_species[[1]]))? 
  ## lapply(names(occ_by_species), f) 
  ## lapply(name, f) 
  
}





#########################################################################################################################
##

# FIT_MAXENT <- function(occ, 
#                        bg, # A Spatial points data frame (SPDF) of candidate background points
#                        sdm.predictors, 
#                        # sdm.predictors is a vector of enviro conditions that you want to include
#                        name, 
#                        outdir, 
#                        template.raster, 
#                        # template.raster is an empty raster with extent, res and projection
#                        # of final output rasters. It is used to reduce
#                        # occurrences to a single point per cell.
#                        min_n                   = 20,
#                        # min_n is the minimum number of records (unique cells)
#                        # required for a model to be fit
#                        max_bg_size             = 100000, 
#                        background_buffer_width = 200000,
#                        shapefiles              = TRUE, 
#                        features, 
#                        replicates, # number of cross-validation replicates
#                        responsecurves          = TRUE, 
#                        rep_args, 
#                        full_args) {
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
#          'l', 'p', 'q', 'h', and 't'.")
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
#       ## Now subset bg to the buffer polygon
#       system.time(o <- over(bg, b))
#       bg <- bg[which(!is.na(o)), ]
#       bg_cells <- cellFromXY(template.raster, bg)
#       
#       ## Clean out duplicates and NAs (including points outside extent of predictor data)
#       bg_not_dupes <- which(!duplicated(bg_cells) & !is.na(bg_cells)) 
#       bg <- bg[bg_not_dupes, ]
#       bg_cells <- bg_cells[bg_not_dupes]
#       
#       ## Reduce background sample if it's larger than max_bg_size
#       if (nrow(bg) > max_bg_size) { 
#         
#         message(nrow(bg), ' target species background records, reduced to random ', 
#                 max_bg_size, '.')
#         
#         bg <- bg[sample(nrow(bg), max_bg_size), ]
#         
#       } else {
#         
#         message(nrow(bg), ' target species background records.')
#         
#       }
#       
#       #####################################################################
#       ## Here is where we need to reduce the predictors to a candidate set
#       sdm.predictors = COR_VARIABLES(occ, bg,
#                                      path,
#                                      species_column  = "species",
#                                      type            = "PI",
#                                      cor_thr         = 0.6,
#                                      quiet           = FALSE)
#       
#       sdm.predictors = sdm.predictors[[1]]  ## Note this must be a character vector, not a list
#       
#       #####################################################################
#       ## Save objects for future reference
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
#       ## Save the background and occurrence points as objects
#       saveRDS(bg,  file.path(outdir_sp, 'bg.rds'))
#       saveRDS(occ, file.path(outdir_sp, 'occ.rds'))
#       
#       #####################################################################
#       ## sdm.predictors: the s_ff object containing sdm.predictors to use in the model
#       ## Sample sdm.predictors at occurrence and background points
#       #####################################################################
#       swd_occ <- occ[, sdm.predictors]
#       saveRDS(swd_occ, file.path(outdir_sp, 'occ_swd.rds'))
#       
#       swd_bg <- bg[, sdm.predictors]
#       saveRDS(swd_bg, file.path(outdir_sp, 'bg_swd.rds'))
#       
#       ## Have a look at the correlation structure
#       chart.Correlation(swd_occ@data, 
#                         histogram = TRUE, pch = 19, 
#                         main = paste0("Uncorrelated predictor subset"))
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
#       ## Combine occurrence and bg SWD data
#       swd <- as.data.frame(rbind(swd_occ@data, swd_bg@data))
#       saveRDS(swd, file.path(outdir_sp, 'swd.rds'))
#       pa <- rep(1:0, c(nrow(swd_occ), nrow(swd_bg)))
#       
#       ## Now check the features and cross validation arguments are correct
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
#         ## This runs the MAXENT for for the cross validation: so 5 replicaes, 0-4
#         ## EG xval = cross validation : "OUT_DIR\Acacia_boormanii\xval\maxent_0.html"
#         me_xval <- maxent(swd, pa, path = file.path(outdir_sp, 'xval'), 
#                           args = c(paste0('replicates=', replicates),
#                                    'responsecurves=true', 
#                                    'outputformat=logistic',
#                                    off, paste(names(rep_args), rep_args, sep = '=')))
#         
#       }
#       
#       ## Not sure why we need all these switches...is this stuff John needs for other purposes? 
#       if(missing(full_args)) full_args <- NULL
#       me_full <- maxent(swd, pa, path = file.path(outdir_sp, 'full'), 
#                         args = c(off, paste(names(full_args), full_args, sep = '='),
#                                  'responsecurves=true',
#                                  'outputformat=logistic'))
#       
#       #####################################################################
#       ## Save fitted model object, and the model-fitting data.
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
#       }
#       
#     }
#     
#   }
#   
# }





#########################################################################################################################
## GET BACKGROUND POINTS AND FIT MAXENT FOR ALL VARIABLES, USE MODEL SELECTION
#########################################################################################################################

# #########################################################################################################################
# ## GET UNCORRELATED VARIABLES
# #########################################################################################################################
# 
# 
# ## This should replace John's code
# COR_VARIABLES = function (occ, 
#                           bg, 
#                           path, 
#                           species_column = "species", 
#                           type = "PI", 
#                           cor_thr,
#                           quiet = TRUE) {
#   
#   if (missing(path)) {
#     
#     save <- FALSE
#     path <- tempdir()
#     
#   }
#   
#   else save <- TRUE
#   occ_by_species <- split(occ, occ[[species_column]])
#   bg_by_species  <- split(bg, bg[[species_column]])
#   
#   ## head(occ)
#   ## head(occ_by_species)
#   ## str(bg_by_species)
#   
#   ## This code breaks on my data, because the background points are taken from points that are not the species 
#   ## background <- subset(SDM.DATA.ALL, searchTaxon != x) ## what species is 
#   
#   if (!identical(sort(names(occ_by_species)), sort(names(bg_by_species)))) {    ##
#     
#     print(paste0("The same set of species names must exist in occ and bg"))
#     
#   }
#   
#   ## This sets the contribution of each variable before dropping each one...
#   type <- switch(type, PI = "permutation.importance", PC = "contribution", 
#                  stop("type must be either \"PI\" or \"PC\".", call. = FALSE))
#   
#   ## Explain
#   lapply(names(occ_by_species), function(name) {
#     if (!quiet) 
#       message("\n\nDoing ", name)
#     
#     ## The problem is in here: conflict between these lines and the main background argument:
#     ## background <- subset(SDM.DATA.ALL, searchTaxon != x)
#     ## Not sure why this dataframe needs to be matched?
#     name_ <- gsub(" ", "_", name)
#     swd <- rbind(occ_by_species[[name]], bg_by_species[[name]])
#     swd <- swd[, -match(species_column, names(swd))]
#     
#     ## Stop if there are an insufficent number of variables
#     if (ncol(swd) < k_thr)
#       stop("Initial number of variables < k_thr")
#     
#     ## What does round do?
#     round(cor(as.data.frame(swd)[2:20], use = "pairwise"), 2)
#     swd.cor = as.data.frame(swd)[2:20]
#     pa <- rep(1:0, c(nrow(occ_by_species[[name]]), nrow(bg_by_species[[name]])))  ## problem here
#     ok <- as.character(usdm::vifcor(swd.cor, maxobservations = nrow(swd), 
#                                     th = cor_thr)@results$Variables)              ## ok
#     
#     ## return the set of predictors which are less correlated than the threshold (e.g. 0.7)
#     return(m)
#     
#     
#     
#     
#   })
#   
# }



#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################