#########################################################################################################################
#############################################  FUNCTIONS FOR RUNNING SDMS ############################################### 
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
      ## Run simplify rmaxent::simplify
      m <- rmaxent::simplify(
        
        swd_occ, 
        swd_bg,
        path            = outdir, 
        species_column  = "species",
        replicates      = replicates,
        response_curves = TRUE, 
        logistic_format = TRUE, 
        cor_thr         = 0.7, 
        pct_thr         = 5, 
        k_thr           = 3, 
        features        ='lpq',  ## change these as necessary (or cor_thr = cor_thr, etc from FIT_MAXENT_SIMP)
        quiet           = FALSE)
      
    }
    
  }
  
}





#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################