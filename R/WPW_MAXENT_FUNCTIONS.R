#########################################################################################################################
############################################# FUNCTIONS FOR RUNNING SDMs AND MAPPING #################################### 
#########################################################################################################################


## flag issues with .........................................................................


#########################################################################################################################
## MAXENT FUNCTIONS
#########################################################################################################################


## Here are the argumetns needed to run the targetted background selection SDMs inside the function itself
# spp                     = GBIF.spp[1]
## This is what is causing the proportional sampling to skip.........................................
# occ <- subset(SDM.SPAT.OCC.BG, searchTaxon == spp)
# occ <- occ[grep(paste(OCC_SOURCE, collapse = '|'), occ$SOURCE, ignore.case = TRUE),]
# message('Using occ records from ', unique(occ$SOURCE))


## Now get the background points. These can come from any species, other than the modelled species.
## However, they should be limited to the same SOURCE as the occ data
# bg <- subset(SDM.SPAT.OCC.BG, searchTaxon != spp)
# bg <- bg[grep(paste(unique(occ$SOURCE), collapse = '|'), bg$SOURCE, ignore.case = TRUE),]
# message('Using occ records from ', unique(bg$SOURCE))
# 
# name                    = spp
# outdir                  = maxent_dir
# bsdir                   = bs_dir
# 
# backwards_sel           = "TRUE"
# cor_thr                 = 0.8      ## The maximum allowable pairwise correlation between predictor variables
# pct_thr                 = 5        ## The minimum allowable percent variable contribution
# k_thr                   = 4
# 
# 
# template.raster         = template.raster.1km   ## 1km, 5km, 10km
# min_n                   = 20
# max_bg_size             = 70000
# Koppen                  = Koppen_1975_1km
# background_buffer_width = 200000
# shapefiles              = TRUE
# features                = 'lpq'
# replicates              = 5
# responsecurves          = TRUE
# shp_path                = "./data/base/CONTEXTUAL/"
# aus_shp                 = "aus_states.rds"# sdm.predictors          = bs.predictors



#########################################################################################################################
## GET TARGETTED BACKGROUND POINTS AND THEN FIT MAXENT WITH STANDARD VARIABLES
#########################################################################################################################

## check the messages regarding sampling
fit_maxent_targ_bg_back_sel <- function(occ,
                                        bg, # A Spatial points data frame (SPDF) of candidate background points
                                        sdm.predictors,
                                        # sdm.predictors is a vector of enviro conditions that you want to include
                                        name,
                                        outdir,
                                        bsdir,
                                        cor_thr,                 
                                        pct_thr, 
                                        k_thr,
                                        
                                        backwards_sel,
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
                                        full_args,
                                        shp_path, 
                                        aus_shp) {
  
  ########################################################################
  ## First, stop if the outdir file exists,
  if(!file.exists(outdir)) stop('outdir does not exist :(', call. = FALSE)
  outdir_sp <- file.path(outdir, gsub(' ', '_', name))
  bsdir_sp  <- file.path(bsdir,  gsub(' ', '_', name))
  
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
  ## h: hinge       ## disabled for this analysis
  ## t: threshold   ## disabled for this analysis
  if(length(setdiff(features, c('l', 'p', 'q', 'h', 't'))) > 1)
    stop("features must be a vector of one or more of ',
         'l', 'p', 'q', 'h', and 't'.")
  
  ## Create a buffer of xkm around the occurrence points
  buffer <- aggregate(gBuffer(occ, width = background_buffer_width, byid = TRUE))
  
  #####################################################################
  ## Get unique cell numbers for species occurrences
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
    ## Subset the background records to the 200km buffered polygon
    message(name, ' creating background cells')
    system.time(o <- over(bg, buffer))
    bg <- bg[which(!is.na(o)), ]
    bg_cells <- cellFromXY(template.raster, bg)
    
    ## Clean out duplicates and NAs (including points outside extent of predictor data)
    bg_not_dupes <- which(!duplicated(bg_cells) & !is.na(bg_cells))
    bg           <- bg[bg_not_dupes, ]
    bg_cells     <- bg_cells[bg_not_dupes]
    
    ## Find which of these cells fall within the Koppen-Geiger zones that the species occupies
    ## Crop the Kopppen raster to the extent of the occurrences, and snap it
    message(name, ' intersecting background cells with Koppen zones')
    Koppen_crop <- crop(Koppen, occ, snap = 'out')
    
    ## Only extract and match those cells that overlap between the ::
    ## 1). cropped koppen zone, 
    ## 2). occurrences and 
    ## 3). background points 
    message(xres(template.raster), ' metre cell size for template raster')
    message(xres(Koppen), ' metre cell size for Koppen raster')
    zones               <- raster::extract(Koppen_crop, occ)
    cells_in_zones_crop <- Which(Koppen_crop %in% zones, cells = TRUE)
    cells_in_zones      <- cellFromXY(Koppen, xyFromCell(Koppen_crop, cells_in_zones_crop))
    bg_cells            <- intersect(bg_cells, cells_in_zones)  ## this is 0 for 5km 
    i                   <- cellFromXY(template.raster, bg)
    bg                  <- bg[which(i %in% bg_cells), ]
    
    ## For some species, we have the problem that the proportion of ALA/INV data is 
    ## very different in the occurrence vs the bg records. 
    ## This should be caused by the 200km / koppen restriction, etc.
    
    #####################################################################
    ## Get the proportion of occurrence records from each source 
    source.prop = round(with(as.data.frame(occ), table(SOURCE)/sum(table(SOURCE))), 3)
    ala.prop    = sum(source.prop["ALA"], source.prop["GBIF"], na.rm = TRUE)
    inv.prop    = source.prop["INVENTORY"]
    
    ## Reduce background sample, if it's larger than max_bg_size
    if (nrow(bg) > max_bg_size) {
      
      message(nrow(bg), ' target species background records for ', name, 
              ', reduced to random ', max_bg_size, ' using random points from :: ', unique(bg$SOURCE))
      bg.samp <- bg[sample(nrow(bg), max_bg_size), ]
      
    } else {
      
      ## If the bg points are smaller that the max_bg_size, just get all the points
      message(nrow(bg), ' target species background records for ', name, 
              ' using all points from :: ', unique(bg$SOURCE))
      bg.samp <- bg
      
    }
    
    ## Only if the occ and bg sources > 2, and include inv, do we need the proportional sampling
    ## otherwise, ditch it.
    if (length(unique(occ$SOURCE)) >= 2 &&
        length(unique(bg$SOURCE))  >= 2 && 
        "INVENTORY" %in% unique(occ$SOURCE) &&
        "INVENTORY" %in% unique(bg$SOURCE)) { 
      
      ## Sample background records from ALA/GBIF and INVENTORY categories in proportion with 
      ## the number of records from each category in the occ data
      message(nrow(bg.samp), ' target species background records for ',name, 
              ', using proportional samples from :: ', unique(bg$SOURCE))
      
      ## Get the ALA
      if(ala.prop*(nrow(bg.samp)) < nrow(subset(bg.samp, SOURCE == "ALA" | SOURCE == "GBIF"))) {
        
        ## These bg records don't overlap with the inventory bg records
        bg.ala = bg.samp[ sample( which(bg.samp$SOURCE == "ALA" | bg.samp$SOURCE == "GBIF"), 
                                  ala.prop * nrow(bg.samp)), ]
        
      } else {
        
        ## Otherwise, get as many inventory records as possible
        bg.ala  = bg.samp[ sample(which(bg.samp$SOURCE == "ALA" | bg.samp$SOURCE == "GBIF")), ]
        
      }
      
      ## If the proportion of inv records in the occ data * the number of background records 
      ## is < than the no. of inventory records in the background data, 
      ## get the proportion of occ inv records * the max bg size
      if(inv.prop*(nrow(bg.samp)) < nrow(subset(bg.samp, SOURCE == "INVENTORY"))) {
        
        ## These bg records don't overlap with the ala bg records
        bg.inv  = bg.samp[ sample( which(bg.samp$SOURCE == "INVENTORY"), inv.prop * nrow(bg.samp)), ]
        
      } else {
        
        ## Otherwise, get as many inventory records as possible from the greater pool of records
        ## This line is why the inventory proportion is sometimes greater. It could be changed to 
        ## bg.samp with Alessandro's data, if there are a lot more records. 
        bg.inv  = bg[ sample(which(bg$SOURCE == "INVENTORY")), ]
        
      }
      
      ## Then, combine the samples from the ALA/GBIF and INV sources.
      ## ALA/GBIF is almost always bigger, so this is the most sensible option
      bg.comb    = rbind(bg.ala, bg.inv)
      
    } else {
      ## If inventory is not in the set, get the background records from any source
      message(nrow(bg.samp), ' target species background records for ', name, 
              ' using random sample from :: ', unique(bg$SOURCE))
      bg.comb = bg.samp
      
    }
    
    #####################################################################
    ## Now save an image of the background points 
    save_name = gsub(' ', '_', name)
    
    aus.mol = readRDS(paste0(shp_path, aus_shp)) %>%
      spTransform(projection(buffer))
    
    aus.kop = crop(Koppen_crop, aus.mol)
    
    occ.mol <- occ %>%
      spTransform(projection(buffer))
    
    ## Print the koppen zones, occurrences and points to screen
    plot(aus.kop, legend = FALSE,
         main = paste0('Occurence SDM records for ', name))
    
    plot(aus.mol, add = TRUE)
    plot(buffer,  add = TRUE, col = "red")
    plot(occ.mol, add = TRUE, col = "blue")
    
    ## Then save the occurrence points
    png(sprintf('%s/%s/%s_%s.png', maxent_path, save_name, save_name, "buffer_occ"),
        16, 10, units = 'in', res = 300)
    
    plot(aus.kop, legend = FALSE,
         main = paste0('Occurence SDM records for ', name))
    
    plot(aus.mol, add = TRUE)
    plot(buffer,  add = TRUE, col = "red")
    plot(occ.mol, add = TRUE, col = "blue")
    
    dev.off()
    
    ## Only save inventory data if it exists in occ
    if ("INVENTORY" %in% unique(occ$SOURCE) == "TRUE") {
      
      inv.mol <- bg.inv  %>%
        spTransform(projection(buffer))  
      
      png(sprintf('%s/%s/%s_%s.png', maxent_path, save_name, save_name, "buffer_inv"),
          16, 10, units = 'in', res = 300)
      
      plot(aus.kop, legend = FALSE, 
           main = paste0('Inventory sdm records for ', save_name))
      plot(aus.mol, add = TRUE)
      plot(buffer,  add = TRUE, col = "red")
      plot(inv.mol, add = TRUE, col = "black")
      
      dev.off()
      
    } else {
      
      message('No inventory data, dont make bg map')
      
    }
    
    #####################################################################
    ## Now save the buffer, the occ and bg points as shapefiles
    if(shapefiles) {
      
      suppressWarnings({
        
        message(name, ' writing occ and bg shapefiles')
        writeOGR(SpatialPolygonsDataFrame(buffer, data.frame(ID = seq_len(length(buffer)))),
                 outdir_sp, paste0(save_name, '_bg_buffer'),          'ESRI Shapefile', overwrite_layer = TRUE)
        writeOGR(bg.comb,  outdir_sp, paste0(save_name, '_bg'),       'ESRI Shapefile', overwrite_layer = TRUE)
        writeOGR(occ,           outdir_sp, paste0(save_name, '_occ'), 'ESRI Shapefile', overwrite_layer = TRUE)
        
      })
      
    }
    
    ## Also save the background and occurrence points as .rds files
    saveRDS(bg.comb,  file.path(outdir_sp, paste0(save_name, '_bg.rds')))
    saveRDS(occ,      file.path(outdir_sp, paste0(save_name, '_occ.rds')))
    
    #####################################################################
    ## SWD = species with data. Now sample the environmental 
    ## variables used in the model at all the occ and bg points
    swd_occ <- occ[, sdm.predictors]
    saveRDS(swd_occ, file.path(outdir_sp, paste0(save_name,'_occ_swd.rds')))
    
    swd_bg <- bg.comb[, sdm.predictors]
    saveRDS(swd_bg, file.path(outdir_sp, paste0(save_name, '_bg_swd.rds')))
    
    ## Save the SWD tables as shapefiles
    if(shapefiles) {
      
      writeOGR(swd_occ, outdir_sp,  paste0(save_name, '_occ_swd'), 'ESRI Shapefile', overwrite_layer = TRUE)
      writeOGR(swd_bg,  outdir_sp,  paste0(save_name, '_bg_swd'),  'ESRI Shapefile', overwrite_layer = TRUE)
      
    }
    
    #####################################################################
    ## Now combine the occurrence and background data
    swd <- as.data.frame(rbind(swd_occ@data, swd_bg@data))
    saveRDS(swd, file.path(outdir_sp, 'swd.rds'))
    pa  <- rep(1:0, c(nrow(swd_occ), nrow(swd_bg)))
    
    ## Now, set the features to be used by maxent ::
    ## Linear, product and quadratic
    off <- setdiff(c('l', 'p', 'q', 't', 'h'), features)
    
    ## This sets threshold and hinge features to "off"
    if(length(off) > 0) {
      
      off <- c(l = 'linear=false',    p = 'product=false', q = 'quadratic=false',
               t = 'threshold=false', h = 'hinge=false')[off]
      
    }
    
    off <- unname(off)
    
    if(replicates > 1) {
      
      if(missing(rep_args)) rep_args <- NULL
      
      ## Run MAXENT for cross validation data splits of swd : so 5 replicaes, 0-4
      ## first argument is the predictors, the second is the occurrence data
      message(name, ' running xval maxent')
      me_xval <- maxent(swd, pa, path = file.path(outdir_sp, 'xval'),
                        args = c(paste0('replicates=', replicates),
                                 'responsecurves=true',
                                 'outputformat=logistic',
                                 off, paste(names(rep_args), rep_args, sep = '=')))
      
    }
    
    ## Run the full maxent model - using all the data in swd
    ## This uses DISMO to output standard files, but the names can't be altered
    if(missing(full_args)) full_args <- NULL
    message(name, ' running full maxent')
    me_full <- maxent(swd, pa, path = file.path(outdir_sp, 'full'),
                      args = c(off, paste(names(full_args), full_args, sep = '='),
                               'responsecurves=true',
                               'outputformat=logistic'))
    
    ## Save the full model. Replicate this line in the backwards selection algortithm
    ## This is needed to project the models.........................................
    ## Also worth checking that the koppen zones can be used at any resolution
    saveRDS(list(me_xval = me_xval, me_full = me_full, swd = swd, pa = pa, 
                 koppen_gridcode=as.character(Koppen_zones$Koppen[match(unique(zones), Koppen_zones$GRIDCODE)])), 
            file.path(outdir_sp, 'full', 'maxent_fitted.rds'))
    
    #####################################################################
    ## Save the chart corrleation file too for the training data set
    # png(sprintf('%s/%s/full/%s_%s.png', outdir,
    #             save_name, save_name, "all_vars_predictor_correlation"),
    #     3236, 2000, units = 'px', res = 300)
    # 
    # ## set margins
    # par(mar   = c(3, 3, 5, 3),
    #     oma   = c(1.5, 1.5, 1.5, 1.5))
    # 
    # ## Add detail to the response plot
    # chart.Correlation(swd_occ@data,
    #                   histogram = TRUE, pch = 19,
    #                   main = paste0('Full variable correlations for ', save_name)) 
    # 
    # ## Finish the device
    # dev.off()
    
    
    if (backwards_sel == "TRUE") {
      
      #####################################################################
      ## Coerce the "species with data" (SWD) files to regular data.frames
      ## This is needed to use the simplify function 
      swd_occ     <- as.data.frame(swd_occ)
      swd_occ$lon <- NULL
      swd_occ$lat <- NULL
      swd_bg      <- as.data.frame(swd_bg)
      swd_bg$lon  <- NULL
      swd_bg$lat  <- NULL
      
      ## Need to create a species column here
      swd_occ$searchTaxon <- name
      swd_bg$searchTaxon  <- name
      
      #####################################################################
      ## Run simplify rmaxent::simplify
      
      # Given a candidate set of predictor variables, this function identifies 
      # a subset that meets specified multicollinearity criteria. Subsequently, 
      # backward stepwise variable selection (VIF) is used to iteratively drop 
      # the variable that contributes least to the model, until the contribution 
      # of each variable meets a specified minimum, or until a predetermined 
      # minimum number of predictors remains. It returns a model object for the 
      # full model, rather than a list of models as does the previous function
      
      ## Using a modified versionof rmaxent::simplify, so that the name of the
      ## maxent model object "maxent_fitted.rds" is the same in both models.
      ## This is needed to run the mapping step over either the full or BS folder
      m <- local_simplify(
        
        swd_occ, 
        swd_bg,
        path            = bsdir, 
        species_column  = "searchTaxon",
        replicates      = replicates,  ## 5 as above
        response_curves = TRUE, 
        logistic_format = TRUE, 
        cor_thr         = cor_thr, 
        pct_thr         = pct_thr, 
        k_thr           = k_thr, 
        features        = features,    ## LPQ as above
        quiet           = FALSE)
      
      ## Save the bg, occ and swd files into the backwards selection folder too
      saveRDS(bg.comb,  file.path(bsdir_sp, paste0(save_name, '_bg.rds')))
      saveRDS(occ,      file.path(bsdir_sp, paste0(save_name, '_occ.rds')))
      saveRDS(swd,      file.path(bsdir_sp, paste0('swd.rds')))
      
      ## Read the model in, because it's tricky to index
      bs.model <- readRDS(sprintf('%s/%s/full/maxent_fitted.rds', bsdir,  save_name))
      identical(length(bs.model@presence$Annual_mean_temp), nrow(occ))
      
      #####################################################################
      ## Save the chart corrleation file too for the training data set
      par(mar   = c(3, 3, 5, 3),
          oma   = c(1.5, 1.5, 1.5, 1.5))
      
      ## Add detail to the response plot
      # chart.Correlation(bs.model@presence,
      #                   histogram = TRUE, pch = 19,
      #                   title = paste0('Reduced variable correlations for ', save_name)) 
      
      png(sprintf('%s/%s/full/%s_%s.png', bsdir,
                  save_name, save_name, "bs_predictor_correlation"),
          3236, 2000, units = 'px', res = 300)
      
      ## set margins
      par(mar   = c(3, 3, 5, 3),
          oma   = c(1.5, 1.5, 1.5, 1.5))
      
      ## Add detail to the response plot
      chart.Correlation(bs.model@presence,
                        histogram = TRUE, pch = 19,
                        title = paste0('Reduced variable correlations for ', save_name))
      
      dev.off()
      
    } else {
      
      message("Don't run backwards selection")
      
    }
    
  }
  
}










#########################################################################################################################
## GET BACKGROUND POINTS AND THEN FIT MAXENT WITH BACKWARDS SELECTION 
#########################################################################################################################


## Variables to run an example within the backwards selection function
# name            = GBIF.spp[1]
# spp             = name
# maxent_path     = './output/maxent/HOLLOW_SPP_PROP_SAMPLE_ALL_VARS/'    ## The directory where files are saved
# maxent_dir      = 'output/maxent/HOLLOW_SPP_PROP_SAMPLE_ALL_VARS'       ## Another version of the path needed to run maxent loop
# bs_dir          = 'output/maxent/HOLLOW_SPP_PROP_SAMPLE_ALL_VARS_BS'
# outdir          = bs_dir
# features        = 'lpq'    ## name
# replicates      = 5
# cor_thr         = 0.8      ## The maximum allowable pairwise correlation between predictor variables
# pct_thr         = 5        ## The minimum allowable percent variable contribution
# k_thr           = 4        ## The minimum number of variables to be kept in the model.
# responsecurves  = TRUE     ## Response curves


## 
fit_maxent_targ_bs <- function(name,
                               maxent_path, ## location of the data
                               outdir, 
                               features, 
                               replicates, # number of cross-validation replicates
                               cor_thr, 
                               pct_thr, 
                               k_thr, 
                               responsecurves) {
  
  ########################################################################
  ## First, stop if the outdir file exists, 
  if(!file.exists(outdir)) stop('outdir does not exist :(', call. = FALSE)
  outdir_sp <- file.path(outdir, gsub(' ', '_', name))
  
  ## Also Skip species which have already been modelled: again, can this be made into an argument?
  if(!file.exists(outdir_sp)) dir.create(outdir_sp)
  features <- unlist(strsplit(features, ''))
  
  #####################################################################
  ## Read in the occ and bg points from the targetted SDM step
  message('Reading previously created occurrence and background data from targetted SDM for ', name)
  save_name = gsub(' ', '_', name)
  
  occ     <- readRDS(sprintf('%s%s/%s_occ.rds',   maxent_path, save_name, save_name))
  bg      <- readRDS(sprintf('%s%s/%s_bg.rds',    maxent_path, save_name, save_name))
  swd     <- readRDS(sprintf('%s%s/swd.rds',      maxent_path, save_name))
  
  ## Then save the occ and bg data to the backwards selection directory
  saveRDS(occ, sprintf('%s/%s/%s_occ.rds',        outdir, save_name, save_name))
  saveRDS(bg,  sprintf('%s/%s/%s_bg.rds',         outdir, save_name, save_name))
  saveRDS(swd, sprintf('%s/%s/swd.rds',           outdir, save_name))
  
  #####################################################################
  ## Coerce the "species with data" (SWD) files to regular data.frames
  ## This is needed to use the simplify function 
  swd_occ     <- as.data.frame(swd_occ)
  swd_occ$lon <- NULL
  swd_occ$lat <- NULL
  swd_bg      <- as.data.frame(swd_bg)
  swd_bg$lon  <- NULL
  swd_bg$lat  <- NULL
  
  ## Need to create a species column here?
  swd_occ$searchTaxon <- name
  swd_bg$searchTaxon  <- name
  
  #####################################################################
  ## Run simplify rmaxent::simplify
  
  # Given a candidate set of predictor variables, this function identifies 
  # a subset that meets specified multicollinearity criteria. Subsequently, 
  # backward stepwise variable selection is used to iteratively drop the variable 
  # that contributes least to the model, until the contribution of each variable 
  # meets a specified minimum, or until a predetermined minimum number of 
  # predictors remains. It returns a model object for the full model, rather 
  # than a list of models as does the previous function
  m <- local_simplify(
    
    swd_occ, 
    swd_bg,
    path            = outdir, 
    species_column  = "searchTaxon",
    replicates      = replicates,  ## 5 as above
    response_curves = TRUE, 
    logistic_format = TRUE, 
    cor_thr         = cor_thr, 
    pct_thr         = pct_thr, 
    k_thr           = k_thr, 
    features        = features,  ## LPQ
    quiet           = FALSE)
  
  ## Read the model in because it's tricky to index
  bs.model <- readRDS(sprintf('%s/%s/full/maxent_fitted.rds', outdir,  save_name)) 
  
  #####################################################################
  ## Save the chart corrleation file too for the training data set
  par(mar   = c(3, 3, 5, 3),
      oma   = c(1.5, 1.5, 1.5, 1.5))
  
  ## Add detail to the response plot
  chart.Correlation(bs.model@presence,
                    histogram = TRUE, pch = 19) 
  
  png(sprintf('%s/%s/full/%s_%s.png', outdir,
              save_name, save_name, "bs_predictor_correlation"),
      3236, 2000, units = 'px', res = 300)
  
  ## set margins
  par(mar   = c(3, 3, 5, 3),
      oma   = c(1.5, 1.5, 1.5, 1.5))
  
  ## Add detail to the response plot
  chart.Correlation(bs.model@presence,
                    histogram = TRUE, pch = 19)
  
  dev.off()
  
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