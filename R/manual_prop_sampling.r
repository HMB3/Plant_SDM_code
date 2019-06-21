## Ale CSIRO variables....................................................................................
## Here are the argumetns needed to run the targetted background selection SDMs inside the function itself
rm(bg);rm(occ);rm(o)

spp                     = GBIF.spp[2]
occ                     = subset(SDM.SPAT.OCC.BG, searchTaxon == spp)
bg                      = subset(SDM.SPAT.OCC.BG, searchTaxon != spp)
sdm.predictors          = sdm.select 
name                    = spp 
outdir                  = maxent_dir
template.raster         = template.raster.1km   ## 1km, 5km, 10km
min_n                   = 20
max_bg_size             = 70000 
Koppen                  = Koppen_1975_1km
background_buffer_width = 200000
shapefiles              = TRUE
features                = 'lpq'
replicates              = 5
responsecurves          = TRUE


## Create directories
outdir <- maxent_dir
dir_name = file.path(maxent_path, gsub(' ', '_', spp))
dir.create(dir_name)
write.csv(data.frame(), file.path(dir_name, "in_progress.txt"))


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

if(length(setdiff(features, c('l', 'p', 'q', 'h', 't'))) > 1)
  stop("features must be a vector of one or more of ','l', 'p', 'q', 'h', and 't'.")


## 
sample_1 <- function(n, ) {
  
  
  
  
}

## Create a buffer of xkm around the occurrence points
buffer <- aggregate(gBuffer(occ, width = background_buffer_width, byid = TRUE))


## Check the points :: could output this as a PNG
aus.mol <- AUS %>%
        spTransform(CRS(sp_epsg54009))
plot(aus.mol)
plot(buffer, add = TRUE, col = "red")
plot(occ, add = TRUE, col = "blue")

#####################################################################
## Get unique cell numbers for species occurrences from the template
## raster at xkm
cells <- cellFromXY(template.raster, occ)
length(cells);dim(occ)

## Clean out duplicate cells and NAs (including points outside extent of predictor data)
## Note this will get rid of a lot of duplicate records not filtered out by GBIF columns, etc.
not_dupes <- which(!duplicated(cells) & !is.na(cells))
occ   <- occ[not_dupes, ]
cells <- cells[not_dupes]
message(nrow(occ), ' occurrence records (unique cells).')

## remove cells that have already been selected? 
message(name, ' creating background cells')
bg_buff  <- sp::over(bg, buffer)
bg       <- bg[which(!is.na(bg_buff )), ]
bg_cells <- cellFromXY(template.raster, bg)
length(bg_cells)


## Clean out duplicates and NAs (including points outside extent of predictor data)
bg_not_dupes <- which(!duplicated(bg_cells) & !is.na(bg_cells))
bg   <- bg[bg_not_dupes, ]
bg_cells <- bg_cells[bg_not_dupes]


## Crop the Kopppen raster to the extent of the occurrences, and snap it
message(name, ' intersecting background cells with Koppen zones')
Koppen_crop        <- crop(Koppen, occ, snap = 'out')

## Only extract and match those cells that overlap between the ::
## 1). cropped koppen zone, 
## 2). occurrences and 
## 3). background points 
## This section is not working at 5km or 10km resolution for template.raster
message(xres(template.raster), ' metre cell size for template raster')
message(xres(Koppen), ' metre cell size for Koppen raster')
zones               <- raster::extract(Koppen_crop, occ)

cells_in_zones_crop <- Which(Koppen_crop %in% zones, cells = TRUE)
cells_in_zones      <- cellFromXY(Koppen, xyFromCell(Koppen_crop, cells_in_zones_crop))

bg_cells            <- intersect(bg_cells, cells_in_zones)  ## this is 0 for 5km 
i                   <- cellFromXY(template.raster, bg)
bg                  <- bg[which(i %in% bg_cells), ]
bg.ind              <- !(i  %in% exclude_list)  #  get indices that are not in some exclusion list
bg                  <- bg[bg.ind, ]


## Exclude the selected list from the INV selection 


## For some species, we have the problem that the proportion of ALA/INV data is 
## very different in the occurrence vs the bg records. 
## This could be caused by the koppen restriction, etc.
message(dim(bg), ' Initial background points')
table(bg$SOURCE)
bg.df = as.data.frame(bg)
round(with(bg.df, table(SOURCE)/sum(table(SOURCE))), 3)


#####################################################################
## Get the proportion of occurrence records from each source 
occ.df      = as.data.frame(occ)
source.prop = round(with(occ.df, table(SOURCE)/sum(table(SOURCE))), 3)
ala.prop    = sum(source.prop["ALA"], source.prop["GBIF"], na.rm = TRUE)
inv.prop    = source.prop["INVENTORY"]


## Reduce background sample, if it's larger than max_bg_size
## same algorithm for all data, if there is a different equation
if (nrow(bg) > max_bg_size) {
  
  ## If there is only one source of occurrence data, use the random sampling
  message(nrow(bg), ' target species background records for ', spp, 
          ', reduced to random ', max_bg_size, ' using random points from :: ', unique(bg$SOURCE))
  bg.samp <- bg[sample(nrow(bg), max_bg_size), ]
  
  ## To be extra thorough, we might need to sample this one proportionally too.....
  bg.samp.df = as.data.frame(bg.samp)
  round(with(bg.samp.df, table(SOURCE)/sum(table(SOURCE))), 3)
  
} else {
  
  message(nrow(bg), ' target species background records for ', spp, 
          ' using all points from :: ', unique(bg$SOURCE))
  bg.samp <- bg
  
}


if (unique(occ$SOURCE) >= 2 && 
    unique(bg$SOURCE)  >= 2 &&
    "INVENTORY" %in% unique(occ$SOURCE) &&
    "INVENTORY" %in% unique(bg$SOURCE) == "TRUE") {
  
  ## Sample background records from ALA/GBIF and INVENTORY categories in proportion with 
  ## the number of records from each category in the occurrence file.........
  message(nrow(bg.samp), ' target species background records for ', spp, 
          ', using proportional samples from :: ', unique(bg$SOURCE))
  
  ## Get the ALA
  bg.ala = bg.samp[ sample( which(bg.samp$SOURCE == "ALA" | bg.samp$SOURCE == "GBIF"), 
                            ala.prop*(nrow(subset(bg.samp, SOURCE == "ALA" | SOURCE == "GBIF")))), ]
  
  ## If the proportion of inv records in the occ data * the number of background records 
  ## is < than the no. of inventory records in the background data, 
  ## get the proportion of occ inv records * the no. of inv bg records
  if(inv.prop*(nrow(bg.samp)) > nrow(subset(bg.samp, SOURCE == "INVENTORY"))) {
    
    bg.inv  = bg.samp[ sample( which(bg.samp$SOURCE == "INVENTORY"), 
                               inv.prop*(nrow(subset(bg.samp, SOURCE == "INVENTORY")))), ]
    
  } else {
    
    ## Otherwise, get the inv.prop * the number of background points
    bg.inv  = bg.samp[ sample( which(bg.samp$SOURCE == "INVENTORY"), 
                               inv.prop*(nrow(bg.samp))), ]
    
  }
  
  ## Then, combine the samples from the ALA/GBIF and INV sources
  bg.samp    = rbind(bg.ala, bg.inv)
  bg.samp.df = as.data.frame(bg.samp)
  round(with(bg.samp.df , table(SOURCE)/sum(table(SOURCE))), 3)
  
} else {
  ## Get the background records from any source
  message(nrow(bg.samp), ' target species background records for ', spp, 
          ' using random sample from :: ', unique(bg$SOURCE))
  
}


save_name = gsub(' ', '_', name)
if(shapefiles) {
  
  suppressWarnings({
    
    message(name, ' writing occ and bg shapefiles')
    writeOGR(SpatialPolygonsDataFrame(b, data.frame(ID = seq_len(length(b)))),
             outdir_sp, paste0(save_name, '_bg_buffer'), 'ESRI Shapefile', overwrite_layer = TRUE)
    writeOGR(bg.samp,  outdir_sp, paste0(save_name, '_bg'),   'ESRI Shapefile', overwrite_layer = TRUE)
    writeOGR(occ,           outdir_sp, paste0(save_name, '_occ'),  'ESRI Shapefile', overwrite_layer = TRUE)
    
  })
  
}


## Also save the background and occurrence points as .rds files
saveRDS(bg.samp,  file.path(outdir_sp, paste0(save_name, '_bg.rds')))
saveRDS(occ,      file.path(outdir_sp, paste0(save_name, '_occ.rds')))

#####################################################################
## SWD = species with data. This line samples the environmental 
## variables used in the model at all the occ and bg points
swd_occ <- occ[, sdm.predictors]
saveRDS(swd_occ, file.path(outdir_sp, paste0(save_name,'_occ_swd.rds')))

swd_bg <- bg.samp[, sdm.predictors]
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

## Now set the features to be used by maxent ::
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
png(sprintf('%s/%s/full/%s_%s.png', outdir,
            save_name, save_name, "predictor_correlation"),
    3236, 2000, units = 'px', res = 300)

## set margins
par(mar   = c(3, 3, 5, 3),
    oma   = c(1.5, 1.5, 1.5, 1.5))

## Add detail to the response plot
chart.Correlation(swd_occ@data,
                  histogram = TRUE, pch = 19) 

## Finish the device
dev.off()





#########################################################################################################################
## 3). MOVE MAXENT FOLDERS TO NEW LOCATION FOR EACH RUN
#########################################################################################################################


# ## Create a list of folders for this run of species:EG hollow bearing species 
run_path         <- "./output/maxent/HOLLOW_SPP"
run_pat          <- map_spp
run_pat          <- paste(run_pat, sep = "", collapse = "|")
maxent_run_list  <- list.files(maxent_path, pattern = run_pat, full.names = TRUE)


## Save the list of maxent directories to file
lapply(maxent_run_list, write, sprintf("%s%s_maxent_folders.txt", maxent_path, save_run), append = TRUE)
maxent_file_list = sprintf("%s_maxent_folders.txt", save_run)


# cd H:/green_cities_sdm/output/maxent/SUA_TREES_ANALYSIS
# cat maxent_file_list | xargs -J % cp % run_path
# cat HOLLOW_SPP_maxent_folders.txt | xargs -J % cp % HOLLOW_SPP


# Copy or move these files to a specific folder
# Then you search for a file pattern in that directory
# This is very slow, would be better done in unix, etc. But
file.copy(from      = maxent_run_list,
          to        = run_path,
          overwrite = FALSE,
          recursive = TRUE,
          copy.mode = TRUE)




