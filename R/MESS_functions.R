#########################################################################################################################
## Load only the packages needed for the analysis
p <- c('ff',      'things',    'raster',        'dismo',             'sp',           'latticeExtra',          'data.table',
       'rgdal',   'rgeos',     'gdalUtils',     'rmaxent',           'readr',        'plyr',                  'dplyr',
       'tidyr',   'readr',     'rnaturalearth', 'rasterVis',         'RColorBrewer', 'latticeExtra',          'parallel',
       'ALA4R',   'stringr',   'Taxonstand',    'CoordinateCleaner', 'gsubfn',       'PerformanceAnalytics',  'utf8',
       'rvest',   'magrittr',  'devtools',      'ggplot2',           'reshape2',     'rmarkdown',             'flexdashboard', 
       'shiny',   'rgbif',     'ENMeval',       'tibble',            'ncdf4',        'Cairo',                 'taxonlookup',  
       'kgc',     'betareg',   'gridExtra',     'grid',              'lattice',      'ConR',                  'writexl',
       'sf',      'ggmap',     'DataCombine',   'exactextractr', 'mgcv')


## Require packages
sapply(p, require, character.only = TRUE)


#########################################################################################################################
## Try to run the mess maps at the same time as the map creation?
project_maxent_grids_mess = function(shp_path, aus_shp, world_shp, scen_list, 
                                     species_list, maxent_path, climate_path,
                                     grid_names, time_slice, current_grids, create_mess, nclust) {
  
  ## Read in the Australian shapefile at the top
  aus_poly = readRDS(paste0(shp_path, aus_shp)) %>%
    spTransform(ALB.CONICAL)
  
  world_poly = readRDS(paste0(shp_path, world_shp)) %>%
    spTransform(CRS.WGS.84)
  
  ###########################################
  ## First, run a loop over each scenario:    
  lapply(scen_list, function(x) {
    
    ## Create a raster stack for each of the 6 GCMs, not for each species
    s <- stack(c(sprintf('%s/20%s/%s/%s%s.tif', climate_path, time_slice, x, x, 1:19)))
    identical(projection(s), projection(aus_poly))
    
    ## Rename both the current and future environmental stack...
    ## critically important that the order of the name.....................................................................
    names(s) <- names(current_grids) <- grid_names 
    
    ########################################################################################################################
    ## Divide the 11 temperature rasters by 10: NA values are the ocean
    ## s[[1:11]] <- s[[1:11]]/10 ## That code doesn't work
    message('First, divide the raster stack for ', x, ' by 10 ')
    for(i in 1:11) {
      ## Simple loop
      message(i)
      s[[i]] <- s[[ i]]/10
      
    }
    
    ## Then apply each GCM to each species.
    ## First, check if the maxent model exists
    ## Then apply each GCM to each species
    ## Define function to then send to one or multiple cores
    maxent_predict_fun <- function(species) {
      
      #############################################
      save_name = gsub(' ', '_', species)
      if(file.exists(sprintf('%s/%s/full/maxent_fitted.rds', maxent_path, species))) {
        message('Then run maxent projections for ', species, ' under ', x, ' scenario')
        
        ## Then, check if the species projection has already been run...
        if(!file.exists(sprintf('%s/%s/full/%s_%s.tif', maxent_path, species, species, x))) {
          
          ########################################################################
          ## Now read in the SDM model, calibrated on current conditions
          ## if it was run with backwards selection, just use the full model
          if (grepl("BS", maxent_path)) {
            
            message('Read in the BS model')
            m   <- readRDS(sprintf('%s/%s/full/maxent_fitted.rds', maxent_path, species)) 
            
          } else {
            
            ## If it was run with targetted selection, index the full model
            message('Read in the full model')
            m   <- readRDS(sprintf('%s/%s/full/maxent_fitted.rds', maxent_path, species))$me_full
            
          }
          
          ## Note that the backwards selection and targetted algorithms output slightly different 
          ## swd objects, df and spdf. best to make these the same
          ## .....................................................................................
          swd <- as.data.frame(readRDS(sprintf('%s%s/swd.rds',    maxent_path, species, species)))
          occ <- readRDS(sprintf('%s%s/%s_occ.rds', maxent_path, species, save_name)) %>%
            spTransform(ALB.CONICAL)  
          
          ## Create a file path for the current raster prediction
          f_current  <- sprintf('%s%s/full/%s_current.tif', maxent_path, species, species)
          
          ## If the current raster prediction has not been run, run it
          if(!file.exists(f_current)) {
            
            ## Report which prediction is in progress :: m$me_full, m$me_full@presence
            message('Running current prediction for ', species) 
            
            pred.current <- rmaxent::project(
              m, current_grids[[colnames(m@presence)]])$prediction_logistic
            writeRaster(pred.current, f_current, overwrite = TRUE)
            
          } else {
            message('Use existing prediction for ', species) 
            pred.current = raster(sprintf('%s/%s/full/%s_current.tif',
                                          maxent_path, species, species))
          }
          
          #####################################################################
          ## Report current mess map in progress
          MESS_dir = sprintf('%s%s/full/%s', 
                             maxent_path, species, 'MESS_output')
          f_mess_current = sprintf('%s/%s%s.tif', MESS_dir, species, "_current_mess")
          
          if(create_mess == "TRUE" && !file.exists(f_mess_current)) {
            message('Running current mess map for ', species)
            
            ## Set the names of the rasters to match the occ data, and subset both
            sdm_vars             = names(m@presence)
            current_grids        = subset(current_grids, sdm_vars)
            swd                  = swd [,sdm_vars]
            identical(names(swd), names(current_grids))
            
            ## Create a map of novel environments for current conditions.
            ## This similarity function only uses variables (e.g. n bioclim), not features
            mess_current  <- similarity(current_grids, swd, full = TRUE)
            novel_current <- mess_current$similarity_min < 0  ##   All novel environments are < 0
            novel_current[novel_current==0] <- NA             ##   0 values are NA
            
            ##################################################################
            ## Write out the current mess maps - 
            ## create a new folder for the mess output - we are going to print it to the maps
            if(!dir.exists(MESS_dir)) {
              message('Creating MESS directory for ', species) 
              dir.create(MESS_dir)
              
            } else {
              
              message(species, ' MESS directory already created') 
              
            }
            
            ## Then write the mess output to a directory inside the 'full' maxent folder
            writeRaster(mess_current$similarity_min, sprintf('%s/%s%s.tif', MESS_dir, species, "_current_mess"), 
                        overwrite = TRUE)
            
            ##################################################################
            ## Create a PNG file of MESS maps for each maxent variable
            ## raster_list  = unstack(mess_current$similarity) :: list of environmental rasters
            ## raster_names = names(mess_current$similarity)   :: names of the rasters
            message('Creating mess maps of each current environmental predictor for ', species)
            mapply(function(raster, raster_name) {
              
              ## Create a level plot of MESS output for each predictor variable, for each species
              p <- levelplot(raster, margin = FALSE, scales = list(draw = FALSE),
                             at = seq(minValue(raster), maxValue(raster), len = 100),
                             colorkey = list(height = 0.6), 
                             main = gsub('_', ' ', sprintf('Current_mess_for_%s (%s)', raster_name, species))) +
                
                latticeExtra::layer(sp.polygons(aus_poly), data = list(aus_poly = aus_poly))  ## need list() for polygon
              
              p <- diverge0(p, 'RdBu')
              f <- sprintf('%s/%s%s%s.png', MESS_dir, species, "_current_mess_", raster_name)
              
              png(f, 8, 8, units = 'in', res = 300, type = 'cairo')
              print(p)
              dev.off()
              
            }, unstack(mess_current$similarity), names(mess_current$similarity))
            
            ## Write the raster of novel environments to the MESS sub-directory
            if(!file.exists(sprintf('%s/%s%s.tif', MESS_dir, species, "_current_novel")))  {
              
              message('Writing currently novel environments to file for ', species) 
              writeRaster(novel_current, sprintf('%s/%s%s.tif', MESS_dir, species, "_current_novel"), 
                          overwrite = TRUE)
              
            } else {
              
              message(species, 'Current MESS directory already created') 
              
            }
            
            ##################################################################
            ## Now mask out novel environments
            ## is.na(novel_current) is a binary layer showing 
            ## not novel [=1] vs novel [=0], 
            ## so multiplying this with hs_current will mask out novel
            hs_current_not_novel <- pred.current * is.na(novel_current) 
            
            ## This layer of currently un-novel environments can be used 
            ## for the next algorithm step, where we combine the models
            plot(hs_current_not_novel, main = save_name)
            
            ## Write out not-novel raster :: this can go to the main directory
            message('Writing currently un-novel environments to file for ', species) 
            writeRaster(hs_current_not_novel, sprintf('%s%s/full/%s%s.tif', maxent_path, 
                                                      species, species, "_current_not_novel"),
                        overwrite = TRUE)
            
          } else {
            message('Dont run current MESS maps for ', species) 
          }
          
          ########################################################################################################################
          ## Create file path for future raster doesn't exist, create it 
          f_future <- sprintf('%s/%s/full/%s_%s.tif', 
                              maxent_path, species, species, x)
          
          if(!file.exists(f_future)) {
            
            ## Report which prediction is in progress
            message('Running future maxent prediction for ', species, ' under ', x) 
            
            ## Create the future raster
            pred.future <- rmaxent::project(
              m, s[[colnames(m@presence)]])$prediction_logistic
            writeRaster(pred.future, f_future, overwrite = TRUE)
            
            #####################################################################
            ## Report future mess map in progress
            f_mess_future = sprintf('%s/%s%s%s.tif', MESS_dir, species, "_future_mess_", x)
            
            if(create_mess == "TRUE" & !file.exists(f_mess_future)) {
              message('Running future mess map for ', species, ' under ', x)
              
              ## Set the names of the rasters to match the occ data, and subset both
              ## Watch the creation of objects in each run
              sdm_vars             = names(m@presence)
              future_grids         = s
              future_grids         = subset(future_grids, sdm_vars)
              swd                  = swd [,sdm_vars]
              identical(names(swd), names(future_grids))
              
              ## Create a map of novel environments for future conditions
              ## This similarity function only uses variables (e.g. n bioclim), not features.
              ## We don't need to repeat the static layer MESS each time - just
              ## include their results here
              mess_future  <- similarity(future_grids, swd, full = TRUE)
              novel_future <- mess_future$similarity_min < 0  ##   All novel environments are < 0
              novel_future[novel_future==0] <- NA             ##   0 values are NA
              
              ##################################################################
              ## Write out the future mess maps, for all variables
              writeRaster(mess_future$similarity_min, sprintf('%s/%s%s%s.tif', MESS_dir, species, "_future_mess_", x), 
                          overwrite = TRUE)
              
              ##################################################################
              ## Create a PNG file of all the future MESS output: 
              ## raster_list  = unstack(mess_current$similarity) :: list of environmental rasters
              ## raster_names = names(mess_current$similarity)   :: names of the rasters
              message('Creating mess maps of each future environmental predictor for ', species, ' under senario ', x)
              mapply(function(raster, raster_name) {
                
                p <- levelplot(raster, margin = FALSE, scales = list(draw = FALSE),
                               at = seq(minValue(raster), maxValue(raster), len = 100),
                               colorkey = list(height = 0.6), 
                               main = gsub('_', ' ', sprintf('Future_mess_for_%s_%s (%s)', raster_name, x, species, x))) +
                  
                  latticeExtra::layer(sp.polygons(aus_poly), data = list(aus_poly = aus_poly))   ## Use this in previous functions
                
                p <- diverge0(p, 'RdBu')
                f <- sprintf('%s/%s%s%s%s%s.png', MESS_dir, species, "_future_mess_", raster_name, "_", x)
                
                png(f, 8, 8, units = 'in', res = 300, type = 'cairo')
                print(p)
                dev.off()
                
              }, unstack(mess_future$similarity), names(mess_future$similarity))
              
              ## Write the raster of novel environments to the maxent directory 
              ## The "full" directory is getting full, could create a sub dir for MESS maps
              message('Writing future novel environments to file for ',    species, ' under senario ', x) 
              writeRaster(novel_future, sprintf('%s/%s%s%s.tif', MESS_dir, species, "_future_novel_",  x), 
                          overwrite = TRUE)
              
              ##################################################################
              # mask out future novel environments 
              # is.na(novel_future) is a binary layer showing 
              # not novel [=1] vs novel [=0], 
              # so multiplying this with hs_future will mask out novel
              hs_future_not_novel <- pred.future * is.na(novel_future)
              
              ## This layer of future un-novel environments can be used 
              ## for the next algorithm step, where we combine the models
              plot(pred.future);plot(hs_future_not_novel)
              summary(pred.future,         maxsamp = 100000);
              summary(hs_future_not_novel, maxsamp = 100000)
              
              ## Write out not-novel raster
              message('Writing un-novel environments to file under ', x, ' scenario for ', species) 
              writeRaster(hs_future_not_novel, sprintf('%s%s/full/%s%s%s.tif', maxent_path, 
                                                       species, species, "_future_not_novel_", x), 
                          overwrite = TRUE)
              
            } else {
              
              message('Dont run future MESS maps for ', species, ' under scenario ',  x ) 
              
            }
            
            ###################################################################
            ## Convert binary rasters of novel climate to polygons
            ## Need to save the polygons to file ::  
            ## this can fail if no ebvironments are novel, e.g the red gum.
            ## Can we add a condtiton in poluygonizer to check if the file has data?..............................
            message('Converting raster MESS maps to polygons under ', x, ' scenario for ', species) 
            novel_current_poly <- polygonizer(sprintf('%s/%s%s.tif',   MESS_dir, species, "_current_novel"))
            novel_future_poly  <- polygonizer(sprintf('%s/%s%s%s.tif', MESS_dir, species, "_future_novel_", x))
            
            ## Re-project the shapefiles
            novel_current_poly = novel_current_poly %>%
              spTransform(ALB.CONICAL)
            
            novel_future_poly = novel_future_poly %>%
              spTransform(ALB.CONICAL)
            
            ## Now save the novel areas as shapefiles
            MESS_shp_path   = sprintf('%s%s/full/%s', 
                                      maxent_path, species, 'MESS_output')
            
            writeOGR(obj    = novel_current_poly, 
                     dsn    = sprintf('%s',  MESS_shp_path), 
                     layer  = paste0(species, "_current_novel_polygon"),
                     driver = "ESRI Shapefile", overwrite_layer = TRUE)
            
            writeOGR(obj    = novel_future_poly, 
                     dsn    = sprintf('%s',  MESS_shp_path), 
                     layer  = paste0(species, "_future_novel_polygon_", x),
                     driver = "ESRI Shapefile", overwrite_layer = TRUE)
            
            ## Create a SpatialLines object that indicates novel areas (this will be overlaid)            
            ## Below, we create a dummy polygon as the first list element (which is the extent
            ## of the raster, expanded by 10%), to plot on panel 1). 50 = approx 50 lines across the polygon
            message('Creating polygon list under ', x, ' scenario for ', species) 
            #  cast the objects into the sf class so we avoid issues with wrong methods being called in hatch()
            novel_hatch <- list(as(extent(pred.current)*1.1, 'SpatialPolygons'),  
                                hatch(as(novel_current_poly, 'sf'), 50),
                                hatch(as(novel_future_poly,  'sf'), 50)) 
            
            ########################################################################################################################
            ## Now create a panel of PNG files for maxent projections and MESS maps
            ## All the projections and extents need to match
            empty_ras <- init(pred.current, function(x) NA) 
            projection(novel_current_poly);projection(occ);projection(empty_ras);projection(poly)
            projection(pred.current);projection(pred.future)
            identical(extent(pred.current), extent(pred.future))
            
            ## Assign the scenario name (to use in the plot below)
            scen_name = eval(parse(text = sprintf('gcms.%s$GCM[gcms.%s$id == x]', time_slice, time_slice)))      
            
            ########################################################################################################################
            ## Use the 'levelplot' function to make a multipanel output: occurrence points, current raster and future raster
            if(create_mess == "TRUE") {
              message('Create MESS panel maps for ', species, ' under ', x, ' scenario')
              
              
              ############################################################
              ## Create level plot of current conditions including MESS                        
              png(sprintf('%s/%s/full/%s_%s.png', maxent_path, species, species, "mess_panel"),      
                  11, 4, units = 'in', res = 300)
              
              print(levelplot(stack(empty_ras,
                                    pred.current, quick = TRUE), margin = FALSE,
                              
                              ## Create a colour scheme using colbrewer: 100 is to make it continuos
                              ## Also, make it a one-directional colour scheme
                              scales      = list(draw = FALSE), 
                              at = seq(0, 1, length = 100),
                              col.regions = colorRampPalette(rev(brewer.pal(11, 'Spectral'))),
                              
                              ## Give each plot a name: the third panel is the GCM
                              names.attr = c('Australian records', 'Current'),
                              colorkey   = list(height = 0.5, width = 3), xlab = '', ylab = '',
                              main       = list(gsub('_', ' ', species), font = 4, cex = 2)) +
                      
                      ## Plot the Aus shapefile with the occurrence points for reference
                      ## Can the current layer be plotted on it's own?
                      ## Add the novel maps as vectors.              
                      latticeExtra::layer(sp.polygons(aus_poly), data = list(aus_poly = aus_poly)) +
                      latticeExtra::layer(sp.points(occ, pch = 19, cex = 0.15, 
                                                    col = c('red', 'transparent', 'transparent')[panel.number()]), data = list(occ = occ)) +
                      latticeExtra::layer(sp.polygons(h[[panel.number()]]), data = list(h = novel_hatch)))
              dev.off()
              
              
              ##################################################################################
              ## Save the global records to PNG :: try to code the colors for ALA/GBIF/INVENTORY
              occ.world <- readRDS(sprintf('%s/%s/%s_occ.rds', maxent_path, species, save_name)) %>%
                spTransform(CRS.WGS.84)
              
              message('writing map of global records for ', species)
              png(sprintf('%s/%s/full/%s_%s.png', maxent_path, save_name, save_name, "global_occ_records"),
                  16, 10, units = 'in', res = 500)
              
              ## Add land
              plot(world_poly, #add = TRUE, 
                   lwd = 0.01, asp = 1, col = 'grey', bg = 'sky blue')
              
              ## Add points
              points(subset(occ.world, SOURCE == "GBIF"), 
                     pch = ".", cex = 3.3, cex.lab = 3, cex.main = 4, cex.axis = 2, 
                     xlab = "", ylab = "", asp = 1,
                     col = "orange", 
                     legend(7,4.3, unique(occ.world$SOURCE), col = "orange", pch = 1))
              
              points(subset(occ.world, SOURCE == "ALA"), 
                     pch = ".", cex = 3.3, cex.lab = 3, cex.main = 4, cex.axis = 2, 
                     xlab = "", ylab = "", asp = 1,
                     col = "blue", 
                     legend(7,4.3, unique(occ.world$SOURCE), col = "blue", pch = 1))
              
              points(subset(occ.world, SOURCE == "INVENTORY"), 
                     pch = ".", cex = 3.3, cex.lab = 3, cex.main = 4, cex.axis = 2, 
                     xlab = "", ylab = "", asp = 1,
                     col = "red", 
                     legend(7,4.3, unique(occ.world$SOURCE), col = "red", pch = 1))
              
              title(main = list(paste0(gsub('_', ' ', species), ' global SDM records'), font = 4, cex = 2),
                    cex.main = 4,   font.main = 4, col.main = "black")
              
              dev.off()
              
              
              ############################################################
              ## Create level plot of scenario x, including MESS                        
              png(sprintf('%s/%s/full/%s_%s.png', maxent_path, species, species, x),      
                  11, 4, units = 'in', res = 300)
              
              ## Create a panel of the Australian occurrences, the current layer and the future layer
              print(levelplot(stack(empty_ras,
                                    pred.current, 
                                    pred.future, quick = TRUE), margin = FALSE,
                              
                              ## Create a colour scheme using colbrewer: 100 is to make it continuos
                              ## Also, make it a one-directional colour scheme
                              scales      = list(draw = FALSE), 
                              at = seq(0, 1, length = 100),
                              col.regions = colorRampPalette(rev(brewer.pal(11, 'Spectral'))),
                              
                              ## Give each plot a name: the third panel is the GCM
                              names.attr = c('Australian records', 'Current', sprintf('%s, 20%s, RCP8.5', scen_name, time_slice)),
                              colorkey   = list(height = 0.5, width = 3), xlab = '', ylab = '',
                              main       = list(gsub('_', ' ', species), font = 4, cex = 2)) +
                      
                      ## Plot the Aus shapefile with the occurrence points for reference
                      ## Can the current layer be plotted on it's own?
                      ## Add the novel maps as vectors.              
                      latticeExtra::layer(sp.polygons(aus_poly), data = list(aus_poly = aus_poly)) +
                      latticeExtra::layer(sp.points(occ, pch = 19, cex = 0.15, 
                                                    col = c('red', 'transparent', 'transparent')[panel.number()]), data = list(occ = occ)) +
                      latticeExtra::layer(sp.polygons(h[[panel.number()]]), data = list(h = novel_hatch)))
              dev.off()
              
            } else {
              
              ############################################################
              ## OR Create level plot without MESS maps                        
              message('Create panel maps for ', species, ' under ', x, ' scenario')
              png(sprintf('%s/%s/full/%s_%s.png', maxent_path, species, species, x),      
                  11, 4, units = 'in', res = 300)
              
              ## Create levelplot without MESS maps
              print(levelplot(stack(empty_ras,
                                    pred.current, 
                                    pred.future, quick = TRUE), margin = FALSE,
                              
                              ## Create a colour scheme using colbrewer: 100 is to make it continuos
                              ## Also, make it a one-directional colour scheme
                              scales      = list(draw = FALSE), 
                              at = seq(0, 1, length = 100),
                              col.regions = colorRampPalette(rev(brewer.pal(11, 'Spectral'))),
                              
                              ## Give each plot a name: the third panel is the GCM
                              names.attr = c('Australian records', 'Current', sprintf('%s, 20%s, RCP8.5', scen_name, time_slice)),
                              colorkey   = list(height = 0.5, width = 3), xlab = '', ylab = '',
                              main       = list(gsub('_', ' ', species), font = 4, cex = 2)) +
                      
                      ## Plot the Aus shapefile with the occurrence points for reference
                      ## Can the points be made more legible for both poorly and well recorded species?
                      ## layer(sp.polygons(aus_albers), data = list(aus_albers = aus_albers))
                      latticeExtra::layer(sp.polygons(aus_poly), data = list(aus_poly = aus_poly)) +
                      latticeExtra::layer(sp.points(occ, pch = 19, cex = 0.15, 
                                                    col = c('red', 'transparent', 'transparent')[panel.number()]), data = list(occ = occ)))
              dev.off()                     
              
              message('Dont run MESS maps for ', species) 
              
            }
            
          }
          
        } else {
          
          message(species, ' ', x, ' skipped - prediction already run') ## Skip species which have already been projected
          
        }
        
      } else {
        
        message(species, ' ', x, ' skipped - SDM not yet run')          ## Skip species with no existing maxent model
        
      }
      
    }
    
    if (nclust==1) {
      lapply(species_list, maxent_predict_fun)  
    } else {
      cl <- makeCluster(nclust)
      clusterExport(cl, c(
        'shp_path',    'aus_shp',       'world_shp',   'scen_list',   'species_list', 
        'maxent_path', 'climate_path',  'grid_names',  'time_slice',  'current_grids',  
        'create_mess', 'hatch',        
        'polygonizer', 'sdm.predictors', 'aus.grids.current', 'nclust',
        'sdm.select',  'diverge0'))
      
      # shp_path, aus_shp, world_shp, scen_list, 
      # species_list, maxent_path, climate_path, 
      # grid_names, time_slice, current_grids, create_mess, nclust
      
      clusterEvalQ(cl, {
        library(rmaxent)
        library(sp)
        library(raster)
        library(rasterVis)
        library(latticeExtra)
      })
      parLapply(cl, species_list, myfun)  
    }
    
    
  })
  
}




#########################################################################################################################
## This function creates hatching on a mess map
hatch <- function(x, density) {
  # x: polygon object (SpatialPolgyons* or sf)
  # density: approx number of lines to plot
  require(sp)
  require(raster)
  e <- extent(x)
  w <- diff(e[1:2])
  x1 <- seq(xmin(e), xmax(e)+w, length.out=floor(density*2))
  x0 <- seq(xmin(e)-w, xmax(e), length.out=floor(density*2))
  y0 <- rep(ymin(e), floor(density*2))
  y1 <- rep(ymax(e), floor(density*2))
  ll <- spLines(mapply(function(x0, y0, x1, y1) {
    rbind(c(x0, y0), c(x1, y1))
  }, x0, y0, x1, y1, 
  SIMPLIFY=FALSE))  
  if(is(x, 'sf')) {
    require(sf)
    ll <- st_as_sf(ll)
    st_crs(ll) <- st_crs(x)
    st_intersection(ll, x)
  } else {
    proj4string(ll) <- proj4string(x)
    intersect(ll, x)
  }
}




#########################################################################################################################
## Plot a rasterVis::levelplot with a colour ramp diverging around zero
diverge0 <- function(p, ramp) {
  # p: a trellis object resulting from rasterVis::levelplot
  # ramp: the name of an RColorBrewer palette (as character), a character 
  #       vector of colour names to interpolate, or a colorRampPalette.
  require(RColorBrewer)
  require(rasterVis)
  if(length(ramp)==1 && is.character(ramp) && ramp %in% 
     row.names(brewer.pal.info)) {
    ramp <- suppressWarnings(colorRampPalette(brewer.pal(11, ramp)))
  } else if(length(ramp) > 1 && is.character(ramp) && all(ramp %in% colors())) {
    ramp <- colorRampPalette(ramp)
  } else if(!is.function(ramp)) 
    stop('ramp should be either the name of a RColorBrewer palette, ', 
         'a vector of colours to be interpolated, or a colorRampPalette.')
  rng <- range(p$legend[[1]]$args$key$at)
  s <- seq(-max(abs(rng)), max(abs(rng)), len=1001)
  i <- findInterval(rng[which.min(abs(rng))], s)
  zlim <- switch(which.min(abs(rng)), `1`=i:(1000+1), `2`=1:(i+1))
  p$legend[[1]]$args$key$at <- s[zlim]
  p[[grep('^legend', names(p))]][[1]]$args$key$col <- ramp(1000)[zlim[-length(zlim)]]
  p$panel.args.common$col.regions <- ramp(1000)[zlim[-length(zlim)]]
  p
}





#########################################################################################################################
## This function turns a raster into a polygon 
polygonizer <- function(x, outshape=NULL, pypath=NULL, readpoly=TRUE, 
                        fillholes=FALSE, aggregate=FALSE, 
                        quietish=TRUE) {
  # x: an R Raster layer, or the file path to a raster file recognised by GDAL 
  # outshape: the path to the output shapefile (if NULL, a temporary file will 
  #           be created) 
  # pypath: the path to gdal_polygonize.py or OSGeo4W.bat (if NULL, the function 
  #         will attempt to determine the location)
  # readpoly: should the polygon shapefile be read back into R, and returned by
  #           this function? (logical) 
  # fillholes: should holes be deleted (i.e., their area added to the containing
  #            polygon)
  # aggregate: should polygons be aggregated by their associated raster value?
  # quietish: should (some) messages be suppressed? (logical)
  if (isTRUE(readpoly) || isTRUE(fillholes)) require(rgdal)
  if (isTRUE(aggregate)) require(rgeos)
  if (is.null(pypath)) {
    cmd <- Sys.which('OSGeo4W.bat')
    pypath <- 'gdal_polygonize'
    if(cmd=='') {
      cmd <- 'python'
      pypath <- Sys.which('gdal_polygonize.py')
      if (!file.exists(pypath)) 
        stop("Could not find gdal_polygonize.py or OSGeo4W on your system.") 
    }
  }
  if (!is.null(outshape)) {
    outshape <- sub('\\.shp$', '', outshape)
    f.exists <- file.exists(paste(outshape, c('shp', 'shx', 'dbf'), sep='.'))
    if (any(f.exists)) 
      stop(sprintf('File already exists: %s', 
                   toString(paste(outshape, c('shp', 'shx', 'dbf'), 
                                  sep='.')[f.exists])), call.=FALSE)
  } else outshape <- tempfile()
  if (is(x, 'Raster')) {
    require(raster)
    writeRaster(x, {f <- tempfile(fileext='.tif')})
    rastpath <- normalizePath(f)
  } else if (is.character(x)) {
    rastpath <- normalizePath(x)
  } else stop('x must be a file path (character string), or a Raster object.')
  
  system2(cmd, args=(
    sprintf('"%s" "%s" %s -f "ESRI Shapefile" "%s.shp"', 
            pypath, rastpath, ifelse(quietish, '-q ', ''), outshape)))
  
  if(isTRUE(aggregate)||isTRUE(readpoly)||isTRUE(fillholes)) {
    shp <- readOGR(dirname(outshape), layer=basename(outshape), 
                   verbose=!quietish)    
  } else return(NULL)
  
  if (isTRUE(fillholes)) {
    poly_noholes <- lapply(shp@polygons, function(x) {
      Filter(function(p) p@ringDir==1, x@Polygons)[[1]]
    })
    pp <- SpatialPolygons(mapply(function(x, id) {
      list(Polygons(list(x), ID=id))
    }, poly_noholes, row.names(shp)), proj4string=CRS(proj4string(shp)))
    shp <- SpatialPolygonsDataFrame(pp, shp@data)
    if(isTRUE(aggregate)) shp <- aggregate(shp, names(shp))
    writeOGR(shp, dirname(outshape), basename(outshape), 
             'ESRI Shapefile', overwrite=TRUE)
  }
  if(isTRUE(aggregate) & !isTRUE(fillholes)) {
    shp <- aggregate(shp, names(shp))
    writeOGR(shp, dirname(outshape), basename(outshape), 
             'ESRI Shapefile', overwrite=TRUE)
  }
  ifelse(isTRUE(readpoly), return(shp), return(NULL))
}





