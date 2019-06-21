project_maxent_grids_mess = function(shp_path, aus_shp, world_shp, scen_list, 
                                     species_list, maxent_path, climate_path, static_path,
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
    ## They need to have exactly the same extent.
    ## Could stack all the rasters, or, keep them separate
    s <- stack(c(sprintf('%s/20%s/%s/%s%s.tif', climate_path, time_slice, x, x, 1:19)))#,  ## Dynamic path
    #list.files(static_path, '\\.tif$', full.names = TRUE)))                  ## static path
    
    identical(projection(s), projection(aus_poly))
    
    ## Rename both the current and future environmental stack...
    ## critically important that the order of the name.....................................................................
    ## So this needs to include all the predicor names :: climate, soil, etc
    
    ## Note this step is only needed if the current grids used in the their original form, rather than being renamed...............
    names(s) <- names(current_grids) <- grid_names 
    
    ####################################################################################
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
      
      ## Why is this checing at the wrong level?                               
      if(file.exists(sprintf('%s/%s/full/maxent_fitted.rds', maxent_path, species))) {
        message('Then run maxent projections for ', species, ' under ', x, ' scenario')
        
        ## Then, check if the species projection has already been run...
        if(!file.exists(sprintf('%s/%s/full/%s_future_not_novel_%s.tif', maxent_path, species, species, x))) {
          
          ########################################################################
          ## Now read in the SDM model, calibrated on current conditions
          ## if it was run with backwards selection, just use the full model
          if (grepl("BS", maxent_path)) {
            
            message('Read in the BS model')
            m   <- readRDS(sprintf('%s/%s/full/maxent_fitted.rds', maxent_path, species)) 
            
          } else {
            ## Otherwise, index the full model
            message('Read in the full model')
            m   <- readRDS(sprintf('%s/%s/full/maxent_fitted.rds', maxent_path, species))$me_full
          }
          
          ## Read in species with data and occurrence files
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
          ## Could optimise this step too......................................
          ## Doesn't need to happen every time for the static variables
          ## Could work out how to the static mess once, before looping through scenarios
          MESS_dir = sprintf('%s%s/full/%s', 
                             maxent_path, species, 'MESS_output')
          f_mess_current = sprintf('%s/%s%s.tif', MESS_dir, species, "_current_mess")
          
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
            message(species, 'Current MESS file already saved') 
          }
          
          ####################################################################################
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
          
          ####################################################################################
          ## Create file path for future raster doesn't exist, create it 
          f_future <- sprintf('%s/%s/full/%s_future_not_novel_%s.tif', 
                              maxent_path, species, species, x)
          
          if(!file.exists(f_future)) {
            
            ## Report which prediction is in progress
            message('Running future maxent prediction for ', species, ' under ', x) 
            
            ## Create the future raster
            pred.future <- rmaxent::project(
              m, s[[colnames(m@presence)]])$prediction_logistic
            writeRaster(pred.future, f_future, overwrite = TRUE)
            
            ####################################################################################
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
              
              ####################################################################################
              ## Write out the future mess maps, for all variables
              writeRaster(mess_future$similarity_min, sprintf('%s/%s%s%s.tif', MESS_dir, species, "_future_mess_", x), 
                          overwrite = TRUE)
              
              ####################################################################################
              ## Create a PNG file of all the future MESS output: 
              ## raster_list  = unstack(mess_current$similarity) :: list of environmental rasters
              ## raster_names = names(mess_current$similarity)   :: names of the rasters
              message('Creating mess maps of each future environmental predictor for ', species, ' under scenario ', x)
              mapply(function(raster, raster_name) {
                
                p <- levelplot(raster, margin = FALSE, scales = list(draw = FALSE),
                               at = seq(minValue(raster), maxValue(raster), len = 100),
                               colorkey = list(height = 0.6), 
                               main = gsub('_', ' ', sprintf('Future_mess_for_%s_%s (%s)', raster_name, x, species, x))) +
                  
                  latticeExtra::layer(sp.polygons(aus_poly), data = list(aus_poly = aus_poly))
                
                p <- diverge0(p, 'RdBu')
                f <- sprintf('%s/%s%s%s%s%s.png', MESS_dir, species, "_future_mess_", raster_name, "_", x)
                
                png(f, 8, 8, units = 'in', res = 300, type = 'cairo')
                print(p)
                dev.off()
                
              }, unstack(mess_future$similarity), names(mess_future$similarity))
              
              ## Write the raster of novel environments to the MESS maps sub-directory
              message('Writing future novel environments to file for ',    species, ' under scenario ', x) 
              writeRaster(novel_future, sprintf('%s/%s%s%s.tif', MESS_dir, species, "_future_novel_",  x), 
                          overwrite = TRUE)
              
              ####################################################################################
              ## mask out future novel environments 
              ## is.na(novel_future) is a binary layer showing 
              ## not novel [=1] vs novel [=0], 
              ## so multiplying this with hs_future will mask out novel
              hs_future_not_novel <- pred.future * is.na(novel_future)
              
              ## This layer of future un-novel environments can be used 
              ## for the next algorithm step, where we combine the models
              plot(pred.future);plot(hs_future_not_novel)
              summary(pred.future,         maxsamp = 100000);
              summary(hs_future_not_novel, maxsamp = 100000)
              
              ## Write out not-novel raster
              ## Try to set up loops so different cores aren't accessing the same files............
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
            
            ## This step was causing lots of species to fail on the HPC
            ## I think this was because _current_novel.tif was being accessed by multiple R sessions on different cores.....                        
            message('Converting raster MESS maps to polygons under ', x, ' scenario for ', species) 
            
            ## If we're on linux, use the standard polygonizer function
            if (!on_windows) {
              
              ## Also check the current MESS map has been made
              novel_current_tif_file <- sprintf('%s/%s%s.tif',   MESS_dir, species, "_current_novel")
              if(!file.exists(novel_current_poly_file)) {
                novel_current_poly <- gBuffer(polygonizer(novel_current_tif_file), width = 0)
              }
              
              novel_future_poly  <- gBuffer(polygonizer(sprintf('%s/%s%s%s.tif', 
                                                                MESS_dir, species, "_future_novel_", x)), width = 0)
            } else {
              
              ## If we're on windows, use the GDAL .bat file
              if(!file.exists(novel_current_tif_file)) {
                
                ## Also check if the future MESS map has been made                   
                novel_current_poly <- gBuffer(polygonizer_windows(novel_current_tif_file), width=0)
              }
              novel_future_poly  <- gBuffer(polygonizer_windows(sprintf('%s/%s%s%s.tif', 
                                                                        MESS_dir, species, "_future_novel_", x)), width = 0)
            }
            
            ###################################################################
            ## Create the MESS path and save shapefiles
            MESS_shp_path   = sprintf('%s%s/full/%s', 
                                      maxent_path, species, 'MESS_output')
            
            ## Check if the current MESS shapefile exists?
            novel_current_shp <- sprintf('%s/%s%s.shp',   MESS_dir, species, "_current_novel_polygon")
            if(!file.exists(novel_current_shp)) {
              
              ## Re-project the shapefiles
              novel_current_poly = novel_current_poly %>%
                spTransform(ALB.CONICAL)
              
              ## Now save the novel areas as shapefiles
              ## There is a problem with accessing the files at the same time
              writeOGR(obj    = novel_current_poly, 
                       dsn    = sprintf('%s',  MESS_shp_path), 
                       layer  = paste0(species, "_current_novel_polygon"),
                       driver = "ESRI Shapefile", overwrite_layer = TRUE)
            }
            
            ## Check if the future MESS shapefile exists? 
            novel_future_shp <- sprintf('%s/%s%s%s.shp',   MESS_dir, species, "_future_novel_polygon_", x)
            if(!file.exists(novel_current_shp)) {
              
              novel_future_poly = novel_future_poly %>%
                spTransform(ALB.CONICAL)
              writeOGR(obj    = novel_future_poly, 
                       dsn    = sprintf('%s',  MESS_shp_path), 
                       layer  = paste0(species, "_future_novel_polygon_", x),
                       driver = "ESRI Shapefile", overwrite_layer = TRUE)
            }
            
            ## Create a SpatialLines object that indicates novel areas (this will be overlaid)            
            ## Below, we create a dummy polygon as the first list element (which is the extent
            ## of the raster, expanded by 10%), to plot on panel 1). 50 = approx 50 lines across the polygon
            message('Creating polygon list under ', x, ' scenario for ', species) 
            
            ## Cast the objects into the sf class so we avoid issues with wrong methods being called in hatch()
            novel_hatch <- list(as(extent(pred.current)*1.1, 'SpatialLines'),
                                hatch(as(novel_current_poly, 'sf'), 50),
                                hatch(as(novel_future_poly,  'sf'), 50))
            
            ####################################################################################
            ## Now create a panel of PNG files for maxent projections and MESS maps
            ## All the projections and extents need to match
            empty_ras <- init(pred.current, function(x) NA) 
            projection(novel_current_poly);projection(occ);projection(empty_ras);projection(poly)
            projection(pred.current);projection(pred.future)
            identical(extent(pred.current), extent(pred.future))
            
            ## Assign the scenario name (to use in the plot below)
            scen_name = eval(parse(text = sprintf('gcms.%s$GCM[gcms.%s$id == x]', time_slice, time_slice)))      
            
            ############################################################
            ## Use the 'levelplot' function to make a multipanel output: occurrence points, current raster and future raster
            message('Create MESS panel maps for ', species, ' under ', x, ' scenario')
            
            current_mess_png = sprintf('%s/%s/full/%s_%s.png', maxent_path, species, species, "mess_panel")
            if(!file.exists(current_mess_png)) {                         
              
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
                      latticeExtra::layer(sp.lines(h[[panel.number()]]), data = list(h = novel_hatch)))
              dev.off()
              
            }
            
            ##################################################################################
            ## Save the global records to PNG :: try to code the colors for ALA/GBIF/INVENTORY
            occ.world <- readRDS(sprintf('%s/%s/%s_occ.rds', maxent_path, species, save_name)) %>%
              spTransform(CRS.WGS.84)
            
            ## If the global map of occurrence points hasn't been created, create it
            global_occ_map = sprintf('%s/%s/full/%s_%s.png', maxent_path, save_name, save_name, "global_occ_records")
            if(!file.exists(novel_current_shp)) {
              
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
              
            }
            
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
                    latticeExtra::layer(sp.lines(h[[panel.number()]]), data = list(h = novel_hatch)))
            dev.off()
            
          }
          
        } else {
          message(species, ' ', x, ' skipped - prediction already run') 
        }
      } else {
        message(species, ' ', x, ' skipped - SDM not yet run')
      }
      
    }
    
    ############################################################
    ## Check this is the best way to run parallel...............
    if (nclust==1) {
      
      lapply(species_list, maxent_predict_fun) 
      
    } else {
      
      ## Export all objects from the function call
      # shp_path, aus_shp, world_shp, scen_list, 
      # species_list, maxent_path, climate_path, 
      # grid_names, time_slice, current_grids, create_mess, nclust
      
      cl <- makeCluster(nclust)
      clusterExport(cl, c(
        'shp_path',    'aus_shp',       'world_shp',   'ALB.CONICAL',  'CRS.WGS.84',   
        'scen_list',   'species_list',  'maxent_path', 'climate_path', 'grid_names',  
        'time_slice',  'current_grids', 'create_mess', 'hatch', 'x',       
        'polygonizer', 'nclust', 'diverge0'),  envir = environment())
      
      clusterEvalQ(cl, {
        
        library(rmaxent)
        library(sp)
        library(raster)
        library(rasterVis)
        library(latticeExtra)
        library(magrittr)
        
      })
      
      message('Running project_maxent_grids_mess for ', length(species_list),
              ' species on ', nclust, ' cores for GCM ', x)
      parLapply(cl, species_list, maxent_predict_fun)  
    }
  })
} 