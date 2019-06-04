#########################################################################################################################
############################################# FUNCTIONS FOR MAPPING SDMS ################################################ 
#########################################################################################################################


#########################################################################################################################
## PROJECT MAXENT MODELS
#########################################################################################################################


## flag issues with.......................................................................................................

## Add if/else for raster divide........................................................................................
## Add argument for just the MESS panel....................................................................................

#########################################################################################################################
## E.G. arguments to run the algorithm inside the function 
# shp_path      = "./data/base/CONTEXTUAL/" ## Path for shapefile
# aus_shp       = "aus_states.rds"          ## Shapefile e.g. Australian states
# world_shp     = "LAND_world.rds"          ## World shapefile
# 
# x             = scen_2030[3]                 ## List of climate scenarios
# species       = map_spp[1]                   ## List of species folders with maxent models
# maxent_path   = bs_path                   ## Output folder
# climate_path  = "./data/base/worldclim/aus/1km/bio" ## climate data
# 
# grid_names    = grid.names
# bs_names      = bs.predictors             ## names of the predictor grids
# time_slice    = 30                        ## Time period
# current_grids = aus.grids.current         ## predictor grids
# create_mess   = "TRUE"
# png_only      = "TRUE"
# nclust        = 1


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
            
            ##  cast the objects into the sf class so we avoid issues with wrong methods being called in hatch()
            novel_hatch <- list(as(extent(pred.current)*1.1, 'SpatialPolygons'),
                                hatch(as(novel_current_poly, 'sf'), 50),
                                hatch(as(novel_future_poly,  'sf'), 50))
            
            # novel_hatch <- list(as(extent(pred.current)*1.1, 'SpatialPolygons'),  
            #                     hatch(novel_current_poly, 50),
            #                     hatch(novel_future_poly, 50))
            
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
                      #latticeExtra::layer(sp.lines(h[[panel.number()]]), data = list(h = novel_hatch)))
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
                      #latticeExtra::layer(sp.polygons(h[[panel.number()]]), data = list(h = novel_hatch)))
                      latticeExtra::layer(sp.lines(h[[panel.number()]]), data = list(h = novel_hatch)))
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
## COMBINE GCM FUNCTIONS
#########################################################################################################################


## Arguments needed to run the function manually
# unit_path     = "./data/base/CONTEXTUAL/SUA/"   ## Data path for the spatial unit of analysis
# unit_shp      = "SUA_2016_AUST.rds"             ## Spatial unit of analysis - E.G. SUAs
# unit_vec      = "SUA_2016_VEC.rds"              ## Vector of rasterized unit cells
# world_shp     = "LAND_world.rds"                ## Polygon for AUS maps
# aus_shp       = "aus_states.rds"                ## Polygon for World maps
# 
# DIR           = SDM.RESULTS.DIR[5]                 ## List of directories with rasters
# species       = map_spp[5]                         ## List of species' directories
# maxent_path   = bs_path                     ## Directory of maxent results
# thresh        = percent.10.log[5]                  ## List of maxent thresholds
# time_slice    = 30                              ## Time period, eg 2030
# write_rasters = TRUE


#########################################################################################################################
## Loop over directories, species and one threshold for each, also taking a time_slice argument
SUA_cell_count = function(unit_path, unit_shp, unit_vec, aus_shp, world_shp,
                          DIR_list, species_list, 
                          maxent_path, thresholds, 
                          time_slice, write_rasters) {
  
  ###################################################################################################################
  ## Read in shapefiles: clunky, but how else will can you read in shapefiles as arguments?  
  areal_unit = readRDS(paste0(unit_path,  unit_shp)) %>%
    spTransform(ALB.CONICAL)
  
  areal_unit      = areal_unit[order(areal_unit$SUA_NAME16),] 
  areal_unit_vec  = readRDS(paste0(unit_path, unit_vec))
  
  aus_poly = readRDS(paste0(unit_path, aus_shp)) %>%
    spTransform(ALB.CONICAL)
  
  world_poly = readRDS(paste0(unit_path, world_shp)) %>%
    spTransform(CRS.WGS.84)
  
  ###################################################################################################################
  ## Loop over each directory
  lapply(DIR_list, function(DIR) { 
    
    ## And each species - although we don't want all possible combinations. use mapply in the function call
    lapply(species_list, function(species) {
      
      ###################################################################################################################
      ## Create a list of the rasters in each species directory for each time period, then take the mean
      message('Running summary of SDM predictions within SUAs for ', species, ' using ', names(areal_unit)[1], " shapefile")
      message('Calcualting mean of 20', time_slice, ' GCMs for ', species)
      
      ## Check if the mean GCM raster exists
      f_mean = sprintf('%s/%s/full/%s_20%s_suitability_mean.tif', maxent_path, species, species, time_slice)
      
      ## The mean of the GCMs doesn't exist, create it
      if(!file.exists(f_mean)) { 
        
        ## Read in all the habitat suitability rasters for each time period which are _not_ novel
        #raster.list       = list.files(as.character(DIR), pattern = sprintf('bi%s.tif$', time_slice), full.names = TRUE)  
        raster.list       = list.files(as.character(DIR), pattern = 'future_not_novel', full.names = TRUE) 
        suit              = stack(raster.list)
        suit.list         = unstack(suit)
        combo_suit_mean   = mean(suit)
        
        writeRaster(combo_suit_mean , sprintf('%s/%s/full/%s_20%s_suitability_mean.tif', 
                                              maxent_path, species, species, time_slice), overwrite = TRUE)
        
      } else {
        
        ## Create another level, without the mean calculation
        raster.list = list.files(as.character(DIR), pattern = sprintf('bi%s.tif$', time_slice), full.names = TRUE)  
        suit        = stack(raster.list)
        suit.list   = unstack(suit)
        
      }
      
      #########################################################################################################################
      ## Then, create rasters that meet habitat suitability criteria thresholds, determined by the rmaxent function
      for (thresh in thresholds) {
        
        ## Check if the SUA summary table exists
        SUA_file =   sprintf('%s/%s/full/%s_20%s_%s%s.tif', maxent_path,
                             species, species, time_slice, "gain_loss_", thresh)
        
        ## If the gain/loss raster already exists, skip the calculation
        if(file.exists(SUA_file)) { 
          
          message(species, ' ', 20, time_slice, ' SUA aggregation already exists, skip')
          next
          
        }
        
        ## Print the species being analysed
        message('doing ', species, ' | Logistic > ', thresh, ' for 20', time_slice)
        
        ## Read in the current suitability raster :: get the current_not_novel raster
        f_current <- raster(sprintf('%s/%s/full/%s_current_not_novel.tif', 
                                    maxent_path, species, species))
        
        ## First, create a simple function to threshold each of the rasters in raster.list,
        ## Then apply this to just the current suitability raster. 
        thresh_greater       = function (x) {x > thresh}
        current_suit_thresh  = thresh_greater(f_current)
        
        #########################################################################################################################
        ## First, calculate the cells which are greater that the: 
        ## Maximum training sensitivity plus specificity Logistic threshold
        message('Running thresholds for ', species, ' | 20', time_slice, ' combined suitability > ', thresh)
        
        suit_ras1_thresh   = thresh_greater(suit.list[[1]])   ## Abbreviate these...
        suit_ras2_thresh   = thresh_greater(suit.list[[2]])
        suit_ras3_thresh   = thresh_greater(suit.list[[3]])
        suit_ras4_thresh   = thresh_greater(suit.list[[4]])   
        suit_ras5_thresh   = thresh_greater(suit.list[[5]])
        suit_ras6_thresh   = thresh_greater(suit.list[[6]])
        
        #########################################################################################################################
        ## Then sum them up: All the threshholds
        combo_suit_thresh   =  Reduce("+", list(suit_ras1_thresh, suit_ras2_thresh, suit_ras3_thresh,
                                                suit_ras4_thresh, suit_ras5_thresh, suit_ras6_thresh))
        
        #########################################################################################################################
        ## For each species, create a binary raster with cells > 4 GCMs above the maxent threshold = 1, and cells with < 4 GCMs = 0. 
        message('Calculating change for ', species, ' | 20', time_slice, ' combined suitability > ', thresh)
        
        ## Functions for thresholding rasters
        band_4           <- function(x) {ifelse(x >=  4, 1, 0) }
        combo_suit_4GCM  <- calc(combo_suit_thresh, fun = band_4)
        
        #########################################################################################################################
        ## Now create a raster of the gain, loss and stable
        ## Create a raster stack of the current and future rasters
        message ("Counting cells lost/gained/stable/never suitable, both across AUS and per SUA")
        
        #########################################################################################################################
        ## Create a table of cell counts using a raster stack of current and future data
        d <- as.data.frame(stack(current_suit_thresh, combo_suit_4GCM)[]) %>% 
          setNames(c('current', 'future')) %>% 
          mutate(SUA_CODE16 = areal_unit_vec,
                 cell_number = seq_len(ncell(current_suit_thresh))) %>% 
          as.tbl
        dim(d);summary(d)
        
        ## Then classify the cells of the raster stack into lost, gained, stable and never suitable
        d2 <- d %>% 
          na.omit %>% 
          
          mutate(lost   = current == 1 & future == 0,
                 gained = current == 0 & future == 1,
                 stable = current == 1 & future == 1,
                 never  = current == 0 & future == 0,
                 nodata = is.na(current) | is.na(future)) 
        d2$class <- apply(select(d2, lost:never), 1, which)
        dim(d2)
        
        ## Then group the cell counts by SUA
        d3 <- d2 %>% 
          group_by(SUA_CODE16) %>%
          
          summarize(CURRENT_SUITABLE = sum(current, na.rm = TRUE),
                    FUTURE_SUITABLE  = sum(future,  na.rm = TRUE),
                    LOST             = sum(lost,    na.rm = TRUE),
                    GAINED           = sum(gained,  na.rm = TRUE),
                    STABLE           = sum(stable,  na.rm = TRUE),
                    NEVER            = sum(never,   na.rm = TRUE),
                    NODAT            = sum(nodata,  na.rm = TRUE),
                    n_cells = n()) %>% 
          
          ## Then calculate change between current and future
          mutate(CHANGE    = FUTURE_SUITABLE - CURRENT_SUITABLE,
                 GAIN_LOSS = ifelse(CHANGE < 0, 'LOSS', ifelse(CHANGE > 0, 'GAIN', 'STABLE')),
                 GAIN_LOSS = ifelse(CURRENT_SUITABLE == 0 & FUTURE_SUITABLE == 0, 'NEVER', GAIN_LOSS))
        dim(d3)
        
        ## Add the species column
        d4 = d3 %>% 
          join(areal_unit@data, .) %>%
          add_column(., SPECIES = species,    .after = "AREASQKM16") %>%
          add_column(., PERIOD  = time_slice, .after = "SPECIES")    %>%
          add_column(., THRESH  = thresh,     .after = "PERIOD")
        
        #########################################################################################################################
        ## Now calculate the number of cells lost/gained/stable across Australia
        message ("Counting cells lost/gained/stable/never suitable across Australia")
        d5 <- stack(current_suit_thresh, combo_suit_4GCM)[]
        r  <- raster(current_suit_thresh)                     ## rename, too cryptic.............................................
        z  <- as.data.frame(d4)                               ## rename, too cryptic.............................................
        
        ## Then classify the raster stack to make each value (i.e. outcome) unique
        r[d5[, 1]==1 & d5[, 2]==0] <- 1  ## 1 in current raster and 0 in future = LOSS
        r[d5[, 1]==0 & d5[, 2]==1] <- 2  ## 0 in current raster and 1 in future = GAIN
        r[d5[, 1]==1 & d5[, 2]==1] <- 3  ## 1 in current raster and 1 in future = STABLE
        r[d5[, 1]==0 & d5[, 2]==0] <- 4  ## 0 in current raster and 0 in future = NEVER_SUIT
        
        ## Create a table of these values, to merge with the levels later. This avoids the problem that not all the categories will
        ## be present for all species
        change_values = data.frame("ID" = 1:4, "CHANGE" = c("LOST", "GAINED", "STABLE", "NEVER"))
        
        
        ## Now convert the raster to a factor and assign lables to the levels
        gain_loss <- as.factor(r)
        levels(gain_loss)[[1]] <- data.frame(ID = 1:4, label = c('Lost', 'Gained', 'Stable', 'Never_Suitable'))
        z <- as.data.frame(d5)
        
        ## Create a table of the gain/loss/stable :: write this to file as well
        gain_loss_table      = table(z[, 1], z[, 2])
        gain_loss_df         = as.data.frame(raster::freq(gain_loss))
        gain_loss_df$SPECIES = species
        gain_loss_df$PERIOD  = time_slice
        
        names(gain_loss_df)  = c("CHANGE", "COUNT", "SPECIES", "PERIOD")
        gain_loss_df         = gain_loss_df[, c("SPECIES", "PERIOD", "CHANGE", "COUNT")]
        
        ## Change values and remove the NA row
        gain_loss_df$CHANGE[gain_loss_df$CHANGE == 1] <- "LOST"
        gain_loss_df$CHANGE[gain_loss_df$CHANGE == 2] <- "GAINED"
        gain_loss_df$CHANGE[gain_loss_df$CHANGE == 3] <- "STABLE"
        gain_loss_df$CHANGE[gain_loss_df$CHANGE == 4] <- "NEVER_SUIT"
        gain_loss_df = head(gain_loss_df, 4)
        head(gain_loss_df)
        
        #########################################################################################################################
        ## Save the gain/loss table for the whole of Australia
        message('Writing ', species, ' gain_loss tables for 20', time_slice)
        write.csv(gain_loss_df, sprintf('%s/%s/full/%s_20%s_%s%s.csv', maxent_path,
                                        species, species, time_slice, "gain_loss_table_", thresh), row.names = FALSE)
        
        ## Save the gain/loss table
        write.csv(d4, sprintf('%s/%s/full/%s_20%s_%s%s.csv', maxent_path,
                              species, species, time_slice, "SUA_cell_count_", thresh), row.names = FALSE)
        
        #########################################################################################################################
        ## Now write the rasters
        ## If the rasters don't exist, write them for each species/threshold
        if(write_rasters == "TRUE") {
          
          ## Write the current suitability raster, thresholded using the Maximum training 
          ## sensitivity plus specificity Logistic threshold
          message('Writing ', species, ' current', ' max train > ', thresh)
          writeRaster(current_suit_thresh, sprintf('%s/%s/full/%s_%s%s.tif', maxent_path,
                                                   species, species, "current_suit_not_novel_above_", thresh), 
                      overwrite = TRUE)
          
          ## Write the combined suitability raster, thresholded using the maximum training value
          message('Writing ', species, ' | 20', time_slice, ' max train > ', thresh)
          writeRaster(combo_suit_thresh, sprintf('%s/%s/full/%s_20%s%s%s.tif', maxent_path,
                                                 species, species, time_slice, "_log_thresh_above_", thresh), 
                      overwrite = TRUE)
          
          ## Write the combined future raster with > 4 GCMs above the maximum training value
          message('Writing ', species, ' | 20', time_slice, ' 4 GCMs > ', thresh)
          writeRaster(combo_suit_4GCM, sprintf('%s/%s/full/%s_20%s%s%s.tif', maxent_path,
                                               species, species, time_slice, "_4GCMs_above_", thresh), 
                      overwrite = TRUE)
          
          ## Write out the gain/loss raster
          writeRaster(gain_loss, sprintf('%s/%s/full/%s_20%s%s%s.tif', maxent_path,
                                         species, species, time_slice, "_gain_loss_", thresh), datatype = 'INT2U', 
                      overwrite = TRUE)
          
          #########################################################################################################################
          #########################################################################################################################
          ## Create a color scheme for the gain_loss plot
          SUA.plot.cols = brewer.pal(12, "Paired")
          ## Create dataframe of colors that match the categories :;
          # 1 <- "LOST"        SUA.plot.cols[8]
          # 2 <- "GAINED"      SUA.plot.cols[4]
          # 3 <- "STABLE"      SUA.plot.cols[1]
          # 4 <- "NEVER_SUIT"  SUA.plot.cols[9]
          ## 1 = STABLE, 4 = GAIN, 8 = LOSS, 9 = NEVER
          lc_id = c(1, 2, 3, 4)
          cover_palette <- c(SUA.plot.cols[8], SUA.plot.cols[4], SUA.plot.cols[1], SUA.plot.cols[11])
          
          colors        <- data.frame(id = lc_id, color = cover_palette, stringsAsFactors = FALSE)
          csort         <- colors[order(colors$id),] 
          
          ## Create Labels for the gain_loss plots 
          gain_plot <- ratify(gain_loss)
          rat <- levels(gain_plot)[[1]]
          
          ## Use an inner join, to accomodate the different categories. EG sometimes, there are 
          ## no cells which are gained for each species, etc.
          rat[["CHANGE"]]   <- join(rat, change_values, type = "inner")$CHANGE 
          levels(gain_plot) <- rat
          
          
          #########################################################################################################################
          ## Save the gain/loss raster to PNG
          save_name = gsub(' ', '_', species)
          identical(projection(aus_poly), projection(gain_plot))
          
          message('writing gain/loss png for ', 'species')
          png(sprintf('%s/%s/full/%s_%s_%s_20%s.png', maxent_path, save_name, save_name, "gain_loss", thresh, time_slice),
              16, 10, units = 'in', res = 500)
          
          ## Could add the SUA polygons as well
          print(levelplot(gain_plot, 
                          col.regions = csort$color, 
                          xlab = NULL, ylab = NULL,
                          main       = list(paste0(gsub('_', ' ', species), ' :: ',  20, 
                                                   time_slice, ' 4GCMs > ',  thresh), font = 4, cex = 2)) +
                  latticeExtra::layer(sp.polygons(aus_poly), data = list(aus_poly = aus_poly)))
          
          dev.off()
          
        } else {
          
          message(' skip raster writing')
          
        }
        
      }
      
    })
    
  })
  
}




#########################################################################################################################
## PLOTTING FUNCTIONS
#########################################################################################################################


#########################################################################################################################
## Create a function to combine the subplots created in Figures 2 and 3
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("top", "right")) {
  
  ## Create a list of plots
  plots <- list(...)
  
  ## Postion stuff
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  
  ## Create a legend
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  
  ## Not sure
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  ## Combine the panels? 
  combined <- switch(position, 
                     "top" = arrangeGrob(do.call(arrangeGrob, gl), legend, ncol = 1,
                                         heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     
                     "right" = arrangeGrob(do.call(arrangeGrob, gl), legend, ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  ## 
  grid.newpage()
  grid.draw(combined)
  
  ## return gtable invisibly
  invisible(combined)
  
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
    raster::intersect(ll, x)
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



## Linux version
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





#########################################################################################################################
## For each scenario, calcualte the mean annual temperature and annual rainfall anomaly
calculate.anomaly = function(scen_list, time_slice, climate_path) {
  
  #######################################################################################################################
  ## Create current rasters :: can the raster creation happen just once?
  ## And read in the current data. Can the specific "current" be removed withouth makin a seconf path argument?
  message('current rasters / 10')
  
  env.grids.current       = stack(file.path(sprintf('%s/current/bio_%02d.tif', climate_path, 1:19)))
  env.grids.current[[1]]  = env.grids.current[[1]]/10
  current.bio1            = env.grids.current[[1]]
  current.bio12           = env.grids.current[[12]]
  
  ## First, run a loop over each scenario :: 
  lapply(scen_list, function(x) {
    
    ## Assign the scenario name (to use later in the plot)
    scen_name = eval(parse(text = sprintf('gcms.%s$GCM[gcms.%s$id == x]', time_slice, time_slice)))
    
    #########################################################################################################################
    ## Create a raster stack for each 2050 GCM - also an empty raster for the final plot
    s <- stack(sprintf('%s/20%s/%s/%s%s.tif', climate_path, time_slice, x, x, 1:19))
    
    ########################################################################################################################
    ## Create future rasters
    message('20', time_slice, ' rasters / 10')
    
    s[[1]]       = s[[1]]/10
    future.bio1  = s[[1]]
    future.bio12 = s[[12]]
    
    ## Now subtract the current from the future...
    temp.anomaly = future.bio1  - current.bio1
    rain.anomaly = future.bio12 - current.bio12
    
    plot(temp.anomaly, main = paste0("BIO1  anomaly ", x))
    plot(rain.anomaly, main = paste0("BIO12 anomaly ", x))
    
    ########################################################################################################################
    ## Write the rasters for each species/threshold
    message('Writing BIO1/12 anomaly rasters for 20', time_slice, ' ', x)
    
    if(!file.exists(sprintf('%s/anomalies/%s_%s.tif', climate_path, "BIO1_anomaly",  x))) {
      
      ## Temp anomalies
      writeRaster(temp.anomaly, sprintf('%s/anomalies/%s_%s.tif', climate_path, "BIO1_anomaly",  x), overwrite = TRUE)
      
    } else {
      
      message('Skipped ', "BIO1 anomaly for ", x, ', already exists')
      
    }
    
    if(!file.exists(sprintf('%s/anomalies/%s_%s.tif', climate_path, "BIO12_anomaly",  x))) {
      
      ## Rain anomalies
      writeRaster(rain.anomaly, sprintf('%s/anomalies/%s_%s.tif', climate_path, "BIO12_anomaly", x), overwrite = TRUE)
      
    } else {
      
      message('Skipped ', "BIO12 anomaly for ", x, ', already exists')
      
    }
    
  })
  
} 





#########################################################################################################################
## For each scenario, calcualte the mean annual temperature and annual rainfall anomaly
gcm.value = function(scen_list, time_slice, climate_path) {
  
  #########################################################################################################################
  ## Create current rasters :: can the raster creation happen just once?
  ## And read in the current data. Can the specific "current" be removed withouth makin a seconf path argument?
  message('current rasters divided by 10')
  
  env.grids.current       = stack(file.path(sprintf('%s/current/bio_%02d.tif', climate_path, 1:19)))
  env.grids.current[[1]]  = env.grids.current[[1]]/10
  current.bio1            = env.grids.current[[1]]
  current.bio12           = env.grids.current[[12]]
  
  ## First, run a loop over each scenario :: 
  lapply(scen_list, function(x) {
    
    ## Assign the scenario name (to use later in the plot)
    scen_name = eval(parse(text = sprintf('gcms.%s$GCM[gcms.%s$id == x]', time_slice, time_slice)))
    
    #########################################################################################################################
    ## Create a raster stack for each 2050 GCM - also an empty raster for the final plot
    s <- stack(sprintf('%s/20%s/%s/%s%s.tif', climate_path, time_slice, x, x, 1:19))
    
    ########################################################################################################################
    ## Create future rasters
    message('20', time_slice, ' rasters divided by 10')
    
    s[[1]]       = s[[1]]/10
    future.bio1  = s[[1]]
    future.bio12 = s[[12]]
    
    ## Now subtract the current from the future...
    temp.anomaly = future.bio1  - current.bio1
    rain.anomaly = future.bio12 - current.bio12
    
    ## Also can we simply indicate the range or rainfall and temperautre values for each ?
    #temp.current = 
    
    plot(temp.anomaly, main = paste0("BIO1  anomaly ", x))
    plot(rain.anomaly, main = paste0("BIO12 anomaly ", x))
    
    ########################################################################################################################
    ## Write the rasters for each species/threshold
    message('Writing BIO1/12 anomaly rasters for 20', time_slice, ' ', x)
    
    if(!file.exists(sprintf('%s/anomalies/%s_%s.tif', climate_path, "BIO1_anomaly",  x))) {
      
      ## Temp anomalies
      writeRaster(temp.anomaly, sprintf('%s/anomalies/%s_%s.tif', climate_path, "BIO1_anomaly",  x), overwrite = TRUE)
      
    } else {
      
      message('Skipped ', "BIO1 anomaly for ", x, ', already exists')
      
    }
    
    if(!file.exists(sprintf('%s/anomalies/%s_%s.tif', climate_path, "BIO12_anomaly",  x))) {
      
      ## Rain anomalies
      writeRaster(rain.anomaly, sprintf('%s/anomalies/%s_%s.tif', climate_path, "BIO12_anomaly", x), overwrite = TRUE)
      
    } else {
      
      message('Skipped ', "BIO12 anomaly for ", x, ', already exists')
      
    }
    
  })
  
} 





#########################################################################################################################
## STACK FUNCTIONS
#########################################################################################################################


# dist_change_binary <- function(from, to) {
#   
#   ## Create a raster stack
#   d <- stack(from, to)[]
#   r <- raster(from)
#   
#   ## Classify the raster stack to make each value (i.e. outcome) unique
#   r[d[, 1]==1 & d[, 2]==0] <- 1  ## 1 in current raster and 0 in future = LOSS
#   r[d[, 1]==0 & d[, 2]==1] <- 2  ## 0 in current raster and 1 in future = GAIN
#   r[d[, 1]==1 & d[, 2]==1] <- 3  ## 1 in current raster and 1 in future = STABLE
#   
#   r_f <- as.factor(r)
#   
#   levels(r_f)[[1]] <- data.frame(ID=1:3, label=c('Lost', 'Gained', 'Stable'))
#   
#   ## Now create a level plot
#   png(sprintf('%s.png', time_slice), 6, 6, units='in', res=300)
#   
#   p <- levelplot(r_f, col.regions=c('purple', 'red', 'green'), scales = list(draw = FALSE)) +
#     layer(sp.polygons(aus))
#   print(p)
#   dev.off()
#   
#   #writeRaster(r_f, sprintf('%s.tif', ), datatype='INT2U', overwrite=T)
#   
#   r_f
#   
# }


# ## Need an empty frame
# print(levelplot(stack(empty,                ## needs to have a different colour scale,
#                       combo_suit_percent,
#                       combo_suit_thresh), margin = FALSE,
#                 
#                 ## Create a colour scheme using colbrewer: 100 is to make it continuos
#                 ## Also, make it a one-directional colour scheme
#                 scales      = list(draw = FALSE), 
#                 at = seq(0, 6, length = 100),
#                 col.regions = colorRampPalette(rev(brewer.pal(11, 'Spectral'))),
#                 
#                 ## Give each plot a name
#                 names.attr = c('Aus occurrences', 
#                                sprintf('20%s GCM 10thp train omission > %s', time_slice, percent), 
#                                sprintf('20%s GCM Max train logis > %s',      time_slice, thresh)),
#                 
#                 colorkey   = list(height = 0.5, width = 3), xlab = '', ylab = '',
#                 main       = list(gsub('_', ' ', species), font = 4, cex = 2)) +
#         
#         ## Plot the Aus shapefile with the occurrence points for reference
#         ## Can we assign different shapefiles to different panels, rather than to them all?
#         
#         layer(sp.polygons(aus)) + ## sp.polygons(aus))
#         layer(sp.points(occ, pch = 20, cex = 0.4, 
#                         col = c('red', 'transparent', 'transparent')[panel.number()]), data = list(occ = occ)))
# 
# ## Finish the device
# dev.off()


#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################
