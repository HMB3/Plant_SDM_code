#########################################################################################################################
################################################# PLOT SPECIES NICHES ###################################################
#########################################################################################################################


## This code creates barplots, histograms, and convex-hull plots for a set of species
## This code can be used to re-do the zonal stats for the SUAs in step 9................................................. 


## Might need to add skip condtions to this..............................................................................
shp_path      = "./data/base/CONTEXTUAL/" ## Path for shapefile
aus_shp       = "aus_states.rds"          ## Shapefile, e.g. Australian states
world_shp     = "LAND_world.rds"          ## World shapefile          
x             = scen_2030[1]            ## List of climate scenarios
species       = map_spp[15]            ## List of species folders with maxent models
maxent_path   = bs_path
time_slice    = 30
nclust        = 2


#########################################################################################################################
## Try to run the mess maps at the same time as the map creation?
create_mess_pngs = function(shp_path, aus_shp, world_shp, scen_list,
                            species_list, maxent_path, time_slice, nclust) {
  
  ## Read in the Australian shapefile at the top
  aus_poly = readRDS(paste0(shp_path, aus_shp)) %>%
    spTransform(ALB.CONICAL)
  
  world_poly = readRDS(paste0(shp_path, world_shp)) %>%
    spTransform(CRS.WGS.84)
  
  ###########################################
  ## First, run a loop over each scenario:    
  lapply(scen_list, function(x) {
    
    ## Create a raster stack for each of the 6 GCMs, not for each species
    ## Define function to then send to one or multiple cores
    maxent_predict_fun <- function(species) {
      
      #############################################
      save_name = gsub(' ', '_', species)
      
      occ <- readRDS(sprintf('%s%s/%s_occ.rds', maxent_path, species, save_name)) %>%
        spTransform(ALB.CONICAL)
      
      MESS_dir = sprintf('%s%s/full/%s', 
                         maxent_path, species, 'MESS_output')
      
      
      ###################################################################
      ## Read in all the objects needed to create the MESS maps for each 
      ## species.
      current = sprintf('%s/%s/full/%s_current.tif',
                        maxent_path, species, species)
      future = sprintf('%s/%s/full/%s_%s.tif', 
                       maxent_path, species, species, x)
      
      
      ## This step assumes the predictions all worked...................
      if(file.exists(current) && file.exists(future)) {
        
        pred.current = raster(sprintf('%s/%s/full/%s_current.tif',
                                      maxent_path, species, species))
        
        pred.future = raster(sprintf('%s/%s/full/%s_%s.tif', 
                                     maxent_path, species, species, x))
        
        ###################################################################
        ## Convert binary rasters of novel climate to polygons
        ## Need to save the polygons to file ::  
        ## this can fail if no ebvironments are novel, e.g the red gum.
        ## Can we add a condtiton in poluygonizer to check if the file has data?..............................
        message('Converting raster MESS maps to polygons under ', x, ' scenario for ', species) 
        novel_current_poly <- polygonizer_windows(sprintf('%s/%s%s.tif',   MESS_dir, species, "_current_novel"))
        novel_future_poly  <- polygonizer_windows(sprintf('%s/%s%s%s.tif', MESS_dir, species, "_future_novel_", x))
        
        ## Create a SpatialLines object that indicates novel areas (this will be overlaid)            
        ## Below, we create a dummy polygon as the first list element (which is the extent
        ## of the raster, expanded by 10%), to plot on panel 1). 50 = approx 50 lines across the polygon
        message('Creating polygon list under ', x, ' scenario for ', species) 
        
        #  cast the objects into the sf class so we avoid issues with wrong methods being called in hatch()
        # novel_hatch <- list(as(extent(pred.current)*1.1, 'SpatialPolygons'),
        #                     hatch(as(novel_current_poly, 'sf'), 50),
        #                     hatch(as(novel_future_poly,  'sf'), 50))
        
        novel_hatch <- list(as(extent(pred.current)*1.1, 'SpatialPolygons'),  
                            hatch(novel_current_poly, 50),
                            hatch(novel_future_poly, 50)) 
        
        
        ########################################################################################################################
        ## Now create a panel of PNG files for maxent projections and MESS maps
        ## All the projections and extents need to match
        empty_ras <- init(pred.current, function(x) NA) 
        projection(novel_current_poly);projection(occ);projection(empty_ras);projection(poly)
        projection(pred.current);projection(pred.future)
        identical(extent(pred.current), extent(pred.future))
        
        ## Assign the scenario name to use in the plot title below
        scen_name = eval(parse(text = sprintf('gcms.%s$GCM[gcms.%s$id == x]', time_slice, time_slice)))      
        
        
        ########################################################################################################################
        ## Use the 'levelplot' function to make a multipanel output: occurrence points, current raster and future raster
        # if(create_mess == "TRUE") {
        message('Create MESS panel maps for ', species, ' under ', x, ' scenario')
        
        
        ############################################################
        ## Create level plot of current conditions including MESS                        
        png(sprintf('%s/%s/full/%s_%s.png', maxent_path, species, species, "local_mess_panel"),      
            11, 4, units = 'in', res = 300)
        
        print(levelplot(stack(empty_ras,
                              pred.current, quick = TRUE), margin = FALSE,
                        
                        ## Create a colour scheme using colbrewer: 100 is to make it continuos
                        ## Also, make it a one-directional colour scheme
                        scales      = list(draw = FALSE), 
                        at = seq(0, 1, length = 100),
                        col.regions = colorRampPalette(rev(brewer.pal(11, 'Spectral'))),
                        
                        ## Give each plot a name: the third panel is the GCM
                        names.attr = c('Australian records', 'Current habitat suitability'),
                        colorkey   = list(height = 0.5, width = 3), xlab = '', ylab = '',
                        main       = list(gsub('_', ' ', species), font = 4, cex = 2)) +
                
                ## Plot the Aus shapefile with the occurrence points for reference
                ## Can the current layer be plotted on it's own?
                ## Add the novel maps as vectors.              
                latticeExtra::layer(sp.polygons(aus_poly), data = list(aus_poly = aus_poly)) +
                latticeExtra::layer(sp.points(occ, pch = 19, cex = 0.15, 
                                              col = c('red', 'transparent', 'transparent')[panel.number()]), data = list(occ = occ)) +
                latticeExtra::layer(sp.polygons(h[[panel.number()]]), data = list(h = novel_hatch)))
        
        ## Finish the device
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
        png(sprintf('%s/%s/full/%s%s_%s.png', maxent_path, species, species, "_local_mess", x),      
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
                        names.attr = c('Australian records', 'Current habitat suitability', 
                                       sprintf('%s, 20%s, RCP8.5', scen_name, time_slice)),
                        colorkey   = list(height = 0.5, width = 3), xlab = '', ylab = '',
                        main       = list(gsub('_', ' ', species), font = 4, cex = 2)) +
                
                ## Plot the Aus shapefile with the occurrence points for reference
                ## Can the current layer be plotted on it's own?
                ## Add the novel maps as vectors.              
                latticeExtra::layer(sp.polygons(aus_poly), data = list(aus_poly = aus_poly)) +
                latticeExtra::layer(sp.points(occ, pch = 19, cex = 0.15, 
                                              col = c('red', 'transparent', 'transparent')[panel.number()]), data = list(occ = occ)) +
                latticeExtra::layer(sp.lines(h[[panel.number()]]), data = list(h = novel_hatch)))
        
        ## Close the device
        dev.off()
        
      } else {
        
        message(species, ' skipped - no current or future predictions')
        
      } 
      
    }
    
    if (nclust==1) {
      
      lapply(species_list, maxent_predict_fun) 
      
    } else {
      
      cl <- makeCluster(nclust)
      clusterExport(cl, c(
        
        ## The number of objects that need exporting is full on - maybe there's a better way
        ## to run a cluster in R?
        'shp_path',    'aus_shp',      'world_shp',   'ALB.CONICAL', 'CRS.WGS.84',
        'scen_list',   'species_list', 'time_slice',  'maxent_path', 'hatch',        
        'polygonizer_windows', 'nclust'), envir = environment())
      
      ## shp_path, aus_shp, world_shp, scen_list, 
      ## species_list, maxent_path, climate_path, 
      ## grid_names, time_slice, current_grids, create_mess, nclust
      
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
      
    }
    
    
  })
  
  
}





#########################################################################################################################
## Create MESS maps for all scenarios and species. To get this done faster, try using the nclust argument
ptm <- proc.time()

tryCatch(
  create_mess_pngs(shp_path      = "./data/base/CONTEXTUAL/", ## Path for shapefile
                   aus_shp       = "aus_states.rds",          ## Shapefile, e.g. Australian states
                   world_shp     = "LAND_world.rds",          ## World shapefile          
                   scen_list     = scen_2030[1:2],            ## List of climate scenarios
                   species_list  = map_spp[15:16],            ## List of species folders with maxent models
                   maxent_path   = bs_path,
                   time_slice    = 30,
                   nclust        = 2),
  
  error = function(cond) {
    
    ## This will write the error message inside the text file, 
    ## but it won't include the species
    file.create(file.path(bs_path, "map_png_failed_2030.txt"))
    cat(cond$message, file=file.path(bs_path, "map_png_failed_2030.txt"))
    warning(cond$message)
    
  })

proc.time() - ptm



