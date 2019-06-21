#########################################################################################################################
## Varibles to run a function over completed models
shp_path      = "./data/base/CONTEXTUAL/" ## Path for shapefile
aus_shp       = "aus_states.rds"          ## Shapefile e.g. Australian states
world_shp     = "LAND_world.rds"          ## World shapefile

x             = scen_2030[3]                 ## List of climate scenarios
species       = map_spp[15:16]                   ## List of species folders with maxent models
maxent_path   = bs_path                   ## Output folder
climate_path  = "./data/base/worldclim/aus/1km/bio" ## climate data
time_slice    = 30                        ## Time period






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
## This function turns a raster into a polygon 
## Windows version
polygonizer_windows <- function(x, outshape=NULL, pypath=NULL, readpoly=TRUE, 
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
    cmd <- Sys.which('C:/OSGeo4W64/OSGeo4W.bat')
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
## Try to run the mess maps at the same time as the map creation?
create_mess_maps = function(shp_path, aus_shp, world_shp, scen_list,
                            species_list, maxent_path, time_slice) {
  
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
      
      swd <- as.data.frame(readRDS(sprintf('%s%s/swd.rds',    maxent_path, species, species)))
      occ <- readRDS(sprintf('%s%s/%s_occ.rds', maxent_path, species, save_name)) %>%
        spTransform(ALB.CONICAL)  
      
      
      ###################################################################
      ## Read in all the objects needed to create the MESS maps for each 
      ## species.
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
      
    } 
    
  })
  
}





#########################################################################################################################
## Call the function
create_mess_maps(shp_path      = "./data/base/CONTEXTUAL/", ## Path for shapefile
                 aus_shp       = "aus_states.rds",          ## Shapefile, e.g. Australian states
                 world_shp     = "LAND_world.rds",          ## World shapefile          
                 
                 scen_list     = scen_2030[3],              ## List of climate scenarios
                 species_list  = map_spp[15:16],            ## List of species folders with maxent models
                 maxent_path   = bs_path,
                 time_slice    = 30)                        ## Output folder)


