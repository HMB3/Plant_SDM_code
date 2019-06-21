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

# devtools::source_gist('c6a1cb61b8b6616143538950e6ec34aa', filename='hatch.R')

# system.time(p1 <- rasterToPolygons(novel_current)) # too slow


polygonizer <- function(x, outshape=NULL, readpoly=TRUE, 
                        fillholes=FALSE, aggregate=FALSE, 
                        quietish=TRUE) {
  # x: an R Raster layer, or the file path to a raster file recognised by GDAL 
  # outshape: the path to the output shapefile (if NULL, a temporary file will 
  #           be created) 
  # readpoly: should the polygon shapefile be read back into R, and returned by
  #           this function? (logical) 
  # fillholes: should holes be deleted (i.e., their area added to the containing
  #            polygon)
  # aggregate: should polygons be aggregated by their associated raster value?
  # quietish: should (some) messages be suppressed? (logical)
  if (isTRUE(readpoly) || isTRUE(fillholes)) require(rgdal)
  if (isTRUE(aggregate)) require(rgeos)
  cmd <- 'c:/OSGeo4W64/OSGeo4W.bat'
  pypath <- 'gdal_polygonize'
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

# Faster to pass file paths to polygonizer
novel_current_poly <- polygonizer('output/maxent/SUA_TREES_ANALYSIS/Abelia_grandiflora/full/Abelia_grandiflora_current_novel_map.tif')
novel_future_poly <- polygonizer('output/maxent/SUA_TREES_ANALYSIS/Abelia_grandiflora/full/Abelia_grandiflora_future_novel_map.tif')
# But can also pass raster objects
novel_current_poly <- polygonizer(novel_current)
novel_future_poly <- polygonizer(novel_future)


# below we create a dummy polygon as the first list element (which is the extent
# of the raster, expanded by 10%), to plot on panel 1
h <- list(as(extent(pred.current)*1.1, 'SpatialPolygons'),  
          hatch(novel_current_poly, 50), hatch(novel_future_poly, 50)) 


(levelplot(stack(empty_ras,
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
    ## Try adding novel maps as vectors.....................................
    latticeExtra::layer(sp.polygons(poly), data = list(poly = poly)) +
    latticeExtra::layer(sp.points(occ, pch = 19, cex = 0.15,
                                  col = c('red', 'transparent', 'transparent')[panel.number()]), data = list(occ = occ)) +
    latticeExtra::layer(sp.polygons(hatch[[panel.number()]], col='#000000dd'), data = list(hatch = h)))

