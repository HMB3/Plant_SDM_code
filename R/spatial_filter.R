## Pascal @ SO - from http://stackoverflow.com/a/22507264
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

# With 10 km spatial filtering
occ_albers_filter10 <- filterByProximity(occ_albers, 10000, mapUnits=TRUE)


d <- read.csv('~/../Downloads/maxent_occurrences_qfly.csv') # ALA and trap, from Linda
r <- raster('~/../Downloads/mask (1)/mask.bil')
b <- buffer(d, width=200000, dissolve=FALSE)
writeOGR(b, tempdir(), 'buffer', 'ESRI Shapefile')

b_rast <- 
  gdal_rasterize(file.path(tempdir(), 'buffer.shp'), 
                 file.path(tempdir(), 'mask.bil'), 
                 of='EHdr',  a_srs='EPSG:3577',
                 burn=1, init=-9999, 
                 a_nodata=-9999, ot='Int16', tr=c(1000, 1000), 
                 te=c(bbox(r)), output_Raster=TRUE)

# 
d <- read.csv('data/occurrence/trap_data_nsw_vic_sa_coords.csv')
coordinates(d) <- ~lon+lat
proj4string(d) <- '+init=epsg:4283'
d <- spTransform(d, CRS('+init=epsg:3577'))
d10000 <- filterByProximity(d, 10000, T)
d25000 <- filterByProximity(d, 25000, T)
d50000 <- filterByProximity(d, 50000, T)
d100000 <- filterByProximity(d, 100000, T)


par(mfrow=c(3, 2), mar=c(0, 0, 3, 0))
plot(aus_albers, main=sprintf('unfiltered (n=%s)', length(d)))
points(d, pch=20, col=2, cex=1)
plot(aus_albers, main=sprintf('10 km spatial filter (n=%s)', length(d10000)))
points(d10000, pch=20, col=2, cex=1)
plot(aus_albers, main=sprintf('25 km spatial filter (n=%s)', length(d25000)))
points(d25000, pch=20, col=2, cex=1)
plot(aus_albers, main=sprintf('50 km spatial filter (n=%s)', length(d50000)))
points(d50000, pch=20, col=2, cex=1)
plot(aus_albers, main=sprintf('100 km spatial filter (n=%s)', length(d100000)))
points(d100000, pch=20, col=2, cex=1)