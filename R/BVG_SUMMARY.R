#########################################################################################################################
########################### WT POINT PAIR ANALYSES : DATA READ AND CHECK  ############################################### 
#########################################################################################################################



## Read in sites and BVGs
BVG      = readOGR("./data/WT_BVG_REM_merge_trim.shp", layer = "WT_BVG_REM_merge_trim")
WT.SITES = read.csv("./data/SITES_TREES.csv", stringsAsFactors = FALSE)
BVG.DESC = read.csv("./data/BVG_SITES_LUT_TABLE.csv", stringsAsFactors = FALSE)
names(BVG)
names(WT.SITES)


## Create current Grids
## Now divide the current environmental grids by 10
world.grids.current <- stack(
  file.path('./data/base/worldclim/world/0.5/bio/current',
            sprintf('bio_%02d.tif', 1:19)))

for(i in 1:11) {
  
  ## simple loop
  message(i)
  world.grids.current[[i]] <- world.grids.current[[i]]/10
  
}


## Convert the WT sites to a shapefile
WT.SITES.SHP = SpatialPointsDataFrame(coords      = WT.SITES[c("LONGITUDE", "LATITUDE")], 
                                      data        = WT.SITES,
                                      proj4string = CRS(projection(BVG)))
class(WT.SITES.SHP)
head(WT.SITES.SHP)


## Read in GPP
GPP.AUS = raster("./data/gpp_av")
GPP.WT  = raster("./data/gpp_av_wt")
plot(GPP.RAS)
plot(GPP.WT)
projection(GPP.RAS)
projection(GPP.WT)


## Projection need to be the same for the spatial overlay 
CRS.new          <- CRS("+proj=longlat +ellps=GRS80 +no_defs")
# CRS.new          <- CRS("+proj=utm +zone=55 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")        ## EPSG:3577
# BVG              = spTransform(BVG, CRS.new)
# WT.SITES.SHP     = spTransform(WT.SITES.SHP, CRS.new)
GPP.WT            = projectRaster(GPP.WT, crs = CRS.new)

## Check they are the same
projection(BVG)
projection(WT.SITES.SHP)
projection(GPP.WT)




#########################################################################################################################
## 1). SUMMARISE GPP BY BVG 1:5M GROUP
#########################################################################################################################


## Can you aggregate the shapefile to the coarser spatial grouping?
length(unique(BVG$DBVG5M))
length(unique(BVG$DBVG1M))


## Yes, but need to check the results
DBVG5M <- unionSpatialPolygons(BVG, BVG$DBVG5M)
BVG.DF =  as.data.frame(BVG)


## Rename IDs of the dataframe to match IDs of the polygons
BVG.DF = BVG.DF[!duplicated(BVG.DF["DBVG5M"]),]
row.names(BVG.DF) = row.names(DBVG5M)
DBVG5M = SpatialPolygonsDataFrame(DBVG5M, BVG.DF)


## Don't need to summarise the rasters for this response...
## Use the zonal stats option? Can calc quantiles somehow?
# bvg.wt.mean   <- spatialEco::zonal.stats(x = DBVG5M, y = GPP.WT, stat = mean,   trace = TRUE, plot = TRUE)
# bvg.wt.max    <- spatialEco::zonal.stats(x = DBVG5M, y = GPP.WT, stat = max,    trace = TRUE, plot = TRUE)
# bvg.wt.min    <- spatialEco::zonal.stats(x = DBVG5M, y = GPP.WT, stat = min,    trace = TRUE, plot = TRUE)
# bvg.wt.median <- spatialEco::zonal.stats(x = DBVG5M, y = GPP.WT, stat = median, trace = TRUE, plot = TRUE)


## How could we summarise quantile(xx[[col]], .95)
str(bvg.wt.mean);str(bvg.wt.max)





#########################################################################################################################
## 2). BOXPLOTS FOR EACH BVG USING WT SITES
#########################################################################################################################


## Intersect BVG with the points
WT.SITES.BVG    = over(WT.SITES.SHP, BVG)
str(WT.SITES.BVG)


## Then join the koppen classification onto the data and add the site descriptions
SITES.BVG     = cbind.data.frame(WT.SITES.SHP, WT.SITES.BVG)
WT.SITES.BVG  = join(SITES.BVG, BVG.DESC)

summary(WT.SITES.BVG$DBVG5M)                                   ## Most sites are in two 1:5M BVGs
summary(WT.SITES.BVG$DBVG2M) 
summary(WT.SITES.BVG$DBVG1M) 


## Boxplots for each BVG
ggplot(data = WT.SITES.BVG, aes(x = DBVG5M, y = TOTAL_AV)) + geom_boxplot(aes(fill = DBVG5M))



########################################################################################################################
## Plot GPP vs. Rainfall for each BVG
for(bvg in bvg.list){
  
  ## Susbet the dataframe by BVG
  DAT  = subset(WT.SITES.BVG, DBVG5M == bvg)
  
  ## Plot GPP vs. RAIN for that subset
  plot(DAT[["RAIN_0112"]], DAT[["TOTAL_AV"]], 
       pch  = 19, col = "blue", #cex = 1.2,
       xlab = "Rainfall (mm)",
       ylab = expression(paste("GPP (gC m"^"-2", "month"^"-1", ")", sep = " ")),
       main = DAT[["BVG_description"]][1])
  
}


## Save as csv
BVG.LUT = SITES.BVG[!duplicated(SITES.BVG["DBVG1M"]),]
write.csv(BVG@data, "./data/BVG_FULL_TABLE.csv", row.names = FALSE)
write.csv(BVG.LUT, "./data/BVG_SITE_TABLE.csv", row.names = FALSE)




#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################
