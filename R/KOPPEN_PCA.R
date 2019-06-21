#########################################################################################################################
########################### CREATE A MATRIX OD KOPPEN BY ENVIRONEMNT  ################################################### 
#########################################################################################################################


## Read in sites and BVGs
source('./R/HIA_LIST_MATCHING.R')
KOPPEN   = readOGR("F:/green_cities_sdm/data/base/CONTEXTUAL/WC05_1975H_Koppen.shp", layer = "WC05_1975H_Koppen")
KOP.BIO  = read.csv("./data/base/HIA_LIST/COMBO/KOPPEN_BIO_MATRIX.csv", stringsAsFactors = FALSE)


## Check dimensions
str(unique(KOPPEN$Koppen))
names(KOP.BIO);dim(KOP.BIO)


## Create current Grids
## Now divide the current environmental grids by 10
world.grids.current <- stack(
  file.path('./data/base/worldclim/world/0.5/bio/current',
            sprintf('bio_%02d', 1:19)))

for(i in 1:11) {
  
  ## simple loop
  message(i)
  world.grids.current[[i]] <- world.grids.current[[i]]/10
  
}

 
# ## Convert the WT sites to a shapefile
# WT.SITES.SHP = SpatialPointsDataFrame(coords      = WT.SITES[c("LONGITUDE", "LATITUDE")], 
#                                       data        = WT.SITES,
#                                       proj4string = CRS(projection(BVG)))
# class(WT.SITES.SHP)
# head(WT.SITES.SHP)


## Plot the variables
class(env.grids.current)
names(env.grids.current)
plot(env.grids.current[[8]]);plot(env.grids.current[[9]])


## Projection need to be the same for the spatial overlay 
# CRS.new          <- CRS("+proj=longlat +ellps=GRS80 +no_defs")
# # CRS.new          <- CRS("+proj=utm +zone=55 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")        ## EPSG:3577
# # BVG              = spTransform(BVG, CRS.new)
# # WT.SITES.SHP     = spTransform(WT.SITES.SHP, CRS.new)
# GPP.WT            = projectRaster(GPP.WT, crs = CRS.new)


## Check they are the same
projection(KOPPEN)
projection(env.grids.current)





#########################################################################################################################
## 1). SUMMARISE GPP BY BIOCLIM VARUABLES BY KOPPEN ZONE
#########################################################################################################################


## Don't need to summarise the rasters for this response...
## Use the zonal stats option? Can calc quantiles somehow?


kop.BIO1.mean   <- spatialEco::zonal.stats(x = KOPPEN, y = env.grids.current[[1]], stat = mean)
kop.BIO1.max    <- spatialEco::zonal.stats(x = DBVG5M, y = GPP.WT, stat = max)
kop.BIO1.min    <- spatialEco::zonal.stats(x = DBVG5M, y = GPP.WT, stat = min)
kop.BIO1.median <- spatialEco::zonal.stats(x = DBVG5M, y = GPP.WT, stat = median)


## How could we summarise quantile(xx[[col]], .95)
str(bvg.wt.mean);str(bvg.wt.max)


#########################################################################################################################
## 2). READ IN ZONAL STATS
#########################################################################################################################


##
KOP.tables = list.files("./data/base/CONTEXTUAL/KOPPEN_MAPS/", pattern = '.csv', full.names = TRUE, recursive = TRUE) 


## Now try combing the tables
KOPPEN.PCA.TABLE <- KOP.tables[c(1:length(KOP.tables))] %>%
  
  ## pipe the list into lapply
  lapply(function(x) {
    
    ## create the character string
    f <- paste0(x)
    
    ## Get name (ugly!!)
    BIO = basename(f)
    BIO =  gsub("KOP_", "", BIO)
    BIO =  gsub(".csv", "", BIO)
    BIO
    
    ## load each .RData file
    d <- read.csv(f)
    
    ## Drop columns
    keeps <- c("KOPPEN", "MIN", "MAX", "RANGE", "MEAN", "STD")
    d = d[ , keeps, drop = FALSE]
    
    ## Add BIO prefix to each colum
    colnames(d) <- paste(BIO, colnames(d), sep = "_")
    names(d)[1] <- "KOPPEN"
    
    ## divide each column by
    # for(i in 2:6) {
    #   
    #   ## Simple loop
    #   message(i)
    #   d[i] <- d[i]/10
    #   
    # }
    d
    
  }) %>%
  
  ## finally, bind columns
  bind_cols


##
dim(KOPPEN.PCA.TABLE)
View(KOPPEN.PCA.TABLE)


##
save(KOPPEN.PCA.TABLE, file = paste("./data/base/CONTEXTUAL/KOPPEN_MAPS/KOPPEN_PCA_TABLE.RData"))
  

  
  
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






#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################
