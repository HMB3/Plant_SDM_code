#########################################################################################################################
################################################# THIN RECORDS ##########################################################
#########################################################################################################################


## This code spatially thins species occurrence records to help address problems associated with spatial sampling biases. 
## Ideally, thinning removes the fewest records necessary to substantially reduce the effects of sampling bias, while 
## simultaneously retaining the greatest amount of useful information. 


## Using the SDM dataset which has already been thinned to 1km, rather than the original dataset

## Use the stratified function: https://www.rdocumentation.org/packages/fifer/versions/1.0/topics/stratified
## spThin didn't really work https://cran.r-project.org/web/packages/spThin/vignettes/spThin_vignette.html


## Create lists
source('./R/HIA_LIST_MATCHING.R')



#########################################################################################################################
## 1). ADD COLUMNS FOR AUS RAINFALL AND STATE CATEGORIES
#########################################################################################################################


#########################################################################################################################
## Load GBIF data and rain shapefile
BIAS.DATA.ALL        = readRDS("./data/base/HIA_LIST/COMBO/SDM_DATA_CLEAN_052018.rds")        
SPP.BIAS             = intersect(SPP.BIAS, SUA.spp)    ## just re-run the models for species on the list


## Project the SDM data into WGS
BIAS.DATA.ALL <- spTransform(BIAS.DATA.ALL, CRS("+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
projection(BIAS.DATA.ALL)


## Get the coordinates
BIAS.COORDS = coordinates(BIAS.DATA.ALL)
class(BIAS.COORDS)
tail(BIAS.COORDS)
tail(BIAS.DATA.ALL)[, c(1:3)]   ## indices match


## Bind the coordinates to the SDM table: can use spCbind for sp data frames
BIAS.COORDS   = as.data.frame(BIAS.COORDS)
BIAS.DATA.DF = cbind(as.data.frame(BIAS.DATA.ALL), BIAS.COORDS)
BIAS.DATA.DF = BIAS.DATA.DF[, c(1:22)]
names(BIAS.DATA.DF)
dim(BIAS.DATA.DF)


#########################################################################################################################
## Intersect BOM with rainfall data
COMBO.RASTER.SP   = SpatialPointsDataFrame(coords      = BIAS.DATA.DF[c("lon", "lat")], 
                                           data        = BIAS.DATA.DF,
                                           proj4string = CRS.WGS.84)


## Project data
AUS.WGS    = spTransform(aus, CRS.WGS.84)
RAIN.WGS   = spTransform(AUS_RAIN, CRS.WGS.84)
AUS.STATE  = AUS.WGS[, c("name")]
names(AUS.STATE)[names(AUS.STATE) == 'name'] <- 'AUS_STATE'
AUS.RAIN   = RAIN.WGS[, c("AUS_RN_ZN")] 


## Save these files out to test in ArcMap
## dim(COMBO.RASTER.SP);names(COMBO.RASTER.SP)
# writeOGR(obj = SDM.DATA.ALL, dsn = "./data/base/HIA_LIST/COMBO", layer = "SDM_DATA_TEST", driver = "ESRI Shapefile")
# writeOGR(obj = RAIN.WGS, dsn = "./data/base/HIA_LIST/COMBO", layer = "RAIN.WGS", driver = "ESRI Shapefile")



## Run join between species records and Australian STATES
## This join is not working...............................................................................................
projection(COMBO.RASTER.SP);projection(AUS.STATE);projection(AUS.RAIN)
plot(COMBO.RASTER.SP)


STATE.JOIN          = over(COMBO.RASTER.SP, AUS.STATE)                   ## =SUA.JOIN      = over(COMBO.RASTER.SP, SUA.WGS) 
COMBO.STATE         = cbind.data.frame(COMBO.RASTER.SP, STATE.JOIN)
names(COMBO.STATE)
unique(COMBO.STATE$AUS_STATE)
saveRDS(COMBO.STATE, file = paste("./data/base/HIA_LIST/COMBO/COMBO_STATE.rds"))

## Recreate sp dataframe
COMBO.STATE.SP   = SpatialPointsDataFrame(coords      = COMBO.STATE[c("lon", "lat")], 
                                          data        = COMBO.STATE,
                                          proj4string = CRS.WGS.84)


## Now join RAIN categories
RAIN.JOIN           = over(COMBO.STATE.SP, AUS.RAIN)
COMBO.STATE.RAIN    = cbind.data.frame(COMBO.STATE.SP, RAIN.JOIN)
names(COMBO.STATE.RAIN)
unique(COMBO.STATE$AUS_STATE)

## Rename and remove
COMBO.STATE.RAIN = subset(COMBO.STATE.RAIN, select = -c(lon.1,  lat.1, optional, optional.1, lon.2, lat.2))
names(COMBO.STATE.RAIN)
unique(COMBO.STATE.RAIN$AUS_STATE)
unique(COMBO.STATE.RAIN$AUS_RN_ZN)


## Save data
saveRDS(COMBO.STATE.RAIN, file = paste("./data/base/HIA_LIST/COMBO/COMBO_STATE_RAIN.rds"))
#COMBO.STATE.RAIN = readRDS("./data/base/HIA_LIST/COMBO/COMBO_STATE_RAIN.rds")


#########################################################################################################################
## Now create the variables needed to access current environmental conditions + their names in the functions
sdm.predictors <- c("Annual_mean_temp",    "Mean_diurnal_range",  "Isothermality",      "Temp_seasonality",  
                    "Max_temp_warm_month", "Min_temp_cold_month", "Temp_annual_range",  
                    "Mean_temp_wet_qu",    "Mean_temp_dry_qu",    
                    "Mean_temp_warm_qu",   "Mean_temp_cold_qu",   
                    
                    "Annual_precip",       "Precip_wet_month",   "Precip_dry_month",    "Precip_seasonality",  
                    "Precip_wet_qu",       "Precip_dry_qu",      
                    "Precip_warm_qu",      "Precip_col_qu")


## A-priori worldclim predictors
sdm.select     <- c("Annual_mean_temp", "Temp_seasonality",    "Max_temp_warm_month", "Min_temp_cold_month",
                    "Annual_precip",    "Precip_seasonality",  
                    "Precip_wet_month", "Precip_dry_month")      


## Create a raster stack of current environmental conditions if needed
i  <- match(sdm.predictors, sdm.predictors)
ff <- file.path('./data/base/worldclim/world/0.5/bio/current',
                sprintf('bio_%02d.tif', i))


## Name the grids :: these should be indentical
env.grids.current = stack(sub('0.5', '1km', ff))
names(env.grids.current) <- sdm.predictors[i]
identical(names(env.grids.current),sdm.predictors)





#########################################################################################################################
## 2). CREATE A RANDOM SUB-SAMPLE OF X% OF SPP RECORDS
#########################################################################################################################


# Shawn :: The Biodiverse approach would snap the coordinates to the centre of each cell.  It would be better to select 
# one point at random within each cell, which should be pretty simple in R.  

# The random sample for a large starting number would probably be reasonable. With either of the above approaches, you also 
# have the problem that you might lose some of your environmental variability, in which case it might be better to thin 
# in environmental space as well as geospatial.  Maybe a stratified thinning would be sensible?


## Use the stratified function: https://www.rdocumentation.org/packages/fifer/versions/1.0/topics/stratified

## For each species in the list :: 


## subset the data to just those species records

## subset the data to just those species records in NSW, and outside NSW

## Subsample within all rainfall bands - take 50% of total records

## Then combine the stratified sample with those records outside NSW




#########################################################################################################################
## Can split sampling by state, before analysis
## First, subset the species data to 
unique(COMBO.STATE.RAIN$AUS_RN_ZN)
unique(COMBO.STATE.RAIN$AUS_STATE)



#########################################################################################################################
## Loop over all the species with bias: spp = SPP.BIAS[1]
lapply(SPP.BIAS, function(spp){ 
  
  ## Skip the species if the directory already exists, before the loop
  ## Create a separate directory for these subsampled species. Then we can compare the folder output with the originals
  outdir <- 'output/maxent/SPP_BIAS'
  if(dir.exists(file.path(outdir, gsub(' ', '_', spp)))) {
    message('Skipping ', spp, ' - already run.')
    invisible(return(NULL))
    
  }
  
  ## Print the taxa being processed to screen
  if(spp %in% COMBO.STATE.RAIN$searchTaxon) {
    message('Doing ', spp) 
    
    ## Do the stratification by rainfall, etc, in here...................................................................
    
    ## First, subset the records to only the taxa being processed
    occurrence <- subset(COMBO.STATE.RAIN, searchTaxon == spp)
    
    ## Now subset to the species data to those records inside and outside NSW
    ## Could change this to a generic state argument
    occurrence.nsw <- subset(occurrence, AUS_STATE == "New South Wales")
    occurrence.out <- subset(occurrence, AUS_STATE != "New South Wales")
    
    ## Now take a 50% random sample from all rainfall groups, just for the points within NSW
    ## Set a large starting number R? How large is large?
    set.seed(123458)
    strat.nsw = stratified(occurrence.nsw, "AUS_RN_ZN", .5)
    
    ## Now combine the two data frames
    occ.strat = rbind(strat.nsw, occurrence.out)
    names(occ.strat)
    dim(occ.strat);dim(occurrence)
    
    ## Now get the background points. These can come from any spp, other than the modelled species.
    background <- subset(COMBO.STATE.RAIN, searchTaxon != spp)
    
    ## Finally fit the models using FIT_MAXENT_TARG_BG. Also use tryCatch to skip any exceptions
    tryCatch(
      FIT_MAXENT_TARG_BG(occ                     = occ.strat, 
                         bg                      = background, 
                         sdm.predictors          = sdm.select, 
                         name                    = spp, 
                         outdir                  = outdir, 
                         template.raster,
                         min_n                   = 20,     ## This should be higher...
                         max_bg_size             = 100000, 
                         Koppen                  = Koppen_1975,
                         background_buffer_width = 200000,
                         shapefiles              = TRUE,
                         features                = 'lpq',
                         replicates              = 5,
                         responsecurves          = TRUE),
      
      ## https://stackoverflow.com/questions/19394886/trycatch-in-r-not-working-properly
      #function(e) message('Species skipped ', spp)) ## skip any species for which the function fails
      error = function(cond) {
        
        message(paste('Species skipped ', spp))
        
      })
    
  } else {
    
    message(spp, ' skipped - no data.')         ## This condition ignores species which have no data...
    
  }  
  
})





#########################################################################################################################
## 3). CREATE A RANDOM SUB-SAMPLE OF X% OF SPP RECORS
#########################################################################################################################


#########################################################################################################################
## OR, split sampling by state ::
SPP.REC.NSW = subset(COMBO.STATE.RAIN, AUS_STATE == "New South Wales") 
SPP.REC.QLD = subset(COMBO.STATE.RAIN, AUS_STATE == "Queensland") 
SPP.REC.VIC = subset(COMBO.STATE.RAIN, AUS_STATE == "Victoria")
SPP.REC.SA  = subset(COMBO.STATE.RAIN, AUS_STATE == "South Australia") 
unique(SPP.REC.NSW$AUS_STATE)
dim(SPP.REC.NSW)


## So 27% of all the clean data from the whole world is in NSW...
dim(SPP.REC.NSW)[1]/dim(COMBO.STATE.RAIN)[1]*100
dim(SPP.REC.QLD)[1]/dim(COMBO.STATE.RAIN)[1]*100
dim(SPP.REC.VIC)[1]/dim(COMBO.STATE.RAIN)[1]*100
dim(SPP.REC.SA)[1]/dim(COMBO.STATE.RAIN)[1]*100


## Of the Austraian data, most is in NSW
unique(COMBO.STATE.RAIN$AUS_STATE)
state.counts <- table(COMBO.STATE.RAIN$AUS_STATE)
barplot(state.counts, main = "All records per state",
        xlab = "State", 
        ylab = "No. records", 
        col = c("lightblue", "pink", "orange", "green", "grey", "beige", "coral", "chocolate4", "aquamarine", "blueviolet"))


## Pie chart
pie(state.counts, main = "AUS records per state")


#########################################################################################################################    
## Now create a random subset of points from NSW, stratified by rainfall category. The stratified function samples from 
## a data.frame in which one of the columns can be used as a "stratification" or "grouping" variable. The result is 
## a new data.frame with the specified number of samples from each group. 


## Using only points in NSW, let's take a 50% sample from all rainfall groups
unique(SPP.REC.NSW$AUS_RN_ZN)
set.seed(1)
RAND.REC.NSW = stratified(SPP.REC.NSW, "AUS_RN_ZN", .5)
dim(RAND.REC.NSW)[1]/dim(SPP.REC.NSW)[1]


## The points are unevenly distributed across the rainfall zones, so we can't sample randomly from these categories 
round(with(SPP.REC.NSW, table(AUS_RN_ZN)/sum(table(AUS_RN_ZN))*100), 1)
round(with(RAND.REC.NSW, table(AUS_RN_ZN)/sum(table(AUS_RN_ZN))*100), 1)

table(SPP.REC.NSW$AUS_RN_ZN)
table(RAND.REC.NSW$AUS_RN_ZN)


## Plot an example species ::
SPP.ORIG    = subset(COMBO.STATE.RAIN, searchTaxon == "Banksia integrifolia") 
unique(SPP.ORIG$searchTaxon)
dim(SPP.TEST)[1]


## Original records are biased to NSW
plot(aus)
points(SPP.ORIG$lon,  SPP.ORIG$lat, 
       cex = 0.8, col = "red", pch = 19, main = "Original data")


## Randomly sampled data
SPP.RAND    = subset(RAND.REC.NSW, searchTaxon == "Banksia integrifolia")
unique(SPP.RAND$searchTaxon)
dim(SPP.RAND)[1]


##
plot(aus)
points(SPP.RAND$lon,  SPP.RAND$lat, 
       cex = 0.8, col = "red", pch = 19, main = "Subsampled data")


## Check these points in ArcMap
SPP.ORIG$SAMPLE = "Original";SPP.RAND$SAMPLE = "Sampled"
SPP.SAMPLE = as.data.frame(rbind(SPP.ORIG, SPP.RAND))
unique(SPP.SAMPLE$SAMPLE)


## They need this data 
SPP.SAMPLE   = SpatialPointsDataFrame(coords      = SPP.SAMPLE[c("lon", "lat")], 
                                      data        = SPP.SAMPLE,
                                      proj4string = CRS.WGS.84)
writeOGR(obj = SPP.SAMPLE, dsn = "./data/base/HIA_LIST/COMBO", layer = "SPP_RAND_SAMPLE", driver = "ESRI Shapefile")


#########################################################################################################################
## Write the points out to check in arc





#########################################################################################################################
#####################################################  TBC ##############################################################
#########################################################################################################################