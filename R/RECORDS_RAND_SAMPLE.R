#########################################################################################################################
################################################# THIN RECORDS ##########################################################
#########################################################################################################################


## This code spatially thins species occurrence records to help address problems associated with spatial sampling biases. 
## Ideally, thinning removes the fewest records necessary to substantially reduce the effects of sampling bias, while 
## simultaneously retaining the greatest amount of useful information. 
## See: https://cran.r-project.org/web/packages/spThin/vignettes/spThin_vignette.html


## Create lists
source('./R/HIA_LIST_MATCHING.R')


#########################################################################################################################
## 1). ADD COLUMN FOR AUS RAINFALL AND STATE CATEGORIES
#########################################################################################################################


## Load GBIF data and rain shapefile
COMBO.RASTER.CONTEXT = readRDS("./data/base/HIA_LIST/COMBO/CLEAN_ONLY_HIA_SPP.rds")
SPP.BIAS             = read.csv("./output/maxent/SPP_BOUNDARY_BIAS.csv", stringsAsFactors = FALSE)
SPP.BIAS             = subset(SPP.BIAS, AUS_BOUND_BIAS == "TRUE")$searchTaxon
View(SPP.BIAS)


#########################################################################################################################
## Intersect rainfall and temperature data
COMBO.RASTER.SP   = SpatialPointsDataFrame(coords      = COMBO.RASTER.CONTEXT[c("lon", "lat")], 
                                           data        = COMBO.RASTER.CONTEXT,
                                           proj4string = CRS.WGS.84)


## Project data
AUS.WGS    = spTransform(aus, CRS.WGS.84)
RAIN.WGS   = spTransform(AUS_RAIN, CRS.WGS.84)
AUS.STATE  = AUS.WGS[, c("name")]
AUS.RAIN   = RAIN.WGS[, c("AUS_RN_ZN")] 


## Run join between species records and Australian STATES
projection(COMBO.RASTER.SP);projection(AUS.STATE);projection(AUS.RAIN)
STATE.JOIN          = over(COMBO.RASTER.SP, AUS.STATE)                   ## =SUA.JOIN      = over(COMBO.RASTER.SP, SUA.WGS) 
COMBO.STATE         = cbind.data.frame(COMBO.RASTER.SP, STATE.JOIN)


## Recreate sp dataframe
COMBO.STATE.SP   = SpatialPointsDataFrame(coords      = COMBO.STATE[c("lon", "lat")], 
                                          data        = COMBO.STATE,
                                          proj4string = CRS.WGS.84)


## Now join RAIN categories
RAIN.JOIN           = over(COMBO.STATE.SP, AUS.RAIN)
COMBO.STATE.RAIN    = cbind.data.frame(COMBO.STATE.SP, RAIN.JOIN)
names(COMBO.STATE.RAIN)


## Reanme and remove
colnames(COMBO.STATE.RAIN)[colnames(COMBO.STATE.RAIN)=="name"] <- "AUS_STATE"
COMBO.STATE.RAIN = subset(COMBO.STATE.RAIN, select = -c(lon.1,  lat.1,  optional.1,  optional.2, lon.2, lat.2))
unique(COMBO.STATE.RAIN$AUS_STATE)
unique(COMBO.STATE.RAIN$AUS_RN_ZN)


## Save data
saveRDS(COMBO.STATE.RAIN, file = paste("./data/base/HIA_LIST/COMBO/COMBO_STATE_RAIN.rds"))
#COMBO.STATE.RAIN = readRDS("./data/base/HIA_LIST/COMBO/COMBO_STATE_RAIN.rds")





#########################################################################################################################
## 2). CREATE A RANDOM SUB-SAMPLE OF X% OF SPP RECORS FROM NSW, THEN STRATIFY BY RAINFALL AMOUNT
#########################################################################################################################


## Try for an example species ::
SPP.TEST            = subset(COMBO.STATE.RAIN, searchTaxon == "Banksia integrifolia") 
unique(SPP.TEST$searchTaxon)
dim(SPP.TEST)[1]


## Looks biased...
plot(aus)
points(SPP.TEST$lon,  SPP.TEST$lat, cex = 0.8, col = "red", pch = 19)


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
#########################################################################################################################
## Now create a random subset of points from NSW, stratified by rainfall category. The stratified function samples from 
## a data.frame in which one of the columns can be used as a "stratification" or "grouping" variable. The result is 
## a new data.frame with the specified number of samples from each group. 


## Let's take a 50% sample from all rainfall groups
unique(SPP.REC.NSW$AUS_RN_ZN)
set.seed(1)
SPP.REC.NSW = stratified(SPP.REC.NSW, "AUS_RN_ZN", .5)








#########################################################################################################################
## Save
saveRDS(TEST.GEO,                'data/base/HIA_LIST/COMBO/CLEAN_FLAGS_HIA_SPP.rds')
saveRDS(CLEAN.TRUE,              'data/base/HIA_LIST/COMBO/CLEAN_ONLY_HIA_SPP.rds')
saveRDS(CLEAN.NICHE.CONTEXT,     'data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT_APRIL_2018_COORD_CLEAN.rds')
write.csv(CLEAN.NICHE.CONTEXT,   "./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT_APRIL_2018_COORD_CLEAN.csv", row.names = FALSE)




#########################################################################################################################
#####################################################  TBC ##############################################################
#########################################################################################################################