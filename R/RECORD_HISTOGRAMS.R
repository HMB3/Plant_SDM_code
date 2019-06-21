#########################################################################################################################
############################################### SUMMARISE GBIF RECORDS ################################################## 
#########################################################################################################################


## This code creates histograms of the number of records per species across the world, and in Australia. The aim is to 
## visualise the justtification for only choosing species with enough records both globally and locally, rather than 
## just picking an arbitray number (e.g. 20). We are underestimating the global niches of the Australian species, and 
## underestimating the Australian niches of the exotics. If we focus the modelling of species in the centre of the 
## records distribution, both globally and in Australia, will the results be more robust than using an arbritary cut off?
## The hypothesis is that having too few records globally and locally leads to poor models, and bad maps.


#########################################################################################################################
## 1). CREATE COUNT OF RECORDS IN AUSTRALIA
#########################################################################################################################


## Load GBIF data
load("./data/base/HIA_LIST/COMBO/COMBO_RASTER_CONTEXT_1601_2018.RData")
load("./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT_1601_2018.RData")
load("./data/base/CONTEXTUAL/urbanareas.rda")
LAND = readRDS("F:/green_cities_sdm/data/base/CONTEXTUAL/LAND_world.rds")
aus  = readRDS("F:/green_cities_sdm/data/base/CONTEXTUAL/aus_states.rds")
names(COMBO.RASTER.CONTEXT)


## Create lists
source('./R/HIA_LIST_MATCHING.R')


## Quickly inspect the niche and raster data
dim(COMBO.NICHE.CONTEXT)
names(COMBO.NICHE.CONTEXT)
summary(COMBO.NICHE.CONTEXT$COMBO.count)
summary(COMBO.NICHE.CONTEXT$LGA.AGG)


## They need to be the same length
dim(COMBO.RASTER.CONTEXT)
length(unique(COMBO.RASTER.CONTEXT$searchTaxon))


#########################################################################################################################
## Create a sp df
COMBO.RASTER.SP   = SpatialPointsDataFrame(coords      = COMBO.RASTER.CONTEXT[c("lon", "lat")], 
                                           data        = COMBO.RASTER.CONTEXT,
                                           #proj4string = CRS('+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'),
                                           proj4string = CRS("+init=epsg:4326"))


#CRS.new  = CRS('+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')
CRS.new  = CRS("+init=epsg:4326")
AUS      = spTransform(aus, CRS.new)
projection(AUS); projection(COMBO.RASTER.SP)


## Put this code into step 4.......................................................................................................
## Create new column for if each record is in Australia
AUS.JOIN  = over(COMBO.RASTER.SP, AUS)
AUS.JOIN  = AUS.JOIN["name"]
names(AUS.JOIN) = c("STATE_NAME")
COMBO.AUS = cbind.data.frame(COMBO.RASTER.SP, AUS.JOIN)


#########################################################################################################################
## AGGREGATE THE NUMBER OF RECORDS EACH SPECIES HAS IN AUSTRALIA 
AUS.AGG   = aggregate(STATE_NAME ~ searchTaxon, data = COMBO.AUS, function(x) {sum(!is.na(x))}, na.action = NULL)
head(AUS.AGG)    ## this number matches with Arcmap
NICHE.AUS = join(COMBO.NICHE.CONTEXT, AUS.AGG)


## Check which species have records in Aus :: 14
dim(subset(NICHE.AUS, STATE_NAME >= 20))[1]
dim(subset(NICHE.AUS, STATE_NAME >= 20 & Top_200 == "TRUE"))[1]
NICHE.AUS.OK = subset(NICHE.AUS, STATE_NAME >= 20)
NICHE.AUS.POP = subset(NICHE.AUS, STATE_NAME >= 20 & Top_200 == "TRUE")


## What is the breakdown here? Need to add phylogenetic
unique(NICHE.AUS.OK$Number.of.States)
summary(NICHE.AUS.OK$STATE_NAME)
with(NICHE.AUS.OK, table(Plant.type)/sum(table(Plant.type))*100)


#########################################################################################################################
## 2). SUMMARISE GLOBAL AND AUSTRALIAN RECORD COUNTS
#########################################################################################################################


#########################################################################################################################
## How the number of records distributed across the world, and Australia
summary(NICHE.AUS$COMBO.count)
summary(NICHE.AUS$STATE_NAME)


## How many species have more global records than the first quartile, and are in the IQR
dim(subset(NICHE.AUS, COMBO.count >= 146.8))[1]
dim(subset(NICHE.AUS, COMBO.count <= 146.8))[1]
dim(subset(NICHE.AUS, COMBO.count > 146.8 | COMBO.count < 3215))[1]


## How many species have less AUS records than the first quartile
dim(subset(NICHE.AUS, STATE_NAME <= 4))[1]
dim(subset(NICHE.AUS, STATE_NAME <= 20))[1]


#########################################################################################################################
## Histograms of RAIN CV vs. TEMP CV at the 510 WT sites


## Some colour parameters
R       = rgb(0,   102,  204,  maxColorValue = 255)
T       = rgb(255, 128,    0,  maxColorValue = 255)


## All sites: rainfall is HEAPS more variable than temperature...
hist(NICHE.AUS$COMBO.count, 
     #xlim = c(0.05, 1.5), 
     breaks = 50, border = "NA", col = R, main = "",
     xlab = "Number of records",
     ylab = "Number of species")

box(lwd = 2)

hist(NICHE.AUS$STATE_NAME, 
     breaks = 50, border = "NA", col = T, add = TRUE,
     main = "",
     xlab = "Number of records",
     ylab = "Number of species")





#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################
