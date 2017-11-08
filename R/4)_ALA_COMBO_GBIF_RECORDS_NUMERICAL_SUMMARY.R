#########################################################################################################################
############################################# ESTIMATE SPECIES NICHES ################################################### 
#########################################################################################################################


#########################################################################################################################
## First, consider the HIA brief again:

# The first module will focus on fifty plant species identified in the project’s Target Species List, and will develop maps 
# that demonstrate each species’ suitability to both current and future climates across Australia.
# 
# These maps will be used to demonstrate how well or poorly a particular species will be able to tolerate future conditions 
# in urban centres across Australia as the climate changes, based on our current understanding of species’ climatic 
# requirements.

# Our research might demonstrate, for example, that a particular species of tree is already at the very limit of its 
# ability to cope with heat, and that the only suitable place to plant this species in the future will be in cool-climate 
# or more temperate locations. This kind of information would be very useful to a council seeking to avoid investing in 
# tree species for street planting that are unlikely to cope with higher temperatures.

# We will also use information from national herbaria and other sources to quantify each species’ climatic limits - the 
# warmest, coldest, driest or wettest conditions they can cope with. This information will then be tested through the 
# Planting Successes and Failures module of the research programme to ensure that the Interactive Plant Features Tool 
# matches the right plant in the right region with an eye on the future.

# We will also be working with growers, nurseries, landscape architects and many others to capture their recordings of 
# major plant traits including:


# Growth rate and form
# height
# canopy density
# Ground cover
# Longevity
# Seasonality
# Water quality
# Allergenicity
# Air and water quality influences and urban temperatures
# Insect resistance
# Ornamental and amenity features
# and biodiversity impacts.


## First step is to estimate the current global realised niche, using the best available data. This will give us 
## a broad indication of the currernt climatic toleance of each species.

## There are two streams here: 

## The broad niches for each species (ie. macroscale env and bio data)
## The bespoke microscale approach (e.g. ecological framework, solar surfaces, etc.)

## The broad niche stream informs the micro stream, based on two tables:

## 1). A table with one row for each species record
## 2). A table with One row for each species, contextual data and species attributes (niches, traits, etc.)





#########################################################################################################################
## WORLDCLIM VARIABLES 
#########################################################################################################################


#########################################################################################################################
## WHICH WORLDCLIM VARIABLES TO ESTIMATE?


## Copy the ones which Rach and Stu have used for the niche finder website for now, ignore edaphic variables
## Use Threshold based traits, use Dave Kendall's approach.

# BIO1  = Annual Mean Temperature                                     ## 
# BIO2  = Mean Diurnal Range (Mean of monthly (max temp - min temp))  ##
# BIO3  = Isothermality (BIO2/BIO7) (* 100)
# BIO4  = Temperature Seasonality (standard deviation *100)           ##
# BIO5  = Max Temperature of Warmest Month                            ## Take out max of max
# BIO6  = Min Temperature of Coldest Month                            ## Take out min of min, use threshold variables 
# BIO7  = Temperature Annual Range (BIO5-BIO6)
# BIO8  = Mean Temperature of Wettest Quarter
# BIO9  = Mean Temperature of Driest Quarter
# BIO10 = Mean Temperature of Warmest Quarter
# BIO11 = Mean Temperature of Coldest Quarter
# BIO12 = Annual Precipitation                                        ##
# BIO13 = Precipitation of Wettest Month                              ##
# BIO14 = Precipitation of Driest Month                               ##
# BIO15 = Precipitation Seasonality (Coefficient of Variation)        ##
# BIO16 = Precipitation of Wettest Quarter
# BIO17 = Precipitation of Driest Quarter
# BIO18 = Precipitation of Warmest Quarter
# BIO19 = Precipitation of Coldest Quarter


## To save time, load in previous data
load("./data/base/HIA_LIST/GBIF/GBIF_LAND_POINTS.RData")
load("./data/base/HIA_LIST/ALA/ALA_LAND_POINTS.RData")


## New names are in the database...
names(GBIF.LAND)
str(GBIF.LAND)
str(ALA.LAND)





#########################################################################################################################
## 1). MERGE GBIF DATA WITH CLEANED AVH/ALA RECORDS
#########################################################################################################################


#########################################################################################################################
## Merge on the ALA data. Consider that GBIF has data for both sources. We are topping up the native ranges with the AVH, 
## so it will be important to get rid of the duplicates. 
names(GBIF.LAND)
names(ALA.LAND)
setdiff(names(GBIF.LAND), names(ALA.LAND))


#########################################################################################################################
## Rename a few fields: Add rename of cultivated field, and make the values the same as mine "CULT" or "UNKNOWN"
## X = CULTIVATED
ALA.LAND     = dplyr::rename(ALA.LAND, 
                             searchTaxon                   = scientificname,
                             coordinateUncertaintyInMeters = uncertainty_m)


## Rename the field values: cultivated and uncultivated?
#ALA.LAND$CULT  = gsub("x", "CULTIVATED",   ALA.LAND$CULT)
#ALA.LAND$CULT  = gsub("x", "UNCULTIVATED", ALA.LAND$CULT)


## Restrict ALA data to just those species on the big list
library(plyr)
HIA.SPP.JOIN     = HIA.SPP
HIA.SPP.JOIN     = dplyr::rename(HIA.SPP.JOIN, searchTaxon = Binomial)


## Set NA to blank, then sort by no. of growers to get them to the top
HIA.SPP.JOIN[is.na(HIA.SPP.JOIN)] <- 0
HIA.SPP.JOIN = HIA.SPP.JOIN[with(HIA.SPP.JOIN, rev(order(Number.of.growers))), ]
head(HIA.SPP.JOIN[, c("searchTaxon", "Number.of.growers")])
View(HIA.SPP.JOIN)


## Get just those ALA species which are on the combined list of HIA and planted/growing
## ALA.LAND.HIA  = ALA.LAND[ALA.LAND$searchTaxon %in% HIA.SPP.JOIN$searchTaxon, ] 
ALA.LAND.HIA  = ALA.LAND[ALA.LAND$searchTaxon %in% all.taxa, ] ## Including Paul's extra species
str(unique(ALA.LAND$searchTaxon))       ##
str(unique(ALA.LAND.HIA$searchTaxon))   ## Reduced from 30k to 4K


## Bind the rows together?
GBIF.ALA.COMBO.LAND = bind_rows(GBIF.LAND, ALA.LAND.HIA)
names(GBIF.ALA.COMBO.LAND)
identical((dim(GBIF.LAND)[1]+dim(ALA.LAND)[1]),dim(GBIF.ALA.COMBO.LAND)[1])   ## only adding the ovelap...
head(GBIF.ALA.COMBO.LAND)


## Check the new data frame has just the species on the HIA list
str(unique(GBIF.ALA.COMBO.LAND$searchTaxon))   ## ok


## What species are unique to each dataset?
length(setdiff(unique(GBIF.LAND$searchTaxon),  unique(ALA.LAND$searchTaxon)))
length(setdiff(unique(ALA.LAND$searchTaxon),   unique(GBIF.LAND$searchTaxon)))
length(intersect(unique(ALA.LAND$searchTaxon), unique(GBIF.LAND$searchTaxon)))


## Create points: consider changing the coordinate system here to a global projected system?
COMBO.POINTS   = SpatialPointsDataFrame(coords = GBIF.ALA.COMBO.LAND[c("lon", "lat")], 
                                        data    = GBIF.ALA.COMBO.LAND[c("lon", "lat")],
                                        proj4string = CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))


## Check
summary(COMBO.POINTS)





#########################################################################################################################
## CREATE NEW COLUMNS FOR ALA CULTIVATED/NOT 
#########################################################################################################################


## Need the ALA data to have the columns for cultivated/not: Need the original data frame. Also we don't need to worry
## that much about the duplicate recrods between GBIF and ALA, given we will just take one records per grid cell. Having
## duplicate records when estimating the niche won't make much difference either, assuming the are close together in space.

## Multiple records of the same specimen from different herbaria could be a problem though. Rachel's criteria of same month,
## year, lat/long could help here though. Ask Stu for this code.



#########################################################################################################################
## 2). EXTRACT RASTER DATA FOR SPECIES RECORDS
#########################################################################################################################


#########################################################################################################################
## Create a stack of rasters to sample: get all the World clim variables just for good measure
env.grids = c("//sci-7910/F/data/worldclim/world/0.5/bio/current/bio_01",
              "//sci-7910/F/data/worldclim/world/0.5/bio/current/bio_02", 
              "//sci-7910/F/data/worldclim/world/0.5/bio/current/bio_03",
              "//sci-7910/F/data/worldclim/world/0.5/bio/current/bio_04",        
              "//sci-7910/F/data/worldclim/world/0.5/bio/current/bio_05",
              "//sci-7910/F/data/worldclim/world/0.5/bio/current/bio_06",
              "//sci-7910/F/data/worldclim/world/0.5/bio/current/bio_07",
              "//sci-7910/F/data/worldclim/world/0.5/bio/current/bio_08",
              "//sci-7910/F/data/worldclim/world/0.5/bio/current/bio_09",
              "//sci-7910/F/data/worldclim/world/0.5/bio/current/bio_10", 
              "//sci-7910/F/data/worldclim/world/0.5/bio/current/bio_11",
              "//sci-7910/F/data/worldclim/world/0.5/bio/current/bio_12",        
              "//sci-7910/F/data/worldclim/world/0.5/bio/current/bio_13",
              "//sci-7910/F/data/worldclim/world/0.5/bio/current/bio_14",
              "//sci-7910/F/data/worldclim/world/0.5/bio/current/bio_15",
              "//sci-7910/F/data/worldclim/world/0.5/bio/current/bio_16",
              "//sci-7910/F/data/worldclim/world/0.5/bio/current/bio_17",
              "//sci-7910/F/data/worldclim/world/0.5/bio/current/bio_18",
              "//sci-7910/F/data/worldclim/world/0.5/bio/current/bio_19")


## Convert all the rasters to a stack
s <- stack(env.grids)


#########################################################################################################################
## Then use the extract function for all the rasters, and finaly bind on the COMBO data to the left of the raster values
## Can we use a cluster to speed this up?

## Best option to speed this up is to use only the unique cells
# COMBO.XY <- cellFromXY(world.temp, COMBO.POINTS) %>% 
#   
#   ## get the unique raster cells
#   unique %>% 
#   
#   ## Get coordinates of the center of raster cells for a row, column, or cell number of WORLDCLIM raster
#   xyFromCell(world.temp, .)


#########################################################################################################################
## Probably best not to use a cluster like this
#beginCluster(n = 8)
COMBO.RASTER <- extract(s, COMBO.POINTS) %>% 
  cbind(GBIF.ALA.COMBO.LAND, .)
#endCluster()

## Multiple rename using dplyr
COMBO.RASTER = dplyr::rename(COMBO.RASTER,
                             Annual_mean_temp     = bio_01, ##
                             Mean_diurnal_range   = bio_02,
                             Isothermality        = bio_03,
                             Temp_seasonality     = bio_04, ##
                             Max_temp_warm_month  = bio_05, ##
                             Min_temp_cold_month  = bio_06, ##
                             Temp_annual_range    = bio_07,
                             Mean_temp_wet_qu     = bio_08,
                             Mean_temp_dry_qu     = bio_09,
                             Mean_temp_warm_qu    = bio_10,
                             Mean_temp_cold_qu    = bio_11,
                             
                             Annual_precip        = bio_12, ##
                             Precip_wet_month     = bio_13, ##
                             Precip_dry_month     = bio_14, ##
                             Precip_seasonality   = bio_15, ##
                             Precip_wet_qu        = bio_16, ##
                             Precip_dry_qu        = bio_17, ##
                             Precip_warm_qu       = bio_18,
                             Precip_col_qu        = bio_19)


## Save/load
save(COMBO.RASTER, file = paste("./data/base/HIA_LIST/GBIF/COMBO_GBIF_ALA_RASTER.RData"))
#load("./data/base/HIA_LIST/GBIF/COMBO_GBIF_ALA_RASTER.RData")


## check
dim(COMBO.RASTER)
names(COMBO.RASTER)




#########################################################################################################################
## 3). INTERSECT SPECIES RECORDS WITH LGA/SUA
#########################################################################################################################


## We want to know the count of species that occur in n LGAs, across a range of climates. Read in LGA and SUA
COMBO.RASTER.SP   = SpatialPointsDataFrame(coords = COMBO.RASTER[c("lon", "lat")], 
                                           data   = COMBO.RASTER,
                                           proj4string = CRS("+init=epsg:4326"))

SUA      = readOGR("./data/base/CONTEXTUAL/SUA_2011_AUST.shp", layer = "SUA_2011_AUST")
LGA      = readOGR("./data/base/CONTEXTUAL/LGA_2016_AUST.shp", layer = "LGA_2016_AUST")

names(SUA)
names(LGA)


## Project
CRS.new  <- CRS("+init=epsg:4326") # EPSG:3577
LGA.WGS  = spTransform(LGA, CRS.new)
SUA.WGS  = spTransform(SUA, CRS.new)

projection(COMBO.RASTER.SP)
projection(LGA.WGS)
projection(SUA.WGS)


#########################################################################################################################
## Run join
LGA.JOIN   = over(COMBO.RASTER.SP, LGA.WGS)              ## [1:300,]
COMBO.LGA  = cbind.data.frame(COMBO.RASTER.SP, LGA.JOIN) ## [1:300,]


#########################################################################################################################
## AGGREGATE THE NUMBER OF LGAs EACH SPECIES IS FOUND IN 
LGA.AGG   = tapply(COMBO.LGA$LGA_NAME16, COMBO.LGA$searchTaxon, function(x) length(unique(x))) ## group LGA by species name
LGA.AGG   = as.data.frame(LGA.AGG)
head(LGA.AGG)


## Save
save(COMBO.LGA, file = paste("./data/base/HIA_LIST/GBIF/COMBO_LGA.RData"))
save(LGA.AGG,   file = paste("./data/base/HIA_LIST/GBIF/LGA_AGG.RData"))





#########################################################################################################################
## 4). CREATE NICHES FOR SELECTED TAXA
#########################################################################################################################


#########################################################################################################################
## Now summarise the niches. But figure out a cleaner way of doing this
env.variables = c("Annual_mean_temp",     
                  "Mean_diurnal_range",   
                  "Isothermality",        
                  "Temp_seasonality",     
                  "Max_temp_warm_month",  
                  "Min_temp_cold_month",  
                  "Temp_annual_range",    
                  "Mean_temp_wet_qu",     
                  "Mean_temp_dry_qu",     
                  "Mean_temp_warm_qu",    
                  "Mean_temp_cold_qu",
                  
                  "Annual_precip",        
                  "Precip_wet_month",     
                  "Precip_dry_month",     
                  "Precip_seasonality",   
                  "Precip_wet_qu",        
                  "Precip_dry_qu",        
                  "Precip_warm_qu",       
                  "Precip_col_qu")


#########################################################################################################################
## Change the raster values here: 
## See http://worldclim.org/formats1 for description of the interger conversion. 
## All temperature variables wer multiplied by 10, so divide by 10 to reverse it.
COMBO.RASTER.CONVERT = as.data.table(COMBO.RASTER)                           ## This is inefficient
COMBO.RASTER.CONVERT[, (env.variables[c(1:11)]) := lapply(.SD, function(x) 
  x / 10 ), .SDcols = env.variables[c(1:11)]]
COMBO.RASTER.CONVERT = as.data.frame(COMBO.RASTER.CONVERT)                   ## Find another method without using data.table


## Check Looks ok?
summary(COMBO.RASTER.CONVERT$Annual_mean_temp)  ## -23 looks too low/. Check where these are ok
summary(COMBO.RASTER$Annual_mean_temp)

summary(COMBO.RASTER.CONVERT$Isothermality)
summary(COMBO.RASTER$Isothermality)


## Plot a few points to see
plot(LAND, col = 'grey', bg = 'sky blue')
points(COMBO.RASTER.CONVERT[ which(COMBO.RASTER.CONVERT$Annual_mean_temp < -5), ][, c("lon", "lat")], 
       pch = ".", col = "red", cex = 3, asp = 1)



#########################################################################################################################
## Create niche summaries for each environmental condition like this...
## Here's what the function will produce :
library(plyr)
head(niche_estimate (DF = COMBO.RASTER.CONVERT, colname = "Annual_mean_temp"))  ## including the q05 and q95


## So lets use lapply on the "Search Taxon". Note additonal flags are needed, and the taxonomic lists need to be managed better...
## test = run_function_concatenate(list, DF, "DF, colname = x") 
COMBO.NICHE <- env.variables[c(1:length(env.variables))] %>% 
  
  ## Pipe the list into lapply
  lapply(function(x) {
    
    ## Now use the niche width function on each colname (so 8 environmental variables)
    ## Also, need to figure out how to make the aggregating column generic (species, genus, etc.)
    ## currently it only works hard-wired
    niche_estimate (DF = COMBO.RASTER.CONVERT, colname = x)
    
    ## would be good to remove the duplicate columns here
    
  }) %>% 
  
  ## finally, create one dataframe for all niches
  as.data.frame


## Remove duplicate Taxon columns and check the output
names(COMBO.NICHE)
drops <- c("searchTaxon.1",  "searchTaxon.2",  "searchTaxon.3",  "searchTaxon.4",
           "searchTaxon.5",  "searchTaxon.6",  "searchTaxon.7",  "searchTaxon.7",
           "searchTaxon.8",  "searchTaxon.9",  "searchTaxon.10", "searchTaxon.11",
           "searchTaxon.12", "searchTaxon.13", "searchTaxon.14", "searchTaxon.15",
           "searchTaxon.16", "searchTaxon.17", "searchTaxon.18")
COMBO.NICHE = COMBO.NICHE[ , !(names(COMBO.NICHE) %in% drops)]


## Add counts for each species, and record the total number of taxa processed
COMBO.count = as.data.frame(table(COMBO.RASTER.CONVERT$searchTaxon))$Freq
Total.taxa.processed = dim(COMBO.NICHE)[1]
COMBO.NICHE  = cbind(COMBO.count, COMBO.NICHE)
names(COMBO.NICHE)





#########################################################################################################################
## 5). CALCULATE AREA OF OCCUPANCY RANGES 
#########################################################################################################################


## These numbers don't look that accurate: Try converting data into a sp data frame and projecting into a projected
## coordinate system
## Create a species list: this should be from the "COMBO.RASTER" file...
spp.geo = as.character(unique(COMBO.RASTER$searchTaxon)) 
data    = COMBO.RASTER


#########################################################################################################################
## AREA OF OCCUPANCY 


## For every species in the list: calculate the AOO
GBIF.AOO <- spp.geo[c(1:length(spp.geo))] %>% 
  
  ## Pipe the list into lapply
  lapply(function(x) {
    
    ## Subset the the data frame 
    DF      = subset(data, searchTaxon == x)[, c("lon", "lat")]
    
    ## Calculate area of occupancy according the the "red" package
    aoo (DF)
    
    ## Warning messages: Ask John if this is a problem
    ## In rgdal::project(longlat, paste("+proj=utm +zone=", zone,  ... :
    ## 3644 projected point(s) not finite
    
  }) %>% 
  
  ## Finally, create one dataframe for all niches
  as.data.frame


# In rgdal::project(longlat, paste("+proj=utm +zone=", zone,  ... :
#                                       6 projected point(s) not finite


#########################################################################################################################
## EXTENT OF OCCURRENCE
## For every species in the list: calculate the EOO
# GBIF.EOO <- spp.geo[c(1:length(spp.geo))] %>% 
#   
#   ## Pipe the list into lapply
#   lapply(function(x) {
#     
#     ## Subset the the data frame 
#     DF      = subset(data, searchTaxon == x)[, c("searchTaxon", "lon", "lat")]
#     DF.GEO  = dplyr::rename(DF, 
#                             identifier = searchTaxon,
#                             XCOOR      = lon,
#                             YCOOR      = lat)
#     
#     ## Calculate area of occupancy according the the "red" package
#     CalcRange (DF.GEO)
#     
#     ## Warning messages: Ask John if this is a problem
#     ## In rgdal::project(longlat, paste("+proj=utm +zone=", zone,  ... :
#     ## 3644 projected point(s) not finite
#     
#   }) %>% 
#   
#   ## Finally, create one dataframe for all niches
#   as.data.frame


## Clean it up. The order of species should be preserved
GBIF.AOO = gather(GBIF.AOO)
str(GBIF.AOO)
str(unique(GBIF.ALA.COMBO.LAND$searchTaxon))   ## same number of species...


## Now join on the GEOGRAPHIC RANGE
COMBO.NICHE$AREA_OCCUPANCY = GBIF.AOO$value    ## vectors same length so don't need to match


## AOO is calculated as the area of all known or predicted cells for the species. The resolution will be 2x2km as 
## required by IUCN. A single value in km2.


#########################################################################################################################
## Add the counts of LGAs for each species in here
COMBO.LGA = cbind.data.frame(COMBO.NICHE, LGA.AGG) ## The tapply needs to go where the niche summaries are
names(COMBO.LGA)





#########################################################################################################################
## 5). JOIN ON CONTEXTUAL DATA
#########################################################################################################################


#########################################################################################################################
## Now join the horticultural contextual data onto one or both tables ()
COMBO.RASTER.CONTEXT = join(COMBO.RASTER.CONVERT, HIA.SPP.JOIN, 
                            by = "searchTaxon", type = "left", match = "all")


## Now join hort context to all the niche
COMBO.NICHE.CONTEXT = join(COMBO.LGA, HIA.SPP.JOIN, 
                           by = "searchTaxon", type = "left", match = "all")


#########################################################################################################################
## For pedantry, reroder columns...
names(COMBO.RASTER.CONTEXT)
names(COMBO.NICHE.CONTEXT)


## 
save(COMBO.RASTER.CONTEXT, file = paste("./data/base/HIA_LIST/COMBO/COMBO_RASTER_CONTEXT.RData", sep = ""))
save(COMBO.NICHE.CONTEXT,  file = paste("./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT.RData",  sep = ""))
write.csv(COMBO.NICHE.CONTEXT, "./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT.csv",       row.names = FALSE)


## Rename...
COMBO.RASTER.CONTEXT = COMBO.RASTER.CONTEXT[, c(1:5,  45:57, 6:43)]                                       ## change order
COMBO.NICHE.CONTEXT  = COMBO.NICHE.CONTEXT[,  c(176, 2, 1, 174:175, 177:189, 3:173)]                      ## change order


## Set NA to blank, then sort by no. of growers.
## The extra species are ones I dowloaded in error
COMBO.NICHE.CONTEXT$Number.of.growers[is.na(COMBO.NICHE.CONTEXT$Number.of.growers)] <- 0
COMBO.NICHE.CONTEXT = COMBO.NICHE.CONTEXT[with(COMBO.NICHE.CONTEXT, rev(order(Number.of.growers))), ]


## View the data
names(COMBO.RASTER.CONTEXT)
names(COMBO.NICHE.CONTEXT)
dim(COMBO.RASTER.CONTEXT)
dim(COMBO.NICHE.CONTEXT)

View(COMBO.RASTER.CONTEXT)
View(COMBO.NICHE.CONTEXT)


#########################################################################################################################
## quickly check how many species match from the original 605. Only 553 are currently there.


## Which species from those with > 25 growers are missing?
missing.25     = setdiff(unique(HIA.SPP.JOIN[ which(HIA.SPP.JOIN$Number.of.growers >= 25), ][["searchTaxon"]]),
                         unique(COMBO.NICHE.CONTEXT[ which(COMBO.NICHE.CONTEXT$Number.of.growers >= 25), ][["searchTaxon"]]))


## Same as for 
missing.all    = setdiff(unique(HIA.SPP.JOIN[["searchTaxon"]]),
                         unique(COMBO.NICHE.CONTEXT[["searchTaxon"]]))

## Which species from the top 200 are missing?
missing.200    = setdiff(unique(spp.200$Binomial),
                         unique(COMBO.NICHE.CONTEXT[ which(COMBO.NICHE.CONTEXT$Number.of.growers >= 25), ][["searchTaxon"]]))


## Which species from Renee's list are missing?
missing.renee  = setdiff(unique(renee.50$Species),
                         unique(COMBO.NICHE.CONTEXT[["searchTaxon"]]))

## Overall list
missing.taxa   = unique(c(missing.25, missing.200, missing.renee))
missing.taxa   = as.data.frame(missing.taxa)
missing.taxa


## The missing species are due to too few records, too many, or taxonomy problems. EG some of the species are varieties, so they 
## only match to the genus. So 605 - 30 = 575. What is the difference?
       

## Save the summary datasets
save(missing.taxa,         file = paste("./data/base/HIA_LIST/COMBO/MISSING_TAXA.RData",         sep = ""))
save(COMBO.RASTER.CONTEXT, file = paste("./data/base/HIA_LIST/COMBO/COMBO_RASTER_CONTEXT.RData", sep = ""))
save(COMBO.NICHE.CONTEXT,  file = paste("./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT.RData",  sep = ""))
write.csv(COMBO.NICHE.CONTEXT, "./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT.csv",       row.names = FALSE)


## Now save .RData file for the next session...
save.image("STEP_4_NICHES.RData")





#########################################################################################################################
## OUTSTANDING NICHE TASKS:
#########################################################################################################################




## Run the niche calculations again for culivated records, as well as native vs. non-native if possible

## Check WORLDCLIM values: some of the numbers don't look right

## Check geographic range: doesn't look right for some species. Calc extent of occurrnece as well

## Improve raster extract: use an index of unique cell values, referring to a second matrix to do the extract 

## Return species EG:                                     -

## Find infrequently sold spp., big environmental & geographic range, but could have similar traits to popular species?

## Find rarest species (are there popular species with not many records?)

## Ask Linda and others about the kind of queries they want to run...


#########################################################################################################################
#####################################################  END ############################################################## 
#########################################################################################################################