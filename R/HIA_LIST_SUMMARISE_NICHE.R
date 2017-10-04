#########################################################################################################################
#######################################  SUMMARY OF SPECIES GLOBAL NICHES ############################################### 
#########################################################################################################################


#########################################################################################################################
## load packages and source functions...should not need any packages to do run this simple code...
source('./R/GREEN_CITIES_FUNCTIONS.R')
#library(plyr)


#########################################################################################################################
## Load two tables: note the GBIF records need further cleaning
load("./data/base/HIA_LIST/COMBO/COMBO_RASTER_CONTEXT.RData")   ## All the environmental data, one row for each record
load("./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT.RData")    ## The niches for each variable, one row for each species
renee.50   = read.csv("./data/base/HIA_LIST/HIA/RENEE_TOP_50.csv", stringsAsFactors = FALSE)  ## Renee's list


## Check
View(COMBO.RASTER.CONTEXT)
View(COMBO.NICHE.CONTEXT)


#########################################################################################################################
## WORLDCLIM VARIABLES: WHICH ARE MOST RELEVANT TO THE EXPERIMENTS?
#########################################################################################################################


## I have summarised all these variables: 

# CODE    NAME                                                          MY NAME
# BIO1  = Annual Mean Temperature                                    ## moslty self explanatory                                
# BIO2  = Mean Diurnal Range (Mean of monthly (max temp - min temp))  
# BIO3  = Isothermality (BIO2/BIO7) (* 100)
# BIO4  = Temperature Seasonality (standard deviation *100)           
# BIO5  = Max Temperature of Warmest Month                            
# BIO6  = Min Temperature of Coldest Month                             
# BIO7  = Temperature Annual Range (BIO5-BIO6)
# BIO8  = Mean Temperature of Wettest Quarter                        ## Mean_temp_wet_qu, etc.
# BIO9  = Mean Temperature of Driest Quarter
# BIO10 = Mean Temperature of Warmest Quarter
# BIO11 = Mean Temperature of Coldest Quarter
# BIO12 = Annual Precipitation                                        
# BIO13 = Precipitation of Wettest Month                              
# BIO14 = Precipitation of Driest Month                               
# BIO15 = Precipitation Seasonality (Coefficient of Variation)        
# BIO16 = Precipitation of Wettest Quarter
# BIO17 = Precipitation of Driest Quarter
# BIO18 = Precipitation of Warmest Quarter
# BIO19 = Precipitation of Coldest Quarter                           ## Precip_col_qu


## For each WORLDCLIM variable, for all GBIF records I have calculated the:
## _min
## _max
## _median
## _median 
## _range
## _q05 : 5th percentile
## _q95 : 95th percentile
## _q95_q05 : 95th percentile - the 5th
## _q98_q02 : 98th percentile - 2nd


## Just let me know if there are other measures you would like. Not all of these make sense (e.g. min of the min).
## But it's easier to just do them all, and ignore the ones we don't need!
names(COMBO.NICHE.CONTEXT)


## Also, note I need to do more cleaning of the GBIF records, and also eliminate the duplicates between the ALA
## and GBIF. It's easy to re-create the niches, but best to treat these niches as an overestimate for now





#########################################################################################################################
## 1). SUMMARISE NICHES
#########################################################################################################################


#########################################################################################################################
## Slice the big and small dataframes to just the ones on Renee's list.
RENEE.SPP          = as.character(unique(renee.50$Species))
GBIF.RASTER.RENEE  = COMBO.RASTER.CONTEXT[COMBO.RASTER.CONTEXT$searchTaxon %in% renee.50$Species, ]
GBIF.NICHE.RENEE   = COMBO.NICHE.CONTEXT[COMBO.NICHE.CONTEXT$searchTaxon %in% renee.50$Species, ]


## 
GBIF.200.RENEE     = GBIF.RASTER.RENEE[ which(GBIF.RASTER.RENEE$Top_200 == "TRUE"), ]
All.spp.map        = unique(GBIF.RASTER.RENEE[["searchTaxon"]])
taxa.n             = "Ficus brachypoda"                  ## EG...


## Have a look 
View(GBIF.NICHE.RENEE)


#########################################################################################################################
## Can use a simple function to slice the niche dataframe to specific columns for one species
GBIF_summary_slice(taxa.list = taxa.n,                       ## Could be a list
                   env.cols  = c("searchTaxon", 
                                 "COMBO.count",
                                 "AREA_OCCUPANCY",
                                 "Number.of.growers",
                                 "Top_200",
                                 "Annual_mean_temp_q05",      
                                 "Annual_mean_temp_q95",     
                                 "Annual_mean_temp_q95_q05"),  
                   GBIF      = COMBO.NICHE.CONTEXT)


## AOO is calculated as the area of all known or predicted cells for the species. The resolution will be 2x2km as required by IUCN.
## A single value in km2.


#########################################################################################################################
## Then plot global map
## note that this still leaves lots of points in the Islands. So we need to decide if those are legitimate.
WORLD <- readOGR("./data/base/CONTEXTUAL/TM_WORLD_BORDERS-0.3.shp", layer = "TM_WORLD_BORDERS-0.3")
LAND  <- readOGR("./data/base/CONTEXTUAL/ne_10m_land.shp", layer = "ne_10m_land")


## Plot global and Australian occurrences for all taxa on the list
## Might need to make plot window bigger
print_occurrence_records(taxa.list = RENEE.SPP, DF = GBIF.RASTER.RENEE)


#########################################################################################################################
## And plot the histograms. Consider how to combine outputs. Might need to make plot window bigger
Print_global_histogram(taxa.list    = RENEE.SPP[1:32], DF = GBIF.RASTER.RENEE,  ## 33 is a problem: Cupianopsis anacardiodes
                       env.var.1    = "Annual_mean_temp",   
                       env.col.1    = "orange",  
                       env.units.1  = "°K",
                       env.var.2    = "Annual_precip",   
                       env.col.2    = "sky blue",     
                       env.units.2  = "mm")


## Do these distributions look sensible? What visual/numerical outputs would be more useful for the other modules?





#########################################################################################################################
## 2). QUERY SPECIES TO FIND NEW ONES
#########################################################################################################################


## Find infrequently sold spp, big environmental & geographic range, but could have similar traits to popular species
summary(COMBO.NICHE.CONTEXT$AREA_OCCUPANCY)
summary(COMBO.NICHE.CONTEXT$Annual_mean_temp_range)
summary(COMBO.NICHE.CONTEXT$Number.of.growers)
summary(COMBO.NICHE.CONTEXT$COMBO.count)


## Rare species we can't model?
rare.spp = subset(COMBO.NICHE.CONTEXT, COMBO.count < 100)[ c("searchTaxon", "Top_200", "Origin")]


## Potential new species
new.spp = subset(COMBO.NICHE.CONTEXT, AREA_OCCUPANCY > 7000 & 
                   Annual_mean_temp_range > 18 & 
                   Number.of.growers < 25)[ c("searchTaxon", "Origin")]


## Other queries?





#########################################################################################################################
## OUTSTANDING NICHE TASKS:
#########################################################################################################################


## Clean the GBIF data and merge with ALA to avoid duplicates, spatial outliers, etc.

## Check on the missing species

## Improve mapping functions to be more useful: Combine maps and historgrams in one

## Correct the worldclim raster values

## Use different rasters too: Manuel?

## Lots more...



#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################