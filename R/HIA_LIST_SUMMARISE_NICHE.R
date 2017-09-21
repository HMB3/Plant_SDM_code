#########################################################################################################################
#######################################  SUMMARY OF SPECIES GLOBAL NICHES ############################################### 
#########################################################################################################################


#########################################################################################################################
## load packages and source functions...should not need any packages to do run this simple code...
source('./R/GREEN_CITIES_FUNCTIONS.R')
library(plyr)


#########################################################################################################################
## Load two tables: note the GBIF records need further cleaning
load("./data/base/HIA_LIST/GBIF/GBIF_RASTER_CONTEXT.RData")   ## All the environmental data, one row for each record
load("./data/base/HIA_LIST/GBIF/GBIF_NICHE_CONTEXT.RData")    ## The niches for each summary, one row for each species
renee.50   = read.csv("./data/base/HIA_LIST/HIA/RENEE_TOP_50.csv", stringsAsFactors = FALSE)  ## Renee's list
## load in LAND map here

## Check
View(GBIF.RASTER.CONTEXT)
View(GBIF.NICHE.CONTEXT)


#########################################################################################################################
## WORLDCLIM VARIABLES: WHICH ARE MOST RELEVANT TO THE EXPERIMENTS?
#########################################################################################################################


## I have summarise them all:

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
## _q95_q05 : 95th percentile - the 5th
## _q98_q02 : 98th percentile - 2nd
## They probabaly wants the top and bottom quintiles themselves: so thats the q95 and q05

## Just let me know if there are other measures you want. Not all of these make sense (e.g. min of the min).
## But it's easier to just do them all and ignore the ones we don't need.
names(GBIF.NICHE.CONTEXT)


#########################################################################################################################
## Now slice the big and small dataframes to just the ones on Renee's list.
RENEE.SPP          = as.character(renee.50$Species)
GBIF.RASTER.RENEE  = GBIF.RASTER.CONTEXT[GBIF.RASTER.CONTEXT$searchTaxon %in% renee.50$Species, ]
GBIF.NICHE.RENEE   = GBIF.NICHE.CONTEXT[GBIF.NICHE.CONTEXT$searchTaxon %in% renee.50$Species, ]

## 
GBIF.200.RENEE     = GBIF.RASTER.RENEE[ which(GBIF.RASTER.RENEE$Top_200 == "TRUE"), ]
All.spp.map        = unique(GBIF.RASTER.RENEE[["searchTaxon"]])
taxa.n             = RENEE.SPP[1]                  ## First species on the list...


## Have a look 
View(GBIF.NICHE.RENEE)


#########################################################################################################################
## Can use a simple function to slice the niche dataframe to specific columns for one species
GBIF_summary_slice(taxa.list = taxa.n,                       ## Could be a list
                   env.cols  = c("searchTaxon", 
                                 "Number.of.growers",
                                 "Top_200",
                                 "Annual_mean_temp_min",     ## Change these using columns above
                                 "Annual_mean_temp_median",  ## Should be q95, q05, etc...
                                 "Annual_mean_temp_max",     
                                 "Annual_mean_temp_range"),  
                   GBIF      = GBIF.NICHE.CONTEXT)


#########################################################################################################################
## Then plot global map
plot(LAND, col = 'grey', bg = 'sky blue')
title(paste0("Global occurrences for ", taxa.n))

points(GBIF.RASTER.RENEE[ which(GBIF.RASTER.RENEE$searchTaxon == taxa.n), ][, c("lon", "lat")],
       pch = ".", col = "red", cex = 3)


## Plot Australian map
plot(GBIF.RASTER.RENEE[ which(GBIF.RASTER.RENEE$searchTaxon == taxa.n
                                & GBIF.RASTER.RENEE$country == "Australia"), ][, c("lon", "lat")],
     pch = ".", cex = 5, col = "red", asp = 1)
title(paste0("Australian occurrences for ", taxa.n))
plot(LAND, add = TRUE, asp = 1)


#########################################################################################################################
## And plot the histograms
## Also consider how to combine outputs?
Print_global_histogram(taxa.list = RENEE.SPP, DF = GBIF.RASTER.RENEE, env.var.1 = "Annual_mean_temp",   env.col.1 = "orange",  env.units.1 = "°K",
                       env.var.2 = "Annual_precip",   env.col.2 = "sky blue",     env.units.2 = "mm")

## the whole list?
histogram_GBIF_records(taxa.list = RENEE.SPP, env.var.1 = "Annual_mean_temp",   env.col.1 = "orange",     env.units.1 = "°K",
                       env.var.2 = "Annual_precip",   env.col.2 = "sky blue",     env.units.2 = "mm")



#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################