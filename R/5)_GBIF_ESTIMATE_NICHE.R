#########################################################################################################################
############################################  VISUALISE SPECIES NICHES ################################################## 
#########################################################################################################################


#########################################################################################################################
## 1). ADD RASTER DATA TO GBIF RECORDS
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
  

# growth rate and form
# height
# canopy density
# ground cover
# longevity
# seasonality
# water quality
# allergenicity
# air and water quality influences and urban temperatures
# insect resistance
# ornamental and amenity features
# and biodiversity impacts.


## First step is to estimate the current niche, using the best available data. This will give us a broad indication of the 
## currernt climatic toleance of each species.

## There are two streams here: 

## the broad niches for each species (ie. macroscale env and bio data)
## the bespoke microscale approach (e.g. ecological framework, solar surfaces, etc.)

## The broad niche stream informs the micro stream, based on two tables:

## 1). All records, contextual data and environmnetal conditions
## 2). One row for each species, contextual data and species attributes (niches, traits, etc.)





#########################################################################################################################
## 1). SUMMARISE NICHES GRAPHICALLY FOR SELECTED TAXA - START WITH THE TOP 200
#########################################################################################################################


## Consider some simple outputs for each species, and each environmental condition, create:

## A map of global and Australian GBIF ocurrences
## Histograms with a bar chart?
## A 3-d representation of the niche for selected variables (e.g. niche volume/hull, etc.)
## Other stuff?


#########################################################################################################################
## Many ways to break this down. Could create a list, then loop/apply over this
Top.200.test = DRAFT.HIA.TAXA.200$searchTaxon[1:10]
taxa.n = "Murraya paniculata"


## example taxa
dim(GBIF.RASTER.CONTEXT[ which(GBIF.RASTER.CONTEXT$searchTaxon == taxa.n), ])
dim(GBIF.RASTER.CONTEXT[ which(GBIF.RASTER.CONTEXT$searchTaxon == taxa.n 
                               & GBIF.RASTER.CONTEXT$country == "Australia"), ][18:17])


## Global map
plot(LAND)
title(paste0("Global occurrences for ", taxa.n))
points(GBIF.RASTER.CONTEXT[ which(GBIF.RASTER.CONTEXT$searchTaxon == taxa.n), ][18:17], 
       pch = ".", col = "red", cex = 1.5)


## Australian map
plot(GBIF.RASTER.CONTEXT[ which(GBIF.RASTER.CONTEXT$searchTaxon == taxa.n 
                                & GBIF.RASTER.CONTEXT$country == "Australia"), ][18:17], 
     pch = ".", cex = 5, col = "red")
title(paste0("Australian occurrences for ", taxa.n))
plot(LAND, add = T)


#########################################################################################################################
## Now try using a mapping function
source('./R/GREEN_CITIES_FUNCTIONS.R')
map_GBIF_records(Top.200.test)

# Warning messages:
#   1: In if (dim(GBIF.RASTER.CONTEXT[which(GBIF.RASTER.CONTEXT$searchTaxon ==  ... :
#                                           the condition has length > 1 and only the first element will be used


## And try using a histogram function
## Is it possible to run multiple environmental variables at once? 
## EG the arguments are just set up for one variable and units, could do multiplese by having env.1, env.2, etc...


##
histogram_GBIF_records(taxa.list = Top.200.test, env.var = "Annual_mean_temp",   env.col = "orange", env.units = "°K")
histogram_GBIF_records(taxa.list = Top.200.test, env.var = "Annual_mean_precip", env.col = "blue",   env.units = "mm") ## etc


## Use lattice histograms




#########################################################################################################################
#####################################################  END ############################################################## 
#########################################################################################################################


