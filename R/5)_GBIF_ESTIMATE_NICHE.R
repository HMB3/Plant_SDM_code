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
  

# Growth rate and form
# Height
# Canopy density
# Ground cover
# Longevity
# Seasonality
# Water quality
# Allergenicity
# Air and water quality influences and urban temperatures
# Insect resistance
# Ornamental and amenity features
# Biodiversity impacts.


## First step is to estimate the current niche, using the best available data. This will give us a broad indication of the 
## currernt climatic toleance of each species.

## There are two streams here: 

## the broad niches for each species (ie. macroscale env and bio data)
## the bespoke microscale approach (e.g. ecological framework, solar surfaces, etc.)

## The broad niche stream informs the micro stream, based on two tables:

## 1). All records, contextual data and environmnetal conditions
## 2). One row for each species, contextual data and species attributes (niches, traits, etc.)





#########################################################################################################################
## 1). SUMMARISE NICHES GRAPHICALLY AND NUMERCIALLY FOR SELECTED TAXA - START WITH THE TOP 200
#########################################################################################################################


## Consider some simple outputs for each species, and each environmental condition, create:

## A map of global and Australian GBIF ocurrences
## Histograms with a bar chart?
## A 3-d representation of the niche for selected variables (e.g. niche volume/hull, etc.)
## Other stuff?


#########################################################################################################################
## Many ways to break this down. Could create a list, then loop/apply over this
## The top 200 list does not have the taxonomic errors corrected
gc()
Top.200.test = GBIF.RASTER.CONTEXT[ which(GBIF.RASTER.CONTEXT$Top_200 == "TRUE"), ]
Top.200.map  = unique(Top.200.test[["searchTaxon"]]) ## [1:10]
All.spp.map  = unique(GBIF.RASTER.CONTEXT[["searchTaxon"]])
taxa.n       = "Magnolia grandiflora"


# ## example taxa
# dim(GBIF.RASTER.CONTEXT[ which(GBIF.RASTER.CONTEXT$searchTaxon == taxa.n), ])
# dim(GBIF.RASTER.CONTEXT[ which(GBIF.RASTER.CONTEXT$searchTaxon == taxa.n 
#                                & GBIF.RASTER.CONTEXT$country == "Australia"), ][18:17])
# 
# 
# ## Global map
# plot(LAND, col = 'grey', bg = 'sky blue')
# title(paste0("Global occurrences for ", taxa.n))
# points(GBIF.RASTER.CONTEXT[ which(GBIF.RASTER.CONTEXT$searchTaxon == taxa.n), ][18:17], 
#        pch = ".", col = "red", cex = 1.5)
# 
# 
# ## Australian map
# plot(GBIF.RASTER.CONTEXT[ which(GBIF.RASTER.CONTEXT$searchTaxon == taxa.n 
#                                 & GBIF.RASTER.CONTEXT$country == "Australia"), ][18:17], 
#      pch = ".", cex = 5, col = "red")
# title(paste0("Australian occurrences for ", taxa.n))
# plot(LAND, add = T)


#########################################################################################################################
## Now try using a mapping function
source('./R/GREEN_CITIES_FUNCTIONS.R')
map_GBIF_records(taxa.list = All.spp.map)
taxa.n = "Magnolia grandiflora"

## And try using a histogram function. Is it possible to run multiple environmental variables at once? EG the arguments 
## are just set up for one variable and units, could do multiples by having env.1, env.2, etc...


## Also consider how to combine outputs?
histogram_GBIF_records(taxa.list = Top.200.map, env.var.1 = "Annual_mean_temp",   env.col.1 = "orange",     env.units.1 = "°K",
                       env.var.2 = "Annual_precip",   env.col.2 = "sky blue",     env.units.2 = "mm")

## the whole list?
histogram_GBIF_records(taxa.list = All.spp.map, env.var.1 = "Annual_mean_temp",   env.col.1 = "orange",     env.units.1 = "°K",
                       env.var.2 = "Annual_precip",   env.col.2 = "sky blue",     env.units.2 = "mm")



#########################################################################################################################
## How bout a a numerical summary for one species? So this is like slicing

## create table for each species?
GBIF_summary_slice(taxa.list = All.spp.map, 
                   env.cols  = c("searchTaxon",             "Annual_mean_temp_min", 
                                 "Annual_mean_temp_median", "Annual_mean_temp_max", 
                                 "Annual_mean_temp_range"), 
                   GBIF      = GBIF.NICHE.CONTEXT)


## ideally, what we want is one summary per species, that has all their particulars.
## Maps, histograms, numerical summaries...





#########################################################################################################################
#####################################################  END ############################################################## 
#########################################################################################################################


