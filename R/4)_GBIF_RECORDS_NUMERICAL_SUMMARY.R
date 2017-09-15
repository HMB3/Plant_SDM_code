#########################################################################################################################
####################################### ESTIMATE SPECIES NICHES AND SUMMARISE RECORDS ################################### 
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

## the broad niches for each species (ie. macroscale env and bio data)
## the bespoke microscale approach (e.g. ecological framework, solar surfaces, etc.)

## The broad niche stream informs the micro stream, based on two tables:

## 1). All records, contextual data and environmnetal conditions
## 2). One row for each species, contextual data and species attributes (niches, traits, etc.)





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


## to save time, load in previous data
load("./data/base/HIA_LIST/GBIF/GBIF_LAND_POINTS.RData")
str(GBIF.LAND)


## create points
GBIF.POINTS   = SpatialPointsDataFrame(coords = GBIF.LAND[c("lon", "lat")], 
                                       data   = GBIF.LAND[c("lon", "lat")],
                                       proj4string = CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))


## check
summary(GBIF.POINTS)


#########################################################################################################################
## Create a stack of rasters to sample
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


## convert all the rasters to a stack
s <- stack(env.grids)


## then use the extract function for all the rasters,
## and finaly bind on the GBIF data to the left of the raster values 
GBIF.RASTER <- extract(s, GBIF.POINTS) %>% 
  cbind(GBIF.LAND, .)


## multiple rename?
GBIF.RASTER = dplyr::rename(GBIF.RASTER,
                            Annual_mean_temp     = bio_01,
                            Mean_diurnal_range   = bio_02,
                            Isothermality        = bio_03,
                            Temp_seasonality     = bio_04,
                            Max_temp_warm_month  = bio_05,
                            Min_temp_cold_month  = bio_06,
                            Temp_annual_range    = bio_07,
                            Mean_temp_wet_qu     = bio_08,
                            Mean_temp_dry_qu     = bio_09,
                            Mean_temp_warm_qu    = bio_10,
                            Mean_temp_cold_qu    = bio_11,
                            
                            Annual_precip        = bio_12,
                            Precip_wet_month     = bio_13,
                            Precip_dry_month     = bio_14,
                            Precip_seasonality   = bio_15,
                            Precip_wet_qu        = bio_16,
                            Precip_dry_qu        = bio_17,
                            Precip_warm_qu       = bio_18,
                            Precip_col_qu        = bio_19)


## Save/load
save(GBIF.RASTER, file = paste("./data/base/HIA_LIST/GBIF/GBIF_RASTER.RData"))
#load("./data/base/HIA_LIST/GBIF/GBIF_RASTER.RData")


## check
dim(GBIF.RASTER)
names(GBIF.RASTER)


## Also consider converting degrees Kelvin to Celsius?
# hist(GBIF.RASTER$Annual_mean_temp)
# hist((GBIF.RASTER$Annual_mean_temp)-273.15)





#########################################################################################################################
## 2). CREATE NICHES FOR SELECTED TAXA
#########################################################################################################################


#########################################################################################################################
## Might not need separate top 200 anymore...
library(plyr)  ## detach plyr again if doing more renaming
HIA.SPP.JOIN     = HIA.SPP
HIA.SPP.JOIN     = dplyr::rename(HIA.SPP.JOIN, searchTaxon = Binomial)
#names(HIA.SPP.JOIN)[names(HIA.SPP.JOIN) == "Binomial"] <- "searchTaxon"


## Set NA to blank, then sort by no. of growers to get them to the top
HIA.SPP.JOIN[is.na(HIA.SPP.JOIN)] <- 0
HIA.SPP.JOIN = HIA.SPP.JOIN[with(HIA.SPP.JOIN, rev(order(Number.of.growers))), ]
HIA.SPP.JOIN.200 = HIA.SPP.JOIN.200[with(HIA.SPP.JOIN.200, rev(order(Number.of.growers))), ]
head(HIA.SPP.JOIN[, c("searchTaxon", "Number.of.growers")])
View(HIA.SPP.JOIN)


## save list to file for Rachel to check Austraits
## Also could join the taxon lookup
## write.csv(HIA.SPP.JOIN, "./data/base/HIA_LIST/HIA/HIA_SPP_JOIN.csv", row.names = FALSE)


## Now summarise the niches. But figure out a cleaner way of doing this?
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
## Create niche summaries for each environmental condition like this...
## Here's what the function will produce :
head(niche_estimate (DF = GBIF.RASTER, colname = "Annual_mean_temp"))
dim(niche_estimate  (DF = GBIF.RASTER, colname = "Annual_mean_temp"))  ## 7 environmetnal summaries so far


## So lets use lapply on the "Search Taxon".
## Note additonal flags are needed, and the taxonomic lists need to be managed better...
GBIF.NICHE <- env.variables[c(1:length(env.variables))] %>% 
  
  ## Pipe the list into lapply
  lapply(function(x) {
    
    ## Now use the niche width function on each colname (so 8 environmental variables)
    ## Also, need to figure out how to make the aggregating column generic (species, genus, etc.)
    ## currently it only works hard-wired
    niche_estimate (DF = GBIF.RASTER, colname = x)
    
    ## would be good to remove the duplicates here, somehow
    
  }) %>% 
  
  ## finally, create one dataframe for all niches
  as.data.frame


## Remove duplicate Taxon columns and check the output
names(GBIF.NICHE)
drops <- c("searchTaxon.1",  "searchTaxon.2",  "searchTaxon.3",  "searchTaxon.4",
           "searchTaxon.5",  "searchTaxon.6",  "searchTaxon.7",  "searchTaxon.7",
           "searchTaxon.8",  "searchTaxon.9",  "searchTaxon.10", "searchTaxon.11",
           "searchTaxon.12", "searchTaxon.13", "searchTaxon.14", "searchTaxon.15",
           "searchTaxon.16", "searchTaxon.17")
GBIF.NICHE = GBIF.NICHE[ , !(names(GBIF.NICHE) %in% drops)]


## Add counts for each species, and record the total number of taxa processed
GBIF.count = as.data.frame(table(GBIF.RASTER$searchTaxon))$Freq
Total.taxa.processed = dim(GBIF.NICHE)[1]
GBIF.NICHE  = cbind(GBIF.count, GBIF.NICHE)
names(GBIF.NICHE)





#########################################################################################################################
## 3). SUMMARISE SPECIES NUMERICALLY: WHICH COULD BE MODELLED?
#########################################################################################################################


#########################################################################################################################
## Now join the horticultural contextual data onto one or both tables ()
GBIF.RASTER.CONTEXT = join(GBIF.RASTER, HIA.SPP.JOIN, 
                           by = "searchTaxon", type = "left", match = "all")


## Now join hort context to all the niche
GBIF.NICHE.CONTEXT = join(GBIF.NICHE, HIA.SPP.JOIN, 
                          by = "searchTaxon", type = "left", match = "all")


#########################################################################################################################
## For pedantry, reroder columns...
## Note that the downloaded species don't all match up to the original list. This is because there are different lists 
## propagating throughout the workflow, I need to watch out for this. 
GBIF.RASTER.CONTEXT = GBIF.RASTER.CONTEXT[, c(38, 2, 1, 3,  39:51, 4:37)]
GBIF.NICHE.CONTEXT  = GBIF.NICHE.CONTEXT[,  c(137, 2, 1, 138:150,  3:133)]


## Set NA to blank, then sort by no. of growers.
## The extra species are ones I dowloaded in error
GBIF.NICHE.CONTEXT$Number.of.growers[is.na(GBIF.NICHE.CONTEXT$Number.of.growers)] <- 0
GBIF.NICHE.CONTEXT = GBIF.NICHE.CONTEXT[with(GBIF.NICHE.CONTEXT, rev(order(Number.of.growers))), ]


## View the data
View(GBIF.RASTER.CONTEXT)
View(GBIF.NICHE.CONTEXT)

names(GBIF.RASTER.CONTEXT)
names(GBIF.NICHE.CONTEXT)
dim(GBIF.NICHE.CONTEXT)


## quickly check how many species match from the original 610 
View(HIA.SPP.JOIN)
missing.25 = setdiff(unique(HIA.SPP.JOIN[ which(HIA.SPP.JOIN$Number.of.growers >= 25), ][["searchTaxon"]]),
                     unique(GBIF.NICHE.CONTEXT[ which(GBIF.NICHE.CONTEXT$Number.of.growers >= 25), ][["searchTaxon"]]))
missing.taxa = unique(c(missing.taxa, missing.25))

renee.extra  = setdiff(unique(renee.50$Species),
                       unique(GBIF.NICHE.CONTEXT[ which(GBIF.NICHE.CONTEXT$Number.of.growers >= 25), ][["searchTaxon"]]))
       

## Save the summary datasets
save(GBIF.RASTER.CONTEXT, file = paste("./data/base/HIA_LIST/GBIF/GBIF_RASTER_CONTEXT.RData", sep = ""))
save(GBIF.NICHE.CONTEXT,  file = paste("./data/base/HIA_LIST/GBIF/GBIF_NICHE_CONTEXT.RData",  sep = ""))
write.csv(GBIF.NICHE.CONTEXT, "./data/base/HIA_LIST/GBIF/GBIF_NICHE_CONTEXT.csv",       row.names = FALSE)



#########################################################################################################################
## Now create a master summary of the recrods so far. What do we need to know?
## These numbers could be variables that update each time code is re-run with changes to the variables
## Also a table could be made for each source (GBIF, ALA, Council, etc.), But ideally it is just one table

## Total number of number of taxa (ie. all variables are dynamic)
## Total number of records (uncleaned, and broken down by source)
## % records retained after filtering 
## 6-figure summary across all taxa (min, max, median)
## number of all taxa with > n records (e.g. +50, +100, +1000, +10,000) - assuming it doesn't matter above certain level
## number of top 200 taxa with > n records
summary(GBIF.NICHE.CONTEXT$GBIF.count)


## How many species where knocked out by using filters?
kable(GBIF.PROBLEMS)


#########################################################################################################################
## Which species have more than n records, are on the top 200, etc?
## first create column headings for final results table  
## these match data frame created by bivariate lm function
n.samples = 1
table.columns = c(  
  "Dataset",
  "Taxa processed",
  "Rank",
  "Total records",
  "Filters applied",
  "Cleaned records",
  "% Records retained",
  "Min",  
  "Max",
  "Median",
  "All taxa count 50+",
  "All taxa count 100+",
  "All taxa count 1000+",
  "All taxa count 10,000+",
  "Top taxa count 50+",
  "Top taxa count 100+",
  "Top taxa count 1000+",
  "Top taxa count 10,000+"
)


## Create big matrices with n rows = n samples
## These will be filled by the linear model functions  
table.length                   = length(table.columns)
GBIF.RECORD.SUMMARY            = matrix(0, n.samples, table.length)
colnames(GBIF.RECORD.SUMMARY)  = table.columns 
GBIF.RECORD.SUMMARY            = as.data.frame(GBIF.RECORD.SUMMARY)


## now put the data in. Could use "melt", etc?
GBIF.RECORD.SUMMARY[1, "Dataset"]                 = "GBIF_occurrence"                # H1.ms[1, "Intercept"]
GBIF.RECORD.SUMMARY[1, "Taxa processed"]          = Total.taxa.processed
GBIF.RECORD.SUMMARY[1, "Rank"]                    = "Species"
GBIF.RECORD.SUMMARY[1, "Total records"]           = Total.count
GBIF.RECORD.SUMMARY[1, "Filters applied"]         = Filters.applied
GBIF.RECORD.SUMMARY[1, "Cleaned records"]         = Remaining.records
GBIF.RECORD.SUMMARY[1, "% Records retained"]      = Remaining.percent

GBIF.RECORD.SUMMARY[1, "Min"]                     = summary(GBIF.NICHE.CONTEXT$GBIF.count)[1]
GBIF.RECORD.SUMMARY[1, "Max"]                     = summary(GBIF.NICHE.CONTEXT$GBIF.count)[6]
GBIF.RECORD.SUMMARY[1, "Median"]                  = summary(GBIF.NICHE.CONTEXT$GBIF.count)[4]

GBIF.RECORD.SUMMARY[1, "All taxa count 50+"]      = sum(GBIF.NICHE.CONTEXT$GBIF.count >= 50,   na.rm = TRUE)
GBIF.RECORD.SUMMARY[1, "All taxa count 100+"]     = sum(GBIF.NICHE.CONTEXT$GBIF.count >= 100,  na.rm = TRUE)
GBIF.RECORD.SUMMARY[1, "All taxa count 1000+"]    = sum(GBIF.NICHE.CONTEXT$GBIF.count >= 1000, na.rm = TRUE)
GBIF.RECORD.SUMMARY[1, "All taxa count 10,000+"]  = sum(GBIF.NICHE.CONTEXT$GBIF.count >= 10000, na.rm = TRUE)

GBIF.RECORD.SUMMARY[1, "Top taxa count 50+"]      = sum(GBIF.NICHE.CONTEXT$GBIF.count >= 50    & GBIF.NICHE.CONTEXT$Top_200 == "TRUE", na.rm = TRUE)
GBIF.RECORD.SUMMARY[1, "Top taxa count 100+"]     = sum(GBIF.NICHE.CONTEXT$GBIF.count >= 100   & GBIF.NICHE.CONTEXT$Top_200 == "TRUE", na.rm = TRUE)
GBIF.RECORD.SUMMARY[1, "Top taxa count 1000+"]    = sum(GBIF.NICHE.CONTEXT$GBIF.count >= 1000  & GBIF.NICHE.CONTEXT$Top_200 == "TRUE", na.rm = TRUE)
GBIF.RECORD.SUMMARY[1, "Top taxa count 10,000+"]  = sum(GBIF.NICHE.CONTEXT$GBIF.count >= 10000 & GBIF.NICHE.CONTEXT$Top_200 == "TRUE", na.rm = TRUE)


## visualise this
str(GBIF.RECORD.SUMMARY)
colnames(GBIF.RECORD.SUMMARY)
rownames(GBIF.RECORD.SUMMARY) = "VALUE"

GBIF.RECORD.SUMMARY = as.data.frame(t(GBIF.RECORD.SUMMARY)) ## do melting, etc later
kable(GBIF.RECORD.SUMMARY)


#########################################################################################################################
## Now visualise the distribution of records
histogram(GBIF.NICHE.CONTEXT$GBIF.count,
          breaks = 50, border = NA, col = "grey",
          xlab = "Number of GBIF records per species", 
          main = "Distribution of GBIF records for HIA species")


## remove older files
rm(GBIF.CLEAN)
rm(GBIF.LAND)
rm(GBIF.TEST)
rm(GBIF.RASTER)
gc()





#########################################################################################################################
#####################################################  END ############################################################## 
#########################################################################################################################