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
## WHICH WORLDCLIM VARIABLES TO ESTIMATE?


## copy the ones which Rach and Stu have used for the niche finder website
## for now, ignore edaphic variables

# BIO1  = Annual Mean Temperature                                     ##
# BIO2  = Mean Diurnal Range (Mean of monthly (max temp - min temp))
# BIO3  = Isothermality (BIO2/BIO7) (* 100)
# BIO4  = Temperature Seasonality (standard deviation *100)           ##
# BIO5  = Max Temperature of Warmest Month                            ##
# BIO6  = Min Temperature of Coldest Month                            ##
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
              "//sci-7910/F/data/worldclim/world/0.5/bio/current/bio_04",
              "//sci-7910/F/data/worldclim/world/0.5/bio/current/bio_05",        
              "//sci-7910/F/data/worldclim/world/0.5/bio/current/bio_06",
              "//sci-7910/F/data/worldclim/world/0.5/bio/current/bio_12",
              "//sci-7910/F/data/worldclim/world/0.5/bio/current/bio_13",
              "//sci-7910/F/data/worldclim/world/0.5/bio/current/bio_14",
              "//sci-7910/F/data/worldclim/world/0.5/bio/current/bio_15")


## convert all the rasters to a stack
s <- stack(env.grids)


## then use the extract function for all the rasters,
## and finaly bind on the GBIF data to the left of the raster values 
GBIF.RASTER <- extract(s, GBIF.POINTS) %>% 
  cbind(GBIF.LAND, .)


## multiple rename?
GBIF.RASTER = rename(GBIF.RASTER,
                     Annual_mean_temp     = bio_01,
                     Temp_seasonality     = bio_04,
                     Max_temp_warm_month  = bio_05,
                     Min_temp_cold_month  = bio_06,
                     Annual_precip        = bio_12,
                     Precip_Wet_month     = bio_13,
                     Precip_dry_month     = bio_14,
                     Precip_seasonality   = bio_15)


## check the data
dim(GBIF.RASTER)
names(GBIF.RASTER) 





#########################################################################################################################
## 2). CREATE NICHES FOR SELECTED TAXA
#########################################################################################################################


#########################################################################################################################
## Focus on the top 200 taxa to begin with: some of the top 200 species have not been downloaded yet
DRAFT.HIA.TAXA.200 = subset(DRAFT.HIA.TAXA, Top_200 == "TRUE")
DRAFT.HIA.TAXA.200 = rename(DRAFT.HIA.TAXA.200, searchTaxon = Species)
DRAFT.HIA.TAXA     = rename(DRAFT.HIA.TAXA,     searchTaxon = Species)


## Set NA to blank, then sort by no. of growers?
DRAFT.HIA.TAXA[is.na(DRAFT.HIA.TAXA)] <- 0
DRAFT.HIA.TAXA = DRAFT.HIA.TAXA[with(DRAFT.HIA.TAXA, rev(order(Number.of.growers))), ]
DRAFT.HIA.TAXA.200 = DRAFT.HIA.TAXA.200[with(DRAFT.HIA.TAXA.200, rev(order(Number.of.growers))), ]

## check
head(DRAFT.HIA.TAXA[, c("searchTaxon", "Number.of.growers")])
View(DRAFT.HIA.TAXA)


## save list to file for Rachel to check Austraits
drops <- c("Count", "Comment")
DRAFT.HIA.TAXA = DRAFT.HIA.TAXA[ , !(names(DRAFT.HIA.TAXA) %in% drops)]
names(DRAFT.HIA.TAXA)
## write.csv(DRAFT.HIA.TAXA, "./data/base/HIA_LIST/HIA/DRAFT_HIA_TAXA.csv", row.names = FALSE)



## now summarise the niches. But figure out a cleaner way of doing this?
env.variables = c("Annual_mean_temp", "Temp_seasonality", "Max_temp_warm_month", "Min_temp_cold_month",
                   "Annual_precip",    "Precip_Wet_month", "Precip_dry_month",    "Precip_seasonality")


#########################################################################################################################
## We want to create niche summaries for each environmental condition like this...
## Here's what the function will produce :
head(niche_estimate (DF = GBIF.RASTER, colname = "Annual_mean_temp"))
dim(niche_estimate  (DF = GBIF.RASTER, colname = "Annual_mean_temp"))  ## 7 environmetnal summaries so far


## So lets use lapply on the "Search Taxon".
## Note additonal flags are needed, and the taxonomic lists need to be managed better...
GBIF.NICHE <- variables[c(1:length(variables))] %>% 
  
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
drops <- c("searchTaxon.1", "searchTaxon.2", "searchTaxon.3", "searchTaxon.4",
           "searchTaxon.5", "searchTaxon.6", "searchTaxon.7")
GBIF.NICHE = GBIF.NICHE[ , !(names(GBIF.NICHE) %in% drops)]
names(GBIF.NICHE)


## Add counts for each species, and record the total number of taxa processed
count = as.data.frame(table(GBIF.RASTER$searchTaxon))$Freq
Total.taxa.processed = dim(GBIF.NICHE)[1]
GBIF.NICHE  = cbind(count, GBIF.NICHE)
names(GBIF.NICHE)





#########################################################################################################################
## 3). SUMMARISE SPECIES NUMERICALLY: WHICH COULD BE MODELLED?
#########################################################################################################################


#########################################################################################################################
## Now join the horticultural contextual data onto one or both tables ()
## Which columns do we need?
names(GBIF.RASTER)
names(GBIF.NICHE)
names(DRAFT.HIA.TAXA.200)
View(DRAFT.HIA.TAXA)


## Now join hort context to all records. Provisional, keep working on it 
GBIF.RASTER.CONTEXT = join(GBIF.RASTER, DRAFT.HIA.TAXA[, c("searchTaxon", "Plant.type", 
                                                           "Number.of.growers", "Origin", "Top_200")], 
                           by = "searchTaxon", type = "left", match = "all")


## Now join hort context to all the niches. Provisional, keep working on it 
GBIF.NICHE.CONTEXT = join(GBIF.NICHE, DRAFT.HIA.TAXA[, c("searchTaxon", "Plant.type", 
                                                         "Number.of.growers", "Origin", "Top_200")], 
                          by = "searchTaxon", type = "left", match = "all")


#########################################################################################################################
## For pedantry, reroder columns...
## Note that the downloaded species don't all match up to the original list. This is because there are different lists 
## propagating throughout the workflow, I need to watch out for this. 
GBIF.RASTER.CONTEXT = GBIF.RASTER.CONTEXT[, c(1:12, 27:30, 13:26)]
GBIF.NICHE.CONTEXT  = GBIF.NICHE.CONTEXT[, c(2,1, 59:62, 3:58)]


## View the data
View(GBIF.RASTER.CONTEXT)
View(GBIF.NICHE.CONTEXT)

names(GBIF.RASTER.CONTEXT)
names(GBIF.NICHE.CONTEXT)


#########################################################################################################################
## now create a master summary of the recrods so far. What do we need to know?
## These numbers could be variables that update each time code is re-run with changes to the variables
## Also a table could be made for each source (GBIF, ALA, Council, etc.), But ideally it is just one table

## Total number of number of taxa (ie. all variables are dynamic)
## Total number of records (uncleaned, and broken down by source)
## % records retained after filtering 
## 6-figure summary across all taxa (min, max, median)
## number of all taxa with > n records (e.g. +50, +100, +1000, +10,000) - assuming it doesn't matter above certain level
## number of taxa with > n records
summary(GBIF.NICHE.CONTEXT$count)


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

GBIF.RECORD.SUMMARY[1, "Min"]                     = summary(GBIF.NICHE.CONTEXT$count)[1]
GBIF.RECORD.SUMMARY[1, "Max"]                     = summary(GBIF.NICHE.CONTEXT$count)[6]
GBIF.RECORD.SUMMARY[1, "Median"]                  = summary(GBIF.NICHE.CONTEXT$count)[4]

GBIF.RECORD.SUMMARY[1, "All taxa count 50+"]      = sum(GBIF.NICHE.CONTEXT$count >= 50,   na.rm = TRUE)
GBIF.RECORD.SUMMARY[1, "All taxa count 100+"]     = sum(GBIF.NICHE.CONTEXT$count >= 100,  na.rm = TRUE)
GBIF.RECORD.SUMMARY[1, "All taxa count 1000+"]    = sum(GBIF.NICHE.CONTEXT$count >= 1000, na.rm = TRUE)
GBIF.RECORD.SUMMARY[1, "All taxa count 10,000+"]  = sum(GBIF.NICHE.CONTEXT$count >= 10000, na.rm = TRUE)

GBIF.RECORD.SUMMARY[1, "Top taxa count 50+"]      = sum(GBIF.NICHE.CONTEXT$count >= 50    & GBIF.NICHE.CONTEXT$Top_200 == "TRUE", na.rm = TRUE)
GBIF.RECORD.SUMMARY[1, "Top taxa count 100+"]     = sum(GBIF.NICHE.CONTEXT$count >= 100   & GBIF.NICHE.CONTEXT$Top_200 == "TRUE", na.rm = TRUE)
GBIF.RECORD.SUMMARY[1, "Top taxa count 1000+"]    = sum(GBIF.NICHE.CONTEXT$count >= 1000  & GBIF.NICHE.CONTEXT$Top_200 == "TRUE", na.rm = TRUE)
GBIF.RECORD.SUMMARY[1, "Top taxa count 10,000+"]  = sum(GBIF.NICHE.CONTEXT$count >= 10000 & GBIF.NICHE.CONTEXT$Top_200 == "TRUE", na.rm = TRUE)


## visualise this
str(GBIF.RECORD.SUMMARY)
colnames(GBIF.RECORD.SUMMARY)
rownames(GBIF.RECORD.SUMMARY) = "VALUE"

GBIF.RECORD.SUMMARY = as.data.frame(t(GBIF.RECORD.SUMMARY)) ## do melting, etc later
kable(GBIF.RECORD.SUMMARY)


#########################################################################################################################
## Now visualise the distribution of records


## Use lattice histograms
histogram(GBIF.NICHE.CONTEXT$count,
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