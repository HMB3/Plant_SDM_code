#########################################################################################################################
#######################################  CREATE SPECIES NICHES FOR HIA LIST ############################################# 
#########################################################################################################################


#########################################################################################################################
## 1). ADD RASTER DATA TO GBIF RECORDS
#########################################################################################################################


## First, consider the HIA brief again:
 
## Our research might demonstrate, for example, that a particular species of tree is already at the very limit of its ability 
## to cope with heat, and that the only suitable place to plant this species in the future will be in cool-climate or more 
## temperate locations.

## First step is to estimate the current niche, using the best available data. This will give us a broad indication of the 
## currernt climatic toleance of each species.


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
## create a stack of rasters to sample
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


## Focus on the top 200 taxa to begin with: some of the top 200 species have not been downloaded yet
head(DRAFT.HIA.TAXA)
DRAFT.HIA.TAXA.200 = subset(DRAFT.HIA.TAXA, Top_200 == "TRUE")


## sort by no. of growers?
DRAFT.HIA.TAXA.200 = DRAFT.HIA.TAXA.200[with(DRAFT.HIA.TAXA.200, rev(order(Number.of.growers))), ] 
dim(DRAFT.HIA.TAXA.200)
head(DRAFT.HIA.TAXA.200)


## now summarise the niches. But figure out a cleaner way of doing this?
variables = c("Annual_mean_temp", "Temp_seasonality", "Max_temp_warm_month", "Min_temp_cold_month",
              "Annual_precip",    "Precip_Wet_month", "Precip_dry_month",    "Precip_seasonality")


#########################################################################################################################
## We want to create niche summaries for each environmental condition like this...
head(niche_estimate (DF = GBIF.RASTER, colname = "Annual_mean_temp"))


## So lets use lapply!
GBIF.NICHE <- variables[c(1:length(variables))] %>% 
  
  ## pipe the list into lapply
  lapply(function(x) {
    
    ## now use the niche width function 
    niche_estimate (DF = GBIF.RASTER, colname = x)
    
    ## would be good to remove the duplicates here
    
  }) %>% 
  
  ## finally, create one dataframe for all niches
  as.data.frame


## check the output
str(GBIF.NICHE)


## Add counts for each species
count = as.data.frame(table(GBIF.RASTER$searchTaxon))$Freq
test = cbind(count, annual.temp.niche.est)
head(annual.temp.niche.est)








