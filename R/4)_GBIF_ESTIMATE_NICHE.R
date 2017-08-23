#########################################################################################################################
#######################################  CREATE SPECIES NICHES FOR HIA LIST ############################################# 
#########################################################################################################################


#########################################################################################################################
## 1). ADD RASTER DATA TO RECORDS
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
# BIO1 = Annual Mean Temperature                                     ##
# BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
# BIO3 = Isothermality (BIO2/BIO7) (* 100)
# BIO4 = Temperature Seasonality (standard deviation *100)           ##
# BIO5 = Max Temperature of Warmest Month                            ##
# BIO6 = Min Temperature of Coldest Month                            ##
# BIO7 = Temperature Annual Range (BIO5-BIO6)
# BIO8 = Mean Temperature of Wettest Quarter
# BIO9 = Mean Temperature of Driest Quarter
# BIO10 = Mean Temperature of Warmest Quarter
# BIO11 = Mean Temperature of Coldest Quarter
# BIO12 = Annual Precipitation                                       ##
# BIO13 = Precipitation of Wettest Month                             ##
# BIO14 = Precipitation of Driest Month                              ##
# BIO15 = Precipitation Seasonality (Coefficient of Variation)       ##
# BIO16 = Precipitation of Wettest Quarter
# BIO17 = Precipitation of Driest Quarter
# BIO18 = Precipitation of Warmest Quarter
# BIO19 = Precipitation of Coldest Quarter


## create points
GBIF.POINTS   = GBIF.LAND[c("lon", "lat")]
GBIF.POINTS   = SpatialPointsDataFrame(coords = GBIF.LAND[c("lon", "lat")], data = GBIF.POINTS,
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


## create a loop to extract all the raster values for the land points
n.env.variables   = length(env.grids)
site.env.mat      = GBIF.POINTS

## this is Karel's

# ## For every environmental variable:
# for(i.env in 1:n.env.variables) {
#   
#   ## Get the column names from GBIF, create the raster, 
#   ## extract the raster values, bind the GBIF points and raster values together,
#   ## then finally get the names 
#   col.heads           = names(site.env.mat)
#   Env1.ras            = raster(env.grids[i.env])
#   point.vals          = extract(Env1.ras, GBIF.POINTS) #, method = 'simple') ## why is this not needed?
#   site.env.mat        = cbind(site.env.mat, point.vals)
#   #names(site.env.mat) = c(col.heads, substr(env.grids[i.env], 1, (nchar(env.grids[i.env])-4)))  ## watch this, probably not needed
#   
# } ## end for i.env


## Can do this all at once in R
s <- stack(env.grids)
point.vals <- extract(s, GBIF.POINTS) %>% 
  cbind(coordinates(GBIF.POINTS), .)




# Error in checkNames(value) : 
#   (converted from warning) attempt to set invalid names: this may lead to problems later on. See ?make.names


#########################################################################################################################
## bind values and rename columns
## check this needs to be update for elevation
GBIF.WORLDCLIM    = cbind(GBIF.LAND, site.env.mat)
colnames(GBIF.WORLDCLIM)[1]  = "bio_01"
colnames(GBIF.WORLDCLIM)[4]  = "bio_04"
colnames(GBIF.WORLDCLIM)[5]  = "bio_05"
colnames(GBIF.WORLDCLIM)[6]  = "bio_06"
colnames(GBIF.WORLDCLIM)[7]  = "bio_12"
colnames(GBIF.WORLDCLIM)[8]  = "bio_13"
colnames(GBIF.WORLDCLIM)[8]  = "bio_14"
colnames(GBIF.WORLDCLIM)[8]  = "bio_15"


## convert back to data frame
SITES.ENV          = as.data.frame(site.env.mat)



#########################################################################################################################
## 2). CREATE NICHES FOR SELECTED TAXA
#########################################################################################################################


## check the HIA taxa again
dim(DRAFT.HIA.TAXA)
head(DRAFT.HIA.TAXA)


## Focus on the top 20 taxa again
DRAFT.HIA.TAXA.200 = subset(DRAFT.HIA.TAXA, Top_200 == "TRUE")
dim(DRAFT.HIA.TAXA.200)
head(DRAFT.HIA.TAXA.200)


## Which species name shall we use to run the analysis? The name from the original file!!!
nicheBreadth

