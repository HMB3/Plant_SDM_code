#########################################################################################################################
############################################# ESTIMATE SPECIES NICHES ################################################### 
#########################################################################################################################


#########################################################################################################################
## This code combines the GBIF records for all species with the ALA data into a single table, extracts environmental 
## values and final adds contextual info for each record (taxonomic and horticultural) 


## It creates two tables:

## 1). A large table with one row for each species record
## 2). A smaller table with one row for each species, including contextual data and species attributes (niches, traits, etc.)
 
## These tables are subsequently used to estimate the current global realised niche/climatic tolerance using the best 
## available data, and susequently model the niches using the maxent algorithm.  


#########################################################################################################################
## Read in all data to run the SDM code :: species lists, shapefile, rasters & tables
#source('./R/HIA_LIST_MATCHING.R')
rasterTmpFile()


## GBIF and ALA data
#GBIF.LAND = readRDS("./data/base/HIA_LIST/GBIF/GBIF_TREES_LAND.rds")
#ALA.TREES.LAND = readRDS("./data/base/HIA_LIST/GBIF/ALA_TREES_LAND.rds")
#load("./data/base/HIA_LIST/ALA/ALA_LAND_POINTS.RData")                      ## save the update in the same place



#########################################################################################################################
## 1). MERGE GBIF DATA WITH CLEANED AVH/ALA RECORDS
#########################################################################################################################


#########################################################################################################################
## Merge on the ALA data. Consider that GBIF has data for both sources. We are topping up the native ranges with the AVH. 
## So there could be duplicates between both sources.
length(unique(GBIF.LAND$searchTaxon)) 
length(unique(ALA.TREES.LAND$searchTaxon)) 
unique(GBIF.LAND$SOURCE)
unique(ALA.TREES.LAND$SOURCE)


dim(GBIF.LAND);dim(ALA.TREES.LAND)
names(GBIF.LAND)
names(ALA.TREES.LAND)
setdiff(names(GBIF.LAND), names(ALA.TREES.LAND))


#########################################################################################################################
## Rename a few fields and and a field for source........................................................................
GBIF.LAND     = dplyr::rename(GBIF.LAND, 
                              coordinateUncertaintyInMetres = coordinateUncertaintyInMeters,
                              rank = taxonRank)


#########################################################################################################################
## Bind the rows together
intersect(names(GBIF.LAND), names(ALA.TREES.LAND))
GBIF.ALA.COMBO = bind_rows(GBIF.LAND, ALA.TREES.LAND)
names(GBIF.ALA.COMBO)
unique(GBIF.ALA.COMBO$SOURCE)
identical((dim(GBIF.LAND)[1]+dim(ALA.TREES.LAND)[1]),dim(GBIF.ALA.COMBO)[1])      ## Only adding the ovelap...
head(GBIF.ALA.COMBO)
length(unique(GBIF.ALA.COMBO$searchTaxon))
length(unique(GBIF.ALA.COMBO$scientificName)) 


#########################################################################################################################
## Check these GBIF/ALA data
summary(GBIF.ALA.COMBO$year)
summary(GBIF.ALA.COMBO$lat)
summary(GBIF.ALA.COMBO$lon)


## Check the dimensions for some example species
dim(subset(GBIF.ALA.COMBO, searchTaxon == "Platanus acerifolia"))
dim(subset(GBIF.ALA.COMBO, searchTaxon == "Eucalyptus cephalocarpa"))



#########################################################################################################################
## Now create table for Alessandro.......................................................................................
COMBO.LUT = as.data.frame(table(GBIF.ALA.COMBO$scientificName))
names(COMBO.LUT) = c("scientificName", "FREQUENCY")
COMBO.LUT = COMBO.LUT[with(COMBO.LUT, rev(order(FREQUENCY))), ] 
head(COMBO.LUT);dim(COMBO.LUT)


## Write out the table
write.csv(COMBO.LUT, "./data/base/HIA_LIST/COMBO/SUA_TREES_GBIF_ALA_LUT.csv", row.names = FALSE)





#########################################################################################################################
## 2). CHECK TAXONOMY
#########################################################################################################################


#########################################################################################################################
## Read the lookup table back in, remap the species...................................................................
## How is the lookup table different to what we already have? Alessandro's species names were not as clean as mine, so searchTaxon 
## might be the same as before.


## How to check the taxonomy?
unique(GBIF.ALA.COMBO$Taxonomic.status)
unique(GBIF.ALA.COMBO$New.Taxonomic.status)


## Is this the correct join?
# GBIF.TRIM.TAXO <- GBIF.TRIM %>%
#   left_join(., GBIF.TAXO, by = c("scientificName" = "Taxon"))

## Create a table
round(with(GBIF.ALA.COMBO, table(Taxonomic.status)/sum(table(Taxonomic.status))*100), 2)


#########################################################################################################################
## Accepted species
Accepted = subset(GBIF.ALA.COMBO, Taxonomic.status == "Accepted")
Accepted = head(Accepted, dim(Accepted)[1])[, c("scientificName",
                                                "searchTaxon", 
                                                "Taxonomic.status")]


## Check the accepted taxa really do match? This doesn't mean some species were not misidentified... 
Match.names = Accepted %>%
  mutate(Match = 
           str_detect(scientificName, searchTaxon))

dim(Match.names)
dim(subset(Match.names, Match == "TRUE"))
mismatch = subset(Match.names, Match == "FALSE")
length(unique(mismatch$scientificName))


## Synonyms
Synonyms = subset(GBIF.ALA.COMBO, Taxonomic.status == "Synonym")
head(Synonyms)[, c("scientificName",
                   "searchTaxon", 
                   "Taxonomic.status")]
length(unique(Synonyms$scientificName))


## Unresolved
Unresolved = subset(GBIF.ALA.COMBO, Taxonomic.status == "Unresolved")
head(Unresolved)[, c("scientificName",
                   "searchTaxon", 
                   "Taxonomic.status")]
length(unique(Unresolved$scientificName))


## Unresolved
Misapplied = subset(GBIF.ALA.COMBO, Taxonomic.status == "Misapplied")
head(Misapplied)[, c("scientificName",
                     "searchTaxon", 
                     "Taxonomic.status"), 10]

## No status
No.status = subset(GBIF.ALA.COMBO, Taxonomic.status == "")
head(No.status, 50)[, c("scientificName",
                    "searchTaxon", 
                    "Taxonomic.status")]
length(unique(No.status$scientificName))



#########################################################################################################################
## Create points: the over function seems to need geographic coordinates for this data...
COMBO.POINTS   = SpatialPointsDataFrame(coords      = GBIF.ALA.COMBO[c("lon", "lat")], 
                                        data        = GBIF.ALA.COMBO[c("lon", "lat")],
                                        proj4string = CRS.WGS.84)


## Check
dim(COMBO.POINTS)
projection(COMBO.POINTS)
names(COMBO.POINTS)





#########################################################################################################################
## 3). PROJECT RASTERS AND EXTRACT ALL WORLDCLIM DATA FOR SPECIES RECORDS
#########################################################################################################################


#########################################################################################################################
## Ignore edaphic variables


# BIO1  = Annual Mean Temperature                                     ## 
# BIO2  = Mean Diurnal Range (Mean of monthly (max temp - min temp))  ##
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


#########################################################################################################################
## Create a stack of rasters to sample: get all the Worldclim variables just for good measure
## Use the Mollweide projection for the points and rasters 
env.grids.current = stack(
  file.path('./data/base/worldclim/world/0.5/bio/current',
            sprintf('bio_%02d', 1:19))) 

## Also get the PET raster
PET               = raster("./data/base/worldclim/world/1km/pet_he_yr1.tif")


#########################################################################################################################
## Is there a way to speed this up?
projection(COMBO.POINTS);projection(env.grids.current)
COMBO.RASTER <- extract(env.grids.current, COMBO.POINTS) %>% 
  cbind(GBIF.ALA.COMBO, .)


## Multiple rename using dplyr
COMBO.RASTER = dplyr::rename(COMBO.RASTER,
                             
                             ## Temperature
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
                             
                             ## Rainfall
                             Annual_precip        = bio_12,
                             Precip_wet_month     = bio_13,
                             Precip_dry_month     = bio_14,
                             Precip_seasonality   = bio_15,
                             Precip_wet_qu        = bio_16,
                             Precip_dry_qu        = bio_17,
                             Precip_warm_qu       = bio_18,
                             Precip_col_qu        = bio_19)


## Save/load
#summary(COMBO.RASTER)
#saveRDS(COMBO.RASTER, file = paste("./data/base/HIA_LIST/GBIF/COMBO_GBIF_ALA_RASTER.rds"))
gc();gc()


#########################################################################################################################
## Extract the raster data for PET
COMBO.POINTS   = SpatialPointsDataFrame(coords      = COMBO.RASTER[c("lon", "lat")], 
                                        data        = COMBO.RASTER[c("lon", "lat")],
                                        proj4string = CRS.WGS.84)


projection(COMBO.POINTS);projection(PET)
POINTS.PET <- extract(PET, COMBO.POINTS) %>% 
  cbind(COMBO.RASTER, .)
#saveRDS(POINTS.PET, file = paste("./data/base/HIA_LIST/GBIF/COMBO_GBIF_ALA_RASTER_PET.rds"))
COMBO.RASTER = POINTS.PET
names(COMBO.RASTER)[names(COMBO.RASTER) == "."] <- 'PET'


## Check 
dim(COMBO.RASTER)
names(COMBO.RASTER)
projection(COMBO.RASTER)


## Check the raster values here..........................................................................................
summary(Annual_mean_temp)





#########################################################################################################################
## 4). CONVERT RASTER VALUES
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
                  "Precip_col_qu",
                  "PET")


#########################################################################################################################
## Change the raster values here: See http://worldclim.org/formats1 for description of the interger conversion. 
## All temperature variables wer multiplied by 10, so divide by 10 to reverse it.
COMBO.RASTER.CONVERT = as.data.table(COMBO.RASTER)                           ## Check this works, also inefficient
COMBO.RASTER.CONVERT[, (env.variables [c(1:11)]) := lapply(.SD, function(x) 
  x / 10 ), .SDcols = env.variables [c(1:11)]]
COMBO.RASTER.CONVERT = as.data.frame(COMBO.RASTER.CONVERT)                   ## Find another method without using data.table


## Check Looks ok?
summary(COMBO.RASTER.CONVERT$Annual_mean_temp)  ## -23 looks too low/. Check where these are ok
summary(COMBO.RASTER$Annual_mean_temp)

summary(COMBO.RASTER.CONVERT$Isothermality)
summary(COMBO.RASTER$Isothermality)


## Print the dataframe dimensions to screen
dim(COMBO.RASTER.CONVERT)
identical(length(unique(COMBO.RASTER.CONVERT$searchTaxon)), length(GBIF.spp))


## Plot a few points to see :: do those look reasonable?
# plot(LAND, col = 'grey', bg = 'sky blue')
# points(COMBO.RASTER.CONVERT[ which(COMBO.RASTER.CONVERT$Annual_mean_temp < -5), ][, c("lon", "lat")], 
#        pch = ".", col = "red", cex = 3, asp = 1, main = "temp records < -5")





#########################################################################################################################
## Combine the existing data with the new species :: use this approach to just process a subset 
## COMBO.RASTER.CONTEXT.UPDATE = COMBO.RASTER.CONTEXT
## saveRDS(COMBO.RASTER.CONTEXT.UPDATE, file = paste("./data/base/HIA_LIST/COMBO/COMBO_RASTER_CONTEXT_UPDATE.rds", sep = ""))
## load("./data/base/HIA_LIST/COMBO/COMBO_RASTER_CONTEXT.rds")
## load("./data/base/HIA_LIST/COMBO/COMBO_RASTER_CONTEXT_UPDATE.rds")
# names(COMBO.RASTER.CONTEXT);names(COMBO.RASTER.CONTEXT.UPDATE)
# COMBO.RASTER.CONTEXT.UPDTATE = bind_rows(COMBO.RASTER.CONTEXT, COMBO.RASTER.CONTEXT.UPDATE)
       

## Save the summary datasets
saveRDS(COMBO.RASTER.CONVERT, file = paste("./data/base/HIA_LIST/COMBO/COMBO_RASTER_CONVERT_AUGUST_2018.rds", sep = ""))



## Now save .RData file for the next session...
#save.image("STEP_4_NICHES.RData")





#########################################################################################################################
## OUTSTANDING NICHE TASKS:
#########################################################################################################################


## Clean species data with Alessandro's method 




#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################