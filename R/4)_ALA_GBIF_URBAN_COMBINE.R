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


## Load GBIF and ALA data
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
## Check these GBIF/ALA data :: No na year or lat/long, all years > 1950
summary(GBIF.ALA.COMBO$year)
summary(GBIF.ALA.COMBO$lat)
summary(GBIF.ALA.COMBO$lon)


#########################################################################################################################
## Now create table for Alessandro.......................................................................................
COMBO.LUT = as.data.frame(table(GBIF.ALA.COMBO$scientificName))
names(COMBO.LUT) = c("scientificName", "FREQUENCY")
COMBO.LUT = COMBO.LUT[with(COMBO.LUT, rev(order(FREQUENCY))), ] 
head(COMBO.LUT);dim(COMBO.LUT)


## Write out the table
write.csv(COMBO.LUT, "./data/base/HIA_LIST/COMBO/SUA_TREES_GBIF_ALA_LUT.csv", row.names = FALSE)





#########################################################################################################################
## 2). CHECK TAXONOMY FOR GBIF, ALA and TREE INVENTORY
#########################################################################################################################

## We Manually create a lookup table of unique species names returned (scientificName) and searchecd (ST_LUT) in GBIF/ALA. 
## We are after the union of the original scientific name, and the lookup table search name ST_LUT. We are taking the
## ST_LUT as the truth for comparisson, and this should inclue valid synonyms.


#########################################################################################################################
## Check taxonomy for the tree inventory data 
head(TI.XY.SPP)
TI.TAXO <- TPL(unique(TI.XY.SPP$searchTaxon), infra = TRUE,
               corr = TRUE, repeats = 100)


## Create a table of unique species names from the TI table
TI.JOIN <- TI.XY.SPP %>%
  left_join(., TI.TAXO, by = c("searchTaxon" = "Taxon"))
TI.LUT = TI.JOIN[, c("searchTaxon",
                     "SOURCE",
                     "Taxonomic.status",
                     "New.Taxonomic.status",
                     "New.Genus",
                     "New.Species")]
Match.TI = TI.LUT[!duplicated(TI.LUT[,c("searchTaxon")]),]
dim(Match.TI)
head(Match.TI)


#########################################################################################################################
## Read the lookup table back in, remap the species...................................................................
SUA.SPP.LUT = read.csv("./data/base/HIA_LIST/COMBO/SUA_SPP_LUT.csv", stringsAsFactors = FALSE)
COMBO.LUT   = join(GBIF.ALA.COMBO[, c("scientificName",
                                      "searchTaxon",
                                      "SOURCE",
                                      "Taxonomic.status",
                                      "New.Taxonomic.status",
                                      "New.Genus",
                                      "New.Species")], SUA.SPP.LUT)

# TI.LUT   = join(TI.LUT[, c("searchTaxon",
#                            "SOURCE",
#                            "Taxonomic.status",
#                            "New.Genus",
#                            "New.Species")], SUA.SPP.LUT)

dim(COMBO.LUT)
dim(TI.LUT)
head(COMBO.LUT)
head(TI.LUT)


## Create a TPL table for just new lookup
# ST.LUT <- TPL(unique(COMBO.LUT$ST_LUT), infra = TRUE,
#               corr = TRUE, repeats = 100)


## See if the Scientific name and searchTaxon match
Match.SN = COMBO.LUT %>%
  mutate(Match.SN.ST = 
           str_detect(scientificName, ST_LUT))    ## str_detect(searchTaxon, ST_LUT))

Match.ST = Match.SN %>%
  mutate(Match.ST.ST = 
           str_detect(ST_LUT, searchTaxon))    ## str_detect(searchTaxon, ST_LUT))



## How many records don't match?
round(with(Match.ST, table(Match.SN.ST)/sum(table(Match.SN.ST))*100), 2)
round(with(Match.ST, table(Match.ST.ST)/sum(table(Match.ST.ST))*100), 2)
with(Match.ST, table(Taxonomic.status))


#########################################################################################################################
## Create a unique of all the cases
#blank.spp = unique$



#########################################################################################################################
## IF the searchTaxon and scientific name match, that is as good as we can do to know they are not misidentified
## Then we need to reassign synonyms to the correct name. 


## So what do the false records mean? If searchTaxon and scientific name do not match, we have searched for one species,
## and GBIF/ALA has either returned a synonym, or the wrong species altogether. To get all possible records, we need to 
## search both GBIF and ALA for all synonyms, which causes mismatches when the data are combined.
## There is no consistent relationship between taxo status and the mismatch.

## 1). We need to check if the SN is the right one
## 2). If it is, update and keep, if not, remove


## Look at the records where the SN and ST don't match
Match.false = subset(Match.ST, Match.SN.ST == "FALSE")
unique(Match.false$Taxonomic.status)
unique(Match.false$SOURCE)
Match.false = Match.false[, c("ST_LUT",
                              "scientificName",
                              "searchTaxon",
                              "SOURCE",
                              "Taxonomic.status",
                              "New.Taxonomic.status",
                              "New.Genus",
                              "New.Species",
                              "FREQUENCY_GBIF",
                              "Match.SN.ST",
                              "Match.SN.ST")]
dim(Match.false)
View(Match.false)


#########################################################################################################################
## If we accept ST_LUT as the taxonomic truth, we can exclude records that are not on this list
# Match.GBIF  = Match.ST[Match.ST$ST_LUT %in% GBIF.spp, ]
# dim(Match.GBIF)
# View(Match.GBIF)


## Take a unique of the updated searchTaxon column (ST_LUT)
length(unique(Match.ST$scientificName))
length(unique(Match.ST$ST_LUT))
length(unique(Match.ST$searchTaxon))
Match.unique.ST = Match.ST[!duplicated(Match.ST[,c("ST_LUT")]),]

length(unique(Match.unique.ST$scientificName))
length(unique(Match.unique.ST$searchTaxon))
length(unique(Match.unique.ST$ST_LUT))


#########################################################################################################################
## We need to create one master table of all the species names used in the analysis
## Use this to re-assign synonyms. Not sure how we want to combine the three sources
length(intersect(TREE.HIA.SPP, Match.unique.ST$searchTaxon))
length(intersect(TREE.HIA.SPP, Match.TI$searchTaxon))
length(intersect(Match.unique.ST$searchTaxon, Match.TI$searchTaxon))


## Not sure what kind of a join we want here
GBIF.ALA.TI.LUT = merge(Match.unique.ST, Match.TI, all = TRUE)
#GBIF.ALA.TI.LUT = merge(Match.unique.ST, Match.TI[,c("searchTaxon", "SOURCE")], by = "searchTaxon")
GBIF.ALA.TI.LUT = GBIF.ALA.TI.LUT[, c("ST_LUT",
                                      "scientificName",
                                      "searchTaxon",
                                      "SOURCE",
                                      "Taxonomic.status",
                                      "New.Taxonomic.status",
                                      "New.Genus",
                                      "New.Species",
                                      "FREQUENCY_GBIF",
                                      "Match.SN.ST",
                                      "Match.SN.ST")]

## Column for SOURCE
#GBIF.ALA.TI.LUT.SPP  = GBIF.ALA.TI.LUT[GBIF.ALA.TI.LUT$ST_LUT %in% TREE.HIA.SPP, ]
#write.csv(GBIF.ALA.TI.LUT, "./data/base/HIA_LIST/COMBO/GBIF_ALA_TI_LUT.csv", row.names = FALSE) 
length(intersect(TREE.HIA.SPP, GBIF.ALA.TI.LUT$ST_LUT))
length(unique(GBIF.ALA.TI.LUT$searchTaxon))
length(unique(GBIF.ALA.TI.LUT$ST_LUT))
table(GBIF.ALA.TI.LUT$SOURCE)


View(GBIF.ALA.TI.LUT)
View(Match.ST)
View(Match.false)


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