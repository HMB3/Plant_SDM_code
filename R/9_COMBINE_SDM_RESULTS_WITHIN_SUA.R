#########################################################################################################################
################################# SUMMARISE SPEICES GAINS BY SUA ######################################################## 
#########################################################################################################################


## This code combines the tables of area occupied by each species in each SUA, and converts it to one big table. 
## This table can be queried/subset to create histograms of species gains and losses inside SUAs, etc.  


#########################################################################################################################
## 1). COMBINE OVERALL TABLES OF SPECIES GAIN/LOSS ACROSS AUSTRALIA
#########################################################################################################################


## Read in the table of species in SUAs and their niches.
SUA.SPP.COUNT       = readRDS(paste0('data/base/HIA_LIST/COMBO/SUA_SPP_COUNT_',        save_run, '.rds'))
COMBO.NICHE.CONTEXT = readRDS(paste0('data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT_',  save_run, '.rds'))
length(unique(SUA.SPP.COUNT$SPECIES))


#########################################################################################################################
## Create arguments for the different settings
## EG ALL_SPP, REC_SPP
## EG ALL_SUA, LARGE_SUA
SUAs       = "ALL_SUAs"   #"LARGE_SUAs"
SUA_SPP    = "ALL_SPP"    #"REC_SPP"
SUA_ORDER  = "CURRENT_MAT"
KOP_ZONE   = "TEMPERATE"
# SUA_ORDER  = "CURRENT_MAP"
# SUA_ORDER  = "CURRENT_MAXT"
# SUA_ORDER  = "CURRENT_AI"
# SUA_ORDER  = "AREASQKM16"
# SUA_ORDER  = "ClimateZ"


#########################################################################################################################
## Create a list of the gain/loss tables. Just run this on the final SUA table 
GAIN.LOSS.list = list.files(maxent_dir, pattern = 'gain_loss_table_', full.names = TRUE, recursive = TRUE) 


## Now combine the SUA tables for each species into one table 
GAIN.LOSS.TABLE <- GAIN.LOSS.list %>%
  
  ## pipe the list into lapply
  lapply(function(x) {
    
    ## create the character string
    f <- paste0(x)
    
    ## load each .RData file
    d <- read.csv(f)
    d
    
  }) %>%
  
  ## finally, bind all the rows together
  bind_rows


## Subset to just the analysis species - some species did not process properly?
GAIN.LOSS.TABLE          =  GAIN.LOSS.TABLE[GAIN.LOSS.TABLE$SPECIES %in% map_spp_list , ] 
GAIN.LOSS.TABLE.COMPLETE =  GAIN.LOSS.TABLE[complete.cases(GAIN.LOSS.TABLE), ]


## Check the exceptions - if all spp models have completed successfully, should be number of species * 3 time periods
length(GAIN.LOSS.TABLE.COMPLETE$SPECIES)/12
summary(GAIN.LOSS.TABLE.COMPLETE)


## How big is the table?
dim(GAIN.LOSS.TABLE)
dim(GAIN.LOSS.TABLE.COMPLETE)
head(GAIN.LOSS.TABLE.COMPLETE, 18)


## How many species are lost?
length(unique(GAIN.LOSS.TABLE$SPECIES))                                   
length(unique(GAIN.LOSS.TABLE.COMPLETE$SPECIES)) 


GAIN.TABLE = subset(GAIN.LOSS.TABLE, CHANGE == "GAINED")
LOSS.TABLE = subset(GAIN.LOSS.TABLE, CHANGE == "LOST")
#View(GAIN.TABLE)
#View(LOSS.TABLE)




#########################################################################################################################
## 2). COMBINE TABLES OF SPECIES PRESENCES IN SUAs
#########################################################################################################################


#########################################################################################################################
## The multiple thresholds could present a problem
SUA.tables = list.files(maxent_dir, pattern = 'SUA_cell_count', full.names = TRUE, recursive = TRUE) 


## Now combine the SUA tables for each species into one table 
SUA.PRESENCE <- SUA.tables %>%
  
  ## pipe the list into lapply
  lapply(function(x) {
    
    ## create the character string
    f <- paste0(x)
    
    ## load each .csv file
    d <- read.csv(f)
    d
    
  }) %>%
  
  ## finally, bind all the rows together
  bind_rows


## This is a summary of maxent output for current conditions
dim(SUA.PRESENCE)
head(SUA.PRESENCE, 10);tail(SUA.PRESENCE, 10)
length(unique(SUA.PRESENCE$SPECIES))
length(unique(SUA.PRESENCE$SUA_NAME16))
length(unique(SUA.PRESENCE$SUA_CODE16))


## Subset to just the analysis species - some species did not process properly?
SUA.PRESENCE  =  SUA.PRESENCE[SUA.PRESENCE$SPECIES %in% map_spp_list , ] 
SUA.PRESENCE  =  SUA.PRESENCE[complete.cases(SUA.PRESENCE), ]
SUA.PRESENCE$SUA_NAME16 = as.character(SUA.PRESENCE$SUA_NAME16)
summary(SUA.PRESENCE)
length(unique(SUA.PRESENCE$SPECIES))


#########################################################################################################################
## Now join on the population. Note that the ABS area shapefile does not have the same SUAs as the ABS table
SUA.DENS = areal_unit@data
SUA.DENS = SUA.DENS[, c("SUA_NAME16", "AREASQKM16")]
names(SUA.DENS) = c("SUA", "AREA_SQKM")
head(SUA.DENS)


## Check the "urban centres" file
length(intersect(ALL.SUA.POP$SUA, SUA.PRESENCE$SUA))
length(intersect(ALL.SUA.POP$SUA, SUA.DENS$SUA))
setdiff(SUA.DENS$SUA, ALL.SUA.POP$SUA)
setdiff(ALL.SUA.POP$SUA, SUA.DENS$SUA)


## Check the "urban centres" file
intersect(SUA.DENS$SUA, URB.POP$SUA)
setdiff(URB.POP$SUA, SUA.DENS$SUA)


## Fill in the missing population fields
ALL.SUA.POP$POP_2017 = as.numeric(ALL.SUA.POP$POP_2017)


## Also, we can calcualte a measure of population density to rank the SUAs. The area's don't look super accurate, but they
## are from the ABS shapefile, so take them to be the standard.
TOP.SUA.POP             = ALL.SUA.POP[, c("SUA", "POP_2017")]
SUA.DEN                 = merge(TOP.SUA.POP, SUA.DENS)
SUA.DEN$POP_DESNITY     = SUA.DEN$POP_2017/SUA.DEN$AREA_SQKM
SUA.DEN                 = SUA.DEN[, c("SUA", "POP_2017", "POP_DESNITY")]


## Try to harmonise the SUAs in both the area and population tables
setdiff(SUA.DEN$SUA, TOP.SUA.POP$SUA)
setdiff(TOP.SUA.POP$SUA, SUA.DENS$SUA)


## Make the SUA field character not factor
setnames(SUA.PRESENCE, old = c("SUA_NAME16", "n_cells"), new = c("SUA", "CELL_COUNT"))
class(SUA.PRESENCE$SUA)
class(SUA.DEN$SUA)


## (94 SUAs overlap between the ABS shapefile and the Population table
str(unique(sort(SUA.DEN$SUA)))
str(unique(sort(SUA.PRESENCE$SUA)))
str(intersect(unique(sort(SUA.DEN$SUA)), unique(sort(SUA.PRESENCE$SUA))))
SUA.PRESENCE = join(SUA.PRESENCE, SUA.DEN)


## Check this combined table
SUA.PRESENCE = SUA.PRESENCE[, c("SUA_CODE16",  "SUA",      "AREASQKM16", "POP_DESNITY",      "POP_2017",
                                "SPECIES",     "PERIOD",   "THRESH",    "CURRENT_SUITABLE",  "FUTURE_SUITABLE",
                                "LOST",        "GAINED",   "STABLE",    "NEVER",             "NODAT",
                                "CELL_COUNT",  "CHANGE",   "GAIN_LOSS")]
SUA.PRESENCE$POP_2017 = as.numeric(SUA.PRESENCE$POP_2017)
URB.POP$POP_2017      = as.numeric(URB.POP$POP_2017)


## Just check there are no NA species
SUA.COMPLETE = completeFun(SUA.PRESENCE, "CURRENT_SUITABLE")
summary(SUA.COMPLETE)


#########################################################################################################################
## Fill in missing population values
SUA.COMPLETE = FillIn(SUA.PRESENCE, URB.POP, "POP_2017", "POP_2017", 
                      KeyVar = c("SUA"), allow.cartesian = FALSE, KeepD2Vars = FALSE)


## Include the maxent rating?
MAXENT.SUMMARY.NICHE = read.csv(paste0('output/maxent/MAXENT_SUMMARY_', save_run, '.csv'), stringsAsFactors = FALSE)
SDM.CHECK = MAXENT.SUMMARY.NICHE[, c("searchTaxon", "Origin", "check.map")]
names(SDM.CHECK) = c("SPECIES", "ORIGIN", "MAXENT_RATING")


## What are the proportions of Origin and complete
table(SDM.CHECK$ORIGIN)
round(with(SDM.CHECK, table(ORIGIN)/sum(table(ORIGIN))*100), 1)
table(SDM.CHECK$MAXENT_RATING)
round(with(SDM.CHECK, table(MAXENT_RATING)/sum(table(MAXENT_RATING))*100), 1)


SUA.COMPLETE$SPECIES = gsub("_", " ", SUA.COMPLETE$SPECIES)
length(intersect(unique(SDM.CHECK$SPECIES), unique(SUA.COMPLETE$SPECIES)))


#########################################################################################################################
## Join on a column for if the species has records in the SUA
SUA.PREDICT = merge(SUA.COMPLETE, SDM.CHECK,  all.x = TRUE)
SUA.PREDICT = join(SUA.PREDICT, SUA.SPP.COUNT, type = "full")


#########################################################################################################################
## If only counting species with records inside the SUA, remove SUA's with < 2 species as their count
if(SUA_SPP == "REC_SPP") {
  
  ## Save basic results and SUA results to file
  message('Analyse only species with records in SUAs') 
  SUA.PREDICT = subset(SUA.PREDICT,  SUA_RECORDS > 0)
  
} else {
  
  message('Analyse all species in SUAs') 
  
}


#SUA.PREDICT = subset(SUA.PREDICT,  SUA_RECORDS > 0)
#SUA.PREDICT = completeFun(SUA.PREDICT, "MAXENT_RATING")
summary(SUA.PREDICT$SUA_RECORDS)
unique(SUA.PREDICT$MAXENT_RATING)
length(unique(SUA.PREDICT$SPECIES))
length(unique(SUA.PREDICT$SUA))
identical(dim(SUA.PREDICT)[1], dim(SUA.COMPLETE)[1])
dim(SUA.PREDICT)





#########################################################################################################################
## 3). ADD CLIMATE DATA TO SUA TABLES
#########################################################################################################################


#########################################################################################################################
## Current climate
SUA.BIO1.current.stats  = read.csv("./data/base/worldclim/aus/1km/bio/SUA_BIO1_current_stats.csv",  stringsAsFactors = FALSE)
SUA.BIO5.current.stats  = read.csv("./data/base/worldclim/aus/1km/bio/SUA_BIO5_current_stats.csv",  stringsAsFactors = FALSE)
SUA.BIO12.current.stats = read.csv("./data/base/worldclim/aus/1km/bio/SUA_BIO12_current_stats.csv", stringsAsFactors = FALSE)
SUA.KOP                 = read.csv("./data/base/worldclim/aus/1km/bio/SUA_KOPPEN.csv",              stringsAsFactors = FALSE)
KOP.LUT                 = read.csv("./data/base/worldclim/aus/1km/bio/KOPPEN_LUT.csv",              stringsAsFactors = FALSE)
SUA.PET                 = read.csv("./data/base/worldclim/aus/1km/bio/SUA_PET_current_stats.csv",   stringsAsFactors = FALSE)
SUA.AI                  = read.csv("./data/base/worldclim/aus/1km/bio/SUA_AI_current_stats.csv",    stringsAsFactors = FALSE)


## Future climate
SUA.BIO1.2030.stats     = read.csv("./data/base/worldclim/aus/1km/bio/SUA_BIO1_2030_stats.csv",     stringsAsFactors = FALSE)
SUA.BIO1.2050.stats     = read.csv("./data/base/worldclim/aus/1km/bio/SUA_BIO1_2050_stats.csv",     stringsAsFactors = FALSE)
SUA.BIO1.2070.stats     = read.csv("./data/base/worldclim/aus/1km/bio/SUA_BIO1_2050_stats.csv",     stringsAsFactors = FALSE)


#########################################################################################################################
## Rename zonal stats
names(SUA.BIO1.current.stats)[names(SUA.BIO1.current.stats)   == "SUA_NAME16"] <- "SUA"
names(SUA.BIO1.current.stats)[names(SUA.BIO1.current.stats)   == "MEAN"] <- "CURRENT_MAT"

names(SUA.BIO5.current.stats)[names(SUA.BIO5.current.stats)   == "SUA_NAME16"] <- "SUA"
names(SUA.BIO5.current.stats)[names(SUA.BIO5.current.stats)   == "MEAN"] <- "CURRENT_MAXT"

names(SUA.BIO12.current.stats)[names(SUA.BIO12.current.stats) == "SUA_NAME16"] <- "SUA"
names(SUA.BIO12.current.stats)[names(SUA.BIO12.current.stats) == "MEAN"] <- "CURRENT_MAP"

names(SUA.PET)[names(SUA.PET) == "SUA_NAME16"] <- "SUA"
names(SUA.PET)[names(SUA.PET) == "MEAN"] <- "CURRENT_PET"

names(SUA.AI)[names(SUA.AI) == "SUA_NAME16"] <- "SUA"
names(SUA.AI)[names(SUA.AI) == "MEAN"] <- "CURRENT_AI"

names(SUA.KOP)[names(SUA.KOP) == "SUA_NAME16"] <- "SUA"
names(SUA.KOP)[names(SUA.KOP) == "MAJORITY"] <- "MAJOR_KOP"


## Divide temperature by 10
SUA.BIO1.2030.stats[["PERIOD"]]  <- 30
colnames(SUA.BIO1.2030.stats)[1] <- "SUA"
SUA.BIO1.current.stats[["CURRENT_MAT"]] = SUA.BIO1.current.stats[["CURRENT_MAT"]]/10


## Divide PET by 10,000
SUA.AI[["CURRENT_AI"]] = SUA.AI[["CURRENT_AI"]]/10000


#########################################################################################################################
## Create Koppen climate zone for centroid of each SUA
writeSpatialShape(SUA.16, "SUA.16")
cents <- coordinates(SUA.16)
cents <- SpatialPointsDataFrame(coords = cents, data = SUA.16@data,
                                proj4string = CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"))

plot(SUA.16)
points(cents, col = "Blue")
writeSpatialShape(cents, "cents")

SUA_centroids <- coordinates(SUA.16)
SUA_location  = data.frame(SUA_centroids)
SUA_location  = cbind(SUA.16$SUA_NAME16, SUA_location)
names(SUA_location)   = c("SUA","rndCoord.lon", "rndCoord.lat")
SUA_location$rndCoord.lon = RoundCoordinates(SUA_location$rndCoord.lon)
SUA_location$rndCoord.lat = RoundCoordinates(SUA_location$rndCoord.lat)

points(centroids, pch = 3, col = "Red")
Kop.loc <- data.frame(SUA_location, ClimateZ = LookupCZ(SUA_location))
Kop.loc = Kop.loc[c("SUA", "ClimateZ")]


#########################################################################################################################
## Merge average climate of SUAs onto the SUA maxent results table
SUA.PREDICT = join(SUA.PREDICT, SUA.BIO1.current.stats[c("SUA", "CURRENT_MAT")])
SUA.PREDICT = join(SUA.PREDICT, SUA.BIO12.current.stats[c("SUA", "CURRENT_MAP")])
SUA.PREDICT = join(SUA.PREDICT, SUA.BIO5.current.stats[c("SUA", "CURRENT_MAXT")])
SUA.PREDICT = join(SUA.PREDICT, SUA.PET[c("SUA", "CURRENT_PET")])
SUA.PREDICT = join(SUA.PREDICT, SUA.AI[c("SUA", "CURRENT_AI")])
SUA.PREDICT = join(SUA.PREDICT, Kop.loc[c("SUA", "ClimateZ")])
summary(SUA.PREDICT)


#########################################################################################################################
## Save tables
if(save_data == "TRUE") {
  
  ## Save basic results and SUA results to file
  message(length(unique(SUA.COMPLETE$SPECIES)), " Species analysed in ", length(unique(SUA.COMPLETE$SUA)), " SUAs")
  
  ## Load GBIF and ALA data
  write.csv(GAIN.LOSS.TABLE,  paste0('output/tables/OVERALL_GAIN_LOSS_TABLE_',   save_run, '.csv'), row.names = FALSE)
  write.csv(SUA.PREDICT,      paste0('output/tables/MAXENT_SUA_PRESENCE_KOPPEN', save_run, '.csv'), row.names = FALSE)
  
} else {
  
  message('Skip file saving, not many species analysed')   ##
  
}


#########################################################################################################################
## We can count of all the species that are being lost, gained or remaining stable in each SUA
## However, this doesn't give us the turner of species, because it ignores the species identities
length(unique(SUA.PREDICT$SPECIES))
SUA.PLOT.GOOD = subset(SUA.PREDICT, MAXENT_RATING < 3) 
unique(SUA.PLOT.GOOD$MAXENT_RATING)
length(unique(SUA.PLOT.GOOD$SPECIES))





#########################################################################################################################
## 4). SUBSET DATA TO MAKE PLOTTING EASIER
#########################################################################################################################


#########################################################################################################################
## Create plots in the format needed for ggplot bar  
SUA.PLOT.GOOD.30 = subset(SUA.PLOT.GOOD, PERIOD == 30)
SUA.PLOT.GOOD.50 = subset(SUA.PLOT.GOOD, PERIOD == 50)
SUA.PLOT.GOOD.70 = subset(SUA.PLOT.GOOD, PERIOD == 70)
unique(SUA.PLOT.GOOD.30$PERIOD);unique(SUA.PLOT.GOOD.50$PERIOD);unique(SUA.PLOT.GOOD.70$PERIOD)
dim(SUA.PLOT.GOOD.30);dim(SUA.PLOT.GOOD.50);dim(SUA.PLOT.GOOD.70)


## Melt the table into the right format : but does this mean the different categories are mutually exclusive?
## Species can only fall in the categories gain, loss, stable, never, in each SUA. So yes, they are exclusive.
SUA.PLOT.30          = table(SUA.PLOT.GOOD.30$SUA, SUA.PLOT.GOOD.30$GAIN_LOSS)
SUA.PLOT.50          = table(SUA.PLOT.GOOD.50$SUA, SUA.PLOT.GOOD.50$GAIN_LOSS)
SUA.PLOT.70          = table(SUA.PLOT.GOOD.70$SUA, SUA.PLOT.GOOD.70$GAIN_LOSS)
SUA.PLOT.30.M        = melt(SUA.PLOT.30)
SUA.PLOT.50.M        = melt(SUA.PLOT.50)
SUA.PLOT.70.M        = melt(SUA.PLOT.70)

names(SUA.PLOT.30.M) = c("SUA", "AREA_CHANGE", "SPECIES_COUNT")
names(SUA.PLOT.50.M) = c("SUA", "AREA_CHANGE", "SPECIES_COUNT")
names(SUA.PLOT.70.M) = c("SUA", "AREA_CHANGE", "SPECIES_COUNT")

SUA.PLOT.30.M        = subset(SUA.PLOT.30.M, AREA_CHANGE != "NEVER")# & AREA_CHANGE != "STABLE")
SUA.PLOT.50.M        = subset(SUA.PLOT.50.M, AREA_CHANGE != "NEVER")# & AREA_CHANGE != "STABLE")
SUA.PLOT.70.M        = subset(SUA.PLOT.70.M, AREA_CHANGE != "NEVER")# & AREA_CHANGE != "STABLE") 

SUA.30.M.LOSS        = subset(SUA.PLOT.30.M, AREA_CHANGE %in% c("LOSS"))
SUA.30.M.GAIN        = subset(SUA.PLOT.30.M, AREA_CHANGE %in% c("GAIN"))
SUA.30.M.STABLE      = subset(SUA.PLOT.30.M, AREA_CHANGE %in% c("STABLE"))

SUA.50.M.LOSS        = subset(SUA.PLOT.50.M, AREA_CHANGE %in% c("LOSS"))
SUA.50.M.GAIN        = subset(SUA.PLOT.50.M, AREA_CHANGE %in% c("GAIN"))
SUA.50.M.STABLE      = subset(SUA.PLOT.50.M, AREA_CHANGE %in% c("STABLE"))

SUA.70.M.LOSS        = subset(SUA.PLOT.70.M, AREA_CHANGE %in% c("LOSS"))
SUA.70.M.GAIN        = subset(SUA.PLOT.70.M, AREA_CHANGE %in% c("GAIN"))
SUA.70.M.STABLE      = subset(SUA.PLOT.70.M, AREA_CHANGE %in% c("STABLE"))

head(SUA.30.M.LOSS)
head(SUA.30.M.GAIN)
head(SUA.30.M.STABLE)


#########################################################################################################################
## Attach the climate
SUA.CLIM      = SUA.PREDICT[!duplicated(SUA.PREDICT[,c('SUA')]),][c("SUA", "CURRENT_MAT", "CURRENT_MAP", 
                                                                    "CURRENT_PET", "CURRENT_AI", "CURRENT_MAXT", 
                                                                    "AREASQKM16", "ClimateZ", "POP_2017")]


## Find a more efficient way to join everything on to the subsets
SUA.PLOT.30.M = join(SUA.PLOT.30.M, SUA.CLIM)
SUA.PLOT.50.M = join(SUA.PLOT.50.M, SUA.CLIM)
SUA.PLOT.70.M = join(SUA.PLOT.70.M, SUA.CLIM)

SUA.30.M.LOSS   = join(SUA.30.M.LOSS,   SUA.CLIM)
SUA.30.M.GAIN   = join(SUA.30.M.GAIN,   SUA.CLIM)
SUA.30.M.STABLE = join(SUA.30.M.STABLE, SUA.CLIM)

SUA.70.M.LOSS   = join(SUA.70.M.LOSS,   SUA.CLIM)
SUA.70.M.GAIN   = join(SUA.70.M.GAIN,   SUA.CLIM)
SUA.70.M.STABLE = join(SUA.70.M.STABLE, SUA.CLIM)



#########################################################################################################################
## Use ggplot to plot the percentage of species inside an LGA, which is being lost or gained
## Calculate the gain/loss for counts of all species?

## Option 1).
## Gain % = species gained / (stable + lost)   * 100
## Lost % = species lost   / (stable + lost)   * 100

## Calculate the gain/loss for for only species that are recorded within each SUA?

## Option 2).
## Gain % = species gained / (stable + lost + gained)   * 100
## Lost % = species lost   / (stable + lost + gained)   * 100



#########################################################################################################################
## Run different variations
#"./R/SUA_GAIN_LOSS_ALL_SUAs.R"       ## 1).
#"./R/SUA_%100_PLOT_LARGER_SUAs.R"    ## 2).
#"./R/SUA_RAIN_PLOT.R"                ## 2), using rain not temp





#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################