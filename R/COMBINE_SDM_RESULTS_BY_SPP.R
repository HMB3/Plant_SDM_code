#########################################################################################################################
################################# SUMMARISE SPEICES GAINS BY SPP ######################################################## 
#########################################################################################################################


## This code tries to summarise the species/losses and gains by species niches, rather than SUA


#########################################################################################################################
## 1). COMBINE OVERALL TABLES OF SPECIES GAIN/LOSS ACROSS AUSTRALIA
#########################################################################################################################


## Read in the table of species in SUAs and their niches
SUA.SPP.COUNT       = readRDS(paste0('data/base/HIA_LIST/COMBO/SUA_SPP_COUNT_',        save_run, '.rds'))
SUA.NICHE           = readRDS(paste0('data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT_',  save_run, '.rds'))
names(SUA.NICHE)[names(SUA.NICHE) == "searchTaxon"] <- "SPECIES"
length(unique(SUA.SPP.COUNT$SPECIES))
length(unique(SUA.NICHE$SPECIES))


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
                                "SPECIES",     "PERIOD",   "THRESH",    "CURRENT_SUITABLE", "FUTURE_SUITABLE",
                                "LOST",        "GAINED",   "STABLE",    "NEVER",            "NODAT",
                                "CELL_COUNT",  "CHANGE",   "GAIN_LOSS")]


## Just check there are no NA species
SUA.COMPLETE = completeFun(SUA.PRESENCE, "CURRENT_SUITABLE")
summary(SUA.COMPLETE)


#########################################################################################################################
## Find the species with the greatest increase/decrease?
summary(SUA.COMPLETE$CHANGE)
unique(SUA.COMPLETE$GAIN_LOSS)


#########################################################################################################################
## Select the major SUAs
top_n   = 8
BIG.SUA = head(TOP.SUA.POP[with(TOP.SUA.POP, rev(order(POP_2017))), ], top_n)
BIG_SUA = BIG.SUA$SUA
MAP_SUA = c("Adelaide", "Brisbane", "Darwin", "Perth", "Sydney", "Hobart", "Melbourne", "Canberra")


## Restrict the big table to just the largest SUAs
SUA.TOP.PRESENCE  = SUA.PRESENCE[SUA.PRESENCE$SUA %in% BIG_SUA, ] 
SUA.TOP.PRESENCE  = SUA.TOP.PRESENCE [with(SUA.TOP.PRESENCE , rev(order(POP_2017))), ]
summary(SUA.TOP.PRESENCE)


#########################################################################################################################
## Save basic results and SUA results to file
message(length(unique(SUA.COMPLETE$SPECIES)), " Species analysed in ", length(unique(SUA.COMPLETE$SUA)), " SUAs")


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


## Create a plot of number of records vs. maxent rating, Boxplot with no. occurrences on y, and maxent rating on x
## No. of global occurrences, No. of Aus records, No. of plantings
## ......................................................................................................................


#########################################################################################################################
## Join on a column for if the species has records in the SUA
SUA.PREDICT = merge(SUA.COMPLETE, SDM.CHECK, all.x = TRUE)
SUA.PREDICT = join(SUA.PREDICT, SUA.SPP.COUNT, type = "full")
SUA.PREDICT = subset(SUA.PREDICT,  SUA_RECORDS > 0)
summary(SUA.PREDICT$SUA_RECORDS)
#SUA.PREDICT = completeFun(SUA.PREDICT, "MAXENT_RATING")
unique(SUA.PREDICT$MAXENT_RATING)
length(unique(SUA.PREDICT$SPECIES))
identical(dim(SUA.PREDICT)[1], dim(SUA.COMPLETE)[1])
dim(SUA.PREDICT)





#########################################################################################################################
## 3). CREATE PLOTS OF SPECIES LOSS AND GAIN
#########################################################################################################################


#########################################################################################################################
## Future rasters
SUA.BIO1.current.stats  = read.csv("./data/base/worldclim/aus/1km/bio/SUA_BIO1_current_stats.csv",  stringsAsFactors = FALSE)
SUA.BIO12.current.stats = read.csv("./data/base/worldclim/aus/1km/bio/SUA_BIO12_current_stats.csv", stringsAsFactors = FALSE)
SUA.KOP                 = read.csv("./data/base/worldclim/aus/1km/bio/SUA_KOPPEN.csv",              stringsAsFactors = FALSE)
KOP.LUT                 = read.csv("./data/base/worldclim/aus/1km/bio/KOPPEN_LUT.csv",              stringsAsFactors = FALSE)
SUA.PET                 = read.csv("./data/base/worldclim/aus/1km/bio/SUA_PET_current_stats.csv",   stringsAsFactors = FALSE)

SUA.BIO1.2030.stats     = read.csv("./data/base/worldclim/aus/1km/bio/SUA_BIO1_2030_stats.csv",     stringsAsFactors = FALSE)
SUA.BIO1.2050.stats     = read.csv("./data/base/worldclim/aus/1km/bio/SUA_BIO1_2050_stats.csv",     stringsAsFactors = FALSE)
SUA.BIO1.2070.stats     = read.csv("./data/base/worldclim/aus/1km/bio/SUA_BIO1_2050_stats.csv",     stringsAsFactors = FALSE)


#########################################################################################################################
## Rename zonal stats
names(SUA.BIO1.current.stats)[names(SUA.BIO1.current.stats) == "SUA_NAME16"] <- "SUA"
names(SUA.BIO1.current.stats)[names(SUA.BIO1.current.stats) == "MEAN"] <- "CURRENT_MAT"

names(SUA.BIO12.current.stats)[names(SUA.BIO12.current.stats) == "SUA_NAME16"] <- "SUA"
names(SUA.BIO12.current.stats)[names(SUA.BIO12.current.stats) == "MEAN"] <- "CURRENT_MAP"

names(SUA.PET)[names(SUA.PET) == "SUA_NAME16"] <- "SUA"
names(SUA.PET)[names(SUA.PET) == "MEAN"] <- "CURRENT_PET"

names(SUA.KOP)[names(SUA.KOP) == "SUA_NAME16"] <- "SUA"
names(SUA.KOP)[names(SUA.KOP) == "MAJORITY"] <- "MAJOR_KOP"


## Divide temperature by 10
SUA.BIO1.2030.stats[["PERIOD"]]  <- 30
colnames(SUA.BIO1.2030.stats)[1] <- "SUA"
SUA.BIO1.current.stats[["CURRENT_MAT"]] = SUA.BIO1.current.stats[["CURRENT_MAT"]]/10


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


## Merge MAP onto the results table
SUA.PREDICT = join(SUA.PREDICT, SUA.BIO1.current.stats[c("SUA", "CURRENT_MAT")])
SUA.PREDICT = join(SUA.PREDICT, SUA.BIO12.current.stats[c("SUA", "CURRENT_MAP")])
SUA.PREDICT = join(SUA.PREDICT, Kop.loc[c("SUA", "ClimateZ")])
SUA.PREDICT = join(SUA.PREDICT, SUA.PET[c("SUA", "CURRENT_PET")])
SUA.PREDICT = join(SUA.PREDICT, SUA.NICHE[c("SPECIES", "Annual_mean_temp_median", "Max_temp_warm_month_median")])
summary(SUA.PREDICT)


#########################################################################################################################
## Save tables
if(save_data == "TRUE") {
  
  ## Load GBIF and ALA data
  write.csv(GAIN.LOSS.TABLE,  paste0('output/tables/OVERALL_GAIN_LOSS_TABLE_',   save_run, '.csv'), row.names = FALSE)
  write.csv(SUA.PREDICT,      paste0('output/tables/MAXENT_SUA_PRESENCE_KOPPEN', save_run, '.csv'), row.names = FALSE)
  write.csv(SUA.TOP.PRESENCE, paste0('output/tables/MAXNET_SUA_TOP_',            save_run, '.csv'), row.names = FALSE)
  
} else {
  
  message('Skip file saving, not many species analysed')   ##
  
}


#########################################################################################################################
## We can count of all the species that are being lost, gained or remaining stable in each SUA
## However, this doesn't give us the turner of species, because it ignores the species identities
length(unique(SUA.PREDICT$SPECIES))
SUA.PLOT.GOOD = subset(SUA.PREDICT, MAXENT_RATING < 3) #& ORIGIN == "Native")
unique(SUA.PLOT.GOOD$MAXENT_RATING)
length(unique(SUA.PLOT.GOOD$SPECIES))


#########################################################################################################################
## Now calculate the change in area for each species, inside each SUA
## Round data to 0 decimal places
SUA.PLOT.GOOD$FUTURE_PROP  <- round(((SUA.PLOT.GOOD$FUTURE_SUITABLE/SUA.PLOT.GOOD$AREASQKM16)*100),0)   # Proportion of SUA that is suitable
SUA.PLOT.GOOD$STABLE_PROP  <- round(((SUA.PLOT.GOOD$STABLE/SUA.PLOT.GOOD$CURRENT_SUITABLE)*100),0)      # Proportion of current SH that remains suitable
SUA.PLOT.GOOD$GAIN_PROP    <- round(((SUA.PLOT.GOOD$GAINED/SUA.PLOT.GOOD$FUTURE_SUITABLE)*100),0)       # Proportion of future SH that is in new cells
SUA.PLOT.GOOD$LOST_PROP    <- round(((SUA.PLOT.GOOD$LOST/SUA.PLOT.GOOD$CURRENT_SUITABLE)*100),0)        # Proportion of current SH that is lost
SUA.PLOT.GOOD$NEVER_PROP   <- round(((SUA.PLOT.GOOD$NEVER/SUA.PLOT.GOOD$AREASQKM16)*100),0)             # Proportion of SUA that is never suitable


## Need a few if/else conditions to catch the different values : If current and future are positive, calculate as a proportion
SUA.PLOT.GOOD$RANGE_CHANGE <- ifelse(SUA.PLOT.GOOD$CURRENT_SUITABLE > 0 & SUA.PLOT.GOOD$FUTURE_SUITABLE > 0, 
                                     round((((SUA.PLOT.GOOD$FUTURE_SUITABLE - 
                                                SUA.PLOT.GOOD$CURRENT_SUITABLE)/SUA.PLOT.GOOD$CURRENT_SUITABLE)*100),0),
                                     
                                     ## If current is positive and future is zero, calculate as a proportion
                                     ifelse(SUA.PLOT.GOOD$CURRENT_SUITABLE > 0 & SUA.PLOT.GOOD$FUTURE_SUITABLE == 0, 
                                            round((((SUA.PLOT.GOOD$FUTURE_SUITABLE - 
                                                       SUA.PLOT.GOOD$CURRENT_SUITABLE)/SUA.PLOT.GOOD$CURRENT_SUITABLE)*100),0),
                                            
                                            ## If current is zero and future is postive, assign 500 to cap the increase, else assign to 0
                                            ifelse(SUA.PLOT.GOOD$CURRENT_SUITABLE == 0 & SUA.PLOT.GOOD$FUTURE_SUITABLE > 0, 100, 0)))

## Should we cap the gains at x %? If not, there are numbers up to 2000% which would look weird to graph
SUA.PLOT.GOOD$RANGE_CHANGE <- ifelse(SUA.PLOT.GOOD$RANGE_CHANGE > 500, 500, SUA.PLOT.GOOD$RANGE_CHANGE)
summary(SUA.PLOT.GOOD$RANGE_CHANGE)
View(SUA.PLOT.GOOD[c("SPECIES", "SUA", "PERIOD", "CURRENT_SUITABLE", "FUTURE_SUITABLE", "GAIN_LOSS", "STABLE_PROP", "RANGE_CHANGE")])



#########################################################################################################################
## Create plots in the format needed for ggplot bar  
SUA.PLOT.GOOD.30 = subset(SUA.PLOT.GOOD, PERIOD == 30)
SUA.PLOT.GOOD.50 = subset(SUA.PLOT.GOOD, PERIOD == 50)
SUA.PLOT.GOOD.70 = subset(SUA.PLOT.GOOD, PERIOD == 70)
unique(SUA.PLOT.GOOD.30$PERIOD);unique(SUA.PLOT.GOOD.50$PERIOD);unique(SUA.PLOT.GOOD.70$PERIOD)
dim(SUA.PLOT.GOOD.30);dim(SUA.PLOT.GOOD.50);dim(SUA.PLOT.GOOD.70)


## Melt the table into the right format : we are now tabulating by species and range change
SUA.PLOT.30          = table(SUA.PLOT.GOOD.30$SPECIES, SUA.PLOT.GOOD.30$RANGE_CHANGE)
SUA.PLOT.50          = table(SUA.PLOT.GOOD.50$SPECIES, SUA.PLOT.GOOD.50$RANGE_CHANGE)
SUA.PLOT.70          = table(SUA.PLOT.GOOD.70$SPECIES, SUA.PLOT.GOOD.70$RANGE_CHANGE)
SUA.PLOT.30.M        = melt(SUA.PLOT.30)
SUA.PLOT.50.M        = melt(SUA.PLOT.50)
SUA.PLOT.70.M        = melt(SUA.PLOT.70)


names(SUA.PLOT.30.M) = c("SPECIES", "RANGE_CHANGE", "SUA_COUNT")
names(SUA.PLOT.50.M) = c("SPECIES", "RANGE_CHANGE", "SUA_COUNT")
names(SUA.PLOT.70.M) = c("SPECIES", "RANGE_CHANGE", "SUA_COUNT")

# SUA.PLOT.30.M        = subset(SUA.PLOT.30.M, AREA_CHANGE != "NEVER")# & AREA_CHANGE != "STABLE")
# SUA.PLOT.50.M        = subset(SUA.PLOT.50.M, AREA_CHANGE != "NEVER")# & AREA_CHANGE != "STABLE")
# SUA.PLOT.70.M        = subset(SUA.PLOT.70.M, AREA_CHANGE != "NEVER")# & AREA_CHANGE != "STABLE") 
# 
# SUA.30.M.LOSS        = subset(SUA.PLOT.30.M, AREA_CHANGE %in% c("LOSS"))
# SUA.30.M.GAIN        = subset(SUA.PLOT.30.M, AREA_CHANGE %in% c("GAIN"))
# SUA.30.M.STABLE      = subset(SUA.PLOT.30.M, AREA_CHANGE %in% c("STABLE"))
# 
# SUA.50.M.LOSS        = subset(SUA.PLOT.50.M, AREA_CHANGE %in% c("LOSS"))
# SUA.50.M.GAIN        = subset(SUA.PLOT.50.M, AREA_CHANGE %in% c("GAIN"))
# SUA.50.M.STABLE      = subset(SUA.PLOT.50.M, AREA_CHANGE %in% c("STABLE"))
# 
# SUA.70.M.LOSS        = subset(SUA.PLOT.70.M, AREA_CHANGE %in% c("LOSS"))
# SUA.70.M.GAIN        = subset(SUA.PLOT.70.M, AREA_CHANGE %in% c("GAIN"))
# SUA.70.M.STABLE      = subset(SUA.PLOT.70.M, AREA_CHANGE %in% c("STABLE"))

head(SUA.PLOT.30.M)
dim(SUA.PLOT.30.M)


#########################################################################################################################
## Attach the climate
SUA.PLOT.30.M = join(SUA.PLOT.30.M, SUA.NICHE[c("SPECIES", "Annual_mean_temp_median", "Max_temp_warm_month_median")])
SUA.PLOT.50.M = join(SUA.PLOT.50.M, SUA.NICHE[c("SPECIES", "Annual_mean_temp_median", "Max_temp_warm_month_median")])
SUA.PLOT.70.M = join(SUA.PLOT.70.M, SUA.NICHE[c("SPECIES", "Annual_mean_temp_median", "Max_temp_warm_month_median")])
head(SUA.PLOT.30.M);dim(SUA.PLOT.30.M)


#########################################################################################################################
## Use ggplot to plot the percentage of species inside an LGA, which is being lost or gained
## Calculate the gain/loss for counts of all species?

## Gain % = gain/(stable + lost)   * 100
## Lost % = lost/(stable + lost)   * 100

## Calculate the gain/loss for for only species that are recorded within each SUA?

## Gain % = gain/(stable + lost + gained)   * 100
## Lost % = lost/(stable + lost + gained)   * 100

#########################################################################################################################
## Create PNG output for all SUAs for 2030, ordered by mean annual temperature
png(sprintf('output/figures/SUA_percent/ALL_SUA_BAR_PLOT_%s_%s.png', 2030, 'temp'),      
    10, 8, units = 'in', res = 500)

## 2030
ggplot(SUA.PLOT.30.M,  aes(x = reorder(SUA, CURRENT_MAT), fill = AREA_CHANGE)) + 
  
  ## Plot species being lost from each SUA
  geom_bar(data        = subset(SUA.PLOT.30.M, AREA_CHANGE %in% c("LOSS")),
           
           ## Create % lost: the species lost, divided by those that were both lost + those always there
           aes(y = -(SPECIES_COUNT/(SUA.30.M.LOSS$SPECIES_COUNT + SUA.30.M.STABLE$SPECIES_COUNT))), 
           position = "stack", stat = "identity") +
  
  ## Plot species being gained inside each SUA
  geom_bar(data  = subset(SUA.PLOT.30.M, AREA_CHANGE %in% c("GAIN")), 
           
           ## Create % gained: the species gained, divided by those that were both lost + those always there
           aes(y = (SPECIES_COUNT/(SUA.30.M.LOSS$SPECIES_COUNT + SUA.30.M.STABLE$SPECIES_COUNT))), 
           position = "stack", stat = "identity") +
  
  scale_y_continuous(labels = scales::percent) +
  
  scale_fill_manual(values = rev(colorRampPalette(c('brown1', 'seagreen3'))(2))) +
  
  labs(#title = "Predicted native species gain/loss (number) within SUAs to 2030", 
    x = "SUA by increasing MAT", y = "% Original species")  +
  
  ## Format axes
  theme(axis.title.x     = element_text(face = "bold", colour = "black", size = 15),
        axis.text.x      = element_text(angle = 90, vjust = 0.5, size = 8),
        axis.title.y     = element_text(face = "bold", colour = "black", size = 15),
        axis.text.y      = element_text(vjust = 0.5, size = 12),
        title            = element_text(face = "bold", colour = "black", size = 15),
        legend.title     = element_text(face = "bold", colour = "black", size = 12),
        legend.text      = element_text(face = "bold", size = 12),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border     = element_rect(colour = "black", fill = NA, size = 2)) +
  
  guides(fill = guide_legend(title = "Current vs. 2030"))

dev.off()


#########################################################################################################################
## Create PNG output for all SUAs for 2070, ordered by mean annual temperature
png(sprintf('output/figures/SUA_percent/ALL_SUA_BAR_PLOT_%s_%s.png', 2070, 'temp'),      
    10, 8, units = 'in', res = 500)

## 2070
ggplot(SUA.PLOT.70.M,  aes(x = reorder(SUA, CURRENT_MAT), fill = AREA_CHANGE)) + 
  
  ## Plot species being lost from each SUA
  geom_bar(data        = subset(SUA.PLOT.70.M, AREA_CHANGE %in% c("LOSS")),
           
           ## Create % lost: the species lost, divided by those that were both lost + those always there
           aes(y = -(SPECIES_COUNT/(SUA.70.M.LOSS$SPECIES_COUNT + SUA.70.M.STABLE$SPECIES_COUNT))), 
           position = "stack", stat = "identity") +
  
  ## Plot species being gained inside each SUA
  geom_bar(data  = subset(SUA.PLOT.70.M, AREA_CHANGE %in% c("GAIN")), 
           
           ## Create % gained: the species gained, divided by those that were both lost + those always there
           aes(y = (SPECIES_COUNT/(SUA.70.M.LOSS$SPECIES_COUNT + SUA.70.M.STABLE$SPECIES_COUNT))), 
           position = "stack", stat = "identity") +
  
  scale_y_continuous(labels = scales::percent) +
  
  scale_fill_manual(values = rev(colorRampPalette(c('brown1', 'seagreen3'))(2))) +
  
  labs(#title = "Predicted native species gain/loss (number) within SUAs to 2070", 
    x = "SUA by increasing MAT", y = "% Original species")  +
  
  ## Format axes
  theme(axis.title.x     = element_text(face = "bold", colour = "black", size = 15),
        axis.text.x      = element_text(angle = 90, vjust = 0.5, size = 8),
        axis.title.y     = element_text(face = "bold", colour = "black", size = 15),
        axis.text.y      = element_text(vjust = 0.5, size = 12),
        title            = element_text(face = "bold", colour = "black", size = 15),
        legend.title     = element_text(face = "bold", colour = "black", size = 12),
        legend.text      = element_text(face = "bold", size = 12),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border     = element_rect(colour = "black", fill = NA, size = 2)) +
  
  guides(fill = guide_legend(title = "Current vs. 2070"))

dev.off()


#########################################################################################################################
## Create PNG output for all SUAs for 2070, ordered by area
png(sprintf('output/figures/SUA_percent/ALL_SUA_BAR_PLOT_%s_%s.png', 2070, 'area'),      
    10, 8, units = 'in', res = 500)

## 2070
ggplot(SUA.PLOT.70.M,  aes(x = reorder(SUA, AREASQKM16), fill = AREA_CHANGE)) + 
  
  ## Plot species being lost from each SUA
  geom_bar(data        = subset(SUA.PLOT.70.M, AREA_CHANGE %in% c("LOSS")),
           
           ## Create % lost: the species lost, divided by those that were both lost + those always there
           aes(y = -(SPECIES_COUNT/(SUA.70.M.LOSS$SPECIES_COUNT + SUA.70.M.STABLE$SPECIES_COUNT))), 
           position = "stack", stat = "identity") +
  
  ## Plot species being gained inside each SUA
  geom_bar(data  = subset(SUA.PLOT.70.M, AREA_CHANGE %in% c("GAIN")), 
           
           ## Create % gained: the species gained, divided by those that were both lost + those always there
           aes(y = (SPECIES_COUNT/(SUA.70.M.LOSS$SPECIES_COUNT + SUA.70.M.STABLE$SPECIES_COUNT))), 
           position = "stack", stat = "identity") +
  
  scale_y_continuous(labels = scales::percent) +
  
  scale_fill_manual(values = rev(colorRampPalette(c('brown1', 'seagreen3'))(2))) +
  
  labs(#title = "Predicted native species gain/loss (number) within SUAs to 2070", 
    x = "SUA by increasing Area", y = "% Original species")  +
  
  ## Format axes
  theme(axis.title.x     = element_text(face = "bold", colour = "black", size = 15),
        axis.text.x      = element_text(angle = 90, vjust = 0.5, size = 8),
        axis.title.y     = element_text(face = "bold", colour = "black", size = 15),
        axis.text.y      = element_text(vjust = 0.5, size = 12),
        title            = element_text(face = "bold", colour = "black", size = 15),
        legend.title     = element_text(face = "bold", colour = "black", size = 12),
        legend.text      = element_text(face = "bold", size = 12),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border     = element_rect(colour = "black", fill = NA, size = 2)) +
  
  guides(fill = guide_legend(title = "Current vs. 2070"))

dev.off()


#########################################################################################################################
## Try plotting a subset of just the mapped SUAs
SUA.PLOT.30.M$SUA = gsub("Canberra - Queanbeyan", "Canberra", SUA.PLOT.30.M$SUA)
SUA.PLOT.50.M$SUA = gsub("Canberra - Queanbeyan", "Canberra", SUA.PLOT.50.M$SUA)
SUA.PLOT.70.M$SUA = gsub("Canberra - Queanbeyan", "Canberra", SUA.PLOT.70.M$SUA)

SUA.PLOT.MAJOR.30 = SUA.PLOT.30.M[SUA.PLOT.30.M$SUA %in% MAP_SUA, ]
SUA.PLOT.MAJOR.50 = SUA.PLOT.50.M[SUA.PLOT.50.M$SUA %in% MAP_SUA, ] 
SUA.PLOT.MAJOR.70 = SUA.PLOT.70.M[SUA.PLOT.70.M$SUA %in% MAP_SUA, ] 
unique(SUA.PLOT.MAJOR.30$SUA);unique(SUA.PLOT.MAJOR.70$SUA)


## Create a subset of just major cities
SUA.30.MJ.LOSS        = subset(SUA.PLOT.MAJOR.30, AREA_CHANGE %in% c("LOSS"))
SUA.30.MJ.GAIN        = subset(SUA.PLOT.MAJOR.30, AREA_CHANGE %in% c("GAIN"))
SUA.30.MJ.STABLE      = subset(SUA.PLOT.MAJOR.30, AREA_CHANGE %in% c("STABLE"))

SUA.50.MJ.LOSS        = subset(SUA.PLOT.MAJOR.50, AREA_CHANGE %in% c("LOSS"))
SUA.50.MJ.GAIN        = subset(SUA.PLOT.MAJOR.50, AREA_CHANGE %in% c("GAIN"))
SUA.50.MJ.STABLE      = subset(SUA.PLOT.MAJOR.50, AREA_CHANGE %in% c("STABLE"))

SUA.70.MJ.LOSS        = subset(SUA.PLOT.MAJOR.70, AREA_CHANGE %in% c("LOSS"))
SUA.70.MJ.GAIN        = subset(SUA.PLOT.MAJOR.70, AREA_CHANGE %in% c("GAIN"))
SUA.70.MJ.STABLE      = subset(SUA.PLOT.MAJOR.70, AREA_CHANGE %in% c("STABLE"))



#########################################################################################################################
## Create PNG output for Major captuials, 2030
png(sprintf('output/figures/SUA_percent/CAPITAL_SUA_BAR_PLOT_%s_%s.png', 2030, 'temp'),      
    10, 8, units = 'in', res = 500)

## 2030
ggplot(SUA.PLOT.MAJOR.30,  aes(x = reorder(SUA, CURRENT_MAT), fill = AREA_CHANGE)) + 
  
  ## Plot species being lost from each SUA
  geom_bar(data        = subset(SUA.PLOT.MAJOR.30, AREA_CHANGE %in% c("LOSS")),
           
           ## Create % lost: the species lost, divided by those that were both lost + those always there
           aes(y = -(SPECIES_COUNT/(SUA.30.MJ.LOSS$SPECIES_COUNT + SUA.30.MJ.STABLE$SPECIES_COUNT))), 
           position = "stack", stat = "identity") +
  
  ## Plot species being gained inside each SUA
  geom_bar(data  = subset(SUA.PLOT.MAJOR.30, AREA_CHANGE %in% c("GAIN")), 
           
           ## Create % gained: the species gained, divided by those that were both lost + those always there
           aes(y = (SPECIES_COUNT/(SUA.30.MJ.LOSS$SPECIES_COUNT + SUA.30.MJ.STABLE$SPECIES_COUNT))), 
           position = "stack", stat = "identity") +
  
  scale_y_continuous(labels = scales::percent) +
  
  scale_fill_manual(values = rev(colorRampPalette(c('brown1', 'seagreen3'))(2))) +
  
  labs(x = "SUA by increasing MAT", y = "% Original species")  +
  
  ## Format axes
  theme(axis.title.x     = element_text(face = "bold", colour = "black", size = 15),
        axis.text.x      = element_text(angle = 45, vjust = 0.5, size = 12),
        axis.title.y     = element_text(face = "bold", colour = "black", size = 15),
        axis.text.y      = element_text(vjust = 0.5, size = 12),
        title            = element_text(face = "bold", colour = "black", size = 15),
        legend.title     = element_text(face = "bold", colour = "black", size = 12),
        legend.text      = element_text(face = "bold", size = 12),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border     = element_rect(colour = "black", fill = NA, size = 2)) +
  
  guides(fill = guide_legend(title = "Current vs. 2030"))

dev.off()



#########################################################################################################################
## Create PNG output for Major captuials, 2070
png(sprintf('output/figures/SUA_percent/CAPITAL_SUA_BAR_PLOT_%s_%s.png', 2070, 'temp'),      
    10, 8, units = 'in', res = 500)

## 2070
ggplot(SUA.PLOT.MAJOR.70,  aes(x = reorder(SUA, CURRENT_MAT), fill = AREA_CHANGE)) + 
  
  ## Plot species being lost from each SUA
  geom_bar(data        = subset(SUA.PLOT.MAJOR.70, AREA_CHANGE %in% c("LOSS")),
           
           ## Create % lost: the species lost, divided by those that were both lost + those always there
           aes(y = -(SPECIES_COUNT/(SUA.70.MJ.LOSS$SPECIES_COUNT + SUA.70.MJ.STABLE$SPECIES_COUNT))), 
           position = "stack", stat = "identity") +
  
  ## Plot species being gained inside each SUA
  geom_bar(data  = subset(SUA.PLOT.MAJOR.70, AREA_CHANGE %in% c("GAIN")), 
           
           ## Create % gained: the species gained, divided by those that were both lost + those always there
           aes(y = (SPECIES_COUNT/(SUA.70.MJ.LOSS$SPECIES_COUNT + SUA.70.MJ.STABLE$SPECIES_COUNT))), 
           position = "stack", stat = "identity") +
  
  scale_y_continuous(labels = scales::percent) +
  
  scale_fill_manual(values = rev(colorRampPalette(c('brown1', 'seagreen3'))(2))) +
  
  labs(x = "SUA by increasing MAT", y = "% Original species")  +
  
  ## Format axes
  theme(axis.title.x     = element_text(face = "bold", colour = "black", size = 15),
        axis.text.x      = element_text(angle = 45, vjust = 0.5, size = 12),
        axis.title.y     = element_text(face = "bold", colour = "black", size = 15),
        axis.text.y      = element_text(vjust = 0.5, size = 12),
        title            = element_text(face = "bold", colour = "black", size = 15),
        legend.title     = element_text(face = "bold", colour = "black", size = 12),
        legend.text      = element_text(face = "bold", size = 12),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border     = element_rect(colour = "black", fill = NA, size = 2)) +
  
  guides(fill = guide_legend(title = "Current vs. 2070"))

dev.off()


#########################################################################################################################
## Create PNG output for Major captuials, 2070
png(sprintf('output/figures/SUA_percent/CAPITAL_SUA_BAR_PLOT_%s_%s.png', 2070, 'area'),      
    10, 8, units = 'in', res = 500)

## 2070
ggplot(SUA.PLOT.MAJOR.70,  aes(x = reorder(SUA, AREASQKM16), fill = AREA_CHANGE)) + 
  
  ## Plot species being lost from each SUA
  geom_bar(data        = subset(SUA.PLOT.MAJOR.70, AREA_CHANGE %in% c("LOSS")),
           
           ## Create % lost: the species lost, divided by those that were both lost + those always there
           aes(y = -(SPECIES_COUNT/(SUA.70.MJ.LOSS$SPECIES_COUNT + SUA.70.MJ.STABLE$SPECIES_COUNT))), 
           position = "stack", stat = "identity") +
  
  ## Plot species being gained inside each SUA
  geom_bar(data  = subset(SUA.PLOT.MAJOR.70, AREA_CHANGE %in% c("GAIN")), 
           
           ## Create % gained: the species gained, divided by those that were both lost + those always there
           aes(y = (SPECIES_COUNT/(SUA.70.MJ.LOSS$SPECIES_COUNT + SUA.70.MJ.STABLE$SPECIES_COUNT))), 
           position = "stack", stat = "identity") +
  
  scale_y_continuous(labels = scales::percent) +
  
  scale_fill_manual(values = rev(colorRampPalette(c('brown1', 'seagreen3'))(2))) +
  
  labs(x = "SUA by increasing Area", y = "% Original species")  +
  
  ## Format axes
  theme(axis.title.x     = element_text(face = "bold", colour = "black", size = 15),
        axis.text.x      = element_text(angle = 45, vjust = 0.5, size = 12),
        axis.title.y     = element_text(face = "bold", colour = "black", size = 15),
        axis.text.y      = element_text(vjust = 0.5, size = 12),
        title            = element_text(face = "bold", colour = "black", size = 15),
        legend.title     = element_text(face = "bold", colour = "black", size = 12),
        legend.text      = element_text(face = "bold", size = 12),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border     = element_rect(colour = "black", fill = NA, size = 2)) +
  
  guides(fill = guide_legend(title = "Current vs. 2070"))

dev.off()


#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################