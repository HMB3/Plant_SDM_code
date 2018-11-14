#########################################################################################################################
################################# SUMMARISE THE HABITAT SUITABILITY RESULTS ############################################# 
#########################################################################################################################


## This code combines the tables of area occupied by each species in each SUA, and converts it to one big table. 
## This table can be queried to create histograms, etc.  


#########################################################################################################################
## 1). COMBINE OVERALL TABLES OF SPECIES GAIN/LOSS ACROSS AUSTRALIA
#########################################################################################################################


## Read in the table of species in SUAs
SUA.SPP.COUNT = readRDS(paste0('data/base/HIA_LIST/COMBO/SUA_SPP_COUNT_', save_run, '.rds'))
length(unique(SUA.SPP.COUNT$SPECIES))


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
length(GAIN.LOSS.TABLE.COMPLETE$SPECIES)/12;length(map_spp)
summary(GAIN.LOSS.TABLE.COMPLETE)


## How big is the table?
dim(GAIN.LOSS.TABLE)
dim(GAIN.LOSS.TABLE.COMPLETE)
head(GAIN.LOSS.TABLE.COMPLETE, 18)


## How many species are lost?
length(unique(GAIN.LOSS.TABLE$SPECIES))                                   
length(unique(GAIN.LOSS.TABLE.COMPLETE$SPECIES)) 





#########################################################################################################################
## 2). COMBINE TABLES OF SPECIES PRESENCES IN SUAs
#########################################################################################################################


#########################################################################################################################
## The multiple thresholds could present a problem
SUA.tables = list.files(maxent_dir, pattern = 'SUA_cell_count', full.names = TRUE, recursive = TRUE) 
length(SUA.tables)/3;length(map_spp) 


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
## ......................................................................................................................


#########################################################################################################################
## Join on a column for if the species has records in the SUA
SUA.PREDICT = merge(SUA.COMPLETE, SDM.CHECK, all.x = TRUE)
SUA.PREDICT = join(SUA.PREDICT, SUA.SPP.COUNT, type = "full")
#SUA.PREDICT = subset(SUA.PREDICT,  SUA_RECORDS > 0)
summary(SUA.PREDICT$SUA_RECORDS)
#SUA.PREDICT = completeFun(SUA.PREDICT, "MAXENT_RATING")
unique(SUA.PREDICT$MAXENT_RATING)
length(unique(SUA.PREDICT$SPECIES))
identical(dim(SUA.PREDICT)[1], dim(SUA.COMPLETE)[1])





#########################################################################################################################
## 3). CREATE PLOTS OF SPECIES LOSS AND GAIN
#########################################################################################################################


#########################################################################################################################
## Future rasters
SUA.BIO1.current.stats  = read.csv("./data/base/worldclim/aus/1km/bio/SUA_BIO1_current_stats.csv",  stringsAsFactors = FALSE)
SUA.BIO12.current.stats = read.csv("./data/base/worldclim/aus/1km/bio/SUA_BIO12_current_stats.csv", stringsAsFactors = FALSE)
SUA.KOP                 = read.csv("./data/base/worldclim/aus/1km/bio/SUA_KOPPEN.csv",              stringsAsFactors = FALSE)
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


## Merge MAP onto the results table
SUA.PREDICT = join(SUA.PREDICT, SUA.BIO1.current.stats[c("SUA", "CURRENT_MAT")])
SUA.PREDICT = join(SUA.PREDICT, SUA.BIO12.current.stats[c("SUA", "CURRENT_MAP")])
SUA.PREDICT = join(SUA.PREDICT, SUA.KOP[c("SUA", "MAJOR_KOP")])
SUA.PREDICT = join(SUA.PREDICT, SUA.PET[c("SUA", "CURRENT_PET")])
summary(SUA.PREDICT)
View(SUA.PREDICT)
#t = SUA.PREDICT[is.na(SUA.PREDICT$MAJOR_KOP),]


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
SUA.CLIM      = SUA.PREDICT[!duplicated(SUA.PREDICT[,c('SUA')]),][c("SUA", "CURRENT_MAT", "CURRENT_MAP", "CURRENT_PET")]
SUA.PLOT.30.M = join(SUA.PLOT.30.M, SUA.CLIM)
SUA.PLOT.50.M = join(SUA.PLOT.50.M, SUA.CLIM)
SUA.PLOT.70.M = join(SUA.PLOT.70.M, SUA.CLIM)
head(SUA.PLOT.30.M);dim(SUA.PLOT.30.M)
   

## Can we use ggplot to plot the percentage of species inside an LGA which is being lost or gained?
## How to calculate the gain/loss? Could do :
## Final - Now / Now *100. OR
## Gain % = gain/(stable + lost)   * 100
## Lost % = lost/(stable + lost) * 100


## Order by current temp
ggplot(SUA.PLOT.30.M,  aes(x = reorder(SUA, CURRENT_MAT), fill = AREA_CHANGE)) + 
  
  ## Subset data to make counting easier
  geom_bar(data        = subset(SUA.PLOT.30.M, AREA_CHANGE %in% c("LOSS")),
           
           ## Create % count
           aes(y = -(SPECIES_COUNT/(SUA.30.M.LOSS$SPECIES_COUNT + SUA.30.M.STABLE$SPECIES_COUNT))), 
           position = "stack", stat = "identity") +
  
  ## Subset data to make counting easier
  geom_bar(data  = subset(SUA.PLOT.30.M, AREA_CHANGE %in% c("GAIN")), 
           
           aes(y = (SPECIES_COUNT/(SUA.30.M.LOSS$SPECIES_COUNT + SUA.30.M.STABLE$SPECIES_COUNT))), 
           position = "stack", stat = "identity") +
  
  scale_y_continuous(labels = scales::percent) +
  
  scale_fill_manual(values = rev(colorRampPalette(c('brown1', 'seagreen3'))(2))) +
  
  labs(title = "Predicted native species gain/loss (number) within SUAs to 2030", 
       x = "SUA by increasing MAT", y = "% Original species gained/lost") 



#########################################################################################################################
## Create PNG output for all SUAs, 2030
png(sprintf('output/figures/SUA_BAR_PLOT/ALL_SUA_BAR_PLOT_%s.png', 2030),      
    10, 8, units = 'in', res = 500)


## 2030
ggplot(SUA.PLOT.30.M,  aes(x = reorder(SUA, CURRENT_MAT), fill = AREA_CHANGE)) + 

  geom_bar(data = subset(SUA.PLOT.30.M, AREA_CHANGE %in% c("LOSS")),
           aes(y = -SPECIES_COUNT), position = "stack", stat = "identity") +
  
  scale_fill_manual(values = rev(colorRampPalette(c('brown1', 'seagreen3'))(2))) +

  geom_bar(data = subset(SUA.PLOT.30.M, !AREA_CHANGE %in% c("LOSS")), 
           aes(y = SPECIES_COUNT), position = "stack", stat = "identity") +
  
  ## The colour scheme
  #scale_fill_manual(values=rev(colorRampPalette(c('seagreen3','brown1', 'grey', 'skyblue3'))(4))) +
  scale_fill_manual(values = rev(colorRampPalette(c('brown1', 'seagreen3'))(2))) +
  
  ## The axes labels
  labs(title = "Predicted native species gain/loss (number) within SUAs to 2030", 
       x = "SUA by increasing MAT", y = "Species Count") +
  
  ## Format axes
  theme(axis.title.x     = element_text(face = "bold", colour = "black", size = 15),
        axis.text.x      = element_text(angle = 90, vjust = 0.5, size = 8),
        axis.title.y     = element_text(face = "bold", colour = "black", size = 15),
        axis.text.y      = element_text(vjust = 0.5, size = 12),
        title            = element_text(face = "bold", colour = "black", size = 15),
        legend.title     = element_blank(),
        legend.text      = element_text(face = "bold", size = 12),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border     = element_rect(colour = "black", fill = NA, size = 2))

dev.off()


#########################################################################################################################
## Create PNG output for all SUAs, 2070
png(sprintf('output/figures/SUA_BAR_PLOT/ALL_SUA_BAR_PLOT_%s.png', 2070),      
    10, 8, units = 'in', res = 500)

## 2070
ggplot(SUA.PLOT.70.M,  aes(x = reorder(SUA, CURRENT_MAT), fill = AREA_CHANGE)) + 
  
  geom_bar(data = subset(SUA.PLOT.70.M, AREA_CHANGE %in% c("LOSS")),
           aes(y = -SPECIES_COUNT), position = "stack", stat = "identity") +
  
  scale_fill_manual(values=rev(colorRampPalette(c('brown1', 'seagreen3'))(2))) +
  
  geom_bar(data = subset(SUA.PLOT.70.M, !AREA_CHANGE %in% c("LOSS")), 
           aes(y = SPECIES_COUNT), position = "stack", stat = "identity") +
  
  ## The colour scheme
  #scale_fill_manual(values=rev(colorRampPalette(c('seagreen3','brown1', 'grey', 'skyblue3'))(4))) +
  scale_fill_manual(values=rev(colorRampPalette(c('brown1', 'seagreen3'))(2))) +
  
  ## The axes labels
  labs(title = "Predicted native species gain/loss (number) within SUAs to 2070", 
       x = "SUA by increasing MAT", y = "Species Count") +
  
  ## Format axes
  theme(axis.title.x     = element_text(face = "bold", colour = "black", size = 15),
        axis.text.x      = element_text(angle = 90, size = 8), # vjust = 0.5, 
        axis.title.y     = element_text(face = "bold", colour = "black", size = 15),
        axis.text.y      = element_text(vjust = 0.5, size = 12),
        title            = element_text(face = "bold", colour = "black", size = 15),
        legend.title     = element_blank(),
        legend.text      = element_text(face = "bold", size = 12),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border     = element_rect(colour = "black", fill = NA, size = 2))

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


# SUA.SPP.M$SUA = gsub("Canberra - Queanbeyan", "Canberra", SUA.SPP.M$SUA)
# SUA.SPP.MAJOR = SUA.SPP.M[SUA.SPP.M$SUA %in% MAP_SUA, ] 
# SUA.SP.1.MAJOR = subset(SUA.SPP.MAJOR, SPECIES == "Eucalyptus camaldulensis")
# unique(SUA.SP.1.MAJOR$SUA)
# unique(SUA.SP.1.MAJOR$SPECIES)


#########################################################################################################################
## Create PNG output for Major captuials, 2030
png(sprintf('output/figures/SUA_BAR_PLOT/SUA_BAR_PLOT_%s.png', 2030),      
    10, 8, units = 'in', res = 500)

## 2030
ggplot(SUA.PLOT.MAJOR.30,  aes(x = reorder(SUA, CURRENT_MAT), fill = AREA_CHANGE)) + 
  
  ## The species being lost
  geom_bar(data = subset(SUA.PLOT.MAJOR.30, AREA_CHANGE %in% c("LOSS")),
           aes(y = -SPECIES_COUNT/sum(SUA.PLOT.MAJOR.30$SPECIES_COUNT)), 
           position = "stack", stat = "identity", colour="black") +
  
  ## The species being gained or remaining stable
  geom_bar(data = subset(SUA.PLOT.MAJOR.30, !AREA_CHANGE %in% c("LOSS")), 
           aes(y = SPECIES_COUNT/sum(SUA.PLOT.MAJOR.30$SPECIES_COUNT)), 
           position = "stack", stat = "identity", colour="black") +
  
  scale_y_continuous(labels=scales::percent) +
  
  ## The colour scheme
  #scale_fill_manual(values=rev(colorRampPalette(c('seagreen3','brown1', 'grey', 'skyblue3'))(4))) +
  scale_fill_manual(values=rev(colorRampPalette(c('brown1', 'seagreen3'))(2))) +
  
  ## The axes labels
  labs(title = "Predicted native species gain/loss (number) within SUAs to 2030", 
       x = "SUA by increasing MAT", y = "Species Count") +

## Format axes
theme(axis.title.x     = element_text(face = "bold", colour = "black", size = 15),
      axis.text.x      = element_text(angle = 45, vjust = 0.5, size = 12),
      axis.title.y     = element_text(face = "bold", colour = "black", size = 15),
      axis.text.y      = element_text(vjust = 0.5, size = 12),
      title            = element_text(face = "bold", colour = "black", size = 15),
      legend.title     = element_blank(),
      legend.text      = element_text(face = "bold", size = 12),
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      panel.border     = element_rect(colour = "black", fill = NA, size = 2))

dev.off()



#########################################################################################################################
## Create PNG output for 2050
png(sprintf('output/figures/SUA_BAR_PLOT/SUA_BAR_PLOT_%s.png', 2050),      
    10, 8, units = 'in', res = 500)

## 2050
ggplot(SUA.PLOT.MAJOR.50,  aes(x = reorder(SUA, CURRENT_MAT), fill = AREA_CHANGE)) + 
  
  ## The species being lost
  geom_bar(data = subset(SUA.PLOT.MAJOR.50, AREA_CHANGE %in% c("LOSS")),
           aes(y = -SPECIES_COUNT), 
           position = "stack", stat = "identity", colour="black") +
  
  ## The species being gained or remaining stable
  geom_bar(data = subset(SUA.PLOT.MAJOR.50, !AREA_CHANGE %in% c("LOSS")), 
           aes(y = SPECIES_COUNT), 
           position = "stack", stat = "identity", colour="black") +
  
  ## The colour scheme
  #scale_fill_manual(values=rev(colorRampPalette(c('seagreen3','brown1', 'grey', 'skyblue3'))(4))) +
  scale_fill_manual(values=rev(colorRampPalette(c('brown1', 'seagreen3'))(2))) +
  
  ## The axes labels
  labs(title = "Predicted native species gain/loss (number) within SUAs to 2050", 
       x = "SUA by increasing MAT", y = "Species Count") +
  
  ## Format axes
  theme(axis.title.x     = element_text(face = "bold", colour = "black", size = 15),
        axis.text.x      = element_text(angle = 45, vjust = 0.5, size = 12),
        axis.title.y     = element_text(face = "bold", colour = "black", size = 15),
        axis.text.y      = element_text(vjust = 0.5, size = 12),
        title            = element_text(face = "bold", colour = "black", size = 15),
        legend.title     = element_blank(),
        legend.text      = element_text(face = "bold", size = 12),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border     = element_rect(colour = "black", fill = NA, size = 2))

dev.off()


#########################################################################################################################
## Create PNG output for 2070
png(sprintf('output/figures/SUA_BAR_PLOT/SUA_BAR_PLOT_%s.png', 2070),      
    10, 8, units = 'in', res = 500)

## 2070
ggplot(SUA.PLOT.MAJOR.70,  aes(x = reorder(SUA, CURRENT_MAT), fill = AREA_CHANGE)) + 
  
  ## The species being lost
  geom_bar(data = subset(SUA.PLOT.MAJOR.70, AREA_CHANGE %in% c("LOSS")),
           aes(y = -SPECIES_COUNT), 
           position = "stack", stat = "identity", colour="black") +
  
  ## The species being gained or remaining stable
  geom_bar(data = subset(SUA.PLOT.MAJOR.70, !AREA_CHANGE %in% c("LOSS")), 
           aes(y = SPECIES_COUNT), 
           position = "stack", stat = "identity", colour="black") +
  
  ## The colour scheme
  #scale_fill_manual(values=rev(colorRampPalette(c('seagreen3','brown1', 'grey', 'skyblue3'))(4))) +
  scale_fill_manual(values=rev(colorRampPalette(c('brown1', 'seagreen3'))(2))) +
  
  ## The axes labels
  labs(title = "Predicted native species gain/loss (number) within SUAs to 2070", 
       x = "SUA by increasing MAT", y = "Species Count") +
  
  ## Format axes
  theme(axis.title.x     = element_text(face = "bold", colour = "black", size = 15),
        axis.text.x      = element_text(angle = 45, vjust = 0.5, size = 12),
        axis.title.y     = element_text(face = "bold", colour = "black", size = 15),
        axis.text.y      = element_text(vjust = 0.5, size = 12),
        title            = element_text(face = "bold", colour = "black", size = 15),
        legend.title     = element_blank(),
        legend.text      = element_text(face = "bold", size = 12),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border     = element_rect(colour = "black", fill = NA, size = 2))

dev.off()



#########################################################################################################################
## Now, lets try excluding any species that has a 0 count in that SUA


#########################################################################################################################
## Could create a loop to plot each species
# png(sprintf('output/figures/SUA_BAR_PLOT/SUA_BAR_PLOT_%s.png', save_run),      
#     15, 8, units = 'in', res = 500)

# ggplot(SUA.SP.1.MAJOR, aes(x = SUA, fill = AREA_CHANGE)) + 
#   
#   ## The species being lost
#   geom_bar(data = subset(SUA.SP.1.MAJOR, AREA_CHANGE %in% c("LOSS")),
#            aes(y = -SPECIES_COUNT), 
#            position = "stack", stat = "identity", colour="black") +
#   
#   ## The species being gained or remaining stable
#   geom_bar(data = subset(SUA.SP.1.MAJOR, !AREA_CHANGE %in% c("LOSS")), 
#            aes(y = SPECIES_COUNT), 
#            position = "stack", stat = "identity", colour="black") +
#   
#   ## The colour scheme
#   scale_fill_manual(values=rev(colorRampPalette(c('seagreen3','brown1', 'grey', 'skyblue3'))(4))) +
#   
#   ## The axes labels
#   labs(title = "Predicted gains and losses within Significant Urban Areas (SUAs)", 
#        x = "SUA", y = "Species Count") +
#   
#   ## Format axes
#   theme(axis.title.x     = element_text(face = "bold", colour = "black", size = 15),
#         axis.text.x      = element_text(angle = 45, vjust = 0.5, size = 12),
#         axis.title.y     = element_text(face = "bold", colour = "black", size = 15),
#         axis.text.y      = element_text(vjust = 0.5, size = 12),
#         title            = element_text(face = "bold", colour = "black", size = 20),
#         legend.title     = element_blank(),
#         legend.text      = element_text(face = "bold", size = 12),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"),
#         panel.border     = element_rect(colour = "black", fill = NA, size = 2))
# 
# dev.off()



#########################################################################################################################
## 4). CREATE A RICHNESS MAP FOR EACH TIME PERIOD
#########################################################################################################################


## Probably don't need this anymore - just use ArcMap to sum the rasters................................................
## Search for


#########################################################################################################################
## Look through the output directory for species where the code didn't finish. This is usually species with < 27 files
## Summing 200 rasters will take ages 
raster.dirs <- list.dirs(path = maxent_path, full.names = FALSE, recursive = FALSE)


# ## Loop over the directories for the best species: current maps
# raster.current <- sapply(raster.dirs, function(x) {
#   
#   ## List the files for that time slice
#   list.files(paste0(maxent_path, x), pattern = 'current_suit_above', full.names = TRUE, recursive = TRUE)
#   
# })
# 
# 
# 
# ## 2030 maps
# raster.2030 <- sapply(raster.dirs, function(x) {
#   
#   ## List the files for that time slice
#   list.files(paste0(maxent_path, x), pattern = '2030_4GCMs_above', full.names = TRUE, recursive = TRUE)
#   
# })
# 
# 
# ## 2050 maps
# raster.2050 <- sapply(raster.dirs, function(x) {
#   
#   ## List the files for that time slice
#   list.files(paste0(maxent_path, x), pattern = '2050_4GCMs_above', full.names = TRUE, recursive = TRUE)
#   
# })
# 
# 
# ## 2070 maps
# raster.2070 <- sapply(raster.dirs, function(x) {
#   
#   ## List the files for that time slice
#   list.files(paste0("./output/maxent/SET_VAR_KOPPEN/", x), pattern = '2070_4GCMs_above', full.names = TRUE, recursive = TRUE)
#   
# })
# 
# 
# #########################################################################################################################
# ## Now unlist so we can create a raster stack
# raster.current  = unlist(raster.current)
# raster.2030  = unlist(raster.2030)
# raster.2050  = unlist(raster.2050)
# raster.2070  = unlist(raster.2070)
# 
# 
# ## Check length
# length(raster.current);length(raster.2030);length(raster.2050);length(raster.2070)
# 
# 
# ## Then create raster stacks and sum
# stack.current   = stack(raster.current, quick = TRUE)
# stack.2030      = stack(raster.2030, quick = TRUE)
# stack.2050      = stack(raster.2050, quick = TRUE)
# stack.2070      = stack(raster.2070, quick = TRUE)
# 
# 
# ## writeRaster(stack.2030, filename = 'output/maxent/multilayer.tif', options = "INTERLEAVE=BAND", overwrite = TRUE)
# ## mystack = stack("multilayer.tif")
# 
# 
# ## Summing takes a long time
# sum.current  = sum(stack.current) 
# sum.2030     = sum(stack.2030) 
# sum.2050     = sum(stack.2050) 
# sum.2070     = sum(stack.2070) 
# 
# 
# #########################################################################################################################
# ## Plot to check
# plot(sum.current)
# plot(sum.2030)
# plot(sum.2050)
# plot(sum.2070)
# 
# 
# #########################################################################################################################
# ## write out rasters
# writeRaster(sum.current, 'output/maxent/checked_spp_current_richness.tif')
# writeRaster(sum.2030,    'output/maxent/checked_spp_2030_richness.tif')
# writeRaster(sum.2050,    'output/maxent/checked_spp_2050_richness.tif')
# writeRaster(sum.2070,    'output/maxent/checked_spp_2070_richness.tif')



#########################################################################################################################
## OUTSTANDING PREDICTION TASKS:
#########################################################################################################################


#########################################################################################################################
## 





#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################