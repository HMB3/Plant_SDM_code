#########################################################################################################################
## Join on a column for if the species has records in the SUA
SUA.PREDICT = merge(SUA.COMPLETE, SDM.CHECK, all.x = TRUE)
SUA.PREDICT = join(SUA.PREDICT, SUA.SPP.COUNT, type = "full")
#SUA.PREDICT = subset(SUA.PREDICT,  SUA_RECORDS > 0)
summary(SUA.PREDICT$SUA_RECORDS)
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
SUA.CLIM      = SUA.PREDICT[!duplicated(SUA.PREDICT[,c('SUA')]),][c("SUA", "CURRENT_MAT", "CURRENT_MAP", "CURRENT_PET", "AREASQKM16")]
SUA.PLOT.30.M = join(SUA.PLOT.30.M, SUA.CLIM)
SUA.PLOT.50.M = join(SUA.PLOT.50.M, SUA.CLIM)
SUA.PLOT.70.M = join(SUA.PLOT.70.M, SUA.CLIM)
head(SUA.PLOT.30.M);dim(SUA.PLOT.30.M)


#########################################################################################################################
## Use ggplot to plot the percentage of species inside an LGA, which is being lost or gained
## How to calculate the gain/loss? Could do :
## Final - Now / Now *100. OR
## Gain % = gain/(stable + lost)   * 100
## Lost % = lost/(stable + lost) * 100


#########################################################################################################################
## Create PNG output for all SUAs for 2030, ordered by mean annual temperature
png(sprintf('output/figures/SUA_percent/ALL_SUA_SUA_BAR_PLOT_SUM_100_%s_%s.png', 2030, 'temp'),      
    10, 8, units = 'in', res = 500)

## 2030
ggplot(SUA.PLOT.30.M,  aes(x = reorder(SUA, CURRENT_MAT), fill = AREA_CHANGE)) + 
  
  ## Plot species being lost from each SUA
  geom_bar(data        = subset(SUA.PLOT.30.M, AREA_CHANGE %in% c("LOSS")),
           
           ## Create % lost: the species lost, divided by those that were both lost + those always there
           aes(y = -(SPECIES_COUNT/(SUA.30.M.LOSS$SPECIES_COUNT + SUA.30.M.STABLE$SPECIES_COUNT + SUA.30.M.GAIN$SPECIES_COUNT))),
           position = "stack", stat = "identity") +
  
  ## Plot species being gained inside each SUA
  geom_bar(data  = subset(SUA.PLOT.30.M, AREA_CHANGE %in% c("GAIN")), 
           
           ## Create % gained: the species gained, divided by those that were both lost + those always there
           aes(y = (SPECIES_COUNT/(SUA.30.M.LOSS$SPECIES_COUNT + SUA.30.M.STABLE$SPECIES_COUNT + SUA.30.M.GAIN$SPECIES_COUNT))),
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
SUA.PLOT.70.M   = subset(SUA.PLOT.70.M,        AREASQKM16 > 200)
SUA.70.M.LOSS   = subset(SUA.70.M.LOSS,   AREASQKM16 > 200)
SUA.70.M.GAIN   = subset(SUA.70.M.GAIN,   AREASQKM16 > 200)
SUA.70.M.STABLE = subset(SUA.70.M.STABLE, AREASQKM16 > 200)
png(sprintf('output/figures/SUA_percent/ALL_SUA_SUA_BAR_PLOT_SUM_100_%s_%s.png', 2070, 'temp'),      
    10, 8, units = 'in', res = 500)

## 2070
ggplot(SUA.PLOT.70.M,  aes(x = reorder(SUA, CURRENT_MAT), fill = AREA_CHANGE)) + 
  
  ## Plot species being lost from each SUA
  geom_bar(data        = subset(SUA.PLOT.70.M, AREA_CHANGE %in% c("LOSS")),
           
           ## Create % lost: the species lost, divided by those that were both lost + those always there
           aes(y = -(SPECIES_COUNT/(SUA.70.M.LOSS$SPECIES_COUNT + SUA.70.M.STABLE$SPECIES_COUNT + SUA.70.M.GAIN$SPECIES_COUNT))),
           position = "stack", stat = "identity") +
  
  ## Plot species being gained inside each SUA
  geom_bar(data  = subset(SUA.PLOT.70.M, AREA_CHANGE %in% c("GAIN")), 
           
           ## Create % gained: the species gained, divided by those that were both lost + those always there
           aes(y = (SPECIES_COUNT/(SUA.70.M.LOSS$SPECIES_COUNT + SUA.70.M.STABLE$SPECIES_COUNT + SUA.70.M.GAIN$SPECIES_COUNT))),
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
png(sprintf('output/figures/SUA_percent/ALL_SUA_SUA_BAR_PLOT_SUM_100_%s_%s.png', 2070, 'area'),      
    10, 8, units = 'in', res = 500)

## 2070
ggplot(SUA.PLOT.70.M,  aes(x = reorder(SUA, AREASQKM16), fill = AREA_CHANGE)) + 
  
  ## Plot species being lost from each SUA
  geom_bar(data        = subset(SUA.PLOT.70.M, AREA_CHANGE %in% c("LOSS")),
           
           ## Create % lost: the species lost, divided by those that were both lost + those always there
           aes(y = -(SPECIES_COUNT/(SUA.70.M.LOSS$SPECIES_COUNT + SUA.70.M.STABLE$SPECIES_COUNT + SUA.70.M.GAIN$SPECIES_COUNT))),
           position = "stack", stat = "identity") +
  
  ## Plot species being gained inside each SUA
  geom_bar(data  = subset(SUA.PLOT.70.M, AREA_CHANGE %in% c("GAIN")), 
           
           ## Create % gained: the species gained, divided by those that were both lost + those always there
           aes(y = (SPECIES_COUNT/(SUA.70.M.LOSS$SPECIES_COUNT + SUA.70.M.STABLE$SPECIES_COUNT + SUA.70.M.GAIN$SPECIES_COUNT))),
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


## Change % cacl to include gains when plotting SUA_RECORDS >0


#########################################################################################################################
## Create PNG output for Major captuials, 2030
png(sprintf('output/figures/SUA_percent/CAPITAL_SUA_SUA_BAR_PLOT_SUM_100_%s_%s.png', 2030, 'temp'),      
    10, 8, units = 'in', res = 500)

## 2030
ggplot(SUA.PLOT.MAJOR.30,  aes(x = reorder(SUA, CURRENT_MAT), fill = AREA_CHANGE)) + 
  
  ## Plot species being lost from each SUA
  geom_bar(data        = subset(SUA.PLOT.MAJOR.30, AREA_CHANGE %in% c("LOSS")),
           
           ## Create % lost: the species lost, divided by those that were both lost + those always there
           aes(y = -(SPECIES_COUNT/(SUA.30.MJ.LOSS$SPECIES_COUNT + SUA.30.MJ.STABLE$SPECIES_COUNT + SUA.30.MJ.GAIN$SPECIES_COUNT))),
           position = "stack", stat = "identity") +
  
  ## Plot species being gained inside each SUA
  geom_bar(data  = subset(SUA.PLOT.MAJOR.30, AREA_CHANGE %in% c("GAIN")), 
           
           ## Create % gained: the species gained, divided by those that were both lost + those always there
           aes(y = (SPECIES_COUNT/(SUA.30.MJ.LOSS$SPECIES_COUNT + SUA.30.MJ.STABLE$SPECIES_COUNT + SUA.30.MJ.GAIN$SPECIES_COUNT))),
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
png(sprintf('output/figures/SUA_percent/CAPITAL_SUA_SUA_BAR_PLOT_SUM_100_%s_%s.png', 2070, 'temp'),      
    10, 8, units = 'in', res = 500)

## 2070
ggplot(SUA.PLOT.MAJOR.70,  aes(x = reorder(SUA, CURRENT_MAT), fill = AREA_CHANGE)) + 
  
  ## Plot species being lost from each SUA
  geom_bar(data        = subset(SUA.PLOT.MAJOR.70, AREA_CHANGE %in% c("LOSS")),
           
           ## Create % lost: the species lost, divided by those that were both lost + those always there
           aes(y = -(SPECIES_COUNT/(SUA.70.MJ.LOSS$SPECIES_COUNT + SUA.70.MJ.STABLE$SPECIES_COUNT + SUA.70.MJ.GAIN$SPECIES_COUNT))),
           position = "stack", stat = "identity") +
  
  ## Plot species being gained inside each SUA
  geom_bar(data  = subset(SUA.PLOT.MAJOR.70, AREA_CHANGE %in% c("GAIN")), 
           
           ## Create % gained: the species gained, divided by those that were both lost + those always there
           aes(y = (SPECIES_COUNT/(SUA.70.MJ.LOSS$SPECIES_COUNT + SUA.70.MJ.STABLE$SPECIES_COUNT + SUA.70.MJ.GAIN$SPECIES_COUNT))),
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
png(sprintf('output/figures/SUA_percent/CAPITAL_SUA_SUA_BAR_PLOT_SUM_100_%s_%s.png', 2070, 'area'),      
    10, 8, units = 'in', res = 500)

## 2070
ggplot(SUA.PLOT.MAJOR.70,  aes(x = reorder(SUA, AREASQKM16), fill = AREA_CHANGE)) + 
  
  ## Plot species being lost from each SUA
  geom_bar(data        = subset(SUA.PLOT.MAJOR.70, AREA_CHANGE %in% c("LOSS")),
           
           ## Create % lost: the species lost, divided by those that were both lost + those always there
           aes(y = -(SPECIES_COUNT/(SUA.70.MJ.LOSS$SPECIES_COUNT + SUA.70.MJ.STABLE$SPECIES_COUNT + SUA.70.MJ.GAIN$SPECIES_COUNT))),
           position = "stack", stat = "identity") +
  
  ## Plot species being gained inside each SUA
  geom_bar(data  = subset(SUA.PLOT.MAJOR.70, AREA_CHANGE %in% c("GAIN")), 
           
           ## Create % gained: the species gained, divided by those that were both lost + those always there
           aes(y = (SPECIES_COUNT/(SUA.70.MJ.LOSS$SPECIES_COUNT + SUA.70.MJ.STABLE$SPECIES_COUNT + SUA.70.MJ.GAIN$SPECIES_COUNT))),
           position = "stack", stat = "identity") +
  
  scale_y_continuous(labels = scales::percent) +
  
  scale_fill_manual(values = rev(colorRampPalette(c('brown1', 'seagreen3'))(2))) +
  
  labs(x = "SUA by increasing Area", y = "% Original species")  +
  
  ## Format axes
  theme(axis.title.x     = element_text(face  = "bold", colour = "black", size = 15),
        axis.text.x      = element_text(angle = 45, vjust = 0.5, size = 12),
        axis.title.y     = element_text(face  = "bold", colour = "black", size = 15),
        axis.text.y      = element_text(vjust = 0.5, size = 12),
        title            = element_text(face  = "bold", colour = "black", size = 15),
        legend.title     = element_text(face  = "bold", colour = "black", size = 12),
        legend.text      = element_text(face  = "bold", size = 12),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border     = element_rect(colour = "black", fill = NA, size = 2)) +
  
  guides(fill = guide_legend(title = "Current vs. 2070"))

dev.off()


#########################################################################################################################
## OUTSTANDING PLOT TASKS:
#########################################################################################################################