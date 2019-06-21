library(ggplot2) 
library(calibrate)
library(maptools)
library(grid)
library(reshape2)
library(dplyr)
library(data.table)
library(plyr)
library(ggrepel)    # to avoid label overplotting
library(rpivotTable)


## NEW CODE 7.11.2018
# Importing Hugh's gain-loss-stable data
data<-read.csv("C:/Users/MQ20174608/Desktop/SUAs_CLIM/MAXENT_SUA_PRESENCE_SUA_ANALYSIS_NATIVE_GOOD_7.11.2018.csv", sep = ",", header=T, na.string = "NA")
data<-subset(data, PERIOD=="30" & MAXENT_RATING < 3)  
names(data)

# preparing data for plotting
SUA.PLOT   = table(data$SUA, data$GAIN_LOSS)
SUA.PLOT.M = melt(SUA.PLOT) 
names(SUA.PLOT.M)
names(SUA.PLOT.M) = c("SUA", "SPECIES_CHANGE", "SPECIES_COUNT")
dim(SUA.PLOT.M)
head(SUA.PLOT.M)

## attach SUAs data
clim.SUA = data[!duplicated(data[,c('SUA')]),] 
SUA.PLOT.M1<-merge(x=SUA.PLOT.M, y=clim.SUA, by="SUA")
head(SUA.PLOT.M1)
dim(SUA.PLOT.M1)

## PLOT SPECIES GAIN-STABLE-LOSS FOR SUAs ordered by increasing MAT not standardised by SUA area without NEVER & STABLE  ######################
SUA.PLOT.M2<-subset(SUA.PLOT.M1, SPECIES_CHANGE != "NEVER" & SPECIES_CHANGE != "STABLE")

dev.new(width=50, height=30)
ggplot(SUA.PLOT.M2, aes(x = reorder(SUA, CURRENT_MAT), fill = SPECIES_CHANGE)) +  

           geom_bar(data = subset(SUA.PLOT.M2, SPECIES_CHANGE %in% c("LOSS")),
           aes(y = -SPECIES_COUNT), position = "stack", stat = "identity") +

           scale_fill_manual(values=rev(colorRampPalette(c('brown1', 'seagreen3'))(2))) +
  
           geom_bar(data = subset(SUA.PLOT.M2, !SPECIES_CHANGE %in% c("LOSS")), 
           aes(y = SPECIES_COUNT), position = "stack", stat = "identity") +

           theme(axis.text.x = element_text(angle = 90, hjust = 1))  + 

           labs(title = "Predicted native species gains and losses (number) within Significant Urban Areas (SUAs) to 2030", 
           x = "SUAs ranked by increasing MAT", y = "Species Count") 
 