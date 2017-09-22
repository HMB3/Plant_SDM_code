#########################################################################################################################
##############################  CALCULATE AREA OF OCCUPANCY RANGES ###################################################### 
#########################################################################################################################


#setwd("D:/backup/D_workstation/GPP_Beta_analysis/git/chapter_3")
## install.packages("dplyr")
## install.packages("plyr")
## install.packages("Hmisc")
##  install.packages("splancs")
##  install.packages("Cairo")
library(plyr)
library(dplyr)
library(Cairo)
library(Hmisc)
library(raster)
library(Cairo)
library(mgcv)
source("BLOCK_SPECIES_RANDOMISATION_FUNCTIONS.R")


#####################################################################################################
## read data
data    = load("./data/./data/base/HIA_LIST/COMBO/COMBO_GBIF_ALA_RASTER.RData")
points  = data.matrix(data["lon", "lat"])
View(data)


#####################################################################################################
## create raster for biodata lat/long bounding box using Aus dimensions, and Rachael Gallagher's 
## 10km * 10km area of occupancy measure
n.x = ceiling(153.64 - 112.92) / 0.1
n.y = ceiling(-10.51-(-43.64)) / 0.1


## create index of all raster rows and columns for AUS at 10km * 10km resolution
index.raster = raster(nrows = n.y,    ncols = n.x, 
                      xmn   = 112.92, xmx   = 153.64, 
                      ymn   = -43.64, ymx   = -10.51, 
                      vals  = c(1:(n.x * n.y)))
plot(index.raster)


## extract raster values
block.index = extract(index.raster, points, method = "simple")
data.index  = cbind(data, block.index) 


## get species IDs as 'unique' on data index: known vascular flora of the WT ~4,300
spp.ids = unique(data.index$searchTaxon)
length(spp.ids)


## create table of area of occupancy for each species 
area.table           = cbind(spp.ids, NA) 
colnames(area.table) = c("searchTaxon", "area_occurrence")





#####################################################################################################
## loop through all the species IDs, counting the number of grid cells each species occurs in
for(i.spp in 1:length(spp.ids)) {
  
  ## create a list of blocks that each species is found in
  spp.blocks = data.index$block.index[data.index$searchTaxon == spp.ids[i.spp]]
  
  ## get unique blocks
  unique.block = unique(spp.blocks)
  
  ## area of occurrence: this needs to change for block size of 10km, 10 * 10 = 100 km squared
  area.table[i.spp, 2] = length(unique.block) * 100
  area.table = as.data.frame(area.table)
  
}


## check area.table
dim(area.table)
area.table = area.table[order(area.table[,"searchTaxon"]),]
head(area.table)


############################################################################################################
## now combine geographic and environmental ranges
## first, read in the species summaries
all       = load("./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT.RData")
dim(all)
head(all)


############################################################################################################
## merge area of occupancy with species niche widths
spp.var.geo.range = merge(all, area.table, by = "searchTaxon", all = FALSE)


## check
dim(spp.var.geo.range)
str(spp.var.geo.range)
spp.var.geo.range$log_area_occurrence = log(spp.var.geo.range$area_occurrence)
summary(spp.var.geo.range$area_occurrence)      ## minimum of 2 blocks a species is found in, = 200 km2
summary(spp.var.geo.range$log_area_occurrence)


## How does the "Freq" column compare? It is not exactly the same, because the count
## of indiviudal records from the table function just counts points, whereas our block
## index approach counts the 10km cells each species occurs in, so should always be less than the "Freq"...?
## so E.G. max Freq = 336600, not 10350 * 100 because obviously the recrods are spatial autocorrelated in the
## same grid cells, etc.
head(area.table)  ## E.G. species 3 had 116 records, that occurred in 9300/100 = 93 unique 10 km *10 km blocks 
head(spp.var.geo.range[c("searchTaxon", "Freq", "area_occurrence", "log_area_occurrence")])


summary(spp.var.geo.range$Freq)
summary(spp.var.geo.range$area_occurrence)
summary(spp.var.geo.range$log_area_occurrence)


## to make the biome map in Fig 2 more useful, perhaps plot Cryptocarya hypospodia over the top?
Cryptocarya_hypospodia = subset(all, binomial == "Cryptocarya hypospodia")
dim(Cryptocarya_hypospodia)
View(Cryptocarya_hypospodia)




#########################################################################################################################
##############################  END CALCULATE AREA OF OCCUPANCY RANGES ################################################## 
#########################################################################################################################