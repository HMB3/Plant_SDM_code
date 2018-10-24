#########################################################################################################################
################################# SUMMARISE THE HABITAT SUITABILITY RESULTS ############################################# 
#########################################################################################################################


## This code combines the tables of area occupied by each species in each SUA, and converts it to one big table. 
## This table can be queried to create histograms, etc.  


#########################################################################################################################
## 1). COMBINE OVERALL TABLES OF SPECIES GAIN/LOSS ACROSS AUSTRALIA
#########################################################################################################################


## Check the match between the SUA shapefile and population fields......................................................


## Read in the table of species in SUAs
SUA.SPP.COUNT = readRDS(paste0('data/base/HIA_LIST/COMBO/SUA_SPP_COUNT_', save_run, '.rds'))
length(unique(SUA.SPP.COUNT$SPECIES))


#########################################################################################################################
## Create a list of the gain/loss tables. Just run this on the final SUA table 
GAIN.LOSS.list = list.files(maxent_dir, pattern = 'gain_loss_table_', full.names = TRUE, recursive = TRUE) 


## Chrck the exceptios - if all spp models have completed successfully, should be number of species * 3 time periods
length(GAIN.LOSS.list)/3;length(map_spp)


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
summary(GAIN.LOSS.TABLE.COMPLETE)


## How big is the table?
dim(GAIN.LOSS.TABLE)
dim(GAIN.LOSS.TABLE.COMPLETE)
head(GAIN.LOSS.TABLE.COMPLETE, 18)


## How many species are lost?
length(unique(GAIN.LOSS.TABLE$SPECIES))                                   
length(unique(GAIN.LOSS.TABLE.COMPLETE$SPECIES)) ## somewhat mysterious?  


#########################################################################################################################
## Now write CSV 
write.csv(GAIN.LOSS.TABLE, paste0('output/tables/OVERALL_GAIN_LOSS_TABLE_', save_run, '.csv'), row.names = FALSE)





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


## Update this to match SUA names between SUA_NAME16 and the population table............................................
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


## Try to harmonise the SUAs in both the area and population tables.......................................................
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


## Which SUA's are missing populations? We should have the Central Coast. An updated table can vbe joined on like this
# SUA.PRESENCE = as.data.frame(rbindlist(list(SUA.PRESENCE, URB.POP), fill = TRUE)[, c('SUA', 'POP_2017') :=
#                                                              lapply(.SD, na.omit) , SUA, .SDcols = SUA:POP_2017][])


## Just check there are no NA species
SUA.COMPLETE = completeFun(SUA.PRESENCE, "CURRENT_SUITABLE")
summary(SUA.COMPLETE)


## Which species are missing?
# intersect(unique(SUA.PRESENCE[is.na(SUA.PRESENCE$POP_2017),]$SUA),
#           intersect(URB.POP$SUA, SUA.COMPLETE$SUA))


#########################################################################################################################
## Find the species with the greatest increase/decrease?
summary(SUA.COMPLETE$CHANGE)
unique(SUA.COMPLETE$GAIN_LOSS)


#########################################################################################################################
## Rank by population
## SUA.PRESENCE = SUA.PRESENCE [with(SUA.PRESENCE , rev(order(POP_2017))), ]


## Use a numerical definition of the most populated areas, OR, just take the capital cities
## MAIN_SUA = c("Sydney", "Melbourne", "Canberra - Queanbeyan", "Brisbane", "Perth", "Adelaide", "Hobart", "Darwin")
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
SDM.CHECK = MXT.CHECK[, c("searchTaxon", "Origin", "Check.map")]
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

SUA.BIO1.2030.stats    = read.csv("./data/base/worldclim/aus/1km/bio/SUA_BIO1_2030_stats.csv",    stringsAsFactors = FALSE)
SUA.BIO1.2050.stats    = read.csv("./data/base/worldclim/aus/1km/bio/SUA_BIO1_2050_stats.csv",    stringsAsFactors = FALSE)
SUA.BIO1.2070.stats    = read.csv("./data/base/worldclim/aus/1km/bio/SUA_BIO1_2050_stats.csv",    stringsAsFactors = FALSE)


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


## Can't join these data
SUA.BIO1.2030.stats[["PERIOD"]]  <- 30
colnames(SUA.BIO1.2030.stats)[1] <- "SUA"

# SUA.BIO1.2050.stats[["PERIOD"]]  <- 50
# colnames(SUA.BIO1.2050.stats)[1] <- "SUA"
# 
# SUA.BIO1.2070.stats[["PERIOD"]]  <- 50
# colnames(SUA.BIO1.2070.stats)[1] <- "SUA"
# 
# SUA.ZONAL.STATS = bind_rows(SUA.BIO1.2030.stats,
#                             SUA.BIO1.2050.stats,
#                             SUA.BIO1.2070.stats)
# SUA.ZONAL.STATS[["CLIMATE"]] = "BIO.01.MAT" 
# 
# 
# ## Divide columns by 10
# SUA.ZONAL.STATS.CONVERT = as.data.table(SUA.ZONAL.STATS)                           ## Check this works, also inefficient
# SUA.ZONAL.STATS.CONVERT[, (names(SUA.ZONAL.STATS[c(5:10)])) := lapply(.SD, function(x) 
#   x / 10 ), .SDcols = names(SUA.ZONAL.STATS[c(5:10)])]
# SUA.ZONAL.STATS.CONVERT = as.data.frame(SUA.ZONAL.STATS.CONVERT) 
# View(SUA.ZONAL.STATS)
SUA.BIO1.current.stats[["CURRENT_MAT"]] = SUA.BIO1.current.stats[["CURRENT_MAT"]]/10

## Merge MAP onto the results table
SUA.PREDICT = join(SUA.PREDICT, SUA.BIO1.current.stats[c("SUA", "CURRENT_MAT")])
SUA.PREDICT = join(SUA.PREDICT, SUA.BIO12.current.stats[c("SUA", "CURRENT_MAP")])
SUA.PREDICT = join(SUA.PREDICT, SUA.KOP[c("SUA", "MAJOR_KOP")])
SUA.PREDICT = join(SUA.PREDICT, SUA.PET[c("SUA", "CURRENT_PET")])
summary(SUA.PREDICT)
View(SUA.PREDICT)
t = SUA.PREDICT[is.na(SUA.PREDICT$MAJOR_KOP),]

#########################################################################################################################
## Save table
write.csv(SUA.PREDICT, paste0('output/tables/MAXENT_SUA_PRESENCE_KOPPEN', save_run, '.csv'), row.names = FALSE)
write.csv(SUA.TOP.PRESENCE, paste0('output/tables/MAXNET_SUA_TOP_',       save_run, '.csv'), row.names = FALSE)



#########################################################################################################################
## We simply need a count of all the species that are being lost, gained or remaining stable in each SUA
length(unique(SUA.PREDICT$SPECIES))
SUA.PLOT.GOOD = subset(SUA.PREDICT, MAXENT_RATING < 3 & ORIGIN == "Native")
unique(SUA.PLOT.GOOD$MAXENT_RATING)
length(unique(SUA.PLOT.GOOD$SPECIES))

       
SUA.PLOT   = table(SUA.PLOT.GOOD$SUA, SUA.PLOT.GOOD$GAIN_LOSS)
SUA.PLOT.M = melt(SUA.PLOT) 
names(SUA.PLOT.M) = c("SUA", "AREA_CHANGE", "SPECIES_COUNT")
head(SUA.PLOT.M)




#########################################################################################################################
## This will plot all the SUAs.
## Could add a column to SUA for what the mean annual temperature is, to see if there is a pattern
ggplot(SUA.PLOT.M, aes(x = SUA, fill = AREA_CHANGE)) + 
  
  geom_bar(data = subset(SUA.PLOT.M, AREA_CHANGE %in% c("LOSS")),
           aes(y = -SPECIES_COUNT), position = "stack", stat = "identity") +
  
  geom_bar(data = subset(SUA.PLOT.M, !AREA_CHANGE %in% c("LOSS")), 
           aes(y = SPECIES_COUNT), position = "stack", stat = "identity")


#########################################################################################################################
## Try plotting a subset of just the mapped SUAs
unique(SUA.PLOT.M$SUA)
SUA.PLOT.M$SUA = gsub("Canberra - Queanbeyan", "Canberra", SUA.PLOT.M$SUA)
SUA.PLOT.MAJOR = SUA.PLOT.M[SUA.PLOT.M$SUA %in% MAP_SUA, ] 
unique(SUA.PLOT.MAJOR$SUA)


## Change PNG output
png(sprintf('output/figures/SUA_BAR_PLOT/SUA_BAR_PLOT_%s.png', save_run),      
    15, 8, units = 'in', res = 500)

ggplot(SUA.PLOT.MAJOR, aes(x = SUA, fill = AREA_CHANGE)) + 
  
  ## The species being lost
  geom_bar(data = subset(SUA.PLOT.MAJOR, AREA_CHANGE %in% c("LOSS")),
           aes(y = -SPECIES_COUNT), 
           position = "stack", stat = "identity", colour="black") +
  
  ## The species being gained or remaining stable
  geom_bar(data = subset(SUA.PLOT.MAJOR, !AREA_CHANGE %in% c("LOSS")), 
           aes(y = SPECIES_COUNT), 
           position = "stack", stat = "identity", colour="black") +
  
  ## The colour scheme
  scale_fill_manual(values=rev(colorRampPalette(c('skyblue3', 'brown1', 'seagreen3'))(3))) +
  
  ## The axes labels
  labs(title = "Predicted gains and losses within Significant Urban Areas (SUAs)", 
       x = "SUA", y = "Species Count") +

## Format axes
theme(axis.title.x     = element_text(face = "bold", colour = "black", size = 15),
      axis.text.x      = element_text(angle = 45, vjust = 0.5, size = 12),
      axis.title.y     = element_text(face = "bold", colour = "black", size = 15),
      axis.text.y      = element_text(vjust = 0.5, size = 12),
      title            = element_text(face = "bold", colour = "black", size = 20),
      legend.title     = element_blank(),
      legend.text      = element_text(face = "bold", size = 12),
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      panel.border     = element_rect(colour = "black", fill = NA, size = 2))

dev.off()




#########################################################################################################################
## 4). CREATE A RICHNESS MAP FOR EACH TIME PERIOD
#########################################################################################################################


## Probably don't need this anymore - just use ArcMap to sum the rasters................................................
## Search for


#########################################################################################################################
## Look through the output directory for species where the code didn't finish. This is usually species with < 27 files
## Summing 200 rasters will take ages 
raster.dirs <- list.dirs(path = maxent_path, full.names = FALSE, recursive = FALSE)


## Loop over the directories for the best species: current maps
raster.current <- sapply(raster.dirs, function(x) {
  
  ## List the files for that time slice
  list.files(paste0(maxent_path, x), pattern = 'current_suit_above', full.names = TRUE, recursive = TRUE)
  
})



## 2030 maps
raster.2030 <- sapply(raster.dirs, function(x) {
  
  ## List the files for that time slice
  list.files(paste0(maxent_path, x), pattern = '2030_4GCMs_above', full.names = TRUE, recursive = TRUE)
  
})


## 2050 maps
raster.2050 <- sapply(raster.dirs, function(x) {
  
  ## List the files for that time slice
  list.files(paste0(maxent_path, x), pattern = '2050_4GCMs_above', full.names = TRUE, recursive = TRUE)
  
})


## 2070 maps
raster.2070 <- sapply(raster.dirs, function(x) {
  
  ## List the files for that time slice
  list.files(paste0("./output/maxent/SET_VAR_KOPPEN/", x), pattern = '2070_4GCMs_above', full.names = TRUE, recursive = TRUE)
  
})


#########################################################################################################################
## Now unlist so we can create a raster stack
raster.current  = unlist(raster.current)
raster.2030  = unlist(raster.2030)
raster.2050  = unlist(raster.2050)
raster.2070  = unlist(raster.2070)


## Check length
length(raster.current);length(raster.2030);length(raster.2050);length(raster.2070)


## Then create raster stacks and sum
stack.current   = stack(raster.current, quick = TRUE)
stack.2030      = stack(raster.2030, quick = TRUE)
stack.2050      = stack(raster.2050, quick = TRUE)
stack.2070      = stack(raster.2070, quick = TRUE)


## writeRaster(stack.2030, filename = 'output/maxent/multilayer.tif', options = "INTERLEAVE=BAND", overwrite = TRUE)
## mystack = stack("multilayer.tif")


## Summing takes a long time
sum.current  = sum(stack.current) 
sum.2030     = sum(stack.2030) 
sum.2050     = sum(stack.2050) 
sum.2070     = sum(stack.2070) 


#########################################################################################################################
## Plot to check
plot(sum.current)
plot(sum.2030)
plot(sum.2050)
plot(sum.2070)


#########################################################################################################################
## write out rasters
writeRaster(sum.current, 'output/maxent/checked_spp_current_richness.tif')
writeRaster(sum.2030,    'output/maxent/checked_spp_2030_richness.tif')
writeRaster(sum.2050,    'output/maxent/checked_spp_2050_richness.tif')
writeRaster(sum.2070,    'output/maxent/checked_spp_2070_richness.tif')



#########################################################################################################################
## OUTSTANDING PREDICTION TASKS:
#########################################################################################################################


#########################################################################################################################
## 





#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################