#########################################################################################################################
################################# SUMMARISE THE HABITAT SUITABILITY RESULTS ############################################# 
#########################################################################################################################


## This code combines the tables of area occupied by each species in each SUA, and converts it to one big table. 
## This table can be queried to create histograms, etc.  


#########################################################################################################################
## Read in all data to run the SDM code :: species lists, shapefile, rasters & tables
#source('./R/HIA_LIST_MATCHING.R')
rasterTmpFile()


## Read in SUAs and create species lists
ALL.SUA.POP       = read.csv("./data/base/CONTEXTUAL/ABS_SUA_POP.csv",    stringsAsFactors = FALSE)
MAXENT.CHECK      = read.csv("./output/maxent/MAXENT_RATING_26_2018.csv", stringsAsFactors = FALSE)

spp.lower.thresh  = subset(MAXENT.CHECK, CHECK_MAP == 2 | CHECK_MAP == 3)$searchTaxon
spp_lower_thresh  = gsub(" ", "_", spp.lower.thresh)

spp.best.thresh   = subset(MAXENT.CHECK, CHECK_MAP == 1 | CHECK_MAP == 2)$searchTaxon
spp_best_thresh   = gsub(" ", "_", spp.best.thresh)





#########################################################################################################################
## 1). COMBINE OVERALL TABLES OF SPECIES GAIN/LOSS ACROSS AUSTRALIA
#########################################################################################################################


## We only want to add the species used for the analysis, not all species folders. Can we restrict the list to just the 
## species on map_list?


## Might not need this code as much, if just summing the rasters in ArcMap..............................................


#########################################################################################################################
## Create a list of the gain/loss tables 
GAIN.LOSS.list = list.files("./output/maxent/SET_VAR_KOPPEN/", pattern = 'gain_loss_table_', full.names = TRUE, recursive = TRUE) 
length(GAIN.LOSS.list)


## Now combine the SUA tables for each species into one table 
GAIN.LOSS.STABLE <- GAIN.LOSS.list %>%
  
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


## This is a summary of maxent output for current conditions :: some of the species did not process properly
summary(GAIN.LOSS.STABLE)
GAIN.LOSS.STABLE[!complete.cases(GAIN.LOSS.STABLE), ]


dim(GAIN.LOSS.STABLE)
head(GAIN.LOSS.STABLE, 30)
length(unique(GAIN.LOSS.STABLE$SPECIES))                                  ## 143 species rated as either 1 or 2 by Linda


#########################################################################################################################
## Get the species with the highest gain
GAIN.OVERALL   = subset(GAIN.LOSS.STABLE, CHANGE == "GAINED") 
GAIN.OVERALL   = GAIN.OVERALL[with(GAIN.OVERALL, rev(order(COUNT))), ]

LOSS.OVERALL   = subset(GAIN.LOSS.STABLE, CHANGE == "LOST")
LOSS.OVERALL   = LOSS.OVERALL[with(LOSS.OVERALL, rev(order(COUNT))), ]


## Have a look
head(GAIN.OVERALL, 50)
head(LOSS.OVERALL, 50)
GAIN.OVERALL[!complete.cases(GAIN.OVERALL), ]
LOSS.OVERALL[!complete.cases(LOSS.OVERALL), ]


## Break it down
summary(GAIN.OVERALL$COUNT)
summary(LOSS.OVERALL$COUNT)


#########################################################################################################################
## Now write CSV 
write.csv(GAIN.LOSS.STABLE, "./output/tables/OVERALL_GAIN_LOSS_STABLE.csv", row.names = FALSE)
write.csv(GAIN.OVERALL,     "./output/tables/OVERALL_GAIN.csv",             row.names = FALSE)
write.csv(LOSS.OVERALL,     "./output/tables/OVERALL_LOSS.csv",             row.names = FALSE)





#########################################################################################################################
## 2). COMBINE TABLES OF SPECIES PRESENCES IN SUAs
#########################################################################################################################


#########################################################################################################################
## The multiple thresholds could present a problem.
SUA.tables = list.files("./output/maxent/SET_VAR_KOPPEN/", pattern = 'area_SUA_summary_', full.names = TRUE, recursive = TRUE) 
length(SUA.tables)


## Now combine the SUA tables for each species into one table 
SUA.PRESENCE <- SUA.tables[c(1:length(SUA.tables))] %>%
  
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


## This is a summary of maxent output for current conditions
dim(SUA.PRESENCE)
head(SUA.PRESENCE, 200)


#########################################################################################################################
## Now join on the population
TOP.SUA.POP      = ALL.SUA.POP[, c("SUA", "POP_2017")]
SUA.PRESENCE$SUA = as.character(SUA.PRESENCE$SUA)
class(SUA.PRESENCE$SUA)
class(TOP.SUA.POP$SUA)


## 87 SUAs overlap between the ABS shapefile and the table
unique(sort(TOP.SUA.POP$SUA))
unique(sort(SUA.PRESENCE$SUA))
intersect(unique(sort(TOP.SUA.POP$SUA)), unique(sort(SUA.PRESENCE$SUA)))
SUA.PRESENCE = join(SUA.PRESENCE, TOP.SUA.POP)


## Join on the Maxent rating
RATING = MAXENT.CHECK
names(RATING)  = c("SPECIES", "RATING")
RATING$SPECIES = gsub(" ", "_", RATING$SPECIES)
length(intersect(RATING$SPECIES, SUA.PRESENCE$SPECIES))
SUA.PRESENCE = join(SUA.PRESENCE, RATING, type = "inner")


## Check this combined table
SUA.PRESENCE = SUA.PRESENCE[, c("SUA",         "POP_2017",     "AREA_SQKM", 
                                "SPECIES",     "PERIOD",       "AREA_THRESH", 
                                "MAX_TRAIN",   "CURRENT_AREA", "FUTURE_AREA", 
                                "AREA_CHANGE", "PRESENT",      "GAIN_LOSS", "RATING")]
names(SUA.PRESENCE)
summary(SUA.PRESENCE)
unique(SUA.PRESENCE$RATING)


## There are still NA results for some taxa
SUA.COMPLETE = completeFun(SUA.PRESENCE, "CURRENT_AREA")
SUA.COMPLETE = completeFun(SUA.COMPLETE, "FUTURE_AREA")
length(unique(SUA.COMPLETE$SPECIES))
names(SUA.COMPLETE)




#########################################################################################################################
## Find the species with the greatest increase
summary(SUA.COMPLETE$AREA_CHANGE)
INCREASE.50 = subset(SUA.COMPLETE, AREA_CHANGE >= 50)
summary(INCREASE.50$AREA_CHANGE)
unique(INCREASE.50$SPECIES)


#########################################################################################################################
## Rank by population
## SUA.PRESENCE = SUA.PRESENCE [with(SUA.PRESENCE , rev(order(POP_2017))), ]


## what are the 20 most populated areas
top_n   = 5
BIG.SUA = head(TOP.SUA.POP[with(TOP.SUA.POP, rev(order(POP_2017))), ], top_n)
BIG_SUA = BIG.SUA$SUA


## Restrict the big table to just these
SUA.TOP.PRESENCE  = SUA.PRESENCE[SUA.PRESENCE$SUA %in% BIG_SUA, ] 
SUA.TOP.PRESENCE  = SUA.TOP.PRESENCE [with(SUA.TOP.PRESENCE , rev(order(POP_2017))), ]
length(unique(SUA.PRESENCE$SPECIES))
length(unique(SUA.TOP.PRESENCE$SPECIES))                            ## So their are 195 species processed
summary(SUA.TOP.PRESENCE)
View(SUA.TOP.PRESENCE)


## Could subset again
SUA.TOP.PRESENCE.2030      = subset(SUA.TOP.PRESENCE, PERIOD == 30)
SUA.TOP.PRESENCE.2050      = subset(SUA.TOP.PRESENCE, PERIOD == 50)
SUA.TOP.PRESENCE.2070      = subset(SUA.TOP.PRESENCE, PERIOD == 70)

SUA.TOP.PRESENCE.2030.MILE = SUA.TOP.PRESENCE.2030[SUA.TOP.PRESENCE.2030$SPECIES %in% spp_mile, ]
SUA.TOP.PRESENCE.2070.MILE = SUA.TOP.PRESENCE.2030[SUA.TOP.PRESENCE.2070$SPECIES %in% spp_mile, ]


## Check
head(SUA.TOP.PRESENCE.2030.MILE)
unique(SUA.TOP.PRESENCE.2030.MILE$SPECIES)


## Just one species
SUA.TOP.PRESENCE.2030.MILE.X.GLAUCA = subset(SUA.TOP.PRESENCE, SPECIES == "Xanthorrhoea_glauca")




#########################################################################################################################
## Save basic results and SUA results to file :: CSV and RData files
write.csv(SUA.COMPLETE,                          "./output/tables/MAXENT_SUA_SUMMARY.csv",     row.names = FALSE)
#write.csv(SUA.TOP.PRESENCE,                     "./output/tables/MAXENT_SUA_SUMMARY.csv",     row.names = FALSE)

write.csv(MAXENT.SUMMARY,       "./output/maxent/MAXENT_STD_VAR_SUMMARY.csv", row.names = FALSE)
write.csv(SUA.TOP.PRESENCE,     "./output/maxent/MAXENT_SUA_SUMMARY.csv",     row.names = FALSE)





#########################################################################################################################
## 3). CREATE A RICHNESS MAP FOR EACH TIME PERIOD
#########################################################################################################################


## Create a list of rasters 
## raster.list[raster.list %like% "Euc"]


#########################################################################################################################
## Look through the output directory for species where the code didn't finish. This is usually species with < 27 files
## in the top "species_name" directory
best.dirs <- list.dirs(path = "./output/maxent/SET_VAR_KOPPEN/", full.names = FALSE, recursive = FALSE)
best.dirs = best.dirs [best.dirs %in% spp_best_thresh]


## Loop over the directories for the best species: current maps
raster.current <- sapply(best.dirs, function(x) {
  
  ## List the files for that time slice
  list.files(paste0("./output/maxent/SET_VAR_KOPPEN/", x), pattern = 'current_suit_above', full.names = TRUE, recursive = TRUE)
  
})



## 2030 maps
raster.2030 <- sapply(best.dirs, function(x) {
  
  ## List the files for that time slice
  list.files(paste0("./output/maxent/SET_VAR_KOPPEN/", x), pattern = '2030_4GCMs_above', full.names = TRUE, recursive = TRUE)
  
})


## 2050 maps
raster.2050 <- sapply(best.dirs, function(x) {
  
  ## List the files for that time slice
  list.files(paste0("./output/maxent/SET_VAR_KOPPEN/", x), pattern = '2050_4GCMs_above', full.names = TRUE, recursive = TRUE)
  
})


## 2070 maps
raster.2070 <- sapply(best.dirs, function(x) {
  
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
## create a list of species which need to be mapped using a different threshold. Either the ::
## 10 percentile training presence Logistic threshold OR
## 10 percentile training presence training omission

## Also Can we calculate these in the table?
## Gain (ie not suitable now but suitable in future)
## Loss (ie suitable now but not in future)
## Stable (suitable in both)

## How should we map each species? 
## Subtract the current raster from the future combined raster





#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################