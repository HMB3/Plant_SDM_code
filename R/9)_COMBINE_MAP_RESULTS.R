#########################################################################################################################
################################# SUMMARISE THE HABITAT SUITABILITY RESULTS ############################################# 
#########################################################################################################################


## This code combines the tables of area occupied by each species in each SUA, and converts it to one big table. 
## This table can be queried to create histograms, etc.  


#########################################################################################################################
## 1). COMBINE OVERALL TABLES OF SPECIES GAIN/LOSS ACROSS AUSTRALIA
#########################################################################################################################


## We only want to add the species used for the analysis, not all species folders. The list can be restricted after 
## raster.list[raster.list %like% "Euc"]


#########################################################################################################################
## Create a list of the gain/loss tables :: this is all species... 
GAIN.LOSS.list = list.files("./output/maxent/SET_VAR_KOPPEN/", pattern = 'gain_loss_table_', full.names = TRUE, recursive = TRUE) 
length(GAIN.LOSS.list)


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
GAIN.LOSS.TABLE  =  GAIN.LOSS.TABLE[GAIN.LOSS.TABLE$SPECIES %in% map_spp_list , ] 


GAIN.LOSS.TABLE[!complete.cases(GAIN.LOSS.TABLE), ]
summary(GAIN.LOSS.TABLE)

dim(GAIN.LOSS.TABLE)
head(GAIN.LOSS.TABLE, 10)
length(unique(GAIN.LOSS.TABLE$SPECIES))                                   ## x species



#########################################################################################################################
## Get the species with the highest gain
GAIN.OVERALL   = subset(GAIN.LOSS.TABLE, CHANGE == "GAINED") 
GAIN.OVERALL   = GAIN.OVERALL[with(GAIN.OVERALL, rev(order(COUNT))), ]

LOSS.OVERALL   = subset(GAIN.LOSS.TABLE, CHANGE == "LOST")
LOSS.OVERALL   = LOSS.OVERALL[with(LOSS.OVERALL, rev(order(COUNT))), ]


## Have a look
head(GAIN.OVERALL)
head(LOSS.OVERALL)
GAIN.OVERALL[!complete.cases(GAIN.OVERALL), ]
LOSS.OVERALL[!complete.cases(LOSS.OVERALL), ]


## Break it down
summary(GAIN.OVERALL$COUNT)
summary(LOSS.OVERALL$COUNT)


#########################################################################################################################
## Now write CSV 
write.csv(GAIN.LOSS.TABLE,  "./output/tables/OVERALL_GAIN_LOSS_TABLE.csv",   row.names = FALSE)
write.csv(GAIN.OVERALL,     "./output/tables/OVERALL_GAIN.csv",              row.names = FALSE)
write.csv(LOSS.OVERALL,     "./output/tables/OVERALL_LOSS.csv",              row.names = FALSE)





#########################################################################################################################
## 2). COMBINE TABLES OF SPECIES PRESENCES IN SUAs
#########################################################################################################################


#########################################################################################################################
## The multiple thresholds could present a problem
SUA.tables = list.files("./output/maxent/SET_VAR_KOPPEN/", pattern = 'area_SUA_summary_', full.names = TRUE, recursive = TRUE) 
length(SUA.tables)


## Now combine the SUA tables for each species into one table 
SUA.PRESENCE <- SUA.tables[c(1:length(SUA.tables))] %>%
  
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
head(SUA.PRESENCE, 50)


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


## Probably don't need the Maxent rating anyore, just use the species which have already been checked.....................
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


## Are there are still NA results for some taxa?........................................................................
SUA.COMPLETE = completeFun(SUA.PRESENCE, "CURRENT_AREA")
SUA.COMPLETE = completeFun(SUA.COMPLETE, "FUTURE_AREA")
length(unique(SUA.COMPLETE$SPECIES))
summary(SUA.COMPLETE$POP_2017)


#########################################################################################################################
## Find the species with the greatest increase?
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


## Restrict the big table to just the largest SUAs
SUA.TOP.PRESENCE  = SUA.PRESENCE[SUA.PRESENCE$SUA %in% BIG_SUA, ] 
SUA.TOP.PRESENCE  = SUA.TOP.PRESENCE [with(SUA.TOP.PRESENCE , rev(order(POP_2017))), ]
length(unique(SUA.PRESENCE$SPECIES))
length(unique(SUA.TOP.PRESENCE$SPECIES))                            ## So their are 195 species processed
summary(SUA.TOP.PRESENCE)
View(SUA.TOP.PRESENCE)


#########################################################################################################################
## Save basic results and SUA results to file
write.csv(SUA.COMPLETE,     "./output/tables/MAXENT_SUA_PRESENCE.csv",    row.names = FALSE)
write.csv(MAXENT.SUMMARY,   "./output/maxent/MAXENT_STD_VAR_SUMMARY.csv", row.names = FALSE)
write.csv(SUA.TOP.PRESENCE, "./output/maxent/MAXNET_SUA_TOP.csv",         row.names = FALSE)





#########################################################################################################################
## 3). CREATE A RICHNESS MAP FOR EACH TIME PERIOD
#########################################################################################################################


## Probably don't need this anymore - just use ArcMap to sum the rasters


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
## TABLE (suitable in both)

## How should we map each species? 
## Subtract the current raster from the future combined raster





#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################