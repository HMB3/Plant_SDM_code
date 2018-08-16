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
## Create a list of the gain/loss tables. Just run this on the final SUA table 
GAIN.LOSS.list = list.files(save_dir, pattern = 'gain_loss_table_', full.names = TRUE, recursive = TRUE) 
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
GAIN.LOSS.TABLE          =  GAIN.LOSS.TABLE[GAIN.LOSS.TABLE$SPECIES %in% map_spp_list , ] 
GAIN.LOSS.TABLE.COMPLETE =  GAIN.LOSS.TABLE[complete.cases(GAIN.LOSS.TABLE), ]
summary(GAIN.LOSS.TABLE.COMPLETE)


## How big is the table?
dim(GAIN.LOSS.TABLE)
dim(GAIN.LOSS.TABLE.COMPLETE)
head(GAIN.LOSS.TABLE.COMPLETE, 10)
head(GAIN.LOSS.TABLE.COMPLETE, 10)


## How many species are lost?
length(unique(GAIN.LOSS.TABLE$SPECIES))                                   
length(unique(GAIN.LOSS.TABLE.COMPLETE$SPECIES))  


#########################################################################################################################
## Now write CSV 
write.csv(GAIN.LOSS.TABLE,  "./output/tables/OVERALL_GAIN_LOSS_TABLE.csv",   row.names = FALSE)





#########################################################################################################################
## 2). COMBINE TABLES OF SPECIES PRESENCES IN SUAs
#########################################################################################################################


#########################################################################################################################
## The multiple thresholds could present a problem
SUA.tables = list.files(save_dir, pattern = 'area_SUA_summary_', full.names = TRUE, recursive = TRUE) 
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
head(SUA.PRESENCE, 50);tail(SUA.PRESENCE, 50)


#########################################################################################################################
## Now join on the population. Note that the area shapefile does not have the same SUAs as the ABS table
SUA.DENS = areal_unit@data
SUA.DENS = SUA.DENS[, c("SUA_NAME11", "AREA_SQKM")]
names(SUA.DENS) = c("SUA", "AREA_SQKM")
head(SUA.DENS)


##
TOP.SUA.POP             = ALL.SUA.POP[, c("SUA", "POP_2017")]
SUA.DEN                 = merge(TOP.SUA.POP, SUA.DENS)
SUA.DEN$POP_DESNITY     = TOP.SUA.DEN$POP_2017/TOP.SUA.DEN$AREA_SQKM
View(SUA.DEN)


## 
SUA.PRESENCE$SUA = as.character(SUA.PRESENCE$SUA)
class(SUA.PRESENCE$SUA)
class(SUA.DEN$SUA)


## 87 SUAs overlap between the ABS shapefile and the table
unique(sort(SUA.DEN$SUA))
unique(sort(SUA.PRESENCE$SUA))
intersect(unique(sort(SUA.DEN$SUA)), unique(sort(SUA.PRESENCE$SUA)))
SUA.PRESENCE = join(SUA.PRESENCE, SUA.DEN)


## Check this combined table
SUA.PRESENCE = SUA.PRESENCE[, c("SUA",         "POP_2017",     "AREA_SQKM", "POP_DESNITY",
                                "SPECIES",     "PERIOD",       "AREA_THRESH", 
                                "MAX_TRAIN",   "CURRENT_AREA", "FUTURE_AREA", 
                                "AREA_CHANGE", "PRESENT",      "GAIN_LOSS")]
names(SUA.PRESENCE)
summary(SUA.PRESENCE)


## Are there are still NA results for some taxa?........................................................................
SUA.COMPLETE = completeFun(SUA.PRESENCE, "CURRENT_AREA")
SUA.COMPLETE = completeFun(SUA.COMPLETE, "FUTURE_AREA")


## Which species are missing?
length(unique(SUA.COMPLETE$SPECIES))
setdiff(unique(SUA.PRESENCE$SPECIES), unique(SUA.COMPLETE$SPECIES))
summary(SUA.COMPLETE$POP_2017)


#########################################################################################################################
## Find the species with the greatest increase/decrease?
summary(SUA.COMPLETE$AREA_CHANGE)


## Increasing
INCREASE.50 = subset(SUA.COMPLETE, AREA_CHANGE >= 50 & PERIOD == 50)
summary(INCREASE.50$AREA_CHANGE)
unique(INCREASE.50$SPECIES)


## Decreasing
DECREASE.50 = subset(SUA.COMPLETE, AREA_CHANGE <= -50 & PERIOD == 50)
summary(DECREASE.50$AREA_CHANGE)
unique(DECREASE.50$SPECIES)


#########################################################################################################################
## Rank by population
## SUA.PRESENCE = SUA.PRESENCE [with(SUA.PRESENCE , rev(order(POP_2017))), ]


## Use a numerical definition of the most populated areas, OR, just take the capital cities
## MAIN_SUA = c("Sydney", "Melbourne", "Canberra - Queanbeyan", "Brisbane", "Perth", "Adelaide", "Hobart", "Darwin")
top_n   = 15
BIG.SUA = head(TOP.SUA.POP[with(TOP.SUA.POP, rev(order(POP_2017))), ], top_n)
BIG_SUA = BIG.SUA$SUA


## Restrict the big table to just the largest SUAs
SUA.TOP.PRESENCE  = SUA.PRESENCE[SUA.PRESENCE$SUA %in% BIG_SUA, ] 
SUA.TOP.PRESENCE  = SUA.TOP.PRESENCE [with(SUA.TOP.PRESENCE , rev(order(POP_2017))), ]
message(length(unique(SUA.PRESENCE$SPECIES)), " Species analysed in ", length(BIG_SUA), " SUAs")
summary(SUA.TOP.PRESENCE)
View(SUA.TOP.PRESENCE)


#########################################################################################################################
## Save basic results and SUA results to file
write.csv(SUA.COMPLETE,     "./output/tables/MAXENT_SUA_PRESENCE.csv",    row.names = FALSE)
write.csv(SUA.TOP.PRESENCE, "./output/maxent/MAXNET_SUA_TOP.csv",         row.names = FALSE)





#########################################################################################################################
## 3). CREATE PLOTS OF SPECIES LOSS AND GAIN
#########################################################################################################################


## How do you convert this table into a format where you can create a bar graph?
length(unique(SUA.COMPLETE$SPECIES))




#########################################################################################################################
## 4). CREATE A RICHNESS MAP FOR EACH TIME PERIOD
#########################################################################################################################


## Probably don't need this anymore - just use ArcMap to sum the rasters................................................
## Search for


#########################################################################################################################
## Look through the output directory for species where the code didn't finish. This is usually species with < 27 files
## Summing 200 rasters will take ages 
raster.dirs <- list.dirs(path = save_dir, full.names = FALSE, recursive = FALSE)


## Loop over the directories for the best species: current maps
raster.current <- sapply(raster.dirs, function(x) {
  
  ## List the files for that time slice
  list.files(paste0(save_dir, x), pattern = 'current_suit_above', full.names = TRUE, recursive = TRUE)
  
})



## 2030 maps
raster.2030 <- sapply(raster.dirs, function(x) {
  
  ## List the files for that time slice
  list.files(paste0(save_dir, x), pattern = '2030_4GCMs_above', full.names = TRUE, recursive = TRUE)
  
})


## 2050 maps
raster.2050 <- sapply(raster.dirs, function(x) {
  
  ## List the files for that time slice
  list.files(paste0(save_dir, x), pattern = '2050_4GCMs_above', full.names = TRUE, recursive = TRUE)
  
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