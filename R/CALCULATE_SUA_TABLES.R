#########################################################################################################################
## 5). SUMARIZE MAXENT RESULTS FOR EACH SPECIES ACROSS MULTIPLE GCMs INSIDE EACH SUA
#########################################################################################################################


#########################################################################################################################
## So we can use the values in these columns to threhsold the rasters of habitat suitability (0-1) when combining them.
## For each species, we will create a binary raster with cell values between 0-6. These cell values represent the number of GCMs 
## where that cell had a suitability value above the threshold determined by maxent (e.g. the max training sensitivity + specif).


## To run the function inside the loop
DIR        = SDM.RESULTS.DIR[59]
species    = map_spp_list[59]
thresh     = percent.10.log[59]
percent    = percent.10.om[59]
time_slice = 50
area_occ   = 10


#########################################################################################################################
## Reverse sort the lists
SDM.DIR.REV     = sort(SDM.RESULTS.DIR, decreasing = TRUE)
map_spp_rev     = sort(map_spp,         decreasing = TRUE) 
percent.log.rev = percent.10.log[sort(order(percent.10.log), decreasing = TRUE)]
percent.om.rev  = percent.10.om[sort(order(percent.10.om),   decreasing = TRUE)]


#load("SDM_COMBINE.RData")
load("SDM_COMBINE..RData")

## Check the length matches - order should be correct, as well as the length
length(SDM.RESULTS.DIR);length(map_spp);length(percent.10.log);length(percent.10.om)
length(SDM.DIR.REV);length(map_spp_rev);length(percent.log.rev);length(percent.om.rev)



## The function is getting stuck on the first list item. What could be causing the problem?
## Changes to the the function/loop
## Problems in the directory structure?



#########################################################################################################################
## Combine output and calculate gain and loss for 2030
suitability.2030 = tryCatch(mapply(combine_gcm_threshold,
                                   DIR_list     = SDM.RESULTS.DIR,
                                   species_list = map_spp,
                                   maxent_path  = maxent_path,
                                   thresholds   = percent.10.log,
                                   percentiles  = percent.10.om,
                                   time_slice   = 30,
                                   area_occ     = 10),
                            
                            error = function(cond) {
                              
                              message(paste('Species skipped - check inputs', spp))
                              
                            })


## Combine GCM output for 2050 
suitability.2050 = tryCatch(mapply(combine_gcm_threshold, 
                                   DIR_list     = SDM.RESULTS.DIR,
                                   species_list = map_spp,
                                   maxent_path  = maxent_path,
                                   thresholds   = percent.10.log,
                                   percentiles  = percent.10.om,
                                   time_slice   = 50,
                                   area_occ     = 10),
                            
                            error = function(cond) {
                              
                              message(paste('Species skipped - check inputs', spp))
                              
                            })


## Combine GCM output for 2070 
suitability.2070 = tryCatch(mapply(combine_gcm_threshold, 
                                   DIR_list     = SDM.RESULTS.DIR,
                                   species_list = map_spp,
                                   maxent_path  = maxent_path,
                                   thresholds   = percent.10.log,
                                   percentiles  = percent.10.om,
                                   time_slice   = 70,
                                   area_occ     = 10),
                            
                            error = function(cond) {
                              
                              message(paste('Species skipped - check inputs', spp))
                              
                            })


#########################################################################################################################
## Reverse the order of the lists to speed up the output
mapply(combine_gcm_threshold,
       DIR_list     = SDM.DIR.REV,
       species_list = map_spp_rev,
       maxent_path  = maxent_path,
       thresholds   = percent.log.rev,
       percentiles  = percent.om.rev,
       time_slice   = 30,  ## 50, 70
       area_occ     = 10)


mapply(SUA_tables,
       DIR_list     = SDM.RESULTS.DIR,
       species_list = map_spp,
       maxent_path  = maxent_path,
       thresholds   = percent.10.log,
       time_slice   = 30,
       area_occ     = 10)


#########################################################################################################################
## OUTSTANDING MAPPING TASKS:
#########################################################################################################################

## intersect(list.files(maxent_path), map_spp_list) 


## 1). Fix species that didn't work on the last run, and run the remaining species

## 2). Then combine species into one big table - should be 238 * 3 time periods - a list of 714 tables

## 3). Create threshold maps for species rated 2 and 3 - could use john's code to summarise the outputted data

## 4). Create mess maps for species rated 2 and 3 - this takes a long time

## 5). Follow up with ALA RE field names 

## 6). Create protcol for species cleaning - initial clean using just the png files, then threshold and mess maps for the 1s
##     and 2s.







#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################