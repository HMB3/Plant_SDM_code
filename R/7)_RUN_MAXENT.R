#########################################################################################################################
#############################################  FIT MAXENT TO GBIF DATA ################################################## 
#########################################################################################################################


#########################################################################################################################
## This code takes a table of all species occurrences (rows) and environmental values (columns), and runs SDMs for a list
## of taxa.


## The HIA brief again:

## The first module will focus on x plant species identified in the project’s Target Species List, and will develop maps 
## that demonstrate each species’ suitability to both current and future climates across Australia.

## These maps will be used to demonstrate how well or poorly a particular species will be able to tolerate future conditions 
## in urban centres across Australia as the climate changes, based on our current understanding of species’ climatic 
## requirements.



#########################################################################################################################
## Load packages, functions and data
SDM.SPAT.ALL = readRDS('data/base/HIA_LIST/COMBO/SDM_SPAT_ALL_OLD_ALA.rds')
SDM.SPAT.TEST = subset(SDM.SPAT.ALL, searchTaxon == "Acacia baileyana" | 
                         searchTaxon == "Acacia boormanii" |
                         searchTaxon == "Acacia cognata" )
unique(SDM.SPAT.TEST$searchTaxon)
unique(SDM.SPAT.TEST$SOURCE)
saveRDS(SDM.SPAT.TEST, 'data/base/HIA_LIST/COMBO/SDM_SPAT_ALL_OLD_ALA_TEST.rds')
dim(subset(SDM.SPAT.ALL, searchTaxon == "Quercus robur")) 
rasterTmpFile()





#########################################################################################################################
## 1). RUN SDMs USING A-PRIORI VARIABLES FOR ALL SPECIES
#########################################################################################################################


#########################################################################################################################
## Check the SDM table
identical(names(env.grids.current),sdm.predictors)
dim(SDM.SPAT.ALL)
length(unique(SDM.SPAT.ALL$searchTaxon))
length(unique(SDM.SPAT.ALL$OBS))
unique(SDM.SPAT.ALL$SOURCE)
## save.image('SDM_ANALYSES.RData')


#########################################################################################################################
## Run Maxent using a random selection of background points. 
length(unique(SDM.SPAT.ALL$searchTaxon))
projection(template.raster);projection(SDM.SPAT.ALL);projection(Koppen_1975)
GBIF.spp[28]


#########################################################################################################################
## Loop over all the species spp = GBIF.spp[1]
lapply(GBIF.spp, function(spp){ 
  
  ## Skip the species if the directory already exists, before the loop
  outdir <- out_dir
  
  if(dir.exists(file.path(outdir, gsub(' ', '_', spp)))) {
    message('Skipping ', spp, ' - already run.')
    invisible(return(NULL))
    
  }
  
  ## Print the taxa being processed to screen
  if(spp %in% SDM.SPAT.ALL$searchTaxon) {
    message('Doing ', spp) 
    
    ## Subset the records to only the taxa being processed
    occurrence <- subset(SDM.SPAT.ALL, searchTaxon == spp)
    
    ## Now get the background points. These can come from any spp, other than the modelled species.
    background <- subset(SDM.SPAT.ALL, searchTaxon != spp)
    
    ## Finally fit the models using FIT_MAXENT_TARG_BG. Also use tryCatch to skip any exceptions
    tryCatch(
      FIT_MAXENT_TARG_BG(occ                     = occurrence, 
                         bg                      = background, 
                         sdm.predictors          = sdm.select, 
                         name                    = spp, 
                         outdir, 
                         template.raster,
                         min_n                   = 20,     ## This should be higher...
                         max_bg_size             = 100000, 
                         Koppen                  = Koppen_1975,
                         background_buffer_width = 200000,
                         shapefiles              = TRUE,
                         features                = 'lpq',
                         replicates              = 5,
                         responsecurves          = TRUE),
      
      ## https://stackoverflow.com/questions/19394886/trycatch-in-r-not-working-properly
      #function(e) message('Species skipped ', spp)) ## skip any species for which the function fails
      error = function(cond) {
        
        message(paste('Species skipped ', spp))
        
      })
    
  } else {
    
    message(spp, ' skipped - no data.')         ## This condition ignores species which have no data...
    
  }  
  
})





#########################################################################################################################
## 2). RUN SDMs USING A-PRIORI VARIABLES FOR THE BIASED SPECIES
#########################################################################################################################


#########################################################################################################################
## Run Maxent using a random selection of background points. 
# SDM.DATA.BIAS = readRDS('./data/base/HIA_LIST/COMBO/BIASED_TREES_RAREFY.rds')
# length(unique(SDM.DATA.BIAS$searchTaxon)) 
# dim(SDM.DATA.BIAS)
# projection(template.raster);projection(SDM.DATA.BIAS);projection(Koppen_1975)
# 
# 
# ## Remove the dodgy columns
# dropList <- c(setdiff(sdm.predictors, sdm.select), "lon", "lat")
# background.bias <- background[, !names(background) %in% dropList]
# SDM.DATA.BIAS   <- SDM.DATA.BIAS[, !names(SDM.DATA.BIAS) %in% dropList]
# identical(names(background.bias), names(SDM.DATA.BIAS))
# 
# 
# ## Now bind on the background points............................................................................................
# intersect(sort(unique(SDM.DATA.BIAS$searchTaxon)), sort(unique(background$searchTaxon)))    ## check the species don't overlap
# SDM.DATA.BIAS = rbind(SDM.DATA.BIAS, background.bias)
# names(SDM.DATA.BIAS)
# 
# 
# ## Change the output directory
# out_dir       = 'output/maxent/SPP_RAREFY_3KM'
# 
# 
# #########################################################################################################################
# ## Loop over all the species spp = test.bias[1]
# lapply(SPP.BIAS, function(spp){
# 
#   ## Skip the species if the directory already exists, before the loop
#   outdir <- out_dir
#   if(dir.exists(file.path(outdir, gsub(' ', '_', spp)))) {
#     message('Skipping ', spp, ' - already run.')
#     invisible(return(NULL))
# 
#   }
# 
#   ## Print the taxa being processed to screen
#   if(spp %in% SDM.DATA.BIAS$searchTaxon) {
#     message('Doing ', spp)
# 
#     ## Subset the records to only the taxa being processed
#     occurrence <- subset(SDM.DATA.BIAS, searchTaxon == spp)
# 
#     ## Now get the background points. These can come from any spp, other than the modelled species.
#     background <- subset(SDM.DATA.BIAS, searchTaxon != spp)
# 
#     ## Finally fit the models using FIT_MAXENT_TARG_BG. Also use tryCatch to skip any exceptions
#     tryCatch(
#       FIT_MAXENT_TARG_BG(occ                     = occurrence,
#                          bg                      = background,
#                          sdm.predictors          = sdm.select,
#                          name                    = spp,
#                          outdir,
#                          template.raster,
#                          min_n                   = 20,     ## This should be higher...
#                          max_bg_size             = 100000,
#                          Koppen                  = Koppen_1975,
#                          background_buffer_width = 200000,
#                          shapefiles              = TRUE,
#                          features                = 'lpq',
#                          replicates              = 5,
#                          responsecurves          = TRUE),
# 
#       ## https://stackoverflow.com/questions/19394886/trycatch-in-r-not-working-properly
#       #function(e) message('Species skipped ', spp)) ## skip any species for which the function fails
#       error = function(cond) {
# 
#         message(paste('Species skipped ', spp))
# 
#       })
# 
#   } else {
# 
#     message(spp, ' skipped - no data.')         ## This condition ignores species which have no data...
# 
#   }
# 
# })





#########################################################################################################################
## OUTSTANDING SDM TASKS:
#########################################################################################################################


## Which of these components can be done now?
## The last steps combining the models would be handy to have finished


## 1). 400 tree spp that are commonly planted and traded, with sufficient spatial data to model robustly: done


## 2). Fix the taxonomy

##     Follow steps in file 3). 

##     Keep the 'source' column in the maxent table :: adjust step 7  

     

## 4). Run modelling and mapping steps through Katana. To do this, we need a loop that processes one species at a time
##     This means processing the data up to step 7, then running steps 7 and 8 through Katana 
     

## 4). Use more forgiving thresholds (10%) for all species: Done
##     "X10.percentile.training.presence.Logistic.threshold"
##     "X10.percentile.training.presence.training.omission"
  

## 5). Combine output: Table
##     "searchTaxon" "No.plantings" "Number.of.growers" "Number.of.States/LGAs/Koppen zones" "Origin" "AUC" "TSS" 

##     Maps: unique raster values for (search in explorer)
##     "current_suit_above"
##     "2030_4GCMs_above"
##     "2050_4GCMs_above"
##     "2070_4GCMs_above"

   
## 6). Check combined maps for all species : make sure step 9 is working. 
##     Create a spreadsheet of all species results. Check as before - 

##     summary statistics, 
##     input points, 
##     current raster,
##     Future suitability


## 7). make sure step 9 is working. Create tables in R to summarise gain and loss in each SUA


## 8). What is the story? Draft results and discussion.


##     If needed, thin records for 28 spp. with boundary bias, using the SDM tool box. Send Alessandro the latest file .shp

##     Settings: Max 5km, min 1km, 5 heterogeneity classes (from a classification of a PCA, using the worldclim layers)
##     Currently using all spp background records: random 100k for every species, sounds ok?

##     Try to thin the background records in the same way - run the BG points through SDM toolbox too. 

##     Compare the models with and without thinning for a 28 species:
##     Acacia_implexa, etc.





#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################