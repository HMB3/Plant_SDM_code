#########################################################################################################################
########################################### CREATE SUMMARY RASTERS ###################################################### 
#########################################################################################################################


#########################################################################################################################
## This code combines the output of the SDM models to estimate species counts under 6 GCMs under current and 2070 climates


## Load packages ::
source('./R/HIA_LIST_MATCHING.R')



#########################################################################################################################
## 1). CREATE SPECIES COUNT RASTERS FOR ALL ANALYSED TAXA
#########################################################################################################################


## Use the species list which has been checked
MAXENT.CHECK = read.csv("./output/maxent/MAXENT_CHECK_RATING.csv", stringsAsFactors = FALSE)
spp.check    = subset(MAXENT.CHECK, CHECK_MAP == 1 | CHECK_MAP == 2)$searchTaxon
spp_check    = gsub(" ", "_", spp.check)
str(spp_check)
  

########################################################################################################################
## Create file list 
DIR               = "./output/maxent/SET_VAR_KOPPEN"
current.list      = list.files(DIR, pattern = 'current_suit_above_', full.names = TRUE, recursive = TRUE)


## restrict the list to only the checked species
current.list.checked <- current.list [c(1:length(current.list))] %>%
  
  ## Pipe the list into lapply
  lapply(function(file) {
    
    ## Create the character string...
    m <- file[file %like% spp_check] 
    #m <-   sprintf('%s%s/full/', path.set.var, species)                ## path.backwards.sel
    m 
    
  }) %>%
  
  ## Bind the list together
  c()


current.list      = current.list [current.list %like% spp_check] 


########################################################################################################################
## Update the
suit              = stack(raster.list)
suit.list         = unstack(suit)
combo_suit_mean   = mean(suit)                            ## plot(mean.suit)

writeRaster(combo_suit_mean , sprintf('%s/%s/full/%s_20%s_suitability_mean.tif', 
                                      maxent_path, species, species, time_slice), overwrite = TRUE)



#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################