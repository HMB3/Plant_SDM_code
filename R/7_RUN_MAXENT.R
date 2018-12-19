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


## Print the species run to the screen
message('Running maxent models for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")


#########################################################################################################################
## Load SDM table
if(read_data == "TRUE") {
  
  ## read in RDS files from previous step
  SDM.SPAT.ALL = readRDS(paste0(DATA_path, 'SDM_SPAT_ALL_', save_run, '.rds'))
  
} else {
  
  message(' skip file reading, not many species analysed')   ##
  
}





#########################################################################################################################
## 1). RUN SDMs USING A-PRIORI VARIABLES FOR ALL SPECIES
#########################################################################################################################


#########################################################################################################################
## Check the SDM table
dim(SDM.SPAT.ALL)
length(intersect(unique(SDM.SPAT.ALL$searchTaxon), GBIF.spp))  ## should be same as the number of species
unique(SDM.SPAT.ALL$SOURCE)


#########################################################################################################################
## Run Maxent using a random selection of background points. 
projection(template.raster);projection(SDM.SPAT.ALL);projection(Koppen_1975)


#########################################################################################################################
## Loop over all the species spp = GBIF.spp[1]
lapply(GBIF.spp, function(spp){ 
  
  ## Skip the species if the directory already exists, before the loop
  outdir <- maxent_dir
  
  dir_name = file.path(maxent_path, gsub(' ', '_', spp))
  if(dir.exists(dir_name)) {
    message('Skipping ', spp, ' - already run.')
    invisible(return(NULL))
    
  }

  #  create the directory so other parallel runs don't try to do it
  dir.create(dir_name)
  write.csv(data.frame(), file.path(dir_name, "in_progress.txt"))
  
  
  ## Print the taxa being processed to screen
  if(spp %in% SDM.SPAT.ALL$searchTaxon) {
    message('Doing ', spp) 
    
    ## Subset the records to only the taxa being processed
    #occurrence <- subset(SDM.SPAT.ALL, searchTaxon == spp)
    occurrence <- subset(SDM.SPAT.ALL, searchTaxon == spp)# & SOURCE != "INVENTORY")
    
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
                         min_n                   = 20,            ## This should be higher...
                         max_bg_size             = 70000,         ## could be 50k or lower, it just depends on the biogeography
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
        write.csv(data.frame(), file.path(dir_name, "failed.txt"))

      })
    
  } else {
    
    message(spp, ' skipped - no data.')         ## This condition ignores species which have no data...
    write.csv(data.frame(), file.path(dir_name, "completed.txt"))
    
  }  

  ## now add a file to the dir to denote that it has completed
  write.csv(data.frame(), file.path(dir_name, "completed.txt"))
  
})





#########################################################################################################################
## OUTSTANDING SDM TASKS:
#########################################################################################################################



## 1). Create list of species with bad models








#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################
