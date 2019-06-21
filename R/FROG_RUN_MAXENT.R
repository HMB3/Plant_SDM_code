#########################################################################################################################
#############################################  FIT MAXENT TO GBIF DATA ################################################## 
#########################################################################################################################


#########################################################################################################################
## This code takes a table of all species occurrences (rows) and environmental values (columns), and runs SDMs for a list
## of taxa.


## For the frogs, we want to impose a few condtions :

## About 6 species for this work, but perhaps more in future. Management applications similar to the Koala study.

## For some taxa you could model at Genus level to bump up the records

## Occurrence data especially gathered, already cleaned - get as CSV

## Could use background points from other frogs, reptiles, etc (I've got occ dat for all the frogs)
                                                             
## Rasters at 30m, e.g. Stream order, aspect, etc. Get as tifs
                                                             
## Just projecting into particular parts of NSW (this will take the resolution from the rasters) - change in project function
                                                             
## Buffer shorter distances to reflect species' lack of mobility - < 50km


## Print the species run to the screen
message('Running maxent models for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")


#########################################################################################################################
## Load SDM table - this will be the data from Tracey
## Could build it from steps 4-6. Need
## Occurrence data as one table
## All the rasters they are using 
if(read_data == "TRUE") {
  
  ## read in RDS files from previous step
  SDM.SPAT.ALL = readRDS(paste0('data/base/HIA_LIST/COMBO/SDM_SPAT_ALL_', save_run, '.rds'))
  
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
## Loop over all the species 
lapply(GBIF.spp, function(spp){ 
  
  ## Skip the species if the directory already exists, before the loop
  outdir <- maxent_dir
  
  dir_name = file.path(maxent_path, gsub(' ', '_', spp))
  if(dir.exists(dir_name)) {
    message('Skipping ', spp, ' - already run.')
    invisible(return(NULL))
    
  }

  ##  Create the directory so other parallel runs don't try to do it
  dir.create(dir_name)
  write.csv(data.frame(), file.path(dir_name, "in_progress.txt"))
  
  
  ## Print the taxa being processed to screen
  if(spp %in% SDM.SPAT.ALL$searchTaxon) {
    message('Doing ', spp) 
    
    ## Subset the records to only the taxa being processed
    occurrence <- subset(SDM.SPAT.ALL, searchTaxon == spp)# & SOURCE != "INVENTORY")
    
    ## Now get the background points. These can come from any spp, other than the modelled species.
    background <- subset(SDM.SPAT.ALL, searchTaxon != spp)
    
    ## Finally fit the models using FIT_MAXENT_TARG_BG. Also use tryCatch to skip any exceptions
    tryCatch(
      FIT_MAXENT_TARG_BG(occ                     = occurrence,    ## Taken from Tracey 
                         bg                      = background,    ## Maybe use other amphibians?
                         sdm.predictors          = sdm.select,    ## Env variables - need these layers from Leroy  
                         name                    = spp,           ## same
                         outdir,                                  ## unchanged                         
                         template.raster,                         ## Empty raster, from sdm.predictors 
                         min_n                   = 20,            ## This should be higher...
                         max_bg_size             = 70000,         ## they won't get there...
                         Koppen                  = Koppen_1975,   ## Try, but prob don't want to use this?
                         background_buffer_width = 50000,         ## 50km, or even shorter?      
                         shapefiles              = TRUE,          ## Unchanged - shapefiles of occ and bg
                         features                = 'lpq',         ## Could change to hinge, etc. 
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

  #  now add a file to the dir to denote that it has completed
  write.csv(data.frame(), file.path(dir_name, "completed.txt"))
  
})





#########################################################################################################################
## OUTSTANDING SDM TASKS:
#########################################################################################################################



## 1). 400 tree spp that are commonly planted and traded, with sufficient spatial data to model robustly: done


##     Update Species table with new ALA data: done 


## 2). Fixed the taxonomy

##     Follow steps in file 3). 

##     Keep the 'source' column in the maxent table :: adjust step 7  

     
##     Recuced BG points to 70k
##     
##     
##     Check mapping function is working for campbelltown species.......................................................
##     setdiff(camp.spp, TPL.SPP)
     

## 4). Use more forgiving thresholds (10%) for all species: Done
##     "X10.percentile.training.presence.Logistic.threshold"
##     "X10.percentile.training.presence.training.omission"
  

## 5). Combine output: 


##     Create table for appendix.....................................................................................

##     Use 2016 SUA shapefile........................................................................................
##     This causes mistmatches in the SUAs - e.g. Central coast is gone..............................................

##      Table
##     "searchTaxon" "No.plantings" "Number.of.growers" "Number.of.States/LGAs/Koppen zones" "Origin" "AUC" "TSS" 

##     Maps: unique raster values for (search in explorer)
##     "current_suit_above"
##     "2030_4GCMs_above"
##     "2050_4GCMs_above"
##     "2070_4GCMs_above"

   
## 6). Check combined maps for all species : make sure step 9 is working. 
##     Create a spreadsheet of all species results. Check as before - 

##     Summary statistics, 
##     input points, 
##     current raster,
##     Future suitability


## 7). Create tables in R to summarise gain and loss in each SUA
##     Update figures and tables - only doing species gain and stable, not loss.............................................


## 8). What is the story? Draft results and discussion. Linda to create story - don't do turnover


##     Fix mapping function to use mapply in function........................................................................

   
##     Compare output for MESS species - which outputs?......................................................................


##     Run modelling and mapping steps through Katana. To do this, we need a loop that processes one species at a time
##     This means processing the data up to step 7, then running steps 7 and 8 through Katana 


##     Check Shawn's link for flattening the data...................................................................


    
##     Where does the mess mask layer get used, in the combine function?
##     If so, we need a separate folder for the MESS species to process these
##     My combine function does everything in one step, so it's not the same process as what John uses





##     If needed, thin records for 28 spp. with boundary bias, using the SDM tool box.

##     Settings: Max 5km, min 1km, 5 heterogeneity classes (from a classification of a PCA, using the worldclim layers)
##     Currently using all spp background records: random 100k for every species, sounds ok?

##     Try to thin the background records in the same way - run the BG points through SDM toolbox too. 

##     Compare the models with and without thinning for a 28 species:
##     Acacia_implexa, etc.





#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################
