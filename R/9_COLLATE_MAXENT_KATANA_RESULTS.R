#########################################################################################################################
############################################# SUMARISE MAXENT RESULTS ################################################### 
#########################################################################################################################


#########################################################################################################################
## This code collates the SDM results for species analyzed on Katana....................................................
## First step is to get the katana models run.


## Read in niche data :: this is needed to combine the context data with the maxent output
if(read_data == "TRUE") {
  
  ## Load GBIF and ALA data :: this would ideally be for all species, not just for those run
  ## the full niche file is :: 'COMBO_NICHE_CONTEXT_ALL_EVREGREEN_MAY_2018.rds'
  COMBO.NICHE.CONTEXT = readRDS(paste0(DATA_path, 'COMBO_NICHE_CONTEXT_',  save_run, '.rds'))
  message('Reading niche data for ', length(unique(COMBO.NICHE.CONTEXT$searchTaxon)), 
          ' species in the set ', "'", save_run, "'")
  
} else {
  message(' skip file reading, running species in parallel')   ##
}





#########################################################################################################################
## 1). TABULATE MAXENT STATISTICS
#########################################################################################################################


## The code that adds niche info is now in './R/COLLATE_MAXENT_RESULTS.R' 
## The code below is just for running on Katana
## Print the species run to the screen
message('Creating summary stats for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")


#########################################################################################################################
## First, make a list of all the species with models, then restrict them to just the models on the GBIF.spp list 
map_spp_list  = gsub(" ", "_", GBIF.spp)
map_spp_patt  = paste0(map_spp_list, collapse = "|")
message ("map_spp_list head:")
message (paste (head(map_spp_list), collapse=","))


#########################################################################################################################
## Now stop R from listing all the maxent files that have completed - this takes a long time
message(results_dir)
maxent.tables = lapply (map_spp_list, FUN = function (x) {paste(results_dir , x, "full/maxent_fitted.rds", sep="/")})


## How many species have been modelled?
message(paste("maxent.tables has this many entries:", length(maxent.tables)))
message(paste(head (maxent.tables), collapse=","))
sdm.exists = lapply(maxent.tables, FUN = function (x) {file.exists (x)})
sdm.exists = unlist(sdm.exists)


## Only list the intersection between the modelled species and 
message(paste(head(sdm.exists), collapse=","))
maxent.tables = maxent.tables[sdm.exists]


message (paste ("maxent.tables has this many entries:", length(maxent.tables)))
maxent.tables = stringr::str_subset(maxent.tables, map_spp_patt)
message (paste ("maxent.tables has this many entries:", length(maxent.tables)))
message (paste (head(maxent.tables), collapse=","))
length(maxent.tables)


#########################################################################################################################
## Now create a table of the results 
## x = maxent.tables[1]
MAXENT.RESULTS <- maxent.tables %>%         
  
  ## Pipe the list into lapply
  lapply(function(x) {
    
    ## We don't need to skip species that haven't been modelled
    x = gsub(paste0(results_dir, "/"), "", x) 
    message (x)
    
    #############################################################
    ## load the backwards selected model
    if (grepl("BS", results_dir)) {
      m = readRDS(paste0(results_dir, '/',  x))
      
    } else {
      ## Get the background records from any source
      m = readRDS(paste0(results_dir, '/',  x))$me_full
      
    }
    
    ## Get the number of Variables
    number.var  = length(m@lambdas) - 4   ## (the last 4 slots of the lambdas file are not variables)
    mxt.records = nrow(m@presence)
    
    ## Get variable importance
    m.vars    = ENMeval::var.importance(m)
    var.pcont = m.vars[rev(order(m.vars[["percent.contribution"]])),][["variable"]][1]
    pcont     = m.vars[rev(order(m.vars[["percent.contribution"]])),][["percent.contribution"]][1]
    
    var.pimp  = m.vars[rev(order(m.vars[["permutation.importance"]])),][["variable"]][1]
    pimp      = m.vars[rev(order(m.vars[["permutation.importance"]])),][["permutation.importance"]][1]
    
    ## Get maxent results columns to be used for model checking
    ## Including the omission rate here
    Training_AUC             = m@results["Training.AUC",]
    Iterations               = m@results["Iterations",]
    Number_background_points = m@results["X.Background.points",]
    Logistic_threshold       = m@results["X10.percentile.training.presence.Logistic.threshold",]
    Omission_rate            = m@results["X10.percentile.training.presence.training.omission",]
    
    ## Now rename the maxent columns that we will use in the results table
    d = data.frame(searchTaxon              = x, 
                   Maxent_records           = mxt.records,
                   Number_var               = number.var, 
                   Var_pcont                = var.pcont,
                   Per_cont                 = pcont,
                   Var_pimp                 = var.pimp,
                   Perm_imp                 = pimp,
                   Iterations,
                   Training_AUC,
                   Number_background_points,  
                   Logistic_threshold,
                   Omission_rate)
    
    ## Remove path gunk, and species
    d$Species     = NULL
    d$searchTaxon = gsub("/full/maxent_fitted.rds", "", d$searchTaxon)
    return(d)
    
  }) %>%
  
  ## Finally, bind all the rows together
  bind_rows


#########################################################################################################################
## Now create a list of the '10th percentile training presence Logistic threshold'. This is used in step 8 to threshold
## the maps to just areas above the threshold.
message ("MAXENT.RESULTS columns") 
message (paste (colnames (MAXENT.RESULTS)))
message (paste (nrow (MAXENT.RESULTS)))
summary(MAXENT.RESULTS["Logistic_threshold"])   
percent.10.log = as.list(MAXENT.RESULTS["Logistic_threshold"])  
percent.10.log = percent.10.log$Logistic_threshold


#########################################################################################################################
## Create a list of the omission files - again, don't do this for all the files, just the intersection
omission.tables = lapply (map_spp_list, FUN = function (x) {paste(results_dir , x, "full/species_omission.csv", sep="/")})
message (head (omission.tables))


## Only process the existing files 
om.exists = lapply (omission.tables, FUN = function (x) {file.exists (x)})
om.exists = unlist(om.exists)


omission.tables = omission.tables[om.exists]
message(head(omission.tables))


## Get the maxium TSS value using the omission data : use _training_ ommission data only
Max_tss <- sapply(omission.tables, function(file) {
  
  ## For eachg species, read in the training data
  d <- read.csv(file)
  i <- which.min(d$Training.omission + d$Fractional.area)
  
  c(Max_tss = 1 - min(d$Training.omission + d$Fractional.area),
    thr     = d$Corresponding.logistic.value[i])
  
})



#########################################################################################################################
## Add a species variable to the TSS results, so we can subset to just the species analysed
Max_tss        = as.data.frame(Max_tss)
MAXENT.RESULTS = cbind(MAXENT.RESULTS, Max_tss)
summary(MAXENT.RESULTS$Max_tss)
summary(MAXENT.RESULTS$Omission_rate)


## Check the results
dim(MAXENT.RESULTS)
head(MAXENT.RESULTS, 10)[1:5]





#########################################################################################################################
## 2). COMBINE SDM RESULTS WITH HIA CONTEXT DATA & VISUALISE RELATIONSHIPS BETWEEN VARIABLES
#########################################################################################################################


#########################################################################################################################
## So we've got a range of quantiative and qualitative measures for each species (e.g. maxent output, number of records,
## model score). Before embarking on a machine learning odessy, can we see some patterns in the 450-odd species that have
## already been rated by hand? Simple categorical visualations in R should be sufficient to enough to get at the patterns.

## So we have a discreate variable - good, fair, poor, and some quantiative and discrete variables. How can they be 
## analysed? 


#########################################################################################################################
## Change the species column so we can join on the contextual data
MAXENT.RESULTS$searchTaxon = gsub("_", " ", MAXENT.RESULTS$searchTaxon)
MAXENT.CONTEXT = join(COMBO.NICHE.CONTEXT, MAXENT.RESULTS)
MAXENT.CONTEXT = dplyr::select(MAXENT.CONTEXT, results.columns)


## 476 species have been rated, but we need to check if the ratings are still applicable, as the data has been updated
length(intersect(MAXENT.RATING.LAT$searchTaxon, WPW.spp))
length(intersect(MAXENT.RATING.LAT$searchTaxon, MAXENT.RESULTS$searchTaxon))


## Add a new text string for good/fair/poor models :: might be easier to use in visualisations
## Join on the new rating
MAXENT.RATING.LAT = mutate(MAXENT.RATING.LAT, 
                           MAXENT_RATING = ifelse(CHECK_MAP %in% 1, "GOOD",
                                                  ifelse(CHECK_MAP %in% 2, "FAIR",
                                                         ifelse(CHECK_MAP %in% 3, "POOR", "UN-RATED"))))
MAXENT.RATING.LAT$CHECK_MAP = NULL
MAXENT.CONTEXT = join(MAXENT.CONTEXT, MAXENT.RATING.LAT)


#########################################################################################################################
if(explore_maps == "TRUE") {
  
  ## save .rds file for the next session
  source('./R/9B_COLLATE_MAXENT_PATTERNS.R')
  
} else {
  message('Dont explore patterns in SDM results')
}


#########################################################################################################################
## Save maxent results 
if(save_data == "TRUE") {
  
  ## save .rds file for the next session
  saveRDS(MAXENT.RESULTS,   paste0(DATA_path, 'MAXENT_RESULTS_', save_run, '.rds'))
  write.csv(MAXENT.RESULTS, paste0(DATA_path, 'MAXENT_RESULTS_', save_run, '.csv'), row.names = FALSE)
  
} else {
  message('Dont save niche summary, only ', length(GBIF.spp), ' species analysed')
}





#########################################################################################################################
## 3). RUN MESS MAPS LOCALLY - BEST TO GET IT WORKING ON LINUX....
#########################################################################################################################


#########################################################################################################################
## Now check the match between the species list, and the results list. 
length(intersect(map_spp_list, MAXENT.RESULTS$searchTaxon)) 
MAXENT.CONTEXT =  MAXENT.RESULTS[MAXENT.RESULTS$searchTaxon %in% map_spp_list , ] 
map_spp         = unique(MAXENT.RESULTS$searchTaxon)
length(map_spp)


## Run MESS maps for all species - this will take awhile!
tryCatch(
  create_mess_pngs(shp_path      = "./data/base/CONTEXTUAL/", ## Path for shapefile
                   aus_shp       = "aus_states.rds",          ## Shapefile, e.g. Australian states
                   world_shp     = "LAND_world.rds",          ## World shapefile          
                   scen_list     = scen_2030,                 ## List of climate scenarios
                   species_list  = map_spp,                   ## rev(map_spp)
                   maxent_path   = bs_path,
                   time_slice    = 30,
                   nclust        = 1),
  
  error = function(cond) {
    
    ## This will write the error message inside the text file, 
    ## but it won't include the species
    file.create(file.path(bs_path, "map_png_failed_2030.txt"))
    cat(cond$message, file=file.path(bs_path, "map_png_failed_2030.txt"))
    warning(cond$message)
    
  })




## Then, we need linux/windows code to find all files with strings, and copy and move them to a new location - the
## same directory where the above table is stored...




#########################################################################################################################
## OUTSTANDING COLLATION TASKS:
#########################################################################################################################


## 1).  Introduce more switches for operating system, etc

## 2).  Tidy up all the code (using piping, etc.)

## 3).  Make the code more modular, less monolithic (e.g. code blocks, markdown, etc.)

## 3).  Iron out some of the points which are not as reproducible (e.g the background points in step 7)

## 4).  Improve the reading of objects and files - e.g. the RData object is cumbersome.... 






#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################