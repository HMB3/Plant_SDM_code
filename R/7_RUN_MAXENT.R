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


## This section needs checking after the HIA_READ_DATA............................................................................
## length(intersect(MAXENT.SUMMARY.NICHE$searchTaxon, APNI$searchTaxon))


## Print the species run to the screen
message('Running maxent models for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")


#########################################################################################################################
## Load SDM table
if(read_data == "TRUE") {
  
  ## read in RDS files from previous step :: use the hollow species as a test
  ## SDM.SPAT.OCC.BG <- readRDS("H:/green_cities_sdm/data/ANALYSIS/SDM_SPAT_OCC_BG_HOLLOW_SPP.rds")
  SDM.SPAT.OCC.BG = readRDS(paste0(DATA_path, 'SDM_SPAT_OCC_BG_',  save_run, '.rds'))
  
} else {
  
  message(' skip file reading, not many species analysed')   ##
  
}





#########################################################################################################################
## 1). RUN SDMs USING ALL A-PRIORI VARIABLES FOR ALL SPECIES
#########################################################################################################################


#########################################################################################################################
## Check the SDM and background tables
dim(SDM.SPAT.OCC.BG)


## The background data is simply the records for all 3.8k species previously downloaded.
## Three souces :: ALA, GBIF, INVENTORY
## So we need to sample bg records from this file in the same proportions as they exist in the occurrence data for each species
## However, some species have high proportions of inventory records in the occurrence data
## But because of the way the background data is selected with 200km of occurrence points and inside the same koppen zone
## That means that we can't choose enough records in proportion with the background?
length(intersect(unique(SDM.SPAT.OCC.BG$searchTaxon), GBIF.spp))  ## This should be same as the number of species


#########################################################################################################################
## Run Maxent using a random selection of background points. Ideally make these projections exactly the same
projection(template.raster.5km);projection(SDM.SPAT.OCC.BG);projection(Koppen_1975)



## Here are the argumetns needed to run the targetted background selection SDMs inside the function itself
spp                     = GBIF.spp[3]
occ                     = subset(SDM.SPAT.OCC.BG, searchTaxon == spp)
bg                      = subset(SDM.SPAT.OCC.BG, searchTaxon != spp)
sdm.predictors          = sdm.select 
name                    = spp 
outdir                  = maxent_dir
template.raster         = template.raster.1km   ## 1km, 5km, 10km
min_n                   = 20            
max_bg_size             = 70000         
Koppen                  = Koppen_1975
background_buffer_width = 200000
shapefiles              = TRUE
features                = 'lpq'
replicates              = 5
responsecurves          = TRUE


## The current code is failing after the background selection, but only for some species.
## why?.................................................................................................................
## Proportionally sampling background data for Corymbia citriodora, OCC + INV records total NA


## First vary the resolution, then run the background selection.........................................................
## back to line 534.........................................


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
  if(spp %in% SDM.SPAT.OCC.BG$searchTaxon) {
    message('Doing ', spp) 
    
    ## Subset the records to only the taxa being processed
    occurrence <- subset(SDM.SPAT.OCC.BG, searchTaxon == spp)
    
    ## Now get the background points. These can come from any spp, other than the modelled species.
    background <- subset(SDM.SPAT.OCC.BG, searchTaxon != spp)
    
    ## Finally fit the models using FIT_MAXENT_TARG_BG. Also use tryCatch to skip any exceptions
    tryCatch(
      fit_maxent_targ_bg_kopp(occ                     = occurrence,    ## name from the .rmd CV doc 
                              bg                      = background,    ## name from the .rmd CV doc  
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
## 2). RUN SDMs USING BACKWARDS SLECTION 
#########################################################################################################################


## Variables to run an example within the backwards selection function
sdm.predictors          = sdm.predictors
name                    = spp
maxent_path             = './output/maxent/SUA_TREES_ANALYSIS/'
outdir                  = outdir
template.raster         = template.raster
min_n                   = 20       ## This should be higher...
shapefiles              = FALSE    ## name this
features                = 'lpq'    ## name 
replicates              = 5
cor_thr                 = 0.7      ## The maximum allowable pairwise correlation between predictor variables 
pct_thr                 = 5        ## The minimum allowable percent variable contribution 
k_thr                   = 3        ## The minimum number of variables to be kept in the model.
responsecurves          = TRUE     ## Response curves


#########################################################################################################################
## Loop over all the species spp = GBIF.spp[1]
lapply(GBIF.spp, function(spp){ 
  
  ## Skip the species if the directory already exists, before the loop
  outdir <- bs_dir
  
  dir_name = file.path(maxent_path, gsub(' ', '_', spp))
  if(dir.exists(dir_name)) {
    message('Skipping ', spp, ' - already run.')
    invisible(return(NULL))
    
  }
  
  #  create the directory so other parallel runs don't try to do it
  dir.create(dir_name)
  write.csv(data.frame(), file.path(dir_name, "in_progress.txt"))
  
  
  ## Print the taxa being processed to screen
  if(spp %in% SDM.SPAT.BG$searchTaxon) {
    message('Doing ', spp) 
    
    ## Subset the records to only the taxa being processed
    ## This should just read in the BG and occ data from the previous function?.........................................
    message('Reading previously created occurrence and background data from targetted SDM for ', spp)
    save_spp = gsub(' ', '_', spp)
    
    # occ     <- readRDS(sprintf('%s%s/%s_occ.rds', maxent_path, save_spp, save_spp))
    # bg      <- readRDS(sprintf('%s%s/%s_bg.rds',  maxent_path, save_spp, save_spp))
    # 
    # swd_occ <- readRDS(sprintf('%s%s/%s_occ_swd.rds', maxent_path, save_spp, save_spp))
    # swd_bg  <- readRDS(sprintf('%s%s/%s_bg_swd.rds',  maxent_path, save_spp, save_spp)) 
    
    ##
    occurrence <- subset(SDM.SPAT.OCC.BG, searchTaxon == spp)# & SOURCE != "INVENTORY")
    
    ## Now get the background points. These can come from any spp, other than the modelled species.
    background <- subset(SDM.SPAT.OCC.BG, searchTaxon != spp)
    
    ## Finally fit the models using FIT_MAXENT_TARG_BG. Also use tryCatch to skip any exceptions
    tryCatch(
      fit_maxent_targ_bs(sdm.predictors          = sdm.predictors, ## List of predictor variables
                         name                    = spp,
                         maxent_path             = './output/maxent/SUA_TREES_ANALYSIS/',
                         outdir,                    ## Change the outdir on the new species 
                         template.raster,
                         min_n                   = 20,       ## This should be higher...
                         shapefiles              = FALSE,    ## name this
                         features                = 'lpq',    ## name 
                         replicates              = 5,
                         cor_thr                 = 0.7,
                         pct_thr                 = 5,
                         k_thr                   = 3,
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
## 2). TABULATE MAXENT STATISTICS
#########################################################################################################################


## This section needs additional tidy up ................................................................................
## Everything that doesn't get used needs to go .........................................................................


## Print the species run to the screen
message('Creating summary stats for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")


## Read in niche data
if(read_data == "TRUE") {
  
  ## Load GBIF and ALA data
  message('Reading niche data for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")
  COMBO.NICHE.CONTEXT = readRDS(paste0(DATA_path, 'COMBO_NICHE_CONTEXT_',  OCC_SOURCE, '_RECORDS_', save_run, '.rds'))
  
} else {
  
  message(' skip file reading, running species in parallel')   ##
  
}


#########################################################################################################################
## Create a file list for each model run: Try crunching this into just the species required
map_spp_list  = gsub(" ", "_", GBIF.spp)
maxent.tables = list.files(maxent_path)                 
maxent.tables = intersect(maxent.tables, map_spp_list)   
maxent_path   = maxent_path                             
length(maxent.tables) 


#########################################################################################################################
## This section illustrates how to index the maxent model object, which could be applied to other methods     #---------#
## Create a table of the results 
## x = maxent.tables[1]
MAXENT.RESULTS <- maxent.tables %>%         
  
  ## Pipe the list into lapply
  lapply(function(x) {
    
    ## Create the character string
    f <- paste0(maxent_path, x, "/full/maxentResults.csv")
    
    ## Load each .RData file
    d <- read.csv(f)
    
    ## load model
    m           = readRDS(sprintf('%s/%s/full/maxent_fitted.rds', maxent_path, x))$me_full
    number.var  = length(m@lambdas) - 4   ## (the last 4 slots of the lambdas file are not variables)
    mxt.records = length(m@presence$Annual_mean_temp)
    
    ## Get variable importance
    m.vars    = ENMeval::var.importance(m)
    var.pcont = m.vars[rev(order(m.vars[["percent.contribution"]])),][["variable"]][1]
    pcont     = m.vars[rev(order(m.vars[["percent.contribution"]])),][["percent.contribution"]][1]
    
    var.pimp  = m.vars[rev(order(m.vars[["permutation.importance"]])),][["variable"]][1]
    pimp      = m.vars[rev(order(m.vars[["permutation.importance"]])),][["permutation.importance"]][1]
    
    ## Get maxent results columns to be used for model checking
    Training_AUC             = m@results["Training.AUC",]
    Number_background_points = m@results["X.Background.points",]
    Logistic_threshold       = m@results["X10.percentile.training.presence.Logistic.threshold",] 
    
    ## Now rename the maxent columns that we will use in the results table
    d = cbind(searchTaxon              = x, 
              Maxent_records           = mxt.records,
              Number_var               = number.var, 
              Var_pcont                = var.pcont,
              Per_cont                 = pcont,
              Var_pimp                 = var.pimp,
              Perm_imp                 = pimp,
              Training_AUC,
              Number_background_points,  
              Logistic_threshold,
              d)  
    dim(d)
    
    ## Remove path gunk, and species
    d$Species    = NULL
    return(d)
    
  }) %>%
  
  ## Finally, bind all the rows together
  bind_rows


#########################################################################################################################
## Create True Skill Statistic values (TSS) for each species
TSS.tables = list.files(maxent_path, pattern = 'species_omission\\.csv$', full.names = TRUE, recursive = TRUE)


## Get the maxium TSS value using the omission data : use _training_ ommision data only
max_tss <- sapply(TSS.tables, function(f) {
  
  ## For eachg species, read in the training data
  d <- read.csv(f)
  i <- which.min(d$Training.omission + d$Fractional.area)
  
  c(max_tss = 1 - min(d$Training.omission + d$Fractional.area),
    thr     = d$Corresponding.logistic.value[i])
  
})


## Add a species variable to the TSS results, so we can subset to just the species analysed
max_tss  = as.data.frame(max_tss)
setDT(max_tss, keep.rownames = TRUE)[]
max_tss  = as.data.frame(max_tss)
names(max_tss)[names(max_tss) == 'rn'] <- 'searchTaxon'


## Remove extra text
max_tss$searchTaxon = gsub("//",         "/", max_tss$searchTaxon)
max_tss$searchTaxon = gsub(maxent_path,   "", max_tss$searchTaxon)
max_tss$searchTaxon = gsub("/full/species_omission.csv.max_tss", "", max_tss$searchTaxon)
head(max_tss)


## Add max TSS to the results table
MAXENT.RESULTS = merge(MAXENT.RESULTS, max_tss, by = "searchTaxon", ALL = FALSE)
head(MAXENT.RESULTS$max_tss)


## This is a summary of maxent output for current conditions
## Also which species have AUC < 0.7?
dim(MAXENT.RESULTS)
head(MAXENT.RESULTS, 20)[1:5]
dim(subset(MAXENT.RESULTS, Training.AUC > 0.7))  ## all models should be above 0.7


#########################################################################################################################
## Plot AUC vs. TSS
if (nrow(MAXENT.RESULTS) > 2) {
  
  lm.auc = lm(MAXENT.RESULTS$max_tss ~ MAXENT.RESULTS$Training.AUC)
  hist(MAXENT.RESULTS$max_tss)
  plot(MAXENT.RESULTS$Training.AUC, MAXENT.RESULTS$max_tss, pch = 19, col  = "blue",
       xlab = "AUC", ylab = "TSS", 
       abline(lm(MAXENT.RESULTS$max_tss ~ MAXENT.RESULTS$Training.AUC)))
  legend("topleft", bty="n", legend=paste("R2 is", format(summary(lm.auc)$adj.r.squared, digits = 4)))
  
  ## If the species has < 2 records, don't plot
} else {
  
  message('Dont plot, only ', length(GBIF.spp), ' species analysed')
  
}


## Now check the match between the species list, and the results list. 
length(intersect(map_spp_list, MAXENT.RESULTS$searchTaxon)) 
MAXENT.RESULTS  =  MAXENT.RESULTS[MAXENT.RESULTS$searchTaxon %in% map_spp_list , ] 
map_spp         = unique(MAXENT.RESULTS$searchTaxon)
length(map_spp);setdiff(sort(map_spp_list), sort(map_spp))


#########################################################################################################################
## Then, make a list of all the directories containing the individual GCM rasters
SDM.RESULTS.DIR <- map_spp %>%
  
  ## Pipe the list into lapply
  lapply(function(species) {
    
    ## Create the character string...
    m <-   sprintf('%s%s/full/', maxent_path, species)                ## path.backwards.sel
    m 
    
  }) %>%
  
  ## Bind the list together
  c()

length(SDM.RESULTS.DIR)
SDM.RESULTS.DIR = unlist(SDM.RESULTS.DIR)


#########################################################################################################################
## Now combine the SDM output with the niche context data 
NICHE.CONTEXT = COMBO.NICHE.CONTEXT[, c("searchTaxon",     "Origin",      "Plant_type", 
                                        "GLOBAL_RECORDS",  "AUS_RECORDS", "Plantings", "SUA_COUNT",
                                        "Total_growers",   "Number_states")]


## Which columns will help with model selection
MAXENT.SUMMARY   = MAXENT.RESULTS[, c("searchTaxon",
                                      "Maxent_records",
                                      "Number_var",
                                      "Var_pcont",
                                      "Per_cont",
                                      "Var_pimp",
                                      "Perm_imp",
                                      "Iterations",                                                                        
                                      "Training_AUC",
                                      "max_tss",
                                      "Number_background_points",  
                                      "Logistic_threshold")]


## Now remove the underscore and join data 
MAXENT.SUMMARY$searchTaxon = gsub("_", " ", MAXENT.SUMMARY$searchTaxon)
MAXENT.SUMMARY.NICHE       = join(MAXENT.SUMMARY,       NICHE.CONTEXT,  type = "left")  ## join does not support the sorting
MAXENT.SUMMARY.NICHE       = MAXENT.SUMMARY.NICHE[order(MAXENT.SUMMARY.NICHE$searchTaxon),]


## Re-order table into a logical order
## What is the difference between SUA_count here, and in the SUA table? Here it is the number of SUA's that species occurs in,
## but in the SUA table, it is the number of records for each species in each SUA
MAXENT.SUMMARY.NICHE = MAXENT.SUMMARY.NICHE[, c("searchTaxon",     ## Columns of evergreen list 
                                                "Origin",
                                                "Plant_type",
                                                "Total_growers",     
                                                "Number_states",
                                                
                                                "Plantings",       ## Columns for counts of records
                                                "GLOBAL_RECORDS",  
                                                "AUS_RECORDS",
                                                "Maxent_records",
                                                "SUA_COUNT",
                                                
                                                "Number_var",      ## Columns for maxent results
                                                "Var_pcont",
                                                "Per_cont",
                                                "Var_pimp",
                                                "Perm_imp",
                                                "Iterations",                                                                        
                                                "Training_AUC",
                                                "max_tss",
                                                "Number_background_points",  
                                                "Logistic_threshold")]


#########################################################################################################################
## Create a taxonomic lookup table
SDM.TAXA <- MAXENT.SUMMARY.NICHE[["searchTaxon"]] %>%  
  lookup_table(., by_species = TRUE) 
SDM.TAXA <- setDT(SDM.TAXA, keep.rownames = TRUE)[]
colnames(SDM.TAXA)[1] <- "searchTaxon"
SDM.TAXA <- join(SDM.TAXA, APNI)                     ## get rid of APNI


## Join on the native data and the APNI
MAXENT.SUMMARY.NICHE <- SDM.TAXA %>% 
  
  join(., MAXENT.SUMMARY.NICHE) #%>% 
  
  #join(., MAXENT.CHECK[c("searchTaxon", "check.map")])    ## get rid of check map, circualar.............................

length(intersect(MAXENT.SUMMARY.NICHE$searchTaxon, GBIF.spp))
View(MAXENT.SUMMARY.NICHE)


#########################################################################################################################
## The supplementary matieral table for the STOTEN MS has :: 
## Species, Records	Plantings	Variable	AUC	TSS	Threshold
## Save these columns as CSV file :: use this as quick check of model performance for different runs
SUA.MS.TABLE              = select(MAXENT.SUMMARY.NICHE, searchTaxon, Maxent_records, Plantings, 
                                   Var_pimp, Training_AUC, max_tss, Logistic_threshold)
SUA.MS.TABLE$Training_AUC       = round(SUA.MS.TABLE$Training_AUC, 3)
SUA.MS.TABLE$max_tss            = round(SUA.MS.TABLE$Training_AUC, 3)
SUA.MS.TABLE$Logistic_threshold = round(SUA.MS.TABLE$Training_AUC, 3)



#########################################################################################################################
## 3). CREATE LISTS OF MAXENT THRESHOLDS FOR MAPPING
#########################################################################################################################


#########################################################################################################################
## Maxent produces a presence threshold for each species (i.e. the columns in MAXENT.RESULTS). 
## John : for AUC, you can report the cross-validated test AUC (if your code currently runs a cross-validated model as well), 
## and for the model threshold (for binarising) you can just use the training value (or the crossval one...there's little 
## guidance about this and you can really get away with either).


#########################################################################################################################
## How do the differnt thresholds compare for the set of species modelled?
summary(MAXENT.RESULTS["Maximum.training.sensitivity.plus.specificity.Logistic.threshold"])    ## The strictest threshold
summary(MAXENT.RESULTS["X10.percentile.training.presence.Logistic.threshold"])                 ## The next strictest
summary(MAXENT.RESULTS["X10.percentile.training.presence.training.omission"])                  ## The most forgiving


#########################################################################################################################
## Now turn the maxent results into lists :: we can use these to generate the consensus layers, whereby habitat suitablity
## is contrained below the threshold. This number is created by the current model, and used for all the future models. 
thresh.max.train  = as.list(MAXENT.RESULTS["Maximum.training.sensitivity.plus.specificity.Logistic.threshold"]) 
thresh.max.train  = thresh.max.train$Maximum.training.sensitivity.plus.specificity.Logistic.threshold

percent.10.log    = as.list(MAXENT.RESULTS["X10.percentile.training.presence.Logistic.threshold"])  
percent.10.log    = percent.10.log$X10.percentile.training.presence.Logistic.threshold

percent.10.om     = as.list(MAXENT.RESULTS["X10.percentile.training.presence.training.omission"])   
percent.10.om     = percent.10.om$X10.percentile.training.presence.training.omission


#########################################################################################################################
## Could create a plot of number of records vs. maxent rating, Boxplot with no. occurrences on y, and maxent rating on x
# plot(MAXENT.SUMMARY.NICHE$GLOBAL_RECORDS, MAXENT.SUMMARY.NICHE$check.map, pch = 19, col  = "blue",
#      xlab = "AUC", ylab = "TSS", 
#      abline(lm(MAXENT.SUMMARY.NICHE$check.map ~ MAXENT.SUMMARY.NICHE$GLOBAL_RECORDS)))
# legend("topleft", bty="n", legend=paste("R2 is", format(summary(lm.auc)$adj.r.squared, digits = 4)))


#########################################################################################################################
## Save maxent results - what happens to these files, if run for one species at a time?
#########################################################################################################################
## save data
if(save_data == "TRUE") {
  
  ## save .rds file for the next session
  saveRDS(MAXENT.SUMMARY.NICHE,   paste0(DATA_path, 'MAXENT_SUMMARY_',         save_run, '.rds'))
  write.csv(MAXENT.SUMMARY.NICHE, paste0(DATA_path, 'MAXENT_SUMMARY_',         save_run, '.csv'), row.names = FALSE)
  write.csv(SUA.MS.TABLE,         paste0(DATA_path, 'MAXENT_SUA_MS_SUPP_MAT_', save_run, '.csv'), row.names = FALSE)
  
} else {
  
  message('Dont save niche summary, only ', length(GBIF.spp), ' species analysed')
  
}


#########################################################################################################################
## OUTSTANDING SDM TASKS:
#########################################################################################################################



## 1). Model species in parallel 

## 2). Harvest folders once they have been run in Katana :: 7-zip

## 3). Clean the joinging of horticultural data so it is clean for the new run of species








#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################