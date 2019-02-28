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
  
  ## read in RDS files from previous step
  SDM.SPAT.OCC.BG = readRDS(paste0(DATA_path, 'SDM_SPAT_OCC_BG',  save_run, '.rds'))
  
} else {
  
  message(' skip file reading, not many species analysed')   ##
  
}





#########################################################################################################################
## 1). RUN SDMs USING ALL A-PRIORI VARIABLES FOR ALL SPECIES
#########################################################################################################################


#########################################################################################################################
## Check the SDM and background tables
dim(SDM.SPAT.OCC.BG)


## The background data has data from about 6.7k species, created earlier in the project. This table is merged with the 
## data for the species being analysed, and the inverse speices are taken so there is no overlap. The GB table hasn't strictly 
## been created with the same processes as the latest SDM data. Re-create it with urban data
length(unique(background.points$searchTaxon))
length(intersect(unique(SDM.SPAT.OCC.BG$searchTaxon), GBIF.spp))  ## This should be same as the number of species


## Sample the backround records in the same proportion as the sources for each species - E.G. 90% from ALA/GBIF, 10% urban
## This would be created as variable from each occurrence file, the source column would need to be included.
## In theory, we can do this with the existing tables for 3.5 species   
unique(SDM.SPAT.OCC.BG$SOURCE)
SDM.SPAT.DF = as.data.frame(subset(SDM.SPAT.OCC.BG, TYPE != "BG"))
round(with(SDM.SPAT.DF, table(SOURCE)/sum(table(SOURCE))*100), 1)
unique(background$SOURCE)                                      ## Could subset by source


#########################################################################################################################
## Run Maxent using a random selection of background points. Ideally make these projections exactly the same
projection(template.raster);projection(SDM.SPAT.BG);projection(Koppen_1975)



occ                     = occurrence    ## name from the .rmd CV doc 
bg                      = background    ## name from the .rmd CV doc  
sdm.predictors          = sdm.select 
name                    = spp 
outdir                  = outdir 
template.raster         = template.raster
min_n                   = 20            ## This should be higher...
max_bg_size             = 70000         ## could be 50k or lower
Koppen                  = Koppen_1975
background_buffer_width = 200000
shapefiles              = TRUE
features                = 'lpq'
replicates              = 5
responsecurves          = TRUE



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
  if(spp %in% SDM.SPAT.BG$searchTaxon) {
    message('Doing ', spp) 
    
    ## Subset the records to only the taxa being processed
    #occurrence <- subset(SDM.SPAT.OCC.BG, searchTaxon == spp)
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


## Variables to run an example within the function
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
    #occurrence <- subset(SDM.SPAT.OCC.BG, searchTaxon == spp)
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


## Create a table of the results 
MAXENT.RESULTS <- maxent.tables[c(1:length(maxent.tables))] %>%         
  
  ## Pipe the list into lapply
  lapply(function(x) {
    
    ## Create the character string
    f <- paste0(maxent_path, x, "/full/maxentResults.csv")
    
    ## Load each .RData file
    d <- read.csv(f)
    
    ## load model
    m <- readRDS(sprintf('%s/%s/full/maxent_fitted.rds', maxent_path, x))
    m <- m$me_full
    number_var = length(m@lambdas) - 4                                   ## (the last 4 slots of lambdas are not variables) 
    
    ## Get variable importance
    m.vars    = ENMeval::var.importance(m)
    var.pcont = m.vars[rev(order(m.vars[["percent.contribution"]])),][["variable"]][1]
    pcont     = m.vars[rev(order(m.vars[["percent.contribution"]])),][["percent.contribution"]][1]
    
    var.pimp  = m.vars[rev(order(m.vars[["permutation.importance"]])),][["variable"]][1]
    pimp      = m.vars[rev(order(m.vars[["permutation.importance"]])),][["permutation.importance"]][1]
    
    ## Now add a model column
    d = cbind(searchTaxon = x, 
              Number_var  = number_var, 
              Var_pcont   = var.pcont,
              Per_cont    = pcont,
              Var_pimp    = var.pimp,
              Perm_imp    = pimp,
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
NICHE.CONTEXT   = COMBO.NICHE.CONTEXT[, c("searchTaxon",     "Origin",      "Plant.type", 
                                          "GLOBAL_RECORDS",  "AUS_RECORDS", "Plantings", "SUA_COUNT",
                                          "Total.growers",   "Number.of.States")]


## Which columns will help with model selection
MAXENT.SUMMARY   = MAXENT.RESULTS[, c("searchTaxon",
                                      "Number_var",
                                      "Var_pcont",
                                      "Per_cont",
                                      "Var_pimp",
                                      "Perm_imp",
                                      "X.Training.samples",                                                                
                                      "Iterations",                                                                        
                                      "Training.AUC",
                                      "max_tss",
                                      "X.Background.points",  
                                      "X10.percentile.training.presence.Logistic.threshold")]


## Now remove the underscore and join data 
MAXENT.SUMMARY$searchTaxon = gsub("_", " ", MAXENT.SUMMARY$searchTaxon)
MAXENT.SUMMARY.NICHE       = join(MAXENT.SUMMARY,       NICHE.CONTEXT,  type = "left")  ## join does not support the sorting
MAXENT.SUMMARY.NICHE       = MAXENT.SUMMARY.NICHE[order(MAXENT.SUMMARY.NICHE$searchTaxon),]


## Re-order table
## What is the difference between SUA_count here and in the SUA table? Here it is the number of SUA's that species occurs in,
## but in the SUA table, it is the number of records for each species in each SUA
MAXENT.SUMMARY.NICHE   = MAXENT.SUMMARY.NICHE[, c("searchTaxon",
                                                  "Origin",
                                                  "Plant.type", 
                                                  "GLOBAL_RECORDS",  
                                                  "AUS_RECORDS", 
                                                  "Plantings", 
                                                  "SUA_COUNT", 
                                                  "Total.growers",     
                                                  "Number.of.States",
                                                  "Number_var",
                                                  "Var_pcont",
                                                  "Per_cont",
                                                  "Var_pimp",
                                                  "Perm_imp",
                                                  "X.Training.samples",  
                                                  "Iterations",                                                                        
                                                  "Training.AUC",
                                                  "max_tss",
                                                  "X.Background.points",  
                                                  "X10.percentile.training.presence.Logistic.threshold")]


#########################################################################################################################
## Create a taxonomic lookup table
## This needs checking after the HIA_READ_DATA............................................................................
length(intersect(MAXENT.SUMMARY.NICHE$searchTaxon, APNI$searchTaxon))


SDM.TAXA <- MAXENT.SUMMARY.NICHE[["searchTaxon"]] %>%  
  lookup_table(., by_species = TRUE) 
SDM.TAXA <- setDT(SDM.TAXA, keep.rownames = TRUE)[]
colnames(SDM.TAXA)[1] <- "searchTaxon"
SDM.TAXA <- join(SDM.TAXA, APNI)                     ## get rid of APNI


## Join on the native data and the APNI
MAXENT.SUMMARY.NICHE <- SDM.TAXA %>% 
  
  join(., MAXENT.SUMMARY.NICHE) %>% 
  
  join(., MAXENT.CHECK[c("searchTaxon", "check.map")])    ## get rid of check map

length(intersect(MAXENT.SUMMARY.NICHE$searchTaxon, GBIF.spp))





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
  saveRDS(MAXENT.SUMMARY.NICHE,   paste0(DATA_path, 'MAXENT_SUMMARY_',  save_run, '.rds'))
  write.csv(MAXENT.SUMMARY.NICHE, paste0(DATA_path, 'MAXENT_SUMMARY_',  save_run, '.csv'), row.names = FALSE)
  
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