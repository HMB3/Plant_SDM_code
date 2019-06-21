#########################################################################################################################
## 2). TABULATE MAXENT STATISTICS
#########################################################################################################################


## How can this be modified to work either on Katana or locally?
## 1). Re-run all the species, 1000 species at a time. 
## 2). Then combine them.
## 3). 

## Print the species run to the screen
message('Creating summary stats for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")


## Read in niche data
if(read_data == "TRUE") {
  
  ## Load GBIF and ALA data :: this would ideally be for all species, not just for those run
  ## On katana, this will be done for just one species. So 
  message('Reading niche data for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")
  COMBO.NICHE.CONTEXT = readRDS(paste0(DATA_path, 'COMBO_NICHE_CONTEXT_', save_run, '.rds'))
  
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
    
    #############################################################
    ## load model
    if (grepl("BS", maxent_path)) {
      m = readRDS(sprintf('%s/%s/full/maxent_fitted.rds', maxent_path, x))
      
    } else {
      ## Get the background records from any source
      m = readRDS(sprintf('%s/%s/full/maxent_fitted.rds', maxent_path, x))$me_full
      
    }
    
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
  
  ## Save this to file
  png(paste0('./output/maxent/', 'Maxent_run_summary_', save_run, '.png'), 16, 12, units = 'in', res = 500)
  
  layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
  
  plot(MAXENT.RESULTS$Training.AUC, MAXENT.RESULTS$max_tss, pch = 19, col  = "blue",
       xlab = "AUC", ylab = "TSS", 
       abline(lm(MAXENT.RESULTS$max_tss ~ MAXENT.RESULTS$Training.AUC)), 
       main = save_run, cex = 3)
  
  legend("topleft", bty = "n", 
         legend = paste("R2 is", format(summary(lm.auc)$adj.r.squared, digits = 4)))
  
  hist(MAXENT.RESULTS$Training.AUC, breaks = 10, col = "blue",   border = FALSE,
       ylab = "Frequency",
       xlab = "Training AUC", main = "AUC", cex = 3)
  hist(MAXENT.RESULTS$max_tss,      breaks = 10, col = "orange", border = FALSE,
       ylab = "Frequency",
       xlab = "Maximum True Skill Statistic", main = "TSS", cex = 3)
  
  ## Finsish the device
  dev.off()
  
  ## If the species list is < 2 records, don't plot
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
MAXENT.SUMMARY.NICHE       = join(MAXENT.SUMMARY,       NICHE.CONTEXT,  type = "inner")  ## join does not support the sorting
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
SUA.MS.TABLE$max_tss            = round(SUA.MS.TABLE$max_tss, 3)
SUA.MS.TABLE$Logistic_threshold = round(SUA.MS.TABLE$Logistic_threshold, 3)


## Write to CSV here....................................................................................................


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
#####################################################  TBC ############################################################## 
#########################################################################################################################