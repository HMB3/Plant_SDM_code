#########################################################################################################################
###################### PREDICT MAXENT TO FUTURE CLIMATES AND SUMARISE THE RESULTS ####################################### 
#########################################################################################################################


#########################################################################################################################
## This code collates the SDM results for species analyzed on Katana....................................................
## First step is to get the katana models run.


## Then, run all the niches locally, and combine the niche data to the maxent results in this code.


## How does this need to be modified to process Katana output?
## All we need to do is splice on the context info from the 


## Print the species run to the screen
message('Creating summary stats for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")


## Read in niche data
if(read_data == "TRUE") {
  
  ## Load GBIF and ALA data :: this would ideally be for all species, not just for those run
  ## On katana, this will be done for just one species. So 
  message('Reading niche data for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")
  COMBO.NICHE.CONTEXT = readRDS(paste0(DATA_path, 'COMBO_NICHE_CONTEXT_',  save_run, '.rds'))
  
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
## First, create a file list for each model run: Try crunching this into just the species required
## GBIF.spp = maxent.tables
## GBIF.spp = gsub("_", " ", GBIF.spp)
map_spp_list  = gsub(" ", "_", GBIF.spp)
maxent.tables = list.files(results_dir)                 
maxent.tables = intersect(maxent.tables, map_spp_list)   
results_dir   = results_dir                             
length(maxent.tables) 


#########################################################################################################################
## This section illustrates how to index the maxent model object, which could be applied to other methods     #---------#
## Create a table of the results 
## x = maxent.tables[1]
MAXENT.RESULTS <- maxent.tables %>%         
  
  ## Pipe the list into lapply
  lapply(function(x) {
    
    ## Create the character string
    f <- paste0(results_dir, "/",  x, "/full/maxentResults.csv")
    
    ## Load each .RData file
    d <- read.csv(f)
    
    #############################################################
    ## load model
    if (grepl("BS", results_dir)) {
      m = readRDS(sprintf('%s/%s/full/maxent_fitted.rds', results_dir, x))
      
    } else {
      ## Get the background records from any source
      m = readRDS(sprintf('%s/%s/full/maxent_fitted.rds', results_dir, x))$me_full
      
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
    Number_background_points = m@results["X.Background.points",]
    Logistic_threshold       = m@results["X10.percentile.training.presence.Logistic.threshold",]
    Omission_rate            = m@results["X10.percentile.training.presence.training.omission",]
    
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
              Omission_rate,
              d)
    
    ## Just get the columns we need
    d = select(d, searchTaxon, Maxent_records, Number_var, Var_pcont,
               Per_cont, Var_pimp, Perm_imp, Training_AUC, Number_background_points,  
               Logistic_threshold, Omission_rate)
    
    dim(d)
    
    ## Remove path gunk, and species
    d$Species     = NULL
    return(d)
    
  }) %>%
  
  ## Finally, bind all the rows together
  bind_rows


#########################################################################################################################
## Now create a list of the '10th percentile training presence Logistic threshold'. This is used in step 8 to threshold
## the maps to just areas above the threshold.
summary(MAXENT.RESULTS["Logistic_threshold"])   
percent.10.log = as.list(MAXENT.RESULTS["Logistic_threshold"])  
percent.10.log = percent.10.log$Logistic_threshold


#########################################################################################################################
## Create a list of the omission files
omission.tables = list.files(results_dir, pattern = 'species_omission\\.csv$', full.names = TRUE, recursive = TRUE)


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
Max_tss  = as.data.frame(Max_tss)
setDT(Max_tss, keep.rownames = TRUE)[]
Max_tss  = as.data.frame(Max_tss)
names(Max_tss)[names(Max_tss) == 'rn'] <- 'searchTaxon'



## Remove extra text
Max_tss$searchTaxon = gsub("//",         "", Max_tss$searchTaxon)
Max_tss$searchTaxon = gsub(results_dir,   "", Max_tss$searchTaxon)
Max_tss$searchTaxon = gsub("/full/species_omission.csv.Max_tss", "", Max_tss$searchTaxon)
Max_tss$searchTaxon = gsub("/",          "", Max_tss$searchTaxon)
head(Max_tss)


## Add max TSS to the results table
MAXENT.RESULTS = join(MAXENT.RESULTS, Max_tss)
summary(MAXENT.RESULTS$Max_tss)
summary(MAXENT.RESULTS$Omission_rate)


## This is a summary of maxent output for current conditions
## Also which species have AUC < 0.7?
dim(MAXENT.RESULTS)
head(MAXENT.RESULTS, nrow(MAXENT.RESULTS))[1:5]





#########################################################################################################################
## 2). PLOT DASHBOARD OF MAXENT STATS
#########################################################################################################################



#########################################################################################################################
## Plot a miny dashboard of the results 
if (nrow(MAXENT.RESULTS) > 2) {
  
  
  ## What other variables are related?
  ## Omission rate vs records, key stats vs. maxent rating, etc. 
  lm.auc = lm(MAXENT.RESULTS$Max_tss ~ MAXENT.RESULTS$Training_AUC)
  lm.auc = lm(MAXENT.RESULTS$Max_tss ~ MAXENT.RESULTS$Training_AUC)
  
  
  ## Save this to file
  png(paste0('./output/maxent/', 'Maxent_run_summary_', save_run, '.png'), 16, 12, units = 'in', res = 500)
  
  layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
  
  plot(MAXENT.RESULTS$Training_AUC, MAXENT.RESULTS$Max_tss, pch = 19, col  = "blue",
       xlab = "AUC", ylab = "TSS", 
       abline(lm(MAXENT.RESULTS$Max_tss ~ MAXENT.RESULTS$Training_AUC)), 
       main = save_run, cex = 3)
  
  legend("topleft", bty = "n", 
         legend = paste("R2 is", format(summary(lm.auc)$adj.r.squared, digits = 4)))
  
  hist(MAXENT.RESULTS$Training_AUC, breaks = 10, col = "blue",   border = FALSE,
       ylab = "Frequency",
       xlab = "Training AUC", main = "AUC", cex = 3)
  hist(MAXENT.RESULTS$Max_tss,      breaks = 10, col = "orange", border = FALSE,
       ylab = "",
       xlab = "Maximum True Skill Statistic", main = "TSS", cex = 3)
  
  ## Finsish the device
  dev.off()
  
  ## If the species list is < 2 records, don't plot
} else {
  
  message('Dont plot, only ', length(GBIF.spp), ' species analysed')
  
}





#########################################################################################################################
## 3). COMBINE SDM RESULTS WITH HIA CONTEXT DATA
#########################################################################################################################


#########################################################################################################################
## Now check the match between the species list, and the results list. 
length(intersect(map_spp_list, MAXENT.RESULTS$searchTaxon)) 
MAXENT.RESULTS  =  MAXENT.RESULTS[MAXENT.RESULTS$searchTaxon %in% map_spp_list , ] 
map_spp         = unique(MAXENT.RESULTS$searchTaxon)
length(map_spp);setdiff(sort(map_spp_list), sort(map_spp))


## Change the species column
MAXENT.RESULTS$searchTaxon = gsub("_", " ", MAXENT.RESULTS$searchTaxon)


## Join the Maxent data to the niche data
MAXENT.RESULTS = join(COMBO.NICHE.CONTEXT, MAXENT.RESULTS, type = "right")


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
## OUTSTANDING SDM TASKS:
#########################################################################################################################


## 1). Check the code can be run through katana


########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################