#########################################################################################################################
## 3). CREATE TABLE OF MAXENT RESULTS FOR CURRENT CONDITIONS
#########################################################################################################################


## Update this with just the native species with good models. This is then joined on at


#########################################################################################################################
## Create a file list for each model run: Try crunching this into just the species required
maxent.tables = list.files(maxent_path)                 
maxent.tables = intersect(maxent.tables, map_spp_list)   
maxent_path   = maxent_path                             
length(maxent.tables)                                                          ## Should match the number of taxa tested


## Create a table of the results 
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
dim(subset(MAXENT.RESULTS, Training.AUC < 0.7))  ## all models should be above 0.7


## Are the TSS values ok?
lm.auc = lm(MAXENT.RESULTS$max_tss ~ MAXENT.RESULTS$Training.AUC)
#hist(MAXENT.RESULTS$max_tss)
plot(MAXENT.RESULTS$Training.AUC, MAXENT.RESULTS$max_tss, pch = 19, col  = "blue",
     xlab = "Area Under the ROC Curve", ylab = "True Skill Statistic", 
     abline(lm(MAXENT.RESULTS$max_tss ~ MAXENT.RESULTS$Training.AUC)))
legend("topleft", bty="n", legend=paste("R2 is", format(summary(lm.auc)$adj.r.squared, digits = 4)))


## Now check the match between the species list, and the results list. These need to match, so we can access
## the right threshold for each species.
length(intersect(map_spp_list, MAXENT.RESULTS$searchTaxon)) ## Accesssing the files from these directories... 
MAXENT.RESULTS  =  MAXENT.RESULTS[MAXENT.RESULTS$searchTaxon %in% map_spp_list , ] 
map_spp = unique(MAXENT.RESULTS$searchTaxon)
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
length(intersect(MAXENT.SUMMARY.NICHE$searchTaxon, APNI$searchTaxon))


SDM.TAXA <- MAXENT.SUMMARY.NICHE[["searchTaxon"]] %>%  
  lookup_table(., by_species = TRUE) 
SDM.TAXA <- setDT(SDM.TAXA, keep.rownames = TRUE)[]
colnames(SDM.TAXA)[1] <- "searchTaxon"
SDM.TAXA <- join(SDM.TAXA, APNI)

# 
# MAXENT.SUMMARY.NICHE = join(SDM.TAXA, MAXENT.SUMMARY.NICHE)
# View(MAXENT.SUMMARY.NICHE)


#########################################################################################################################
## Join on the native data and the APNI
MAXENT.SUMMARY.NICHE <- SDM.TAXA %>% 
  
  join(., MAXENT.SUMMARY.NICHE) %>% 
  
  join(., MAXENT.CHECK[c("searchTaxon", "check.map")]) 

View(MAXENT.SUMMARY.NICHE)


#View(MAXENT.SUMMARY.TABLE)
length(intersect(MAXENT.SUMMARY.NICHE$searchTaxon, GBIF.spp))


#########################################################################################################################
## Save - could add date as a sprintf variable to save multiple versions?
## write.csv(MAXENT.SUMMARY.NICHE, paste0('output/maxent/MAXENT_SUMMARY_', save_run, '.csv'), row.names = FALSE)
## MAXENT.SUMMARY.NICHE = read.csv(paste0('output/maxent/MAXENT_SUMMARY_', save_run, '.csv'), stringsAsFactors = FALSE)





#########################################################################################################################
## 4). CREATE LISTS OF HABITAT SUITABILITY THRESHOLDS
#########################################################################################################################


#########################################################################################################################
## Maxent produces a presence threshold for each species (i.e. the columns in MAXENT.RESULTS). 
## The trouble here is that we might need to change the threshold for different species, rather than using the same one 
## for all of them. That changes the order of lists, which is a problem for looping over them.


## John : for AUC you can report the cross-validated test AUC (if your code currently runs a cross-validated model as well), 
## and for the model threshold (for binarising) you can just use the training value (or the crossval one...there's little 
## guidance about this and you can really get away with either).


## How do the differnt thresholds compare for the set of species modelled?
summary(MAXENT.RESULTS["Maximum.training.sensitivity.plus.specificity.Logistic.threshold"])    ## The strictest threshold
summary(MAXENT.RESULTS["X10.percentile.training.presence.Logistic.threshold"])                 ## The next strictest
summary(MAXENT.RESULTS["X10.percentile.training.presence.training.omission"])                  ## The most forgiving


#########################################################################################################################
## Turn the maxent results into lists :: we can use these to generate the consensus layers 
thresh.max.train       = as.list(MAXENT.RESULTS["Maximum.training.sensitivity.plus.specificity.Logistic.threshold"]) 
thresh.max.train       = thresh.max.train$Maximum.training.sensitivity.plus.specificity.Logistic.threshold

percent.10.log         = as.list(MAXENT.RESULTS["X10.percentile.training.presence.Logistic.threshold"])  ## for the twos 
percent.10.log         = percent.10.log$X10.percentile.training.presence.Logistic.threshold


percent.10.om          = as.list(MAXENT.RESULTS["X10.percentile.training.presence.training.omission"])   ## discount
percent.10.om          = percent.10.om$X10.percentile.training.presence.training.omission


## Check that the order of the species names directories and threshold values match before final run
head(SDM.RESULTS.DIR, 20);head(map_spp, 20); head(MAXENT.RESULTS, 20)[, c("searchTaxon",
                                                                          "Maximum.training.sensitivity.plus.specificity.Logistic.threshold", 
                                                                          "X10.percentile.training.presence.Logistic.threshold")]

tail(SDM.RESULTS.DIR, 20);tail(map_spp, 20); tail(MAXENT.RESULTS, 20)[, c("searchTaxon",
                                                                          "Maximum.training.sensitivity.plus.specificity.Logistic.threshold", 
                                                                          "X10.percentile.training.presence.Logistic.threshold")]



## Save this as an object
#save.image("MAPPING_DATA_NEW_ALA.RData")





#########################################################################################################################
## 5). SUMARIZE MAXENT RESULTS FOR EACH SPECIES ACROSS MULTIPLE GCMs
#########################################################################################################################


#########################################################################################################################
## So we can use the values in these columns to threhsold the rasters of habitat suitability (0-1) when combining them.
## For each species, we will create a binary raster with cell values between 0-6. These cell values represent the number of GCMs 
## where that cell had a suitability value above the threshold determined by maxent (e.g. the max training sensitivity + specif).


#########################################################################################################################
## Then use the more forgiving threshold
#########################################################################################################################


## Test problematic species by entering the values and running the function manually
## Common errors are due to corrupt rasters in the previous step
DIR        = SDM.RESULTS.DIR[1]
species    = map_spp[1]
thresh     = percent.10.log[1]
percent    = percent.10.om[1]
time_slice = 30



#########################################################################################################################
## Reverse sort the lists
SDM.DIR.REV     = sort(SDM.RESULTS.DIR, decreasing = TRUE)
map_spp_rev     = sort(map_spp,         decreasing = TRUE) 
percent.log.rev = percent.10.log[sort(order(percent.10.log), decreasing = TRUE)]
percent.om.rev  = percent.10.om[sort(order(percent.10.om),   decreasing = TRUE)]


## Check the length matches - order should be correct, as well as the length
length(SDM.RESULTS.DIR);length(map_spp);length(percent.10.log);length(percent.10.om)
length(SDM.DIR.REV);length(map_spp_rev);length(percent.log.rev);length(percent.om.rev)