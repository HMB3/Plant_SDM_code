#########################################################################################################################
###################### PREDICT MAXENT TO FUTURE CLIMATES AND SUMARISE THE RESULTS ####################################### 
#########################################################################################################################


#########################################################################################################################
## This code takes the output of the Maxent species distribution models (i.e. using current conditions), and generates 
## a prediction of habitat suitability for current and future environmental conditions. The input data table is in the 
## format of all species occurrences (rows) and environmental variables (columns).


## The predictions from all 6 GCMs are then combined into a single habitat suitability layer.
## And the total area of habitat gained, lost or remaining stable is calculated (i.e. for AUS). 


## Using this combined layer, the % of area occupied by species within areal units (significant urban areas or SUAs), 
## under each projection (2030, 2050 and 2070) is also be calculated. 


## Print the species run to the screen
message('Creating habitat suitability maps for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")
COMBO.NICHE.CONTEXT = readRDS(paste0('data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT_',  save_run, '.rds'))


#########################################################################################################################
## Read in niche data for table creation
# if(read_data == "TRUE") {
#   
#   ## Load GBIF and ALA data
#   COMBO.NICHE.CONTEXT = readRDS(paste0('data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT_',  save_run, '.rds'))
#   
# } else {
#   
#   message(' skip file reading, not many species analysed')   ##
#   
# }


#########################################################################################################################
## 1). READ IN THE CURRENT AND FUTURE ENVIRONMENTAL DATA FOR SIX GCMs
#########################################################################################################################


#########################################################################################################################
## Use a list of GCM scenarios to create maps of habitat suitability 
## Eight of the 40 CMIP5 models assessed in this project have been selected for use in provision of application-ready data. 
## This facilitates efficient exploration of climate projections for Australia.
## https://www.climatechangeinaustralia.gov.au/en/support-and-guidance/faqs/eight-climate-models-data/
head(gcms.50) ; head(gcms.70) ; head(gcms.30)



#########################################################################################################################
## 2). PROJECT MAXENT MODELS FOR MULTIPLE CLIMATE SCEANARIOS AT 2030, 2050 AND 2070
#########################################################################################################################


#########################################################################################################################
## For each species, use a function to create raster files and maps under all six GCMs at each time step


#########################################################################################################################
## Create 2030 maps
env.grids.2030 = tryCatch(project_maxent_grids(scen_list     = scen_2030,
                                               species_list  = map_spp_list,
                                               maxent_path   = maxent_path, 
                                               climate_path  = "./data/base/worldclim/aus/1km/bio",
                                               grid_names    = grid.names,
                                               time_slice    = 30,
                                               current_grids = aus.grids.current),
                          
                          ## Skip species
                          error = function(cond) {
                            
                            message(paste('Species skipped - check', spp))
                            
                          })


#########################################################################################################################
## Create 2050 maps
env.grids.2050 = tryCatch(project_maxent_grids(scen_list     = scen_2050,
                                               species_list  = map_spp_list,
                                               time_slice    = 50,
                                               maxent_path   = maxent_path,
                                               climate_path  = "./data/base/worldclim/aus/1km/bio",
                                               grid_names    = grid.names,
                                               current_grids = aus.grids.current),
                          
                          error = function(cond) {
                            
                            message(paste('Species skipped - check', spp))
                            
                          })


#########################################################################################################################
## Create 2070 maps
env.grids.2070 = tryCatch(project_maxent_grids(scen_list     = scen_2070,
                                               species_list  = map_spp_list,
                                               time_slice    = 70,
                                               maxent_path   = maxent_path,
                                               climate_path  = "./data/base/worldclim/aus/1km/bio",
                                               grid_names    = grid.names,
                                               current_grids = aus.grids.current),
                          
                          error = function(cond) {
                            
                            message(paste('Species skipped - check inputs', spp))
                            
                          })





#########################################################################################################################
## 3). CREATE TABLE OF MAXENT RESULTS FOR CURRENT CONDITIONS
#########################################################################################################################


#########################################################################################################################
## Create a file list for each model run: Try crunching this into just the species required
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
dim(subset(MAXENT.RESULTS, Training.AUC < 0.7))  ## all models should be above 0.7


## Are the TSS values ok?
lm.auc = lm(MAXENT.RESULTS$max_tss ~ MAXENT.RESULTS$Training.AUC)
hist(MAXENT.RESULTS$max_tss)
plot(MAXENT.RESULTS$Training.AUC, MAXENT.RESULTS$max_tss, pch = 19, col  = "blue",
    xlab = "AUC", ylab = "TSS", 
    abline(lm(MAXENT.RESULTS$max_tss ~ MAXENT.RESULTS$Training.AUC)))
legend("topleft", bty="n", legend=paste("R2 is", format(summary(lm.auc)$adj.r.squared, digits = 4)))


## Now check the match between the species list, and the results list. 
length(intersect(map_spp_list, MAXENT.RESULTS$searchTaxon)) 
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


#########################################################################################################################
## Join on the native data and the APNI
MAXENT.SUMMARY.NICHE <- SDM.TAXA %>% 
  
  join(., MAXENT.SUMMARY.NICHE) %>% 
  
  join(., MAXENT.CHECK[c("searchTaxon", "check.map")]) 

#View(MAXENT.SUMMARY.NICHE)

length(intersect(MAXENT.SUMMARY.NICHE$searchTaxon, GBIF.spp))


#########################################################################################################################
## Create a plot of number of records vs. maxent rating, Boxplot with no. occurrences on y, and maxent rating on
# plot(MAXENT.SUMMARY.NICHE$GLOBAL_RECORDS, MAXENT.SUMMARY.NICHE$check.map, pch = 19, col  = "blue",
#      xlab = "AUC", ylab = "TSS", 
#      abline(lm(MAXENT.SUMMARY.NICHE$check.map ~ MAXENT.SUMMARY.NICHE$GLOBAL_RECORDS)))
# legend("topleft", bty="n", legend=paste("R2 is", format(summary(lm.auc)$adj.r.squared, digits = 4)))


#########################################################################################################################
## Save maxent results
write.csv(MAXENT.SUMMARY.NICHE, paste0('output/maxent/MAXENT_SUMMARY_', save_run, '.csv'), row.names = FALSE)





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

percent.10.log         = as.list(MAXENT.RESULTS["X10.percentile.training.presence.Logistic.threshold"])  
percent.10.log         = percent.10.log$X10.percentile.training.presence.Logistic.threshold


percent.10.om          = as.list(MAXENT.RESULTS["X10.percentile.training.presence.training.omission"])   
percent.10.om          = percent.10.om$X10.percentile.training.presence.training.omission





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


# #load("SDM_COMBINE.RData")
# #save.image("STEP_8_MAPPING.RData")
# load("STEP_8_MAPPING.RData")


## Test problematic species by entering the values and running the function manually
## Common errors are due to corrupt rasters in the previous step
# DIR        = SDM.RESULTS.DIR[1]
# species    = map_spp[1]
# thresh     = percent.10.log[1]
# percent    = percent.10.om[1]
# time_slice = 30



#########################################################################################################################
## Reverse sort the lists
SDM.DIR.REV     = sort(SDM.RESULTS.DIR, decreasing = TRUE)
map_spp_rev     = sort(map_spp,         decreasing = TRUE) 
percent.log.rev = percent.10.log[sort(order(percent.10.log), decreasing = TRUE)]
percent.om.rev  = percent.10.om[sort(order(percent.10.om),   decreasing = TRUE)]


## Check the length matches - order should be correct, as well as the length
length(SDM.RESULTS.DIR);length(map_spp);length(percent.10.log);length(percent.10.om)
length(SDM.DIR.REV);length(map_spp_rev);length(percent.log.rev);length(percent.om.rev)

  
#########################################################################################################################
## Combine output and calculate gain and loss for 2030
suitability.2030 = tryCatch(mapply(SUA_cell_count,                                               ## combine_gcm_threshold
                                   DIR_list     = SDM.RESULTS.DIR,
                                   species_list = map_spp,
                                   maxent_path  = maxent_path,
                                   thresholds   = percent.10.log,
                                   percentiles  = percent.10.om,
                                   time_slice   = 30,
                                   write_rasters = FALSE),
                            
                            error = function(cond) {
                              
                              message(paste('Species skipped - check inputs', spp))
                              
                            })


## Combine GCM output for 2050 
suitability.2050 = tryCatch(mapply(SUA_cell_count,                                               
                                   DIR_list     = SDM.RESULTS.DIR,
                                   species_list = map_spp,
                                   maxent_path  = maxent_path,
                                   thresholds   = percent.10.log,
                                   percentiles  = percent.10.om,
                                   time_slice   = 50,
                                   write_rasters = FALSE),
                            
                            error = function(cond) {
                              
                              message(paste('Species skipped - check inputs', spp))
                              
                            })


## Combine GCM output for 2070 
suitability.2070 = tryCatch(mapply(SUA_cell_count,                                               
                                   DIR_list     = SDM.RESULTS.DIR,
                                   species_list = map_spp,
                                   maxent_path  = maxent_path,
                                   thresholds   = percent.10.log,
                                   percentiles  = percent.10.om,
                                   time_slice   = 70,
                                   write_rasters = FALSE),
                            
                            error = function(cond) {
                              
                              message(paste('Species skipped - check inputs', spp))
                              
                            })


#########################################################################################################################
## Reverse the order of the lists to speed up the output
# mapply(SUA_cell_count,
#        DIR_list     = SDM.DIR.REV,
#        species_list = map_spp_rev,
#        maxent_path  = maxent_path,
#        thresholds   = percent.log.rev,
#        percentiles  = percent.om.rev,
#        time_slice   = 70,  ## 50, 70
#        write_rasters = FALSE)





#########################################################################################################################
## 6). MOVE MAXENT FOLDERS TO NEW LOCATION FOR EACH RUN
#########################################################################################################################


## Create a list of folders for this run of species:EG hollow bearing species 
run_path         <- "./output/maxent/HOLLOW_SPP"
run_pat          <- map_spp
run_pat          <- paste(run_pat, sep = "", collapse = "|")
maxent_run_list  <- list.files(maxent_path, pattern = run_pat, full.names = TRUE)


## Copy or move these files to a specific folder
## Then you search for a file pattern in that directory
## This is very slow, would be better done in unix, etc. But 
file.copy(from      = maxent_run_list, 
          to        = run_path, 
          overwrite = FALSE, 
          recursive = TRUE, 
          copy.mode = TRUE)


#########################################################################################################################
## OUTSTANDING MAPPING TASKS:
#########################################################################################################################


## 1). 


## 2). Fix the species mapping loop over two lists at once : mapply




#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################