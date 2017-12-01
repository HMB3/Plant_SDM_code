#########################################################################################################################
################################################ COLLATE SDM RESULTS  ###################################################
#########################################################################################################################


## This code takes all the randomised results - derrived from 1000 realisations of the site * spp matrix - and calculates 
## the mean, standard deviation and variance. These stats are needed to calculate the staistic Shawn was talking about
## the other day.


#########################################################################################################################
## 1). Read in observed results, and get contextual columns to join on
#########################################################################################################################


## Read in the list of files for the standard variables
std.var.all = list.files("./output/maxent/STD_VAR_ALL/")
table.list = std.var.all
path = "./output/maxent/STD_VAR_ALL/"


## All variables...
sel.var.all   = list.files("./output/maxent/SEL_VAR_ALL/")





#########################################################################################################################
## 2). CREATE SUMMARY TABLE OF MAXENT RESULTS
#########################################################################################################################


## Loop over a list of subfolders...Could turn this into a function
MAXENT.STD.VAR.SUMMARY <- table.list[c(1:length(table.list))] %>%
  
  ## pipe the list into lapply
  lapply(function(x) {
    
    ## create the character string
    f <- paste0(path, x, "/full/maxentResults.csv")
    
    ## load each .RData file
    d <- read.csv(f)
    
    ## now add a model column
    d = cbind(GBIF_Taxon = x, Model_run  = path, d) 
    dim(d)
    
    ## Remove path gunk, and species
    d$GBIF_Taxon = gsub("_", " ", d$GBIF_Taxon)
    d$Model_run  = gsub("./output/maxent/", "", d$Model_run)
    d$Model_run  = gsub("/", "", d$Model_run)
    d$Species    = NULL
    d
  }) %>%
  
  ## finally, bind all the rows together
  bind_rows


## Not sure why it is not returning only one row per species...
dim(MAXENT.STD.VAR.SUMMARY)
names(MAXENT.STD.VAR.SUMMARY)
head(MAXENT.STD.VAR.SUMMARY)[1:8]
View(MAXENT.STD.VAR.SUMMARY)

## Run as a function?
# MAXENT.SUM.STD.VAR.ALL = read_bind_tables(table.list = sel.var.all[1:2], 
#                                           path = "./output/maxent/STD_VAR_ALL/")


#########################################################################################################################
## Save results
write.csv(MAXENT.STD.VAR.SUMMARY, "./output/maxent/MAXENT_STD_VAR_SUMMARY.csv", row.names = FALSE)


#########################################################################################################################
## OUTSTANDING SDM RESULTS TASKS:
#########################################################################################################################


## Can we compare the model runs? 






#########################################################################################################################
################################################# TBC  ################################################################## 
#########################################################################################################################