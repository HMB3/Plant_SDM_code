#########################################################################################################################
################################################ COLLATE SDM RESULTS  ###################################################
#########################################################################################################################


## This code takes all the randomised results - derrived from 1000 realisations of the site * spp matrix - and calculates 
## the mean, standard deviation and variance. These stats are needed to calculate the staistic Shawn was talking about
## the other day.


#########################################################################################################################
## 1). Read in observed results, and get contextual columns to join on
#########################################################################################################################


## Look at the output
RESULTS.EG = read.csv("./output/maxent/SEL_VAR_ALL/Yucca_elephantipes/full/maxentResults.csv", stringsAsFactors = FALSE)
dim(RESULTS.EG)
names(RESULTS.EG)


## Read in the list of 
std.var.all   = list.files("./output/maxent/STD_VAR_ALL/")
sel.var.all   = list.files("./output/maxent/SEL_VAR_ALL/")




#########################################################################################################################
## 2). CREATE SUMMARY TABLE OF MAXENT RESULTS
#########################################################################################################################


## Loop over a list of subfolders
table.list = std.var.all [1:2]
path = "./output/maxent/STD_VAR_ALL/"


## Could turn this into a function
MAXENT.SUMMARY <- table.list[c(1:length(table.list))] %>%
  
  ## pipe the list into lapply
  lapply(function(x) {
    
    ## create the character string
    f <- paste0(path, x, "/full/maxentResults.csv")
    
    ## load each .RData file
    d <- read.csv(f)
    
    ## now add a model column
    d = cbind(GBIF_Taxon = table.list, Model_run  = path, d) 
    
    ## Remove path gunk, and species
    d$GBIF_Taxon = gsub("_", " ", d$GBIF_Taxon)
    d$Model_run  = gsub("./output/maxent/", "", d$Model_run)
    d$Model_run  = gsub("/", "", d$Model_run)
    head(d)[1:8]
    
    ## Not sure why it is not returning only one row per species...
    d = d[!duplicated(d$GBIF_Taxon), ]
    
    ## need to return d from the function
    return(d)

  
  }) %>%
  
  ## finally, bind all the rows together
  bind_rows


## Not sure why it is not returning only one row per species...
MAXENT.SUMMARY = MAXENT.SUMMARY[!duplicated(MAXENT.SUMMARY$GBIF_Taxon), ]
head(MAXENT.SUMMARY)[1:8]
View(MAXENT.SUMMARY)

## Run as a function?
# MAXENT.SUM.STD.VAR.ALL = read_bind_tables(table.list = sel.var.all[1:2], 
#                                           path = "./output/maxent/STD_VAR_ALL/")




#########################################################################################################################
################################################# TBC  ################################################################## 
#########################################################################################################################