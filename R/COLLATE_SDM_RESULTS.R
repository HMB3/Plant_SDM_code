#########################################################################################################################
################################################ COLLATE SDM RESULTS  ###################################################
#########################################################################################################################


## This code takes all the randomised results - derrived from 1000 realisations of the site * spp matrix - and calculates 
## the mean, standard deviation and variance. These stats are needed to calculate the staistic Shawn was talking about
## the other day.


#########################################################################################################################
## 1). CREATE SUMMARY TABLE OF MAXENT RESULTS
#########################################################################################################################


#########################################################################################################################
## Read in the list of files for the standard variables, and specify the file path
table.list = list.files("./output/maxent/STD_VAR_ALL/")
path = "./output/maxent/STD_VAR_ALL/"


## Could turn this into a function, and loop over a list of subfolders...
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


## If there are multiple model runs, how would we compare them? Could have a row for each species and model run? 
dim(MAXENT.STD.VAR.SUMMARY)
head(MAXENT.STD.VAR.SUMMARY)[1:8]


#########################################################################################################################
## Save results
write.csv(MAXENT.STD.VAR.SUMMARY, "./output/maxent/MAXENT_STD_VAR_SUMMARY.csv", row.names = FALSE)





#########################################################################################################################
## OUTSTANDING RESULTS TASKS:
#########################################################################################################################


## What are the most sensible ways to compare model fits (i.e. what do good and bad models look like)? 

## Can we calculate TSS from the existing Maxent code?

## We should record global settings of the SDMs too - e.g. variable correlation, background points, cross vaildation settings, etc.

## What is the best way to summarise the results into a final table - which results columns, and which format?






#########################################################################################################################
################################################# TBC  ################################################################## 
#########################################################################################################################