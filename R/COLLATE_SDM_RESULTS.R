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


## Read in the list of 
std.var.all   = list.files("./output/maxent/STD_VAR_ALL/")
sel.var.all   = list.files("./output/maxent/SEL_VAR_ALL/")




#########################################################################################################################
## 3). RE-CREATE FIGURES 4 and 5
#########################################################################################################################


## Loop over a list of subfolders
read_bind_maxent = function (table.list, path) {
  
  READ.BIND.TABLE <- table.list[c(1:length(table.list))] %>% 
    
    ## pipe the list into lapply
    lapply(function(x) {
      
      ## create the character string
      f <- paste0(path, x, "/full/maxentResults.csv")
      
      ## read each .csv file
      d <- read.csv(f)
      
      ## now add a model column
      cbind(GBIF_Taxon = x,
            Model_run  = path, 
            d)
      
      ## Remove path gunk, and species
      d$GBIF_Taxon = gsub("_", " ", d$GBIF_Taxon)
      d$Model_run  = gsub("./output/maxent/", "", d$Model_run)
      d$Model_run  = gsub("/", "", d$Model_run)
      
      
    }) %>% 
    
    ## finally, bind all the rows together
    bind_rows
  
}


## Plot the observed coeficient and t-value over the random ones? First, get the values for the histogram
## T2 random histos 
MAXENT.SUM.STD.VAR.ALL = read_bind_tables(table.list = sel.var.all, 
                                          path = "./output/maxent/STD_VAR_ALL/")





## Check: 
dim(T2.RANDOM.HISTO)   ## 18 results rows for each table *1000 randomisations = 18000 rows







## Then, what is the easiest way to re-create the figures? Probably avoid all the other code, and just plot the randomised
## results on the right panel: 

## 1). Take the existing figure in Arcmap and just sub in the new code. Code to make the panels starts on line 400 of:
## FIGURE_4_MEDIAN_NICHES_GPP_DELTA_20_SUBSETS_SEPARATE_NICHE_PANELS.R

## 2). Or sub th new code into



#########################################################################################################################
################################################# END  ################################################################## 
#########################################################################################################################