#########################################################################################################################
## 3). SDM FOR ALL ENV VARIABLES, ALL RECORDS
#########################################################################################################################


## Fit maxent function needs to be updated to do the variable selection


########################################################################################################################
## We can run Maxent from a cluster of cores on the local computer. Here we send (i.e. export) all the necessary ingredients 
## to the cluster. So that's the:


## template.raster raster, 
## data frame and 
## maxent function


## the simplify function rmaxent::simplify


## 100 species takes about 4 hours...
cl <- makeCluster(6)
clusterExport(cl, c('template.raster', 'SDM.DATA.ALL', 'FIT_MAXENT'))
clusterEvalQ(cl, {
  
  require(ff)
  require(rgeos)
  require(sp)
  require(raster)
  require(rJava)
  require(dismo)
  require(things)
  
})


## Now use 'lapply' to run maxent for multiple species
## Note that running the code in parallel causes problems
lapply(test.spp[1:length(test.spp)], function(x) { # for serial, parLapply(cl, species[1:8], function(x) { # for parallel 
  
  ## Print the taxa being processed to screen
  message('Doing ', x)
  
  ## Subset the records to only the taxa being processed
  occurrence <- subset(SDM.DATA.ALL, searchTaxon == x)
  
  ## Now get the background points. These can come from anywhere in the whole dataset,
  ## other than the species used.
  background <- subset(SDM.DATA.ALL, searchTaxon != x)
  
  ## The create a vector of the sdm.predictors used. 
  ## This should be based on an ecological framework! 
  sdm.predictors <- sdm.predictors.all # vector of used sdm.predictors 
  
  ## Finally fit the models using FIT_MAXENT
  ## There is no switch in the function to skip outputs that exist.
  ## Given all the changes likely to be made to the models, this could be wise...
  FIT_MAXENT(occ                     = occurrence, 
             bg                      = background, 
             sdm.predictors          = sdm.predictors.all, 
             name                    = x, 
             outdir                  = 'output/maxent/SEL_VAR_ALL', 
             template.raster,
             min_n                   = 20,   ## This should be higher...
             max_bg_size             = 100000,
             background_buffer_width = 200000,
             shapefiles              = TRUE,
             features                = 'lpq',
             replicates              = 5,
             responsecurves          = TRUE)
  
})


stopCluster(cl)



#########################################################################################################################
## 4). SDM FOR ALL VARIABLES USING CULTIVATED RECORDS
#########################################################################################################################


## Fit maxent function needs to be update to do the variable selection


########################################################################################################################
## We can run Maxent from a cluster of cores on the local computer. Here we send (i.e. export) all the necessary ingredients 
## to the cluster. So that's the:


## template.raster raster, 
## data frame and 
## maxent function


## 100 species takes about 4 hours...
cl <- makeCluster(6)
clusterExport(cl, c('template.raster', 'SDM.DATA.ALL', 'FIT_MAXENT'))
clusterEvalQ(cl, {
  
  require(ff)
  require(rgeos)
  require(sp)
  require(raster)
  require(rJava)
  require(dismo)
  require(things)
  
})


## Now use 'lapply' to run maxent for multiple species


## This code would be exactly the same, except for the cultivated/non cultivated output
lapply(test.spp[1:length(test.spp)], function(x) { # for serial, parLapply(cl, species[1:8], function(x) { # for parallel 
  
  ## Print the taxa being processed to screen
  message('Doing ', x)
  
  ## Subset the records to only the taxa being processed
  ## Here, add condition for cultivated
  ## occurrence <- subset(SDM.DATA.ALL, searchTaxon == x & CULTIVATED == "CULTIVATED")
  occurrence <- subset(SDM.DATA.ALL, searchTaxon == x & CULTIVATED == "CULTIVATED")
  
  ## Now get the background points. These can come from anywhere in the whole dataset,
  ## other than the species used.
  background <- subset(SDM.DATA.ALL, searchTaxon != x & CULTIVATED == "CULTIVATED")
  
  ## The create a vector of the sdm.predictors used. 
  ## This should be based on an ecological framework! 
  sdm.predictors <- sdm.predictors.all # vector of used sdm.predictors 
  
  ## Finally fit the models using FIT_MAXENT
  ## There is no switch in the function to skip outputs that exist.
  ## Given all the changes likely to be made to the models, this could be wise...
  FIT_MAXENT(occ                     = occurrence, 
             bg                      = background, 
             sdm.predictors          = sdm.predictors.all, 
             name                    = x, 
             outdir                  = 'output/maxent/SEL_VAR_CULT', 
             template.raster,
             min_n                   = 20,   ## This should be higher...
             max_bg_size             = 100000,
             background_buffer_width = 200000,
             shapefiles              = TRUE,
             features                = 'lpq',
             replicates              = 5,
             responsecurves          = TRUE)
  
})


stopCluster(cl)





#########################################################################################################################
## 5). SDM FOR ALL VARIABLES USING UN-CULTIVATED RECORDS
#########################################################################################################################


## Fit maxent function needs to be updated to do the variable selection


########################################################################################################################
## We can run Maxent from a cluster of cores on the local computer. Here we send (i.e. export) all the necessary ingredients 
## to the cluster. So that's the:


## template.raster raster, 
## data frame and 
## maxent function


## 100 species takes about 4 hours...
cl <- makeCluster(6)
clusterExport(cl, c('template.raster', 'SDM.DATA.ALL', 'FIT_MAXENT'))
clusterEvalQ(cl, {
  
  require(ff)
  require(rgeos)
  require(sp)
  require(raster)
  require(rJava)
  require(dismo)
  require(things)
  
})


## Now use 'lapply' to run maxent for multiple species


## This code would be exactly the same, except for the cultivated/non cultivated output
lapply(test.spp[1:length(test.spp)], function(x) { # for serial, parLapply(cl, species[1:8], function(x) { # for parallel 
  
  ## Print the taxa being processed to screen
  message('Doing ', x)
  
  ## Subset the records to only the taxa being processed
  ## Here, add condition for cultivated
  ## occurrence <- subset(SDM.DATA.ALL, searchTaxon == x & CULTIVATED == "CULTIVATED")
  occurrence <- subset(SDM.DATA.ALL, searchTaxon == x & CULTIVATED == "UNKNOWN")
  
  ## Now get the background points. These can come from anywhere in the whole dataset,
  ## other than the species used.
  background <- subset(SDM.DATA.ALL, searchTaxon != x & CULTIVATED == "UNKNOWN")
  
  ## The create a vector of the sdm.predictors used. 
  ## This should be based on an ecological framework! 
  sdm.predictors <- sdm.predictors.all # vector of used sdm.predictors 
  
  ## Finally fit the models using FIT_MAXENT
  ## There is no switch in the function to skip outputs that exist.
  ## Given all the changes likely to be made to the models, this could be wise...
  FIT_MAXENT(occ                     = occurrence, 
             bg                      = background, 
             sdm.predictors          = sdm.predictors.all, 
             name                    = x, 
             outdir                  = 'output/maxent/SEL_VAR_CULT', 
             template.raster,
             min_n                   = 20,   ## This should be higher...
             max_bg_size             = 100000,
             background_buffer_width = 200000,
             shapefiles              = TRUE,
             features                = 'lpq',
             replicates              = 5,
             responsecurves          = TRUE)
  
})


stopCluster(cl)





#########################################################################################################################
## 6). SDM FOR SELECTED VARIABLES USING CULTIVATED RECORDS
#########################################################################################################################


## Fit maxent function needs to be update to do the variable selection


########################################################################################################################
## We can run Maxent from a cluster of cores on the local computer. Here we send (i.e. export) all the necessary ingredients 
## to the cluster. So that's the:


## template.raster raster, 
## data frame and 
## maxent function


## 100 species takes about 4 hours...
cl <- makeCluster(6)
clusterExport(cl, c('template.raster', 'SDM.DATA', 'FIT_MAXENT'))
clusterEvalQ(cl, {
  
  require(ff)
  require(rgeos)
  require(sp)
  require(raster)
  require(rJava)
  require(dismo)
  require(things)
  
})


## Now use 'lapply' to run maxent for multiple species


## This code would be exactly the same, except for the cultivated/non cultivated output
lapply(test.spp[1:length(test.spp)], function(x) { # for serial, parLapply(cl, species[1:8], function(x) { # for parallel 
  
  ## Print the taxa being processed to screen
  message('Doing ', x)
  
  ## Subset the records to only the taxa being processed
  ## Here, add condition for cultivated
  ## occurrence <- subset(SDM.DATA, searchTaxon == x & CULTIVATED == "CULTIVATED")
  occurrence <- subset(SDM.DATA, searchTaxon == x & CULTIVATED == "CULTIVATED")
  
  ## Now get the background points. These can come from anywhere in the whole dataset,
  ## other than the species used.
  background <- subset(SDM.DATA, searchTaxon != x & CULTIVATED == "CULTIVATED")
  
  ## The create a vector of the sdm.predictors used. 
  ## This should be based on an ecological framework! 
  sdm.predictors <- sdm.predictors # vector of used sdm.predictors 
  
  ## Finally fit the models using FIT_MAXENT
  ## There is no switch in the function to skip outputs that exist.
  ## Given all the changes likely to be made to the models, this could be wise...
  FIT_MAXENT(occ                     = occurrence, 
             bg                      = background, 
             sdm.predictors          = sdm.predictors, 
             name                    = x, 
             outdir                  = 'output/maxent/STD_VAR_CULT', 
             template.raster,
             min_n                   = 20,   ## This should be higher...
             max_bg_size             = 100000,
             background_buffer_width = 200000,
             shapefiles              = TRUE,
             features                = 'lpq',
             replicates              = 5,
             responsecurves          = TRUE)
  
})


stopCluster(cl)





#########################################################################################################################
## 7). RUN SDM FOR ALL VARIABLES USING UN-CULTIVATED RECORDS
#########################################################################################################################


## Fit maxent function needs to be updated to do the variable selection


########################################################################################################################
## We can run Maxent from a cluster of cores on the local computer. Here we send (i.e. export) all the necessary ingredients 
## to the cluster. So that's the:


## template.raster raster, 
## data frame and 
## maxent function


## 100 species takes about 4 hours...
cl <- makeCluster(6)
clusterExport(cl, c('template.raster', 'SDM.DATA', 'FIT_MAXENT'))
clusterEvalQ(cl, {
  
  require(ff)
  require(rgeos)
  require(sp)
  require(raster)
  require(rJava)
  require(dismo)
  require(things)
  
})


## Now use 'lapply' to run maxent for multiple species


## This code would be exactly the same, except for the cultivated/non cultivated output
lapply(test.spp[1:length(test.spp)], function(x) { # for serial, parLapply(cl, species[1:8], function(x) { # for parallel 
  
  ## Print the taxa being processed to screen
  message('Doing ', x)
  
  ## Subset the records to only the taxa being processed
  ## Here, add condition for cultivated
  ## occurrence <- subset(SDM.DATA, searchTaxon == x & CULTIVATED == "CULTIVATED")
  occurrence <- subset(SDM.DATA, searchTaxon == x & CULTIVATED == "UNKNOWN")
  
  ## Now get the background points. These can come from anywhere in the whole dataset,
  ## other than the species used.
  background <- subset(SDM.DATA, searchTaxon != x & CULTIVATED == "UNKNOWN")
  
  ## The create a vector of the sdm.predictors used. 
  ## This should be based on an ecological framework! 
  sdm.predictors <- sdm.predictors # vector of used sdm.predictors 
  
  ## Finally fit the models using FIT_MAXENT
  ## There is no switch in the function to skip outputs that exist.
  ## Given all the changes likely to be made to the models, this could be wise...
  FIT_MAXENT(occ                     = occurrence, 
             bg                      = background, 
             sdm.predictors          = sdm.predictors, 
             name                    = x, 
             outdir                  = 'output/maxent/STD_VAR_UNCULT', 
             template.raster,
             min_n                   = 20,   ## This should be higher...
             max_bg_size             = 100000,
             background_buffer_width = 200000,
             shapefiles              = TRUE,
             features                = 'lpq',
             replicates              = 5,
             responsecurves          = TRUE)
  
})


stopCluster(cl)




#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################