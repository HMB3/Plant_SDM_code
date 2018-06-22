#########################################################################################################################
################################################# THIN RECORDS ##########################################################
#########################################################################################################################


## This code spatially thins species occurrence records to help address problems associated with spatial sampling biases. 
## Ideally, thinning removes the fewest records necessary to substantially reduce the effects of sampling bias, while 
## simultaneously retaining the greatest amount of useful information. See: https://cran.r-project.org/web/packages/spThin/vignettes/spThin_vignette.html


## Create lists
source('./R/HIA_LIST_MATCHING.R')


#########################################################################################################################
## 1). CHECK DATA FOR AN EXAMPLE SPECIES...
#########################################################################################################################


## Load GBIF data
COMBO.RASTER.CONTEXT = readRDS("./data/base/HIA_LIST/COMBO/CLEAN_ONLY_HIA_SPP.rds")
SPP.BIAS             = read.csv("./output/maxent/SPP_BOUNDARY_BIAS.csv", stringsAsFactors = FALSE)
SPP.BIAS             = subset(SPP.BIAS, AUS_BOUND_BIAS == "TRUE")$searchTaxon


## Try for an example species
SPP.TEST            = subset(COMBO.RASTER.CONTEXT, searchTaxon == "Acacia floribunda") 
dim(SPP.TEST)


## HB hard to see - but we'd love to remove the point in Madagascar 
plot(LAND)
points(SPP.TEST$lon,  SPP.TEST$lat, cex = 0.2, col = "red", pch = 19)


#########################################################################################################################
## We can see the boundary problem in QLD
plot(aus)
points(SPP.TEST$lon,  SPP.TEST$lat, cex = 0.8, col = "red", pch = 19)


## Run spThin::thin on the full dataset. thin involves multiple settings. This allows for extensive flexibility in how 
## the user spatially thins a dataset. However, many have default values. See ?thin for further information.
## This is really slow...
thinned_dataset_full <-
  thin(loc.data                 = SPP.TEST, 
       lat.col                  = "lat", 
       long.col                 = "lon", 
       spec.col                 = "searchTaxon", 
       thin.par                 = 100, 
       reps                     = 1, 
       locs.thinned.list.return = TRUE, 
       write.files              = TRUE, 
       max.files                = 5, 
       out.dir                  = "data/base/HIA_LIST/COMBO/ficus_thinned_full/", 
       out.base                 = "ficus_thinned", 
       write.log.file           = TRUE,
       log.file                 = "hanomalus_thinned_full_log_file.txt",
       verbose                  = TRUE)



## In the case above, we found that 10 repetitions were sufficient to return spatially thinned datasets with the 
## optimal number of occurrence records (124). Because this is a random process, it is possible that a similarly 
## repeated run would not return any datasets with the optimal number of occurrence records. To visually assess 
## whether we are using enough reps to approach the optimal number we use the function plotThin, This function 
## produces three plots: 1) the cumulative number of records retained versus the number of repetitions, 2) the 
## log cumulative number of records retained versus the log number of repetitions, and 3) a histogram of the 
## maximum number of records retained for each thinned dataset.





#########################################################################################################################
## FLAG GBIF AND SPATIAL OUTLIERS FOR ONE SPECIES
#########################################################################################################################








#########################################################################################################################
## 2). NOW, CLEAN WHOLE DATASET
#########################################################################################################################


#########################################################################################################################
## FLAG GBIF OUTLIERS
#########################################################################################################################



#########################################################################################################################
## 3). FLAG SPATIAL OUTLIERS
#########################################################################################################################





#########################################################################################################################
## Join data :: exclude the decimal lat/long, check the length 


## Check the values for each flag
# summary(TEST.GEO$validity)
# summary(TEST.GEO$equal)
# summary(TEST.GEO$zeros)
# summary(TEST.GEO$capitals)
# summary(TEST.GEO$centroids)
# summary(TEST.GEO$gbif)
# summary(TEST.GEO$institution)
# summary(TEST.GEO$gbif)
# summary(TEST.GEO$summary)
#summary(TEST.GEO$GBIF.SPAT.OUT)





#########################################################################################################################
## 4). SUBSET FOR THE NICHES : JUST USE COORDCLEAN SUMMARY
#########################################################################################################################


## So ~2.6% of the data is dodgy according to the GBIF fields or spatial outliers
## This seems ok as a median figure across the data set?
#dim(subset(TEST.GEO, summary == "FALSE" | GBIF.SPAT.OUT == "FALSE"))[1]/dim(TEST.GEO)[1]*100


#########################################################################################################################
## RE-CREATE NICHES
#########################################################################################################################


#########################################################################################################################
## 5). INTERSECT SPECIES RECORDS WITH LOCAL GOV AREAS AND SIGNIFICANT URBAN AREAS
#########################################################################################################################


#########################################################################################################################
## See the ABS for details :: there are 563 LGAs
## http://www.abs.gov.au/ausstats/abs@.nsf/Lookup/by%20Subject/1270.0.55.003~July%202016~Main%20Features~Local%20Government%20Areas%20(LGA)~7


## Remove the columns we don't need
LGA.WGS = LGA.WGS[, c("LGA_CODE16", "LGA_NAME16")] 



#########################################################################################################################
## Save
saveRDS(TEST.GEO,                'data/base/HIA_LIST/COMBO/CLEAN_FLAGS_HIA_SPP.rds')
saveRDS(CLEAN.TRUE,              'data/base/HIA_LIST/COMBO/CLEAN_ONLY_HIA_SPP.rds')
saveRDS(CLEAN.NICHE.CONTEXT,     'data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT_APRIL_2018_COORD_CLEAN.rds')
write.csv(CLEAN.NICHE.CONTEXT,   "./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT_APRIL_2018_COORD_CLEAN.csv", row.names = FALSE)




#########################################################################################################################
#####################################################  TBC ##############################################################
#########################################################################################################################