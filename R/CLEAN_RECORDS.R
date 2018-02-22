#########################################################################################################################
################################################ CLEAN RECORDS ########################################################## 
#########################################################################################################################


## This code creates individual shapefiles for each species, so we can create a list of the spatial outliers manually


#########################################################################################################################
## 1). CREATE UNIQUE ID FIELD FOR EACH RECORD
#########################################################################################################################


## This column can be applied across the project  
dim(GBIF.ALA.COMBO.HIA)
GBIF.ALA.COMBO.HIA$OBS <- 1:nrow(GBIF.ALA.COMBO.HIA)
dim(GBIF.ALA.COMBO.HIA)[1];length(GBIF.ALA.COMBO.HIA$OBS)                                          ## This matches step 6


#########################################################################################################################
## Check the field names :: are ALA and GBIF ID fields mutually exclusive?
test_exotic = subset(GBIF.ALA.COMBO.HIA, searchTaxon == "Cycas revoluta")
test_native = subset(GBIF.ALA.COMBO.HIA, searchTaxon == "Syzygium floribundum")

View(head(test_exotic, 162)[, c("searchTaxon",
                                "lon",
                                "lat",
                                "gbifID", 
                                "id")])

View(head(test_native, 521)[, c("searchTaxon",
                                "lon",
                                "lat",
                                "gbifID", 
                                "id")])


## So, excluding the outliers is not as simple to just exclude the GBIF records for native species...
## This will work for some taxa, but for others in could lead to equal bias in the model by further underestimating
## the fundamental niche.

## So the preferable option is to either use a density-based approach, and/or clean them spatially using manual inspection.
## No approach will get everything, and the scaling problems with occurrence vs spatial resolution remain. The key is to
## get sensible looking maps at the end :]


## Also, we should consider that some species (e.g. Agave attenuata) are simply so poorly sampled, that we can't model them
## well, or estimate their distribution. We should discuss how to present those, if we indeed do.


#########################################################################################################################
## 2). CREATE INDIVIDUAL SHAPEFILES WITH UNIQUE ID FOR MANUAL CLEANING
#########################################################################################################################


## The whole table of records is too big, so split it up into each species...
GBIF.ALA.CLEAN               = GBIF.ALA.COMBO.HIA
coordinates(GBIF.ALA.CLEAN)  <- ~lon+lat
proj4string(GBIF.ALA.CLEAN)  <- '+init=epsg:4326'
GBIF.ALA.CLEAN               <- spTransform(GBIF.ALA.CLEAN, CRS('+init=ESRI:54009'))


## Could also remove the species with insufficent records to run maxent?
#COMBO.RASTER.CLEAN   = COMBO.RASTER.CLEAN[!COMBO.RASTER.CLEAN$searchTaxon %in% unique(MISSING$searchTaxon), ]
#FINAL.MAXENT.LIST =  unique(COMBO.RASTER.CLEAN$searchTaxon[!COMBO.RASTER.CLEAN$searchTaxon %in% unique(MISSING$searchTaxon) ])


## The OBS field is truncated in ArcMap by default
View(head(subset(GBIF.ALA.CLEAN, searchTaxon == "Araucaria bidwillii"), 368))


## Reorder the columns to see in ArcMap
GBIF.ALA.CLEAN  = GBIF.ALA.CLEAN[,c("searchTaxon",
                                    "OBS",
                                    "scientificName",
                                    "commonname",  
                                    "taxonRank", 
                                    "genus", 
                                    "family", 
                                    "TPL_binomial",                
                                    "taxo_agree", 
                                    "gbifID",
                                    "id",
                                    "cloc",                          
                                    "basisOfRecord",                 
                                    "locality",                      
                                    "establishmentMeans",           
                                    "institutionCode",               
                                    "datasetName",                   
                                    "habitat",                       
                                    "catalogNumber",                 
                                    "country",                       
                                    "coordinateUncertaintyInMeters",
                                    "geodeticDatum",                 
                                    "year",                          
                                    "month",                         
                                    "day",                           
                                    "eventDate",                    
                                    "eventID",                       
                                    "CULTIVATED",                    
                                    "POLY_ID",                       
                                    "SUB_CODE_7",                    
                                    "REG_CODE_7")]


## Then, loop over the species list and create a shapefile for each 
for (i in 1:length(searchTaxon)) {
  
  ## Need to check the OBS column matches up - or do we not need this again?
  tmp <- GBIF.ALA.CLEAN[GBIF.ALA.CLEAN$searchTaxon == searchTaxon[i], ] 
  writeOGR(tmp, dsn = "./data/base/HIA_LIST/COMBO/CLEAN_GBIF", searchTaxon[i], driver = "ESRI Shapefile", overwrite_layer = TRUE)
  
}





#########################################################################################################################
## 3). MANUALLY INSPECT THE SHAPEFILES
#########################################################################################################################


#########################################################################################################################
## Need a criteria for removing records :: key ecological information is in the density of records 

## Spatial outliers     - long way from where they are most dense, e.g. in a garden of central Australian homestead
## Logistical outliers  - Herbarium specimens, etc (we know the location of Australian herbaria, but not overseas)







#########################################################################################################################
## Then Read back in a list of suspect points and remove them ................................................................





