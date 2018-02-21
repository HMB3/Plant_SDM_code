#########################################################################################################################
################################################ CLEAN RECORDS ########################################################## 
#########################################################################################################################


#########################################################################################################################
## REMOVE DODGY RECORDS FROM NICHE CALCULATION
#########################################################################################################################
dim(GBIF.ALA.COMBO.HIA)
GBIF.ALA.COMBO.HIA$OBS <- 1:nrow(GBIF.ALA.COMBO.HIA)
dim(GBIF.ALA.COMBO.HIA)[1];length(GBIF.ALA.COMBO.HIA$OBS)                                          ## This matches step 6


#########################################################################################################################
## Now check the field names :: ALA and GBIF fields are uniquely identified?
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


## So, not simple to just exclude the GBIF records for ALA species...


#########################################################################################################################
## CREATE INDIVIDUAL SHAPEFILES FOR MANUAL CLEANING
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



## Read in a list of bad points and remove...................................................................
## Check the OBS feature works on some points