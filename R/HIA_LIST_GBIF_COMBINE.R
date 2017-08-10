#########################################################################################################################
##############################  DOWNLOAD GLOBAL SPECIES RECORDS FOR HORT AUS LIST ####################################### 
#########################################################################################################################


## 
library(tidyverse)
library(stringr)
library(ALA4R)
library(magrittr)


## set download variables
rm(list = ls())
#gbif_config(download_reason_id = 7)




#########################################################################################################################
## 1). READ IN DRAFT HIA LIST
#########################################################################################################################


## get final taxon list (i.e. the cleaned list)
str(draft.taxa) 


## for all draft.taxa in a list...
for (i in 1:nrow(draft.taxa)) {
  
  ## print the draft.taxa to screen
  ## not "i" is created inside the loop ()
  taxon <- draft.taxa[i, ]$Species
  message(taxon)
  
  ## get occurrence data for taxon from ALA, with extra fields specified
  occ <- occurrences(taxon, extra = c("data_hub_uid", "data_provider", "data_resource"))
  
  ## pipe the draft.taxa into a .CSV file... 
  taxon %>%
    str_replace(" ", "_") %>%
    paste0("output/occurrence/ala/original/", ., "_data.csv") %>%
    write_csv(occ$data, .)
  
}





#########################################################################################################################
## 2). COMBINE
#########################################################################################################################


## what are the fields?
load("./data/base/HIA_LIST/GBIF/Viburnum suspensum_GBIF_records.RData")
names(GBIF)


#########################################################################################################################
## a bunch of GBIF fields we don't need
#########################################################################################################################
gbifColsToDrop <- c("cloc",
                    "crawlId",
                    "disposition",
                    "dynamicProperties",
                    "elevationAccuracy",
                    "endDayOfYear",
                    "startDayOfYear",
                    "eventRemarks",
                    "eventTime",
                    "eventID",
                    "month",
                    
                    "familyKey",
                    "genusKey",
                    "kingdomKey",
                    "phylumKey",
                    "publishingOrgKey",
                    
                    "georeferenceVerificationStatus",
                    "eventID",
                    "lastCrawled",
                    "higherGeography",
                    "municipality",
                    
                    "informationWithheld",
                    "language",
                    "protocol",
                    "rightsHolder")
                    
                    # "kingdom", 
                    # "phylum", 
                    # "class", 
                    # "order", 
                    # "IBRA7Regions", 
                    # "IMCRA4Regions", 
                    # "localGovernmentAreas",
                    # "country", 
                    # "basisOfRecordOriginal", 
                    # "sex", 
                    # "altitudeInFeet", 
                    # "altitudeNonNumeric",
                    # "minimumElevationInMetres", 
                    # "maximumElevationInMetres", 
                    # "minimumDepthInMeters", 
                    # "maximumDepthInMeters",
                    # "latitudeOriginal", 
                    # "longitudeOriginal", 
                    # "geodeticDatum")


#########################################################################################################################
## read data frames to list
#########################################################################################################################


## this seems to be a way to combine thousands of data frames without it getting progressively slower and s l o w e r
browser()
ALL.GBIF.HIA.RECORDS <- list()


## for all draft.taxa
for (i in seq_len(nrow(draft.taxa))) {
  
  data <-
    
    draft.taxa[i, ]$Species %>%
    #str_replace(" ", " ")   %>%
    paste0("./data/base/HIA_LIST/GBIF/", ., "_GBIF_records.RData") %>%
    
    load %>%
    
    ## here we are creating a column which is the draft.taxa searched for in GBIF
    ## match this against what was returned, to check for other errors...
    mutate(searchTaxon   = draft.taxa[i, ]$Species,
           catalogNumber = as.character(catalogNumber)) %>%
    
    ## this is the function which combines dataframes with different columns
    dplyr::select(-one_of(gbifColsToDrop))
  
  ## the ith element of data becomes big GBIF dataframe
  ALL.GBIF.HIA.RECORDS[[i]] <- data
  
}


## create the combined GBIF object
GBIF <-
  
  ## check the searched species matches the found species...
  bind_rows(ALL.GBIF.HIA.RECORDS) %>%
  filter(scientificName == searchTaxon)

GBIF %>% write_csv("./data/base/HIA_LIST/GBIF/ALL_GBIF_HIA_SPP_RECORDS.csv")

## get sample of 100 taxa
# GBIF %>%
#   
#   filter(scientificName %in% taxa[sample(nrow(taxa), 100), ]$taxonName) %>%
#   write_csv("output/occurrence/occurrence_combined_sample.csv")





#########################################################################################################################
#############################################  FUNCTIONS FOR HORT AUS LIST ############################################## 
#########################################################################################################################