#########################################################################################################################
#######################################  COMBINE ALL SPECIES DATA FRAMES INTO ONE ####################################### 
#########################################################################################################################


#########################################################################################################################
## 1). COMBINE ALL SPECIES DATA FRAMES USING SHAWN'S APPROACH
#########################################################################################################################


## check GBIF field names for a key species...
Magnolia.grandiflora = gbif('Magnolia grandiflora', download = TRUE)
GBIF.names = sort(names(Magnolia.grandiflora))
GBIF.names


## create empty dataframes to store all the data?
GBIF.HIA.SPP.RECORDS = data.frame()
GBIF.HIA.GEN.RECORDS = data.frame()


## the species list doesn't match the downloaded species, so create a list from the downloaded files
spp.download = list.files("./data/base/HIA_LIST/GBIF/SPECIES/", pattern = ".RData")
spp.download = gsub("_GBIF_records.RData", "", spp.download)

gen.download = list.files("./data/base/HIA_LIST/GBIF/GENERA/", pattern = ".RData")
gen.download = gsub("_GBIF_records.RData", "", gen.download )

str(spp.download)
str(gen.download)





#########################################################################################################################
## make a list of all the taxa, then combine the species
#########################################################################################################################

## take the list of taxa
# Start the clock!
ptm <- proc.time()
GBIF.HIA.SPP.RECORDS <- sprintf("./data/base/HIA_LIST/GBIF/SPECIES/%s_GBIF_records.RData", spp.download[c(1:length(spp.download))]) %>% ## [c(1:4)])
  
  ## pipe it into the lapply function to load all the .Rdata files
  lapply(function(x) get(load(x))) %>% 
  
  ## create a column for the searched taxon
  # mutate(searchTaxon   = spp.download[i], #taxa[i, ]$taxonName,
  #        catalogNumber = as.character(catalogNumber)) %>%

  ## then bind the rows together
  bind_rows

# Start the clock!
proc.time() - ptm


## Error: cannot allocate vector of size 32.0 Mb


## have a look at the output of bind rows
class(GBIF.HIA.SPP.RECORDS)
names(GBIF.HIA.SPP.RECORDS)
dim(GBIF.HIA.SPP.RECORDS)
unique(GBIF.HIA.SPP.RECORDS$scientificName)
unique(GBIF.HIA.SPP.RECORDS$searchTaxon)


ala <-
  bind_rows(df) %>%
  filter(scientificName == searchTaxon)

ala %>% write_csv("output/occurrence/occurrence_combined.csv")



## this is slow, because every file was saved as GBIF, so creating the list of data frames doesn't work...
for(sp.n in spp.download[c(1:2)]) {
  
  ## Load GBIF records for each species as a .CSV file
  filename = paste0("./data/base/HIA_LIST/GBIF/SPECIES/", sp.n, "_GBIF_records.RData")
  load(filename)

  ## print
  print(sp.n)
  
  ## then bind the .RData files together
  GBIF                  = GBIF
  GBIF.HIA.SPP.RECORDS  = bind_rows(GBIF.HIA.SPP.RECORDS, GBIF)
  
}


## now save the combined species GBIF data as one .RData file
class(GBIF.HIA.SPP.RECORDS)
dim(GBIF.HIA.SPP.RECORDS)
str(GBIF.HIA.SPP.RECORDS)

save(GBIF.HIA.SPP.RECORDS,  file = paste("./data/base/HIA_LIST/GBIF/SPECIES/", sep = ""))





#########################################################################################################################
## combine genera
## for all the species in the HIA list
for(gen.n in gen.download) {
  
  # Load GBIF records for each species as a .CSV file
  filename = paste0("./data/base/HIA_LIST/GBIF/GENERA/", sp.n, "_GBIF_records.RData")
  load(filename)
  
  ## then bind the .RData files together
  GBIF                  = GBIF
  GBIF.HIA.GEN.RECORDS  = bind_rows(GBIF.HIA.GEN.RECORDS, GBIF)
  
}


## now save the combined species GBIF data as one .RData file
save(GBIF.HIA.SPP.RECORDS,  file = paste("./data/base/HIA_LIST/GBIF/GENERA/", sep = ""))


## now check the combined data frames
dim(GBIF.HIA.SPP.RECORDS)
dim(GBIF.HIA.GEN.RECORDS)

head(GBIF.HIA.SPP.RECORDS, 20)
head(GBIF.HIA.GEN.RECORDS, 20) 





#########################################################################################################################
## 2). COMBINE ALL SPECIES DATA FRAMES INTO ONE WITH STU'S APPROACH
#########################################################################################################################

# Stu: for binding multiple data frames that may not all have the same columns: dplyr::bind_rows() will do this. I do this once 
# all the individual files are downloaded for each species. This can get very slow for a large number of data frames if you 
# do them one by one, I found a quicker way was to read all data frames to a list, then use a single bind_rows() to do them 
# all in one hit. The code below is also here (https://gist.github.com/snubian/b0819f1ad9861d7663161f49bc91663a)


## set download variables
rm(list = ls())
#gbif_config(download_reason_id = 7)


## get final taxon list (i.e. the cleaned list)
str(draft.taxa) 


## for all draft.taxa in a list...
for (i in 1:nrow(draft.taxa)) {
  
  ## print the draft.taxa to screen
  ## not "i" is created inside the loop ()
  taxon <- draft.taxa[i, ]$Species
  message(taxon)
  
  ## get occurrence data for taxon from ALA, with extra fields specified
  #occ <- occurrences(taxon, extra = c("data_hub_uid", "data_provider", "data_resource"))
  occ <- gbif(taxon, download = TRUE)
  
  ## pipe the draft.taxa into a .CSV file... 
  taxon %>%
    str_replace(" ", "_") %>%
    #paste0("output/occurrence/ala/original/", ., "_data.csv") %>%
    paste0("./data/base/HIA_LIST/GBIF/CSV/", ., "_data.csv") %>%
  write_csv(occ$data, .)
  
}





#########################################################################################################################
## 2). COMBINE
#########################################################################################################################


## what are the fields?
load("./data/base/HIA_LIST/GBIF/Viburnum suspensum_GBIF_records.RData")
names(GBIF)


## check GBIF field names for a key species...
Magnolia.grandiflora = gbif('Magnolia grandiflora', download = TRUE)
GBIF.names = sort(names(Magnolia.grandiflora))
GBIF.names


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
## Temporary fix
#########################################################################################################################


##
setwd("./data/base/HIA_LIST/GBIF/SPECIES/")

file_list <- list.files()

ala <- bind_rows(file_list)

for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    
    dataset <- read.table(file, header = TRUE, sep = "\t")
  }
  
  # if the merged dataset does exist, append to it
  if (exists("dataset")){
    
    temp_dataset <- read.table(file, header = TRUE, sep = "\t")
    dataset      <- bind_rows(dataset, temp_dataset)
    
    rm(temp_dataset)
  }
  
}




#########################################################################################################################
## this seems to be a way to combine thousands of data frames without it getting progressively slower and s l o w e r
ALL.GBIF.HIA.RECORDS <- list()


## for all draft.taxa
for (i in seq_len(nrow(draft.taxa))) {
  
  data <-
    
    #browser()
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


list.taxa <- sprintf("./data/base/HIA_LIST/GBIF/%s_GBIF_records.RData", draft.taxa$Species) %>%
  
  lapply(function(x) get(load(x))) %>% 
  
  bind_rows





#########################################################################################################################
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


df <- list()

for (i in seq_len(nrow(taxa))) {
  
  data <-
    
    taxa[i, ]$taxonName %>%
    
    str_replace(" ", "_") %>%
    
    paste0("output/occurrence/ala/original/", ., "_data.csv") %>%
    
    read_csv %>%
    
    mutate(searchTaxon = taxa[i, ]$taxonName,
           catalogNumber = as.character(catalogNumber)) %>%
    
    dplyr::select(-one_of(alaColsToDrop))
  
  df[[i]] <- data
  
}

ala <- bind_rows(df)





#########################################################################################################################
#####################################################  END ############################################################## 
#########################################################################################################################