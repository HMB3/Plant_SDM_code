#########################################################################################################################
#######################################  COMBINE ALL SPECIES DATA FRAMES INTO ONE ####################################### 
#########################################################################################################################


#########################################################################################################################
## 1). CREATE LIST OF TAXA
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
## 2). LOAD TAXA, ADD SEARCH TAXA, DROP COLUMNS AND COMBINE SPECIES DATAFRAMES INTO ONE
#########################################################################################################################


## take the list of taxa
GBIF.HIA.SPP.RECORDS.300 <- spp.download[c(1:300)] %>% ##length(spp.download)
  
  ## pipe the list into lapply
  lapply(function(x) {
    
    ## create the character string
    f <- sprintf("./data/base/HIA_LIST/GBIF/SPECIES/%s_GBIF_records.RData", x)
    
    ## load each .RData file
    d <- get(load(f))
    
    ## now drop the columns which we don't need
    data.frame(searchTaxon = x, d[, !colnames(d) %in% gbifColsToDrop], 
               stringsAsFactors = FALSE)
    
  }) %>% 

  ## finally, bind all the rows together
  bind_rows





## take the list of taxa
GBIF.HIA.SPP.RECORDS.600 <- spp.download[c(301:600)] %>% ##length(spp.download)
  
  ## pipe the list into lapply
  lapply(function(x) {
    
    ## create the character string
    f <- sprintf("./data/base/HIA_LIST/GBIF/SPECIES/%s_GBIF_records.RData", x)
    
    ## load each .RData file
    d <- get(load(f))
    
    ## now drop the columns which we don't need
    data.frame(searchTaxon = x, d[, !colnames(d) %in% gbifColsToDrop], 
               stringsAsFactors = FALSE)
    
  }) %>% 
  
  ## finally, bind all the rows together
  bind_rows





## take the list of taxa
GBIF.HIA.SPP.RECORDS.800 <- spp.download[c(601:length(spp.download))] %>% ##length(spp.download)
  
  ## pipe the list into lapply
  lapply(function(x) {
    
    ## create the character string
    f <- sprintf("./data/base/HIA_LIST/GBIF/SPECIES/%s_GBIF_records.RData", x)
    
    ## load each .RData file
    d <- get(load(f))
    
    ## now drop the columns which we don't need
    data.frame(searchTaxon = x, d[, !colnames(d) %in% gbifColsToDrop], 
               stringsAsFactors = FALSE)
    
  }) %>% 
  
  ## finally, bind all the rows together
  bind_rows





#########################################################################################################################
## bind all 3 together
dim(GBIF.HIA.SPP.RECORDS.300)[1]+dim(GBIF.HIA.SPP.RECORDS.600)[1]+dim(GBIF.HIA.SPP.RECORDS.800)[1]
GBIF.HIA.SPP.RECORDS.ALL = bind_rows(GBIF.HIA.SPP.RECORDS.300,
                                     GBIF.HIA.SPP.RECORDS.600,
                                     GBIF.HIA.SPP.RECORDS.800)


## 
dim(GBIF.HIA.SPP.RECORDS.ALL)
GBIF.HIA.SPP.RECORDS.TEST = 
TEST =   GBIF.HIA.SPP.RECORDS.ALL[, !colnames(GBIF.HIA.SPP.RECORDS.ALL) %in% gbifColsToDrop]



#########################################################################################################################
## 2). COMBINE ALL SPECIES DATA FRAMES INTO ONE WITH STU'S APPROACH
#########################################################################################################################

# Stu: for binding multiple data frames that may not all have the same columns: dplyr::bind_rows() will do this. I do this once 
# all the individual files are downloaded for each species. This can get very slow for a large number of data frames if you 
# do them one by one, I found a quicker way was to read all data frames to a list, then use a single bind_rows() to do them 
# all in one hit. The code below is also here (https://gist.github.com/snubian/b0819f1ad9861d7663161f49bc91663a)


# ## set download variables
# rm(list = ls())
# #gbif_config(download_reason_id = 7)
# 
# 
# ## get final taxon list (i.e. the cleaned list)
# str(draft.taxa) 
# 
# 
# ## for all draft.taxa in a list...
# for (i in 1:nrow(draft.taxa)) {
#   
#   ## print the draft.taxa to screen
#   ## not "i" is created inside the loop ()
#   taxon <- draft.taxa[i, ]$Species
#   message(taxon)
#   
#   ## get occurrence data for taxon from ALA, with extra fields specified
#   #occ <- occurrences(taxon, extra = c("data_hub_uid", "data_provider", "data_resource"))
#   occ <- gbif(taxon, download = TRUE)
#   
#   ## pipe the draft.taxa into a .CSV file... 
#   taxon %>%
#     str_replace(" ", "_") %>%
#     #paste0("output/occurrence/ala/original/", ., "_data.csv") %>%
#     paste0("./data/base/HIA_LIST/GBIF/CSV/", ., "_data.csv") %>%
#     write_csv(occ$data, .)
#   
# }





#########################################################################################################################
## 2). COMBINE
#########################################################################################################################


# #########################################################################################################################
# ## this seems to be a way to combine thousands of data frames without it getting progressively slower and s l o w e r
# ALL.GBIF.HIA.RECORDS <- list()
# 
# 
# ## for all draft.taxa
# for (i in seq_len(nrow(draft.taxa))) {
#   
#   data <-
#     
#     #browser()
#     draft.taxa[i, ]$Species %>%
#     #str_replace(" ", " ")   %>%
#     paste0("./data/base/HIA_LIST/GBIF/", ., "_GBIF_records.RData") %>%
#     
#     load %>%
#     
#     ## here we are creating a column which is the draft.taxa searched for in GBIF
#     ## match this against what was returned, to check for other errors...
#     mutate(searchTaxon   = draft.taxa[i, ]$Species,
#            catalogNumber = as.character(catalogNumber)) %>%
#     
#     ## this is the function which combines dataframes with different columns
#     dplyr::select(-one_of(gbifColsToDrop))
#   
#   ## the ith element of data becomes big GBIF dataframe
#   ALL.GBIF.HIA.RECORDS[[i]] <- data
#   
# }





#########################################################################################################################
## create the combined GBIF object
# GBIF <-
#   
#   ## check the searched species matches the found species...
#   bind_rows(ALL.GBIF.HIA.RECORDS) %>%
#   filter(scientificName == searchTaxon)
# 
# GBIF %>% write_csv("./data/base/HIA_LIST/GBIF/ALL_GBIF_HIA_SPP_RECORDS.csv")

## get sample of 100 taxa
# GBIF %>%
#   
#   filter(scientificName %in% taxa[sample(nrow(taxa), 100), ]$taxonName) %>%
#   write_csv("output/occurrence/occurrence_combined_sample.csv")





#########################################################################################################################
#####################################################  END ############################################################## 
#########################################################################################################################