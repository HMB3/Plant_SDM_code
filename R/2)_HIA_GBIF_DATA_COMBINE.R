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


## memory is a problem. So we need more RAM
memory.limit()


## Take the first 300 taxa
GBIF.300 <- spp.download[c(1:300)] %>% 
  
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


## Take the second 300 taxa
GBIF.600 <- spp.download[c(301:600)] %>% 
  
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



## take the last 300 taxa
GBIF.800 <- spp.download[c(601:length(spp.download))] %>% 
  
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
## 3). CHECK SEARCHED AND RETURNED TAXONOMIC NAMES FROM GBIF
#########################################################################################################################


## store the total no.of records as a varible
total.records = dim(GBIF.300)[1]+dim(GBIF.600)[1]+dim(GBIF.800)[1]
GBIF.ALL = bind_rows(GBIF.300,
                     GBIF.600,
                     GBIF.800)


## now get just the columns we want to keep. Note gc() frees up RAM
GBIF.TRIM <- GBIF.ALL %>% 
  select(one_of(gbif.keep))


## remove working dataframes from memory
rm(GBIF.ALL)
rm(GBIF.300)
rm(GBIF.600)
rm(GBIF.800)


## now check if this term has worked
head(GBIF.TRIM)


## remove the varieties, etc., from the scientific name returned by GBIF.
Returned.binomial <- unlist(lapply(GBIF.TRIM$scientificName, string_fun_first_two_words))
GBIF.TRIM = cbind(Returned.binomial, GBIF.TRIM)


## check how many records match the search term?
head(GBIF.TRIM)
head(GBIF.TRIM$scientificName, 10)
head(GBIF.TRIM$searchTaxon, 10)
head(GBIF.TRIM$scientificName, 10)
gc()


# ## remove columns again
# GBIF.ALL = GBIF.ALL[, !colnames(GBIF.ALL) %in% gbifColsToDrop]
# dim(GBIF.ALL)





#########################################################################################################################
#####################################################  END ############################################################## 
#########################################################################################################################