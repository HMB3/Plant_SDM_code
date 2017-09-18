#########################################################################################################################
#######################################  COMBINE ALL SPECIES DATA FRAMES INTO ONE ####################################### 
#########################################################################################################################


#########################################################################################################################
## 1). CREATE LIST OF TAXA
#########################################################################################################################


## check ALA field names for a key species...
Magnolia.grandiflora = gbif('Magnolia grandiflora', download = TRUE)
GBIF.names = sort(names(Magnolia.grandiflora))

Syzygium.floribundum = occurrences(taxon = 'Syzygium floribundum', download_reason_id = 7)[["data"]]
ALA.names = sort(names(Syzygium.floribundum))
ALA.names
length(ALA.names)
length(GBIF.names)
setdiff(GBIF.names, ALA.names)

## the species list doesn't match the downloaded species, so create a list from the downloaded files
spp.download = list.files("./data/base/HIA_LIST/ALA/SPECIES/", pattern = ".RData")
spp.download = gsub("_ALA_records.RData", "", spp.download)





#########################################################################################################################
## 2). LOAD TAXA, ADD SEARCH TAXA, DROP COLUMNS AND COMBINE SPECIES DATAFRAMES INTO ONE
#########################################################################################################################


## memory is a problem. So we need more RAM
memory.limit()
gc()


#########################################################################################################################
## Take the first 300 taxa
ALA.300 <- spp.download[c(1:300)] %>% 
  
  ## Pipe the list into lapply
  lapply(function(x) {
    
    ## Create the character string
    f <- sprintf("./data/base/HIA_LIST/ALA/SPECIES/%s_ALA_records.RData", x)
    
    ## Load each .RData file
    d <- get(load(f))
    
    ## Now return the searched bionomial, but don't need to drop the columns which we don't need
    data.frame(searchTaxon = x, d[, !colnames(d) %in% ColsToDrop], 
               stringsAsFactors = FALSE)
    
  }) %>% 

  ## finally, bind all the rows together
  bind_rows


#########################################################################################################################
## Take the second 300 taxa
ALA.600 <- spp.download[c(301:length(spp.download))] %>% 
  
  ## pipe the list into lapply
  lapply(function(x) {
    
    ## create the character string
    f <- sprintf("./data/base/HIA_LIST/ALA/SPECIES/%s_ALA_records.RData", x)
    
    ## load each .RData file
    d <- get(load(f))
    
    ## now drop the columns which we don't need
    data.frame(searchTaxon = x, #d[, !colnames(d) %in% gbifColsToDrop], 
               stringsAsFactors = FALSE)
    
  }) %>% 
  
  ## finally, bind all the rows together
  bind_rows





#########################################################################################################################
## 3). CHECK SEARCHED AND RETURNED TAXONOMIC NAMES FROM ALA
#########################################################################################################################


## store the total no.of records as a variable
ALA.ALL = bind_rows(ALA.300,
                    ALA.600)
total.records = dim(ALA.ALL)[1]


## now get just the columns we want to keep. Note gc() frees up RAM
ALA.TRIM <- ALA.ALL #%>% 
  #select(one_of(gbif.keep))


## remove the varieties, etc., from the scientific name returned by ALA: probably get rid of this!
Returned.binomial <- unlist(lapply(ALA.TRIM$scientificName, string_fun_first_two_words))
ALA.TRIM = cbind(Returned.binomial, ALA.TRIM)


## Check how many records match the search term?
head(ALA.TRIM)
head(ALA.TRIM$scientificName, 10)
head(ALA.TRIM$searchTaxon, 10)
head(ALA.TRIM$scientificName, 10)


## Remove working dataframes from memory
rm(ALA.ALL)
rm(ALA.300)
rm(ALA.600)
rm(ALA.800)
gc()


# ## remove columns again
# ALA.ALL = ALA.ALL[, !colnames(ALA.ALL) %in% gbifColsToDrop]
# dim(ALA.ALL)





#########################################################################################################################
#####################################################  END ############################################################## 
#########################################################################################################################