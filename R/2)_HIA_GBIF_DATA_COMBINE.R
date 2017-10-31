#########################################################################################################################
#######################################  COMBINE ALL SPECIES DATA FRAMES INTO ONE ####################################### 
#########################################################################################################################


#########################################################################################################################
## 1). CREATE LIST OF TAXA FROM DOWLOADED FILES
#########################################################################################################################


#########################################################################################################################
## Create a list of species from the downloaded files
spp.download = list.files("./data/base/HIA_LIST/GBIF/SPECIES/", pattern = ".RData")
spp.download = gsub("_GBIF_records.RData", "", spp.download)





#########################################################################################################################
## 2). LOAD TAXA, ADD SEARCH TAXA, DROP COLUMNS AND COMBINE SPECIES DATAFRAMES INTO ONE
#########################################################################################################################


## memory is a problem. So we need more RAM. Can these sorts of operations be run in parallel?
memory.limit()
gc()


## Check GBIF column names:
sort(unique(gbifColsToDrop))
sort(unique(gbif.keep))
intersect(unique(gbifColsToDrop), unique(gbif.keep))     ## don't get rid of any of the columns we want to keep
 
 
#########################################################################################################################
## Combine all the taxa at once
GBIF.ALL <- spp.download[c(1:length(spp.download))] %>%

  ## pipe the list into lapply
  lapply(function(x) {

    ## create the character string
    f <- sprintf("./data/base/HIA_LIST/GBIF/SPECIES/%s_GBIF_records.RData", x)

    ## load each .RData file
    d <- get(load(f))

    ## now drop the columns which we don't need
    dat <- data.frame(searchTaxon = x, d[, !colnames(d) %in% gbifColsToDrop],
                      stringsAsFactors = FALSE)
    
    if(!is.character(dat$gbifID)) {
      
      dat$gbifID <- as.character(dat$gbifID)
      
    }
    
    dat
    
  }) %>%

  ## finally, bind all the rows together
  bind_rows





#########################################################################################################################
## 3). CHECK SEARCHED AND RETURNED TAXONOMIC NAMES FROM GBIF
#########################################################################################################################


## Store the total no.of records as a variable
# GBIF.ALL = bind_rows(GBIF.1000, GBIF.2000, GBIF.3000, GBIF.4000, GBIF.5000, GBIF.6000)
total.records = dim(GBIF.ALL)[1]


## What are the names? Can't get date for GBIF across 
sort(names(GBIF.ALL))
str(GBIF.ALL$recordedBy)
str(GBIF.ALL$dateIdentified)


## Now get just the columns we want to keep. Note gc() frees up RAM
GBIF.TRIM <- GBIF.ALL %>% 
  select(one_of(gbif.keep))


## Check names
dim(GBIF.TRIM)
names(GBIF.TRIM)


## Have a quick look at date values
# sort(unique(GBIF.TRIM$year))
# sort(unique(GBIF.TRIM$month))
# sort(unique(GBIF.TRIM$eventDate))


## Remove the varieties, etc., from the scientific name returned by GBIF: probably get rid of this!
## Also, replace with John's example?
#Returned.binomial <- unlist(lapply(GBIF.TRIM$scientificName, string_fun_first_two_words))
# GBIF.TRIM = cbind(Returned.binomial, GBIF.TRIM)


#########################################################################################################################
## How can we check if the searched and returned taxa match?
## filter(scientificName == searchTaxon)
# GBIF.TRIM.CHECK = GBIF.TRIM
# GBIF.TRIM$CHECK = pmatch(GBIF.TRIM.CHECK$Returned.binomial, GBIF.TRIM.CHECK$scientificName)
# GBIF.TRIM.NO    = subset(GBIF.TRIM, CHECK == NA)


## Check how many records match the search term?
gc()

## Remove working dataframes from memory
save(GBIF.TRIM, file = paste("./data/base/HIA_LIST/COMBO/GBIF_TRIM.RData"))


## Now save .RData file for the next session
save.image("STEP_2_GBIF_RAW.RData")





#########################################################################################################################
## OUTSTANDING DOWNLOADING TASKS:
#########################################################################################################################


## Get the species with > 200k records
## Keep checking the taxonomy





#########################################################################################################################
#####################################################  END ############################################################## 
#########################################################################################################################