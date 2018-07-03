#########################################################################################################################
#######################################  COMBINE ALL SPECIES DATA FRAMES INTO ONE ####################################### 
#########################################################################################################################


#########################################################################################################################
## This code combines all downloaded records for each species into a single database.


#########################################################################################################################
## 1). CREATE LIST OF TAXA FROM DOWNLOADED FILES
#########################################################################################################################


#########################################################################################################################
## Read in all data to run the SDM code :: species lists, shapefile, rasters & tables
#source('./R/HIA_LIST_MATCHING.R')
spp.download = list.files("./data/base/HIA_LIST/GBIF/SPECIES/", pattern = ".RData")
spp.download = gsub("_GBIF_records.RData", "", spp.download)
spp.download = trimws(spp.download)
bind.spp     = trimws(GBIF.spp)
spp.download = intersect(bind.spp, spp.download)


## Check the difference between the risk list and the downloaded list
#outstanding.spp = setdiff(spp.download, unique(COMBO.RASTER.CONTEXT$searchTaxon))
#intersect(NURSE.MATCH$species, spp.download)


#########################################################################################################################
## 2). LOAD TAXA, ADD SEARCH TAXA, DROP COLUMNS AND COMBINE SPECIES DATAFRAMES INTO ONE
#########################################################################################################################


## memory is a problem. So we need more RAM. Can these sorts of operations be run in parallel?
memory.limit()
gc()


## Check GBIF column names:
sort(unique(gbifColsToDrop))
sort(unique(gbif.keep))
intersect(unique(gbifColsToDrop), unique(gbif.keep))     ## Check we don't get rid of any of the columns we want to keep
 
 
#########################################################################################################################
## Combine all the taxa into a single dataframe at once
GBIF.ALL <- spp.download[c(1:length(spp.download))] %>%   ## spp.download[c(1:length(spp.download))] 

  ## Pipe the list into lapply
  lapply(function(x) {

    ## Create a character string of each .RData file
    f <- sprintf("./data/base/HIA_LIST/GBIF/SPECIES/%s_GBIF_records.RData", x)

    ## Load each file
    d <- get(load(f))

    ## Now drop the columns which we don't need
    dat <- data.frame(searchTaxon = x, d[, !colnames(d) %in% gbifColsToDrop],
                      stringsAsFactors = FALSE)
    
    if(!is.character(dat$gbifID)) {
      
      dat$gbifID <- as.character(dat$gbifID)
      
    }
    
    ## Need to print the object within the loop
    dat
    
  }) %>%

  ## Finally, bind all the rows together
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
setdiff(names(GBIF.TRIM), gbif.keep)
intersect(names(GBIF.TRIM), gbif.keep)


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
dim(GBIF.TRIM)
setdiff(RISK.BINOMIAL, GBIF.TRIM$searchTaxon) 


## Remove working dataframes from memory
#saveRDS(GBIF.TRIM, file = paste("./data/base/HIA_LIST/COMBO/GBIF_TRIM_LATEST.rds")) ##


## Now save .RData file for the next session
#save.image("STEP_2_GBIF_RAW_NEW_SPP.RData")





#########################################################################################################################
## OUTSTANDING DOWNLOADING TASKS:
#########################################################################################################################


## Get the species with > 200k records
## Keep checking the taxonomy





#########################################################################################################################
#####################################################  END ############################################################## 
#########################################################################################################################