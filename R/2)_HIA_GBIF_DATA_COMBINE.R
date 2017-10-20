#########################################################################################################################
#######################################  COMBINE ALL SPECIES DATA FRAMES INTO ONE ####################################### 
#########################################################################################################################


#########################################################################################################################
## 1). CREATE LIST OF TAXA FROM DOWLOADED FILES
#########################################################################################################################


## quickly check the GBIF names:
sp.n = "Magnolia grandiflora"
Magnolia.grandiflora = gbif(sp.n, download = TRUE)
GBIF.names = sort(names(Magnolia.grandiflora))


## Rachel uses month and identifier. We could add a final check inside the Maxent code that will look for these.
## But the taxonomic problems are more important
unique(Magnolia.grandiflora$month)
unique(Magnolia.grandiflora$identifiedBy)


## The species list doesn't match the downloaded species, so create a list from the downloaded files
spp.download = list.files("./data/base/HIA_LIST/GBIF/SPECIES/", pattern = ".RData")
spp.download = gsub("_GBIF_records.RData", "", spp.download)


# ## genera
# gen.download = list.files("./data/base/HIA_LIST/GBIF/GENERA/", pattern = ".RData")
# gen.download = gsub("_GBIF_records.RData", "", gen.download )
# 
# str(spp.download)
# str(gen.download)





#########################################################################################################################
## 2). LOAD TAXA, ADD SEARCH TAXA, DROP COLUMNS AND COMBINE SPECIES DATAFRAMES INTO ONE
#########################################################################################################################


# ## memory is a problem. So we need more RAM
# memory.limit()
# gc()
# 
# 
# #########################################################################################################################
# ## Take the first 300 taxa
# GBIF.1000 <- spp.download[c(1:1000)] %>%
# 
#   ## pipe the list into lapply
#   lapply(function(x) {
# 
#     ## create the character string
#     f <- sprintf("./data/base/HIA_LIST/GBIF/SPECIES/%s_GBIF_records.RData", x)
# 
#     ## load each .RData file
#     d <- get(load(f))
# 
#     ## now drop the columns which we don't need
#     dat <- data.frame(searchTaxon = x, d[, !colnames(d) %in% gbifColsToDrop],
#                       stringsAsFactors = FALSE)
#     if(!is.character(dat$gbifID)) {
#       dat$gbifID <- as.character(dat$gbifID)
#     }
#     dat
#   }) %>%
# 
#   ## finally, bind all the rows together
#   bind_rows
# 
# 
# #########################################################################################################################
# ## Take the second 300 taxa
# GBIF.2000 <- spp.download[c(1001:2000)] %>%
# 
#   ## pipe the list into lapply
#   lapply(function(x) {
# 
#     ## create the character string
#     f <- sprintf("./data/base/HIA_LIST/GBIF/SPECIES/%s_GBIF_records.RData", x)
# 
#     ## load each .RData file
#     d <- get(load(f))
# 
#     ## now drop the columns which we don't need
#     dat <- data.frame(searchTaxon = x, d[, !colnames(d) %in% gbifColsToDrop],
#                       stringsAsFactors = FALSE)
#     if(!is.character(dat$gbifID)) {
#       dat$gbifID <- as.character(dat$gbifID)
#     }
#     dat
#   }) %>%
# 
#   ## finally, bind all the rows together
#   bind_rows
# 
# 
# #########################################################################################################################
# ## take the last 300 taxa
# GBIF.3000 <- spp.download[c(2001:3000)] %>%
# 
#   ## pipe the list into lapply
#   lapply(function(x) {
# 
#     ## create the character string
#     f <- sprintf("./data/base/HIA_LIST/GBIF/SPECIES/%s_GBIF_records.RData", x)
# 
#     ## load each .RData file
#     d <- get(load(f))
# 
#     ## now drop the columns which we don't need
#     dat <- data.frame(searchTaxon = x, d[, !colnames(d) %in% gbifColsToDrop],
#                       stringsAsFactors = FALSE)
#     if(!is.character(dat$gbifID)) {
#       dat$gbifID <- as.character(dat$gbifID)
#     }
#     dat
#   }) %>%
# 
#   ## finally, bind all the rows together
#   bind_rows
# 
# 
# #########################################################################################################################
# ## take the last 300 taxa
# GBIF.4000 <- spp.download[c(3001:4000)] %>%
#   
#   ## pipe the list into lapply
#   lapply(function(x) {
#     
#     ## create the character string
#     f <- sprintf("./data/base/HIA_LIST/GBIF/SPECIES/%s_GBIF_records.RData", x)
#     
#     ## load each .RData file
#     d <- get(load(f))
#     
#     ## now drop the columns which we don't need
#     dat <- data.frame(searchTaxon = x, d[, !colnames(d) %in% gbifColsToDrop],
#                       stringsAsFactors = FALSE)
#     if(!is.character(dat$gbifID)) {
#       dat$gbifID <- as.character(dat$gbifID)
#     }
#     dat
#   }) %>%
#   
#   ## finally, bind all the rows together
#   bind_rows
# 
# 
# #########################################################################################################################
# ## take the last 300 taxa
# GBIF.5000 <- spp.download[c(4001:5000)] %>%
#   
#   ## pipe the list into lapply
#   lapply(function(x) {
#     
#     ## create the character string
#     f <- sprintf("./data/base/HIA_LIST/GBIF/SPECIES/%s_GBIF_records.RData", x)
#     
#     ## load each .RData file
#     d <- get(load(f))
#     
#     ## now drop the columns which we don't need
#     dat <- data.frame(searchTaxon = x, d[, !colnames(d) %in% gbifColsToDrop],
#                       stringsAsFactors = FALSE)
#     if(!is.character(dat$gbifID)) {
#       dat$gbifID <- as.character(dat$gbifID)
#     }
#     dat
#   }) %>%
#   
#   ## finally, bind all the rows together
#   bind_rows
# 
# 
# #########################################################################################################################
# ## take the last 300 taxa
# GBIF.6000 <- spp.download[c(5001:length(spp.download))] %>%
#   
#   ## pipe the list into lapply
#   lapply(function(x) {
#     
#     ## create the character string
#     f <- sprintf("./data/base/HIA_LIST/GBIF/SPECIES/%s_GBIF_records.RData", x)
#     
#     ## load each .RData file
#     d <- get(load(f))
#     
#     ## now drop the columns which we don't need
#     dat <- data.frame(searchTaxon = x, d[, !colnames(d) %in% gbifColsToDrop],
#                       stringsAsFactors = FALSE)
#     if(!is.character(dat$gbifID)) {
#       dat$gbifID <- as.character(dat$gbifID)
#     }
#     dat
#   }) %>%
#   
#   ## finally, bind all the rows together
#   bind_rows


#########################################################################################################################
## Take the last 300 taxa
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


## store the total no.of records as a variable
# GBIF.ALL = bind_rows(GBIF.1000, GBIF.2000, GBIF.3000, GBIF.4000, GBIF.5000, GBIF.6000)
total.records = dim(GBIF.ALL)[1]


## now get just the columns we want to keep. Note gc() frees up RAM
GBIF.TRIM <- GBIF.ALL %>% 
  select(one_of(gbif.keep))


## Remove the varieties, etc., from the scientific name returned by GBIF: probably get rid of this!
## Also, replace with John's example?
#Returned.binomial <- unlist(lapply(GBIF.TRIM$scientificName, string_fun_first_two_words))
# GBIF.TRIM = cbind(Returned.binomial, GBIF.TRIM)


#########################################################################################################################
## Here we could check if the searched and returned taxa match
## filter(scientificName == searchTaxon)
# GBIF.TRIM.CHECK = GBIF.TRIM
# GBIF.TRIM$CHECK = pmatch(GBIF.TRIM.CHECK$Returned.binomial, GBIF.TRIM.CHECK$scientificName)
# GBIF.TRIM.NO    = subset(GBIF.TRIM, CHECK == NA)


## Check how many records match the search term?
# head(GBIF.TRIM)
# head(GBIF.TRIM$scientificName, 10)
# head(GBIF.TRIM$searchTaxon, 10)
# head(GBIF.TRIM$scientificName, 10)
# rm(GBIF.ALL)
# rm(GBIF.300)
# rm(GBIF.600)
# rm(GBIF.800)
gc()

## Remove working dataframes from memory
save(GBIF.TRIM, file = paste("./data/base/HIA_LIST/GBIF/GBIF_TRIM.RData"))





#########################################################################################################################
## OUTSTANDING DOWNLOADING TASKS:
#########################################################################################################################


## Get the species with > 200k records

## Keep in touch with Renee RE the list





#########################################################################################################################
#####################################################  END ############################################################## 
#########################################################################################################################