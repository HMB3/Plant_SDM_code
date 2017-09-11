#########################################################################################################################
##############################  DOWNLOAD GLOBAL SPECIES RECORDS FOR HORT AUS LIST ####################################### 
#########################################################################################################################


## This code downloads all occurrence records for species on the Horticulture Australia list, from the GBIF database. 
## The trick here will be substituing the functions from the ALA4R package for the rgbif functions...
## Matching the core fields could be tricky for some taxa. Also, 


## setup 
library(gtools)
library(GISTools)
library(devtools)
library(Rtools)
#devtools::install_github("ropensci/rgbif")
library(data.table)
library(RCurl)
library(Rcpp)
library(raster)
library(rgdal)
#library(plyr)
library(dplyr)
library(sfsmisc)

library(SDMTools)
library(rmaxent)
library(dismo)
library(AdaptR)
library(ALA4R)
library(rgbif)
library(speciesgeocodeR)
library(raster)

library(knitr)
library(htmltools)
library(yaml)
library(caTools)
library(bitops)
library(rmarkdown)
library(gsubfn)
library(functional)
library(splitstackshape)

library(tidyverse)
library(stringr)
library(maptools)
library(rgeos)
library(magrittr)
library(datastorr)
library(baad.data)

library(Cairo)
library(lattice)
library(latticeExtra)

library("biglm")
library("bigmemory")
library("biganalytics")
library("bigtabulate")
library(taxonlookup)

## source functions
source('./R/GREEN_CITIES_FUNCTIONS.R')
source('./R/SDM_FUNCTIONS.R')




#########################################################################################################################
## 1). READ IN DRAFT HIA LIST AND CLEAN
#########################################################################################################################


# ## this list derives from all species and varieties sold anywhere in Australia in the last 5 years. Anthony Manea cleaned 
# ## up the data and cross-linked to growth form and exotic/native status and derived a list of ~1000 species that are the 
# ## most commonly sold, covering the right ratio of growth forms, regional representation and native/exotic
# spp.list = read.csv("./data/base/HIA_LIST/HIA/GREEN_CITIES_DRAFT_LIST.csv", stringsAsFactors = FALSE)
# top.200  = read.csv("./data/base/HIA_LIST/HIA/HIA_TOP_200.csv", stringsAsFactors = FALSE)
# 
# 
# ## have a look
# dim(spp.list)
# str(spp.list)
# head(spp.list)
# 
# 
# ## also, add the "Top 200 species in here"
# spp.200          = top.200[c("Species")]
# spp.200$Top_200  = "TRUE"
# 
# spp.list = merge(spp.list, spp.200, by = "Species", all.x = TRUE) 
# spp.list$Top_200[is.na(spp.list$Top_200)] <- "FALSE"
# spp.list$Origin <- gsub(" ",  "", spp.list$Origin)
# 
# 
# ## check
# str(spp.list)
# head(spp.list)
# unique(spp.list$Top_200)
# unique(spp.list$Origin)
# 
# 
# #########################################################################################################################
# ## WHAT ARE SOME USEFUL MEASURES OF THE DATASET?
# #########################################################################################################################
# 
# 
# ## ~ HALF the total species are native, ~47% of the top 200
# dim(subset(spp.list, Origin == "Native"))[1]/dim(spp.list)[1]*100
# dim(subset(top.200, Origin == "Native"))[1]/dim(top.200)[1]*100
# 
# 
# 
# 
# 
# #########################################################################################################################
# ## 2). DRAFT LIST PREP: THIS NEEDS TO CHANGE IN CONSULTATION WITH RACH, PAUL, ETC
# #########################################################################################################################
# 
# 
# ## for now, remove all the subspecies, varieties etc.
# ## but check if some of them work on GBIF, e.g. a separate list of varities
# spp.list$Species <- gsub("spp.", "", spp.list$Species)
# spp.list$Species <- gsub(" x ",  " ", spp.list$Species)
# spp.list$Species <- gsub("  ",  " ", spp.list$Species)
# ## spp.list$Species = gsub(" $","",     spp.list$Species, perl = TRUE) ## also trailing space?
# 
# 
# ## then remove the varieties and subspecies
# spp.list$Species = vapply(lapply(strsplit(spp.list$Species, " "),
#                                  unique), paste, character(1L), collapse = " ")
# 
# 
# ## then get just the first two words (again cleaning up the subspecies, and single genera)
# spp.list$Species = vapply(lapply(strsplit(spp.list$Species, " "), 
#                                  string_fun_first_two_words), paste, character(1L), collapse = " ")
# 
# 
# ## if exclude NA, spp.list = spp.list[!grepl("NA", spp.list$Species),]
# spp.list = spp.list[with(spp.list, order(Species)), ] 
# 
# 
# ## check
# str(spp.list)
# View(spp.list)
# 
# 
# ## Now create list of HIA taxa. 768 unique species, mius the corrections, etc. 
# spp.list = spp.list[with(spp.list, order(Species)), ]
# spp = unique(as.character(spp.list$Species))
# str(spp)   ## why 680? later check on what happens with the different queries
# head(spp, 50)
# tail(spp, 50)
# 
# 
# ## also make a genus list?
# genera = unique(vapply(lapply(strsplit(spp.list$Species, " "), 
#                               string_fun_first_word), paste, character(1L), collapse = " "))
# 
# 
# ## check
# str(genera)
# head(genera, 50)
# tail(genera, 50)
# 
# 
# ## also, for Stuarts code EG, I need a df not a list. Get just the rows of spp.list which have
# ## unique species names. 
# DRAFT.HIA.TAXA = spp.list[!duplicated(spp.list$Species), ]
# dim(DRAFT.HIA.TAXA)
# head(DRAFT.HIA.TAXA)
# 
# 
# ########################################################################################################################
# ## Try using taxonlookup to check the taxonomy
# DRAFT.TAXA.LOOKUP = lookup_table(DRAFT.HIA.TAXA[["Species"]], by_species = TRUE) ## convert rows to column and merge
# DRAFT.TAXA.LOOKUP = setDT(DRAFT.TAXA.LOOKUP , keep.rownames = TRUE)[]
# DRAFT.TAXA.LOOKUP = rename(DRAFT.TAXA.LOOKUP, searchTaxon = rn)
# head(DRAFT.TAXA.LOOKUP)



#########################################################################################################################
## 3). DOWNLOAD RECORDS FROM GBIF USING HIA LIST
#########################################################################################################################


## now run lOOPs to dowload species in the "spp" list from GBIF
## not including any data quality checks here, just downloading everything ##
## 


## set a few global variables to be used inside the functions...
#GBIF.download.limit = 200000
# skip.spp.list       = list()
# skip.gen.list       = list()


## run the download function on the species and genera lists
## these functions need to download at least one file, or they will return NULL
skipped.species = download_GBIF_all_species(spp)    ## saves each spp as .Rdata file, returning list of skipped spp 
skipped.genera  = download_GBIF_all_genera(genera)  ## saves each gen as .Rdata file, returning list of skipped genera 


## get the extra species
setdiff.species = download_GBIF_setdiff_species(missing.taxa)
setdiff.species = download_GBIF_Renee_species(spp.renee)


## check an eg file...not sure why it was working with RData files, but not .csv files...
load("./data/base/HIA_LIST/GBIF/Viburnum suspensum_GBIF_records.RData")
str(GBIF)





#########################################################################################################################
## 4). CHECK THE SPECIES WHICH WERE SKIPPED? 
#########################################################################################################################


## now convert the lists of skipped species and genera into a dataframe
length(skipped.species)   
length(skipped.genera)


## converting the lists of skipped species and genera into a dataframe
skipped.species.df <- data.frame(matrix(unlist(skipped.species), nrow = length(skipped.species), byrow = TRUE))
skipped.genera.df  <- data.frame(matrix(unlist(skipped.genera),  nrow = length(skipped.genera),  byrow = TRUE))
skipped.setdiff.df <- data.frame(matrix(unlist(setdiff.species),  nrow = length(setdiff.species),  byrow = TRUE))

## split the reason and the species into separate columns
skipped.species.df  <- cSplit(skipped.species.df, 1:ncol(skipped.species.df), sep = "|", stripWhite = TRUE, type.convert = FALSE)
skipped.genera.df   <- cSplit(skipped.genera.df,  1:ncol(skipped.genera.df),  sep = "|", stripWhite = TRUE, type.convert = FALSE)
skipped.setdiff.df  <- cSplit(skipped.setdiff.df,  1:ncol(skipped.setdiff.df),  sep = "|", stripWhite = TRUE, type.convert = FALSE)


## update names
colnames(skipped.species.df)[1] <- "Reason_skipped"
colnames(skipped.species.df)[2] <- "Species"
colnames(skipped.genera.df)[1]  <- "Reason_skipped"
colnames(skipped.genera.df)[2]  <- "Genus"  ## head(skipped.species.df), head(skipped.genera.df)
colnames(skipped.setdiff.df)[1]  <- "Reason_skipped"
colnames(skipped.setdiff.df)[2]  <- "Genus" 


## get subset for each type
max.records.spp  <- skipped.species.df[ which(skipped.species.df$Reason_skipped == "Number of records > 200,000"), ]
max.records.gen  <- skipped.genera.df[ which(skipped.genera.df$Reason_skipped   == "Number of records > 200,000"), ]

name.records.spp <- skipped.species.df[ which(skipped.species.df$Reason_skipped == "Possible incorrect nomenclature"), ]
name.records.gen <- skipped.genera.df[ which(skipped.genera.df$Reason_skipped   == "Possible incorrect nomenclature"), ]
no.records.spp   <- skipped.species.df[ which(skipped.species.df$Reason_skipped == "No GBIF records"), ]


## create lists for each category
max.records.spp.list  = unique(as.character(max.records.spp$Taxa))
name.records.spp.list = unique(as.character(name.records.spp$Taxa))
no.records.spp.list   = unique(as.character(no.records.spp$Taxa))

max.records.gen.list  = unique(as.character(max.records.gen$Taxa))
name.records.gen.list = unique(as.character(name.records.gen$Taxa))


## save lists just in case
## save(skipped.species.df, file = paste("./data/base/HIA_LIST/GBIF/skipped_species.RData", sep = ""))
## save(skipped.genera.df,  file = paste("./data/base/HIA_LIST/GBIF/skipped_genera.RData",  sep = ""))


## have a look at the list of skipped species
View(skipped.species.df)


## which of these species are on the top 200?
skipped.200.spp = merge(spp.200, skipped.species.df, by = "Species", all = FALSE)
kable(skipped.200.spp)





#########################################################################################################################
## 5). DOWNLOAD SPECIES WITH >200k RECORDS USING THE GBIF API
#########################################################################################################################


## create a list of taxa with too many records, to pass to the curl code..

##
# function gbifapi { 
#   
#   curl -i –user popple_1500:Popple1500 -H "Content-Type: 
#   application/json" -H "Accept: application/json" -X POST -d "{\"creator\":
#   \”yourGbifUserName\", \"notification_address\": [\”hugh.burley@mq.edu.au\"], 
#   \"predicate\": {\"type\":\"and\", \"predicates\": [{\"type\":\"equals\",\"key\":
#   \"HAS_COORDINATE\",\"value\":\"true\"}, {\"type\":\"equals\", \"key\":
#   \"TAX O N _KEY\", \"value\":\"$1\"}] }}" 
#   http://api.gbif.org/v1/occurrence/download/request  >> log.txt 
#   echo -e "\r\n$1 $2\r\n\r\n----------------\r\n\r\n" >> log.txt
#   
# }  
# 
# $ gbifapi 4140730 "Betula pendula" 
# $ gbifapi 2882316 "Fagus sylvatica" 
# $ gbifapi 3172358 "Fraxinus excelsior" 
# $ gbifapi 8351737 "Hedera helix" 
# $ gbifapi 2878688 "Quercus robur"





#########################################################################################################################
############################################  END OF GBIF DOWLOAD CODE ################################################## 
#########################################################################################################################