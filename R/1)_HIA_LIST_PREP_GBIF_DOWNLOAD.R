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
## 1). DOWNLOAD RECORDS FROM GBIF USING HIA LIST
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
setdiff.25      = download_GBIF_setdiff_species(missing.25)
setdiff.Renee   = download_GBIF_Renee_species(spp.renee)


## check an eg file...not sure why it was working with RData files, but not .csv files...
load("./data/base/HIA_LIST/GBIF/Viburnum suspensum_GBIF_records.RData")
str(GBIF)

source("./R/2)_HIA_GBIF_DATA_COMBINE.R")
source("./R/3)_GBIF_DATA_FILTER.R")
source("./R/4)_GBIF_RECORDS_NUMERICAL_SUMMARY.R")
source("./R/4)_GBIF_ESTIMATE_NICHE.R")



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