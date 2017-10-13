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
#library(Rtools)
#devtools::install_github("ropensci/rgbif")
library(data.table)
library(Rcpp)
library(raster)
library(rgdal)
#library(plyr)
library(dplyr)
library(sfsmisc)
library(spatstat)

library(SDMTools)
library(rmaxent)
library(dismo)
library(AdaptR)
library(red)
library(ConR)

library(ff)
library(rgeos)
library(sp)
library(raster)
library(rJava)
library(things)


library(ALA4R)
library(rgbif)
library(scrubr)
library(RCurl)
library(httr)
library(taxonlookup)
library(speciesgeocodeR)
library(raster)
library(raster)
library(rnaturalearth)
library(gdalUtils)

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


## source functions
source('./R/GREEN_CITIES_FUNCTIONS.R')
#source('./R/SDM_FUNCTIONS.R')



###########################################################################################
## remote directoty for raster data
## login from workstation: ## mqauth.uni.mq.edu.au\MQ20174568 ## Popple2016
## \\SCI-7910
## \\sci-7910\F\data\worldclim\world\0.5\bio




#########################################################################################################################
## 1). DOWNLOAD RECORDS FROM GBIF USING HIA LIST
#########################################################################################################################


## Now run loops to dowload species in the "spp" list from GBIF. Not including any data quality checks here, just 
## downloading everything
source('./R/HIA_LIST_MATCHING.R')


## Run the download function on the species and genera lists these functions need to download at least one file, or they will return NULL
skipped.species = download_GBIF_all_species(spp)    ## saves each spp as .Rdata file, returning list of skipped spp 
skipped.ALA     = download_ALA_all_species(spp)     ## saves each spp as .Rdata file, returning list of skipped spp 
skipped.genera  = download_GBIF_all_genera(genera)  ## saves each gen as .Rdata file, returning list of skipped genera 


#########################################################################################################################
## Run all the GBIF code 
source("./R/2)_HIA_GBIF_DATA_COMBINE.R")
source("./R/3)_GBIF_DATA_FILTER.R")
source("./R/4)_ALA_COMBO_GBIF_RECORDS_NUMERICAL_SUMMARY.R")





#########################################################################################################################
## 4). CHECK THE SPECIES WHICH WERE SKIPPED? 
#########################################################################################################################


## Now convert the lists of skipped species and genera into a dataframe
length(skipped.species)   
length(skipped.genera)


## Converting the lists of skipped species and genera into a dataframe
skipped.species.df <- data.frame(matrix(unlist(skipped.species), nrow = length(skipped.species), byrow = TRUE))
skipped.grow.df    <- data.frame(matrix(unlist(skipped.grow),    nrow = length(skipped.grow),    byrow = TRUE))
skipped.clean.df   <- data.frame(matrix(unlist(skipped.clean),   nrow = length(skipped.clean),   byrow = TRUE))


## Split the reason and the species into separate columns
skipped.species.df <- cSplit(skipped.species.df, 1:ncol(skipped.species.df), sep = "|", stripWhite = TRUE, type.convert = FALSE)
skipped.grow.df    <- cSplit(skipped.grow.df,    1:ncol(skipped.grow.df),    sep = "|", stripWhite = TRUE, type.convert = FALSE)
skipped.clean.df   <- cSplit(skipped.clean.df,   1:ncol(skipped.clean.df),   sep = "|", stripWhite = TRUE, type.convert = FALSE)


## Update names
colnames(skipped.species.df)[1] <- "Reason_skipped"
colnames(skipped.species.df)[2] <- "Species"
colnames(skipped.grow.df)[1]    <- "Reason_skipped"
colnames(skipped.grow.df)[2]    <- "Genus"  ## head(skipped.species.df), head(skipped.grow.df)
colnames(skipped.clean.df)[1]   <- "Reason_skipped"
colnames(skipped.clean.df)[2]   <- "Genus" 


## Get subset for each type
max.records.spp    <- skipped.species.df[ which(skipped.species.df$Reason_skipped == "Number of records > 200,000"), ]
max.records.grow   <- skipped.grow.df[ which(skipped.grow.df$Reason_skipped       == "Number of records > 200,000"), ]
max.records.clean  <- skipped.clean.df[ which(skipped.clean.df$Reason_skipped      == "Number of records > 200,000"), ]


name.records.spp    <- skipped.species.df[ which(skipped.species.df$Reason_skipped == "Possible incorrect nomenclature"), ]
name.records.grow   <- skipped.grow.df[ which(skipped.grow.df$Reason_skipped       == "Possible incorrect nomenclature"), ]
name.records.clean  <- skipped.clean.df[ which(skipped.clean.df$Reason_skipped       == "Possible incorrect nomenclature"), ]


no.records.spp     <- skipped.species.df[ which(skipped.species.df$Reason_skipped == "No GBIF records"), ]
no.records.grow    <- skipped.grow.df[ which(skipped.grow.df$Reason_skipped == "No GBIF records"), ]
no.records.clean   <- skipped.clean.df[ which(skipped.clean.df$Reason_skipped == "No GBIF records"), ]


## Create lists for each category
max.records.spp.list  = unique(as.character(max.records.spp$Taxa))
name.records.spp.list = unique(as.character(name.records.spp$Taxa))
no.records.spp.list   = unique(as.character(no.records.spp$Taxa))

max.records.grow.list  = unique(as.character(max.records.grow$Taxa))
name.records.grow.list = unique(as.character(name.records.grow$Taxa))
no.records.grow.list   = unique(as.character(no.records.grow$Taxa))

max.records.clean.list  = unique(as.character(max.records.clean$Taxa))
name.records.clean.list = unique(as.character(name.records.clean$Taxa))
no.records.clean.list   = unique(as.character(no.records.clean$Taxa))


## save lists just in case
## save(skipped.species.df, file = paste("./data/base/HIA_LIST/GBIF/skipped_species.RData", sep = ""))
## save(skipped.genera.df,  file = paste("./data/base/HIA_LIST/GBIF/skipped_genera.RData",  sep = ""))
## save(skipped.ALA.df,     file = paste("./data/base/HIA_LIST/GBIF/skipped_ALA.RData",     sep = ""))


## have a look at the list of skipped species
View(skipped.species.df)


## which of these species are on the top 200?
skipped.200.spp = merge(spp.200, skipped.species.df, by = "Species", all = FALSE)
kable(skipped.200.spp)





#########################################################################################################################
## 5). DOWNLOAD SPECIES WITH >200k RECORDS USING THE GBIF API
#########################################################################################################################



## Read in each file as text?
# Magnolia.g = gbif('Magnolia grandiflora', download = TRUE)
# Magnolia.g.names = sort(names(Magnolia.g))
# 
# 
# Betula.p = read.csv("./data/base/HIA_LIST/GBIF/SPECIES/Betula_pendula.csv", stringsAsFactors = FALSE)
# Betula.p.names = sort(names(Betula.p))
# 
# 
# 
# ## How different are the manually downloaded species?
# setdiff(Magnolia.g.names, Betula.p.names)
# setdiff(Betula.p.names, Magnolia.g.names)
# setdiff(gbif.keep, Betula.p.names)
# 
# 
# ##
# unique(Magnolia.g$country)
# unique(Betula.p$countryCode)
# unique(Betula.p$coordinateUncertaintyInMeters)
# unique(Betula.p$continent)


## Save all files as .Rdata, then try concatenating them.



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
## OUTSTANDING DOWNLOADING TASKS:
#########################################################################################################################


## Got the species with > 200k records

## Keep in touch with Renee RE the list





#########################################################################################################################
############################################  END OF GBIF DOWLOAD CODE ################################################## 
#########################################################################################################################