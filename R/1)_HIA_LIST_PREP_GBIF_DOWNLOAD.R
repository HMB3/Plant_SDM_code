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
library(Rcpp)
library(raster)
library(rgdal)
library(plyr)
library(dplyr)
library(sfsmisc)
library(spatstat)
library(data.table)

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

p <- c('ff',    'things',         'raster',    'dismo',        'sp',           'latticeExtra', 'data.table', 
       'rgdal', 'rgeos',          'gdalUtils', 'rmaxent',      'readr',        'dplyr',        'tidyr',
       'readr', 'rnaturalearth',  'rasterVis', 'RColorBrewer', 'latticeExtra', 'parallel')


## Require packages
sapply(p, require, character.only = TRUE)

## source functions
source('./R/GREEN_CITIES_FUNCTIONS.R')
source('./R/SDM_FUNCTIONS.R')



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
source('./R/HIA_CLEAN_MATCHING.R')


## Run the download function on the species and genera lists these functions need to download at least one file, or they will return NULL
all.taxa        = unique(c(spp, spp.grow, spp.clean))
str(all.taxa)
class(all.taxa)
skipped.taxa    = download_GBIF_all_species(all.taxa)        ## saves each spp as .Rdata file, returning list of skipped spp 





#########################################################################################################################
## 2). CHECK THE SPECIES WHICH WERE SKIPPED? 
#########################################################################################################################


## Convert the list of skipped species and genera into a dataframe
skipped.taxa.df    <- data.frame(matrix(unlist(skipped.taxa), nrow = length(skipped.taxa), byrow = TRUE))

## Split the reason for skipping and the species into separate columns
skipped.taxa.df    <- cSplit(skipped.taxa.df,    1:ncol(skipped.taxa.df),    sep = "|", stripWhite = TRUE, type.convert = FALSE)


## Update the column names
colnames(skipped.taxa.df)[1] <- "Reason_skipped"
colnames(skipped.taxa.df)[2] <- "Species"

 
## Get subset for each type
max.records.taxa    <- skipped.taxa.df[ which(skipped.taxa.df$Reason_skipped == "Number of records > 200,000"), ]
name.records.taxa   <- skipped.taxa.df[ which(skipped.taxa.df$Reason_skipped == "Possible incorrect nomenclature"), ]
no.records.taxa     <- skipped.taxa.df[ which(skipped.taxa.df$Reason_skipped == "No GBIF records"), ]

## Create lists for each category
max.records.taxa.list  = unique(as.character(max.records.taxa$Taxa))


## Save list just in case
save(skipped.taxa.df, file = paste("./data/base/HIA_LIST/GBIF/skipped_taxa.RData", sep = ""))


## Have a look at the list of skipped species
View(skipped.taxa.df)
View(max.records.taxa)


## Which of these species are on the top 200? These species are all those which cannot be resolved to the species 
## level by the GBIF taxonomy, only the geunus level. So we can't feasibly model them...
intersect(skipped.taxa.df$Species, spp.200$Binomial)
# skipped.200.spp = merge(spp.200, skipped.taxa.df, by = "Species", all = FALSE)
# kable(skipped.200.spp)





#########################################################################################################################
## 3). RUN GBIF PROCESSING CODE
#########################################################################################################################


#########################################################################################################################
## Run all the GBIF code 
source("./R/2)_HIA_GBIF_DATA_COMBINE.R")
source("./R/3)_GBIF_DATA_FILTER.R")
source("./R/4)_ALA_COMBO_GBIF_RECORDS_NUMERICAL_SUMMARY.R")





#########################################################################################################################
## DOWNLOAD SPECIES WITH >200k RECORDS USING THE GBIF API
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


## Get the species with > 200k records

## Keep in touch with Paul and Renee RE the list





#########################################################################################################################
############################################  END OF GBIF DOWLOAD CODE ################################################## 
#########################################################################################################################