#########################################################################################################################
##############################  DOWNLOAD GLOBAL SPECIES RECORDS FOR HORT AUS LIST ####################################### 
#########################################################################################################################


## This code downloads all occurrence records for species on the Horticulture Australia list, from the GBIF database. 
## The trick here will be substituing the functions from the ALA4R package for the rgbif functions...
## Matching the core fields could be tricky for some taxa. Also, 


## setup 
library(gtools)
library(devtools)
#devtools::install_github("ropensci/rgbif")
library(data.table)
library(RCurl)
library(Rcpp)
library(raster)
library(rgdal)
library(plyr)

library(SDMTools)
library(rmaxent)
library(dismo)
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


## source functions
source('./R/GREEN_CITIES_FUNCTIONS.R')





#########################################################################################################################
## 1). READ IN DRAFT HIA LIST AND CLEAN
#########################################################################################################################


## this list derives from all species and varieties sold anywhere in Australia in the last 5 years. Anthony Maneahas cleaned 
## up the data and cross-linked to growth form and exotic/native status and derived a list of ~1000 species that are the most 
## commonly sold, covering the right ratio of growth forms, regional representation and native/exotic
spp.list = read.csv("./data/base/HIA_LIST/HIA/GREEN_CITIES_DRAFT_LIST.csv", stringsAsFactors = FALSE)
top.200  = read.csv("./data/base/HIA_LIST/HIA/HIA_TOP_200.csv", stringsAsFactors = FALSE)


## have a look
dim(spp.list)
str(spp.list)
head(spp.list)


## also, add the "Top 200 species in here"
spp.200          = top.200[c("Species")]
spp.200$Top_200  = "TRUE"

spp.list = merge(spp.list, spp.200, by = "Species", all.x = TRUE) 
spp.list$Top_200[is.na(spp.list$Top_200)] <- "FALSE"
spp.list$Origin <- gsub(" ",  "", spp.list$Origin)


## check
str(spp.list)
head(spp.list)
unique(spp.list$Top_200)
unique(spp.list$Origin)


#########################################################################################################################
## WHAT ARE SOME USEFUL MEASURES OF THE DATASET?
#########################################################################################################################


## ~ HALF the total species are native, ~47% of the top 200
dim(subset(spp.list, Origin == "Native"))[1]/dim(spp.list)[1]*100
dim(subset(top.200, Origin == "Native"))[1]/dim(top.200)[1]*100


#########################################################################################################################
## DRAFT CLEAN: THIS NEEDS TO CHANGE IN CONSULTATION WITH RACH, PAUL, ETC
#########################################################################################################################


## for now, remove all the subspecies, varieties etc.
## but check if some of them work on GBIF, e.g. a separate list of varities
spp.list$Species <- gsub("spp.", "", spp.list$Species)
spp.list$Species <- gsub(" x ",  " ", spp.list$Species)
spp.list$Species <- gsub("  ",  " ", spp.list$Species)


## then remove the varieties and subspecies
spp.list$Species = vapply(lapply(strsplit(spp.list$Species, " "),
                                                 unique), paste, character(1L), collapse = " ")


## then get just the first two words (again cleaning up the subspecies, and single genera)
spp.list$Species = vapply(lapply(strsplit(spp.list$Species, " "), 
                                 string_fun_first_two_words), paste, character(1L), collapse = " ")


## if exclude NA, spp.list = spp.list[!grepl("NA", spp.list$Species),]
spp.list = spp.list[with(spp.list, order(Species)), ] 


## check
str(spp.list)
View(spp.list)


## Now create list of HIA taxa. 768 unique species, mius the corrections, etc. 
spp.list = spp.list[with(spp.list, order(Species)), ]
spp = unique(as.character(spp.list$Species))
str(spp)   ## why 680? later check on what happens with the different queries
head(spp, 50)
tail(spp, 50)


## also make a genus list?
genera = unique(vapply(lapply(strsplit(spp.list$Species, " "), 
                              string_fun_first_word), paste, character(1L), collapse = " "))


## check
str(genera)
head(genera, 50)
tail(genera, 50)


## also, for Stuarts code EG, I need a df not a list. Get just the rows of spp.list which have
## unique species names. 
draft.taxa = spp.list[!duplicated(spp.list$Species), ]
dim(draft.taxa)
head(draft.taxa)





#########################################################################################################################
## 2). DOWNLOAD RECORDS FROM GBIF USING HIA LIST
#########################################################################################################################


## now run lOOPs to dowload species in the "spp" list from GBIF
## not including any data quality checks here, just downloading everything ##
## 


## set a few global variables to be used inside the functions...
#GBIF.download.limit = 200000
skip.spp.list       = list()
skip.gen.list       = list()


## run the download function on the species and genera lists
## why is this loop producing NULL?
skipped.species = download_GBIF_all_species(spp)    ## saves each spp as .Rdata file, returning list of skipped spp 
skipped.genera  = download_GBIF_all_genera(genera)  ## saves each gen as .Rdata file, returning list of skipped genera 


## now convert the lists of skipped species and genera into a dataframe
str(skipped.species)
str(skipped.genera)


## converting the lists of skipped species and genera into a dataframe
skipped.species <- data.frame(matrix(unlist(skipped.species), nrow = length(skipped.species), byrow = TRUE))
skipped.genera  <- data.frame(matrix(unlist(skipped.genera),  nrow = length(skipped.spp), byrow = TRUE))


## split the reason and the species into separate columns
skipped.species <- cSplit(skipped.species, 1:ncol(skipped.species), sep = "|", stripWhite = TRUE, type.convert = FALSE)


## update names
colnames(skipped.species)[1] <- "Reason_skipped"
colnames(skipped.species)[2] <- "Taxa"


## check
str(skipped.species)
head(skipped.species)
tail(skipped.species)
View(skipped.species)


## get subset for each type
max.records  <- skipped.species[ which(skipped.species$Reason_skipped == "Number of records > 200,000"), ]
name.records <- skipped.species[ which(skipped.species$Reason_skipped == "Possible incorrect nomenclature"), ]
no.records   <- skipped.species[ which(skipped.species$Reason_skipped == "No GBIF records"), ]


## create lists for each category
max.records.spp  = unique(as.character(max.records$Taxa))
name.records.spp = unique(as.character(name.records$Taxa))
no.records.spp   = unique(as.character(no.records$Taxa))





#########################################################################################################################
############################################  END OF GBIF DOWLOAD CODE ################################################## 
#########################################################################################################################