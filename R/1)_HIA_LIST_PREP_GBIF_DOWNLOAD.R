#########################################################################################################################
##############################  DOWNLOAD GLOBAL SPECIES RECORDS FOR HORT AUS LIST ####################################### 
#########################################################################################################################


## This code downloads all occurrence records for species on the Horticulture Australia list, from the GBIF database. 
## The trick here will be substituing the functions from the ALA4R package for the rgbif functions...
## Matching the core fields could be tricky for some taxa.


###########################################################################################
## remote directoty for raster data
## login from workstation: ## mqauth.uni.mq.edu.au\MQ20174568 ## Popple2016
## \\SCI-7910
## \\sci-7910\F\data\worldclim\world\0.5\bio


#########################################################################################################################
## 1). DOWNLOAD RECORDS FROM GBIF USING HIA LIST
#########################################################################################################################


## Now run loops to dowload species in the "spp" list from GBIF. Not including any data quality checks here, just 
## downloading everything...
load("./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT.RData")
source('./R/HIA_LIST_MATCHING.R')
source('./R/HIA_CLEAN_MATCHING.R')


## Create one big list of all the taxa
all.taxa = unique(c(spp, spp.grow, spp.clean))
length(all.taxa)   


########################################################################################################################
## Use taxonlookup to check the taxonomy
HIA.SPP.LOOKUP = lookup_table(all.taxa, by_species = TRUE, missing_action = "NA")    ## convert rows to column and merge
HIA.SPP.LOOKUP = setDT(HIA.SPP.LOOKUP , keep.rownames = TRUE)[]
HIA.SPP.LOOKUP = dplyr::rename(HIA.SPP.LOOKUP, Binomial = rn)


## So we have searched for every species on the list
length(all.taxa) - dim(HIA.SPP.LOOKUP)[1]
head(HIA.SPP.LOOKUP)                                                                 ## Can merge on the bilogical data here..
View(HIA.SPP.LOOKUP)


## But, just get the species that don't match (i.e. the NA rows...)
HIA.SPP.LOOKUP.MATCH  = na.omit(HIA.SPP.LOOKUP)
HIA.SPP.TAXO.ERRORS  <- HIA.SPP.LOOKUP[rowSums(is.na(HIA.SPP.LOOKUP)) > 0,]
dim(HIA.SPP.TAXO.ERRORS)
head(HIA.SPP.TAXO.ERRORS)


## Write each table out to file:
write.csv(HIA.SPP.LOOKUP,       "./data/base/TRAITS/HIA_SPP_LOOKUP.csv",       row.names = FALSE)
write.csv(HIA.SPP.LOOKUP.MATCH, "./data/base/TRAITS/HIA_SPP_LOOKUP_MATCH.csv", row.names = FALSE)
write.csv(HIA.SPP.TAXO.ERRORS,  "./data/base/TRAITS/HIA_SPP_TAXO.ERRORS.csv",  row.names = FALSE)


########################################################################################################################
## Finally, check the taxonomy for the data already downloaded against this list
## How can we confirm the taxonomy is ok? Load big dataset in and check
taxon.match          = intersect(HIA.SPP.LOOKUP.MATCH$Binomial, COMBO.NICHE.CONTEXT$searchTaxon)
popular.spp.match    = intersect(spp, HIA.SPP.LOOKUP.MATCH$Binomial)                             ## Check this with Rach
Top.200.match        = intersect(top.200$Species, HIA.SPP.LOOKUP.MATCH$Binomial) 


## So all the species on the downloadedlist are on the taxonomically matched list
length(spp %in% HIA.SPP.LOOKUP.MATCH$Binomial) - length(spp)


## Record the differences
taxon.difference = setdiff(HIA.SPP.LOOKUP.MATCH$Binomial, COMBO.NICHE.CONTEXT$searchTaxon)


## In future, use the checked taxa as the list for downloading
taxon.search = as.list(HIA.SPP.LOOKUP.MATCH$Binomial)



#########################################################################################################################
## Run the download function on the species and genera lists these functions need to download at least one file, or they 
## will return NULL
skipped.taxa    = download_GBIF_all_species(all.taxa)     ## saves each spp as .Rdata file, returning list of skipped spp





#########################################################################################################################
## 2). CHECK THE SPECIES WHICH WERE SKIPPED? 
#########################################################################################################################


## Converting the lists of skipped species and genera into a dataframe
skipped.taxa.df <- data.frame(matrix(unlist(skipped.taxa), nrow = length(skipped.taxa), byrow = TRUE))


## Split the reason and the species into separate columns
skipped.taxa.df    <- cSplit(skipped.taxa.df,    1:ncol(skipped.taxa.df),    sep = "|", stripWhite = TRUE, type.convert = FALSE)


## Update names
colnames(skipped.taxa.df)[1] <- "Reason_skipped"
colnames(skipped.taxa.df)[2] <- "Species"


## Get subset for each type
max.records.taxa    <- skipped.taxa.df[ which(skipped.taxa.df$Reason_skipped == "Number of records > 200,000"), ]
name.records.taxa   <- skipped.taxa.df[ which(skipped.taxa.df$Reason_skipped == "Possible incorrect nomenclature"), ]
no.records.taxa     <- skipped.taxa.df[ which(skipped.taxa.df$Reason_skipped == "No GBIF records"), ]


## Create lists for each category
# max.records.spp.list  = unique(as.character(max.records.spp$Taxa))
# name.records.spp.list = unique(as.character(name.records.spp$Taxa))
# no.records.spp.list   = unique(as.character(no.records.spp$Taxa))

## save lists just in case
## save(skipped.species.df, file = paste("./data/base/HIA_LIST/GBIF/skipped_species.RData", sep = ""))
## save(skipped.grow.df,  file = paste("./data/base/HIA_LIST/GBIF/skipped_growing.RData",  sep = ""))
## save(skipped.clean.df,     file = paste("./data/base/HIA_LIST/GBIF/skipped_clean.RData",     sep = ""))


## have a look at the list of skipped species
View(skipped.species.df)


## which of these species are on the top 200?
skipped.200.spp = merge(spp.200, skipped.species.df, by = "Species", all = FALSE)
kable(skipped.200.spp)





#########################################################################################################################
## Run all the GBIF code 
# source("./R/2)_HIA_GBIF_DATA_COMBINE.R")
# source("./R/3)_GBIF_DATA_FILTER.R")
# source("./R/4)_ALA_COMBO_GBIF_RECORDS_NUMERICAL_SUMMARY.R")





#########################################################################################################################
## 3). DOWNLOAD SPECIES WITH >200k RECORDS USING THE GBIF API
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


## Get the extra species with > 200k records

## Get the missing species on Manuel's list.





#########################################################################################################################
############################################  END OF GBIF DOWLOAD CODE ################################################## 
#########################################################################################################################