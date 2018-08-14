#########################################################################################################################
##############################  DOWNLOAD GLOBAL SPECIES RECORDS FOR HORT AUS LIST ####################################### 
#########################################################################################################################


#########################################################################################################################
## This code downloads all occurrence records for species on the Evergreen list using the GBIF database. .


#########################################################################################################################
## source packages and functions
source('./R/HIA_LIST_MATCHING.R')


#########################################################################################################################
## 1). DOWNLOAD RECORDS FROM GBIF USING LIST
#########################################################################################################################


## Create one big list of all the taxa
all.taxa     =  GBIF.spp
all.taxa.rev =  all.taxa[rev(order(all.taxa))]


########################################################################################################################
## Use lookup to get the order, group, etc.
SPP.LOOKUP = lookup_table(all.taxa, by_species = TRUE, missing_action = "NA")    ## convert rows to column and merge
SPP.LOOKUP = setDT(SPP.LOOKUP , keep.rownames = TRUE)[]
SPP.LOOKUP = dplyr::rename(SPP.LOOKUP, Binomial = rn)

head(SPP.LOOKUP)
View(SPP.LOOKUP)


## But, just get the species that don't match (i.e. the NA rows...)
SPP.LOOKUP.MATCH  = na.omit(SPP.LOOKUP)
SPP.TAXO.ERRORS  <- SPP.LOOKUP[rowSums(is.na(SPP.LOOKUP)) > 0,]
dim(SPP.TAXO.ERRORS)
head(SPP.TAXO.ERRORS)


########################################################################################################################
## Use Taxonstand to check the taxonomy
SPP.TAXO <- TPL(GBIF.spp, infra = TRUE,
                corr = TRUE, repeats = 100)  ## to stop it timing out...


## How do we use this table to clean the taxonomy? Best to do it in concert with Alessandro
names(SPP.TAXO)
names(SPP.TAXO)[names(SPP.TAXO) == 'Taxon'] <- 'Binomial'
SPP.LOOKUP = join(SPP.LOOKUP, SPP.TAXO)
View(SPP.TAXO)


## Write each table out to file:
write.csv(SPP.LOOKUP,       "./data/base/HIA_LIST/COMBO/SPP_TPL_LOOKUP.csv",       row.names = FALSE)


########################################################################################################################
## To get synonyms of a species  you have to use its taxonconceptkey (GBIF parlance) in the `occurrencelist` function. 
# So lets get the taxonconceptkey first using taxonsearch function. You have to set dataresourcekey=1 to get the sort 
#   of master taxonconceptkey for the species
palustris <- occ_data(scientificName = 'Quercus palustris')
coccinea  <- occ_data(scientificName = 'Quercus coccinea')


##
dim(palustris$data)
dim(coccinea$data)
sort(names(palustris$data))
unique(palustris$data$scientificName)


##
palustris.key <- name_backbone(name = 'Quercus palustris', rank = 'species')$usageKey
coccinea.key  <- name_backbone(name = 'Quercus coccinea',  rank = 'species')$usageKey

palustris.occ = occ_data(taxonKey = palustris.key, limit = 200000)
coccinea.occ  = occ_data(taxonKey = coccinea.key,  limit = 200000)


## Check the dimensions
dim(palustris.occ$data)
dim(coccinea.occ$data)


## 
unique(palustris.occ$data$scientificName)
unique(coccinea.occ$data$scientificName)



#########################################################################################################################
## Run the download function on the species lists. This function needs to download at least one file, or they 
## will return NULL. Saves each spp as .Rdata file, returning list of skipped spp
skipped.taxa    = download_GBIF_all_species(species_list = all.taxa, 
                                            path         = "./data/base/HIA_LIST/GBIF/OCC_SEARCH/") ## insert path 


ALA.taxa    = download_ALA_all_species(species_list = all.taxa, 
                                       path         = "./data/base/HIA_LIST/ALA/TREE_SPECIES/")


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


## have a look at the list of skipped species
View(skipped.species.df)




#########################################################################################################################
## OUTSTANDING DOWNLOADING TASKS:
#########################################################################################################################


## Get the extra species with > 200k records






#########################################################################################################################
############################################  END OF GBIF DOWLOAD CODE ################################################## 
#########################################################################################################################