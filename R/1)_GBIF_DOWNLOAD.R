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
all.taxa = read.csv("./data/base/HIA_LIST/HIA/RISK_BINOMIAL_DF.csv", stringsAsFactors = FALSE)  ## RISK.BINOMIAL
all.taxa = as.list(all.taxa[1])


########################################################################################################################
## Use taxonlookup to check the taxonomy: can use this code to check the planted and growing list, etc
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


## Write each table out to file:
write.csv(SPP.LOOKUP,       "SPP_LOOKUP.csv",       row.names = FALSE)
write.csv(SPP.LOOKUP.MATCH, "SPP_LOOKUP_MATCH.csv", row.names = FALSE)
write.csv(SPP.TAXO.ERRORS,  "SPP_TAXO.ERRORS.csv",  row.names = FALSE)


#########################################################################################################################
## Run the download function on the species lists. This function needs to download at least one file, or they 
## will return NULL. Saves each spp as .Rdata file, returning list of skipped spp
skipped.taxa    = download_GBIF_all_species(species_list = all.taxa, 
                                            path         = "./data/base/HIA_LIST/GBIF/SPECIES/") ## insert path 





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