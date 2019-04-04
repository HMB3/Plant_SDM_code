#########################################################################################################################
##############################  DOWNLOAD GLOBAL SPECIES RECORDS FOR HORT AUS LIST ####################################### 
#########################################################################################################################


#########################################################################################################################
## This code downloads all occurrence records for species on the Evergreen list using the GBIF database..................


#########################################################################################################################
## 1). DOWNLOAD RECORDS FROM GBIF USING LIST
#########################################################################################################################


## Create a list of all the taxa to download.
all.taxa     =  GBIF.spp
all.taxa.rev =  all.taxa[rev(order(all.taxa))]
#message (GBIF_path)
#
#

# 
# 5651 Records returned for Zamia loddigesii
# Error in check_status_code(status_code, on_redirect = on_redirect, on_client_error = on_client_error,  : 
#                              HTTP status code 504 received.
#                            Either there was an error with the request, or the servers may be down (try again later). 
#                            If this problem persists please notify the ALA4R maintainers by lodging an issue at 
#                            https://github.com/AtlasOfLivingAustralia/ALA4R/issues/ or emailing support@ala.org.au


# 11111 Records returned for Acacia mangium
# Downloading GBIF records for Acacia mearnsii using rgbif :: occ_data
# Error in curl::curl_fetch_memory(x$url$url, handle = x$url$handle) : 
#   transfer closed with outstanding read data remaining



#########################################################################################################################
## Run the download function on the species lists. This function needs to download at least one file, or they 
## will return NULL. Saves each spp as .Rdata file, returning list of skipped spp
# TEMP DISABLE for katana - we have all the records so no need to check again
GBIF.taxa    = download_GBIF_all_species(species_list = all.taxa, 
                                         path         = GBIF_path) ## insert path 


ALA.taxa     = download_ALA_all_species(species_list = all.taxa.rev, 
                                        path         = ALA_path)





#########################################################################################################################
## 2). READ IN DODGY SPECIES
#########################################################################################################################


##
# Eremophila.bignoniiflora = read.csv("./data/base/HIA_LIST/ALA/Eremophila_bignoniiflora .csv", stringsAsFactors = FALSE, check.names=FALSE)
# Eucalyptus.populnea      = read.csv("./data/base/HIA_LIST/ALA/Eucalyptus_populnea.csv",       stringsAsFactors = FALSE, check.names=FALSE)
# 
# Eremophila.bignoniiflora = convert_ala_cols_no_space(Eremophila.bignoniiflora)
# Eucalyptus.populnea      = convert_ala_cols_no_space(Eucalyptus.populnea)
# 
# 
# save(Eremophila.bignoniiflora,  file = paste0(ALA_path, 'Eremophila bignoniiflora_ALA_records.RData'))
# save(Eucalyptus.populnea,       file = paste0(ALA_path, 'Eucalyptus populnea_ALA_records.RData'))


#########################################################################################################################
## 2). CHECK THE SPECIES WHICH WERE SKIPPED? 
#########################################################################################################################


# ## Converting the lists of skipped species and genera into a dataframe
# skipped.taxa.df <- data.frame(matrix(unlist(skipped.taxa), nrow = length(skipped.taxa), byrow = TRUE))
# 
# 
# ## Split the reason and the species into separate columns
# skipped.taxa.df    <- cSplit(skipped.taxa.df,    1:ncol(skipped.taxa.df),    sep = "|", stripWhite = TRUE, type.convert = FALSE)
# 
# 
# ## Update names
# colnames(skipped.taxa.df)[1] <- "Reason_skipped"
# colnames(skipped.taxa.df)[2] <- "Species"
# 
# 
# ## Get subset for each type
# max.records.taxa    <- skipped.taxa.df[ which(skipped.taxa.df$Reason_skipped == "Number of records > 200,000"), ]
# name.records.taxa   <- skipped.taxa.df[ which(skipped.taxa.df$Reason_skipped == "Possible incorrect nomenclature"), ]
# no.records.taxa     <- skipped.taxa.df[ which(skipped.taxa.df$Reason_skipped == "No GBIF records"), ]
# 
# 
# ## have a look at the list of skipped species
# View(skipped.species.df)




#########################################################################################################################
## OUTSTANDING DOWNLOADING TASKS:
#########################################################################################################################


## Get the extra species with > 200k records






#########################################################################################################################
############################################  END OF GBIF DOWLOAD CODE ################################################## 
#########################################################################################################################
