#########################################################################################################################
##############################  DOWNLOAD GLOBAL SPECIES RECORDS FOR HORT AUS LIST ####################################### 
#########################################################################################################################


#########################################################################################################################
## This code downloads all occurrence records for species on the Evergreen list using the GBIF database..................


#########################################################################################################################
## 1). DOWNLOAD RECORDS FROM GBIF USING LIST
#########################################################################################################################


## Change the ALA biocache to the old one
# ala_config(caching = "off")
# ala_config(verbose = TRUE)
# 
# 
# server_config = getOption("ALA4R_server_config")
# 
# server_config$base_url_biocache = "https://biocache-lb.ala.org.au/ws"
# server_config$base_url_biocache_download = "https://biocache-lb.ala.org.au/biocache-download/"
# server_config


## Create a list of all the taxa to download.
all.taxa     =  GBIF.spp
all.taxa.rev =  all.taxa[rev(order(all.taxa))]

#message (GBIF_path)

#########################################################################################################################
## Run the download function on the species lists. This function needs to download at least one file, or they 
## will return NULL. Saves each spp as .Rdata file, returning list of skipped spp
# TEMP DISABLE for katana - we have all the records so no need to check again
if (FALSE) {
    GBIF.taxa    = download_GBIF_all_species(species_list = all.taxa, 
                                            path         = GBIF_path) ## insert path 


    ALA.taxa        = download_ALA_all_species(species_list = all.taxa.rev, 
                                           path         = ALA_path)

}



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
