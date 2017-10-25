#########################################################################################################################
############################################ CHECK TAXONOMY OF RETURNED DATA ############################################ 
#########################################################################################################################


#########################################################################################################################
## 1). RUN TAXONLOOKUP
#########################################################################################################################

## Load data
load("./data/base/HIA_LIST/GBIF/GBIF_LAND_POINTS.RData")
returned.taxa = unique(GBIF_LAND_POINTS$scientificName)


########################################################################################################################
## Use taxonlookup to check the taxonomy
HIA.SPP.LOOKUP = lookup_table(returned.taxa, by_species = TRUE, missing_action = "NA")    ## convert rows to column and merge
HIA.SPP.LOOKUP = setDT(HIA.SPP.LOOKUP , keep.rownames = TRUE)[]
HIA.SPP.LOOKUP = dplyr::rename(HIA.SPP.LOOKUP, Binomial = rn)


## So we have searched for every species on the list
length(returned.taxa) - dim(HIA.SPP.LOOKUP)[1]
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
## 5). RUN TAXONLOOKUP
#########################################################################################################################




