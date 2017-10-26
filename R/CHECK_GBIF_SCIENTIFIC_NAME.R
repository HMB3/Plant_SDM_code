#########################################################################################################################
############################################ CHECK TAXONOMY OF RETURNED DATA ############################################ 
#########################################################################################################################


#########################################################################################################################
## 1). RUN TAXONLOOKUP
#########################################################################################################################


## Load data
load("./data/base/HIA_LIST/GBIF/GBIF_LAND_POINTS.RData")
dim(GBIF.LAND)


## Get species list
returned.taxa = unique(GBIF.LAND$scientificName)
str(returned.taxa)


########################################################################################################################
## Use taxonlookup to check the taxonomy
GBIF.LOOKUP = lookup_table(returned.taxa, by_species = TRUE, missing_action = "NA")    
GBIF.LOOKUP = setDT(GBIF.LOOKUP , keep.rownames = TRUE)[]
GBIF.LOOKUP = dplyr::rename(GBIF.LOOKUP, Binomial = rn)


## So we have searched for every species on the list
head(GBIF.LOOKUP)                                                               
View(GBIF.LOOKUP)


## But, just get the species that don't match (i.e. the NA rows...)
GBIF.LOOKUP.MATCH  = na.omit(GBIF.LOOKUP)
GBIF.TAXO.ERRORS  <- GBIF.LOOKUP[rowSums(is.na(GBIF.LOOKUP)) > 0,]
dim(GBIF.TAXO.ERRORS)
head(GBIF.TAXO.ERRORS)


## Write each table out to file:
write.csv(GBIF.LOOKUP,       "./data/base/TRAITS/HIA_SPP_LOOKUP.csv",       row.names = FALSE)
write.csv(GBIF.LOOKUP.MATCH, "./data/base/TRAITS/HIA_SPP_LOOKUP_MATCH.csv", row.names = FALSE)
write.csv(GBIF.TAXO.ERRORS,  "./data/base/TRAITS/HIA_SPP_TAXO.ERRORS.csv",  row.names = FALSE)


########################################################################################################################
## Now remove the dodgy species from the existing data...do this in the step 3 code


## So all the species on the downloaded list are on the taxonomically matched list
GBIF.LAND.TAXO = GBIF.LAND[GBIF.LAND$scientificName %in% GBIF.LOOKUP.MATCH$Binomial, ]
dim(GBIF.LAND)[1] - dim(GBIF.LAND.TAXO)[1] ## A difference of 200k records



#########################################################################################################################
#####################################################  END ############################################################## 
#########################################################################################################################