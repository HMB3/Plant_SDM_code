#########################################################################################################################
## Load Manuels' list
AUS.GROW = read.csv("./data/base/HIA_LIST/HIA/DATABSE_AUS_SPP_GROW.csv", stringsAsFactors = FALSE) ## change to directory
HIA.SPP  = read.csv("./data/base/HIA_LIST/HIA/HIA_SPP_LOOKUP.csv", stringsAsFactors = FALSE)
names(AUS.GROW)
names(HIA.SPP)
unique(AUS.GROW$scientific_name)


## Extracing only unique Rows based on only 1 Column: "scientific_name"
## The last bit gets the columns we want and re-arranges them
dim(AUS.GROW)
AUS.GROW.UNIQUE = subset(AUS.GROW, !duplicated(AUS.GROW$scientific_name))[, c("scientific_name", "state", "SUA", 
                                                                              "LGA", "origin", "common_name", "ref")]

## Check the dimensions
dim(AUS.GROW)
dim(AUS.GROW.UNIQUE)
length(unique(AUS.GROW$scientific_name))  ## looks right!


## First what are the differences between the lists?
setdiff(AUS.GROW.UNIQUE$scientific_name, HIA.SPP$searchTaxon)  ## only on the grown list
setdiff(HIA.SPP$searchTaxon, AUS.GROW.UNIQUE$scientific_name)  ## only on the HIA list


## Now merge them:  change to "all = TRUE" to get all rows in X
HIA.JOIN = dplyr::rename(HIA.SPP, scientific_name = searchTaxon)[, c(1, 3:19)]
GROW.HIA = merge(AUS.GROW.UNIQUE, HIA.JOIN, by = "scientific_name", all = FALSE)[, c(1, 6, 21:24, 2:4, 7:20)]


## Check the data
names(GROW.HIA)
dim(GROW.HIA)
head(GROW.HIA)


## get fi


## Can save it out
write.csv(GROW.HIA, "./data/base/HIA_LIST/GBIF/GROW_HIA_OVERLAP.csv", row.names = FALSE)




