#########################################################################################################################
## Load Manuels' list
HIA.SPP      = read.csv("./data/base/HIA_LIST/HIA/HIA_SPP_LOOKUP.csv", stringsAsFactors = FALSE)
RAW.HIA      = read.csv("./data/base/HIA_LIST/HIA/GREEN_CITIES_DRAFT_LIST_2709_2017.csv", stringsAsFactors = FALSE)
GROWING      = read.csv("./data/base/HIA_LIST/HIA/database_aus_sp_growing.csv", stringsAsFactors = FALSE)


## can also use a "trim" in R
trim <- function (x) gsub("^\\s+|\\s+$", "", x)
HIA.RAW$Species = trim(HIA.RAW$Species)
GROWING$scientific_name = trim(GROWING$scientific_name)
RAW.HIA.SPP = RAW.HIA.SPP$Species


## check the names
names(GROWING)
names(HIA.SPP)
length(RAW.HIA.SPP)
unique(GROWING$scientific_name)


## Extracing only unique Rows based on only 1 Column: "scientific_name"
## The last bit gets the columns we want and re-arranges them
dim(GROWING)
GROWING.UNIQUE = subset(GROWING, !duplicated(GROWING$scientific_name))[, c("scientific_name", "state", "SUA", 
                                                                           "LGA", "origin", "common_name", "ref")]

## Check the dimensions
dim(GROWING)
dim(GROWING.UNIQUE)
length(unique(GROWING$scientific_name))  ## looks right!


## First what are the differences between the lists?
species.on.grow.not.hia =  setdiff(GROWING.UNIQUE$scientific_name, RAW.HIA.SPP)  ## only on the grown list
species.on.hia.not.grow =  setdiff(RAW.HIA.SPP, GROWING.UNIQUE$scientific_name)  ## only on the HIA list

## intersect
intersect(RAW.HIA.SPP, GROWING.UNIQUE$scientific_name)


## First what are the differences between the lists?
# setdiff(GROWING.UNIQUE$scientific_name, HIA.SPP$searchTaxon)  ## only on the grown list
# setdiff(HIA.SPP$searchTaxon, GROWING.UNIQUE$scientific_name)  ## only on the HIA list
# intersect(HIA.SPP$searchTaxon, GROWING.UNIQUE$scientific_name) 

## Now merge them:  change to "all = TRUE" to get all rows in X
HIA.JOIN = dplyr::rename(HIA.RAW, scientific_name = Species) #[, c(1, 3:19)]
GROW.HIA = merge(GROWING.UNIQUE, HIA.JOIN, by = "scientific_name", all = FALSE)#[, c(1, 6, 21:24, 2:4, 7:20)]


## Check the data
names(GROW.HIA)
dim(GROW.HIA)
head(GROW.HIA)


## get fi


## Can save it out
write.csv(GROW.HIA, "./data/base/HIA_LIST/GBIF/GROW_HIA_OVERLAP.csv", row.names = FALSE)




