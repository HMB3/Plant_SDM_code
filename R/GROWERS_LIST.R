#########################################################################################################################
############################################  MATCH GROWING AND ALA LISTS ###############################################
#########################################################################################################################


#########################################################################################################################
## Load Manuels' growing list, the list of ~600 binomials from the HIA data and the full list of ~1100 taxa
GROWING      = read.csv("./data/base/HIA_LIST/HIA/database_aus_sp_growing.csv", stringsAsFactors = FALSE)   ## change directory
#HIA.SPP      = read.csv("./data/base/HIA_LIST/HIA/HIA_SPP_LOOKUP.csv", stringsAsFactors = FALSE)
RAW.HIA      = read.csv("./data/base/HIA_LIST/HIA/GREEN_CITIES_DRAFT_LIST_2709_2017.csv", stringsAsFactors = FALSE)


## Can also use a "trim" function in R
trim <- function (x) gsub("^\\s+|\\s+$", "", x)
RAW.HIA$Species = trim(RAW.HIA$Species)
GROWING$scientific_name = trim(GROWING$scientific_name)
RAW.HIA.SPP = RAW.HIA$Species


## check the names
names(GROWING)
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


## First, what are the differences between the growers list and the full HIA list?
species.on.grow.not.hia =  setdiff(GROWING.UNIQUE$scientific_name, RAW.HIA.SPP)  ## only on the grown list
species.on.hia.not.grow =  setdiff(RAW.HIA.SPP, GROWING.UNIQUE$scientific_name)  ## only on the HIA list


## Most of these differences are genuine, rathter than being due to character errors. I would suggest telling 
## Paul that the overlap is not that big...I guess there could still be taxonomic and typing errors?


## Now what is the overlap between the lists?
overlap.hia.grow = intersect(RAW.HIA.SPP, GROWING.UNIQUE$scientific_name)


## This is the old code to match the growing list with my list of 605 unique HIA binomials
# setdiff(GROWING.UNIQUE$scientific_name, HIA.SPP$searchTaxon)  ## only on the grown list
# setdiff(HIA.SPP$searchTaxon, GROWING.UNIQUE$scientific_name)  ## only on the HIA list
# intersect(HIA.SPP$searchTaxon, GROWING.UNIQUE$scientific_name) 

## Now merge them:  change to "all = TRUE" to get all rows in the first data frame
HIA.JOIN = dplyr::rename(RAW.HIA, scientific_name = Species) #[, c(1, 3:19)]
GROW.HIA = merge(GROWING.UNIQUE, HIA.JOIN, by = "scientific_name", all = FALSE)#[, c(1, 6, 21:24, 2:4, 7:20)]


## Check the data
names(GROW.HIA)
dim(GROW.HIA)
head(GROW.HIA)


## Save?
# write.csv(GROW.HIA,                "GROW_HIA_JOIN.csv",           row.names = FALSE)
# write.csv(species.on.grow.not.hia, "species_on_grow_not_hia.csv", row.names = FALSE)
# write.csv(species.on.hia.not.grow, "species_on_hia_not_grow.csv", row.names = FALSE)


## Maybe try matching the RAW.HIA file in excel, after trimming?


#########################################################################################################################
##
#########################################################################################################################


## Now try creating a list for each LGA?
load("")





#########################################################################################################################
#######################################################  TBC ############################################################
#########################################################################################################################