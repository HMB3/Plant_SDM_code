#########################################################################################################################
##############################  DOWNLOAD GLOBAL SPECIES RECORDS FOR HORT AUS LIST ####################################### 
#########################################################################################################################


#########################################################################################################################
## 1). READ IN DRAFT HIA LIST AND CLEAN
#########################################################################################################################


## this list derives from all species and varieties sold anywhere in Australia in the last 5 years. Anthony Manea cleaned 
## up the data and cross-linked to growth form and exotic/native status and derived a list of ~1000 species that are the 
## most commonly sold, covering the right ratio of growth forms, regional representation and native/exotic
spp.list = read.csv("./data/base/HIA_LIST/HIA/GREEN_CITIES_DRAFT_LIST.csv", stringsAsFactors = FALSE)
top.200  = read.csv("./data/base/HIA_LIST/HIA/HIA_TOP_200.csv", stringsAsFactors = FALSE)


## have a look
dim(spp.list)
str(spp.list)
head(spp.list)


## also, add the "Top 200 species in here"
spp.200          = top.200[c("Species")]
spp.200$Top_200  = "TRUE"

spp.list = merge(spp.list, spp.200, by = "Species", all.x = TRUE) 
spp.list$Top_200[is.na(spp.list$Top_200)] <- "FALSE"
spp.list$Origin <- gsub(" ",  "", spp.list$Origin)


## check
str(spp.list)
head(spp.list)
unique(spp.list$Top_200)
unique(spp.list$Origin)





#########################################################################################################################
## WHAT ARE SOME USEFUL MEASURES OF THE DATASET?
#########################################################################################################################


## ~ HALF the total species are native, ~47% of the top 200
dim(subset(spp.list, Origin == "Native"))[1]/dim(spp.list)[1]*100
dim(subset(top.200, Origin == "Native"))[1]/dim(top.200)[1]*100


#########################################################################################################################
## 2). DRAFT LIST PREP: THIS NEEDS TO CHANGE IN CONSULTATION WITH RACH, PAUL, ETC
#########################################################################################################################


## for now, remove all the subspecies, varieties etc.
## but check if some of them work on GBIF, e.g. a separate list of varities
spp.list$Species <- gsub("spp.", "", spp.list$Species)
spp.list$Species <- gsub(" x ",  " ", spp.list$Species)
spp.list$Species <- gsub("  ",  " ", spp.list$Species)


## then remove the varieties and subspecies
spp.list$Species = vapply(lapply(strsplit(spp.list$Species, " "),
                                 unique), paste, character(1L), collapse = " ")


## then get just the first two words (again cleaning up the subspecies, and single genera)
spp.list$Species = vapply(lapply(strsplit(spp.list$Species, " "), 
                                 string_fun_first_two_words), paste, character(1L), collapse = " ")


## if exclude NA, spp.list = spp.list[!grepl("NA", spp.list$Species),]
spp.list = spp.list[with(spp.list, order(Species)), ] 


## check
str(spp.list)
View(spp.list)


## Now create list of HIA taxa. 768 unique species, mius the corrections, etc. 
spp.list = spp.list[with(spp.list, order(Species)), ]
spp = unique(as.character(spp.list$Species))                         ## 
str(spp)   ## why 680? later check on what happens with the different queries
head(spp, 50)
tail(spp, 50)


## also make a genus list?
genera = unique(vapply(lapply(strsplit(spp.list$Species, " "), 
                              string_fun_first_word), paste, character(1L), collapse = " "))


## check
str(genera)
head(genera, 50)
tail(genera, 50)


## also, for Stuarts code EG, I need a df not a list. Get just the rows of spp.list which have
## unique species names. 
DRAFT.HIA.TAXA = spp.list[!duplicated(spp.list$Species), ]
dim(DRAFT.HIA.TAXA)
head(DRAFT.HIA.TAXA)


########################################################################################################################
## Try using taxonlookup to check the taxonomy
DRAFT.TAXA.LOOKUP = lookup_table(DRAFT.HIA.TAXA[["searchTaxon"]], by_species = TRUE) ## convert rows to column and merge
DRAFT.TAXA.LOOKUP = setDT(DRAFT.TAXA.LOOKUP , keep.rownames = TRUE)[]
DRAFT.TAXA.LOOKUP = rename(DRAFT.TAXA.LOOKUP, searchTaxon = rn)
head(DRAFT.TAXA.LOOKUP)


