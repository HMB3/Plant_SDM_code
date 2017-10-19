#########################################################################################################################
##################################### MATCH GROWING LIST TO APNI ######################################################## 
#########################################################################################################################


#########################################################################################################################
## 1). READ IN AUS PLANT NAME INDEX
#########################################################################################################################


## Load data and packages
load("./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT.RData")
source('./R/HIA_LIST_MATCHING.R')
source('./R/HIA_CLEAN_MATCHING.R')


## Read in data
APNI.list = read.csv("./data/base/HIA_LIST/HIA/APNI_LIST.csv", stringsAsFactors = FALSE)
dim(APNI.list)


## Check the dimensions of what we are joining on
dim(COMBO.NICHE.CONTEXT)
str(spp.grow)
str(COMBO.NICHE.CONTEXT$searchTaxon)


## These are the weirdos
missing.grow = setdiff(spp.grow, COMBO.NICHE.CONTEXT$searchTaxon) 
## These names could be fixed in the orginal database Manuel sent, then re-run the download (280 extra species.
# APNI.SPECIES = gsub("(L.)",   "", APNI.list$Species)





#########################################################################################################################
## 2). DRAFT 'CLEAN' LIST PREP: REMOVE SUBSPECIES, VARIETIES, ETC
#########################################################################################################################


#########################################################################################################################
## Get binomials
APNI.list$Binomial <- sub('(^\\S+ \\S+).*', '\\1', APNI.list$Species)


## First, remove taxa with "spp", etc: see setdiff at the end for the ones which are removed...
DRAFT.APNI.TAXA         = APNI.list
dim(DRAFT.APNI.TAXA)    ## 

## Don't need to remove weird characters...
## APNI.SPECIES = gsub("(L.)",   "", APNI.list$Species)


#######################################################################################################################
## Which taxa have the second word capitalised, and are not captured by eliminating the "spp."?
EXLUDE.APNI     = grep('^\\S+ [A-Z]', DRAFT.APNI.TAXA$Species, val = TRUE)
DRAFT.APNI.TAXA = DRAFT.APNI.TAXA[!DRAFT.APNI.TAXA$Species %in% EXLUDE.APNI, ]


## Now for GBIF, just get the unique species.........................
APNI.SPP = DRAFT.APNI.TAXA[!duplicated(DRAFT.APNI.TAXA["Binomial"]),]
APNI.SPP = APNI.SPP[, c("Species", "Binomial", "APNI")] 


## Reorder by species
APNI.SPP = APNI.SPP[with(APNI.SPP, order(Binomial)), ]


## Save with both taxonomic columns
save(APNI.SPP, file = paste("./data/base/HIA_LIST/APNI_BINOMIALS.RData", sep = ""))


## Now create list of APNI species. 44752 seems too many, but we will only match to the spatial data 
spp.apni            = unique(as.character(APNI.SPP$Binomial))
length(spp.apni)   





#########################################################################################################################
## 3). JOIN THE APNI NAMES ONTO THE NICHE DATA SET
#########################################################################################################################


#########################################################################################################################
## Just get the columns we need and rename
APNI = APNI.SPP[, c("Binomial", "APNI")]
APNI = dplyr::rename(APNI, searchTaxon = Binomial)
## The species that are different are all exotic, good! setdiff(COMBO.NICHE.CONTEXT$searchTaxon, APNI$searchTaxon)


## Add growing species as a column
SPP.GROW = as.data.frame(spp.grow)
SPP.GROW$PLANTED_GROWING = TRUE
SPP.GROW = dplyr::rename(SPP.GROW, searchTaxon = spp.grow)


#########################################################################################################################
## Now merge the data
names(APNI)
names(SPP.GROW)


## Merge APNI
COMBO.APNI = merge(COMBO.NICHE.CONTEXT, APNI, by = "searchTaxon", all.x = TRUE) 
COMBO.APNI$APNI[is.na(COMBO.APNI$APNI)] <- "FALSE"
unique(COMBO.APNI$APNI)


## Merge combo
COMBO.APNI = merge(COMBO.APNI,      SPP.GROW, by = "searchTaxon", all.x = TRUE) 
COMBO.APNI$PLANTED_GROWING[is.na(COMBO.APNI$PLANTED_GROWING)] <- "FALSE"
unique(COMBO.APNI$APNI)


## Set NA for no.of growers to blank, then sort by no. of growers.The extra species are ones I dowloaded in error
COMBO.APNI$Number.of.growers[is.na(COMBO.APNI$Number.of.growers)] <- 0
COMBO.APNI = COMBO.APNI[with(COMBO.APNI, rev(order(Number.of.growers))), ]


## Reorder
names(COMBO.APNI)
COMBO.APNI = COMBO.APNI[, c(1:18, 190, 191, 19:189)]


## Check
View(COMBO.APNI[, c(1:20)])




#########################################################################################################################
## 4). QUERY TABLE TO CREATE GROWING LISTS FOR PAUL
#########################################################################################################################


########################################
## Create lists for the growing species:
# a) native and restricted distribution; 
# b) native and broad distribution; 
# c) exotic and restricted distribution; 
# d) exotic and broad distribution. 


#########################################################################################################################
## Use the six figure summary to decide what is "broad" and "restricted". Currently using the 3rd quartile
summary(COMBO.APNI$AREA_OCCUPANCY)
summary(COMBO.APNI$LGA.AGG)
summary(COMBO.APNI$Annual_mean_temp_range)    ## I am suspicious of this...
summary(COMBO.APNI$Number.of.growers)
summary(COMBO.APNI$COMBO.count)               ## Looks weird


#########################################################################################################################
## a). Native and restricted distribution; 
NATIVE.RESTRICTED = subset(COMBO.APNI, AREA_OCCUPANCY < 200 &  
                             APNI == "TRUE")[ c("searchTaxon", "COMBO.count", "AREA_OCCUPANCY", "LGA.AGG",
                                                "Annual_mean_temp_range", "APNI", "PLANTED_GROWING", "Top_200")]

## Can check the maps for these species
dim(NATIVE.RESTRICTED)


#########################################################################################################################
## b). Native and restricted distribution; 
NATIVE.BROAD = subset(COMBO.APNI, AREA_OCCUPANCY > 4000 &  
                             APNI == "TRUE")[ c("searchTaxon", "COMBO.count", "AREA_OCCUPANCY", "LGA.AGG",
                                                "Annual_mean_temp_range", "APNI", "PLANTED_GROWING", "Top_200")]

## Can check the maps for these species
dim(NATIVE.BROAD)


#########################################################################################################################
## c). Native and restricted distribution; 
EXOTIC.RESTRICTED = subset(COMBO.APNI, AREA_OCCUPANCY < 200 &  
                             APNI == "FALSE")[ c("searchTaxon", "COMBO.count", "AREA_OCCUPANCY", "LGA.AGG",
                                                "Annual_mean_temp_range", "APNI", "PLANTED_GROWING", "Top_200")]

## Can check the maps for these species
dim(EXOTIC.RESTRICTED)


#########################################################################################################################
## d). Native and restricted distribution; 
EXOTIC.BROAD = subset(COMBO.APNI, AREA_OCCUPANCY > 4000 &  
                        APNI == "FALSE")[ c("searchTaxon", "COMBO.count", "AREA_OCCUPANCY", "LGA.AGG",
                                           "Annual_mean_temp_range", "APNI", "PLANTED_GROWING", "Top_200")]

## Can check the maps for these species
dim(EXOTIC.BROAD)





#########################################################################################################################
## 5). QUERY TABLE TO CREATE HIA LISTS FOR PAUL
#########################################################################################################################


########################################
## Create lists for the growing species:
# e) native highly produced; 
# f) native lowly produced; 
# g) exotic highly produced; 
# h) exotic lowly produced.
summary(COMBO.APNI$Number.of.growers)



#########################################################################################################################
## e). Native highly produced 
NATIVE.POPULAR = subset(COMBO.APNI, Number.of.growers > 200 &  PLANTED_GROWING == "TRUE" & 
                          APNI == "TRUE")[ c("searchTaxon", "COMBO.count", "AREA_OCCUPANCY", "LGA.AGG",
                                             "Annual_mean_temp_range", "APNI", "PLANTED_GROWING", "Top_200")]

## Can check the maps for these species
dim(NATIVE.POPULAR)


#########################################################################################################################
## f). Native lowly produced 
NATIVE.UNPOPULAR = subset(COMBO.APNI, Number.of.growers < 25 &  PLANTED_GROWING == "TRUE" &
                            APNI == "TRUE")[ c("searchTaxon", "COMBO.count", "AREA_OCCUPANCY", "LGA.AGG",
                                               "Annual_mean_temp_range", "APNI", "PLANTED_GROWING", "Top_200")]

## Can check the maps for these species
dim(NATIVE.UNPOPULAR)


#########################################################################################################################
## g). Exotic highly produced 
EXOTIC.POPULAR = subset(COMBO.APNI, Number.of.growers > 200 &  PLANTED_GROWING == "TRUE" & 
                          APNI == "TRUE")[ c("searchTaxon", "COMBO.count", "AREA_OCCUPANCY", "LGA.AGG",
                                             "Annual_mean_temp_range", "APNI", "PLANTED_GROWING", "Top_200")]

## Can check the maps for these species
dim(EXOTIC.POPULAR)


#########################################################################################################################
## h). Exotic lowly produced 
EXOTIC.UNPOPULAR = subset(COMBO.APNI, Number.of.growers < 25 &  PLANTED_GROWING == "TRUE" &
                            APNI == "FALSE")[ c("searchTaxon", "COMBO.count", "AREA_OCCUPANCY", "LGA.AGG",
                                               "Annual_mean_temp_range", "APNI", "PLANTED_GROWING", "Top_200")]

## Can check the maps for these species
dim(EXOTIC.UNPOPULAR)





#########################################################################################################################
## 6). QUERY TABLE TO CREATE LISTS OF NEW SPP FOR SALLY
#########################################################################################################################


#########################################################################################################################
## Rare species we can't model?
RARE.SPP = subset(COMBO.APNI, COMBO.count < 50)[ c("searchTaxon", "COMBO.count", "AREA_OCCUPANCY", "LGA.AGG",
                                                            "Top_200", "Origin", "Annual_mean_temp_range")]


## Reorder DF by species
RARE.SPP = RARE.SPP[with(RARE.SPP, order(searchTaxon)), ] 
dim(RARE.SPP)


#########################################################################################################################
## Potential new species: For some species, could add in traits here too
## Find infrequently sold spp, big environmental & geographic range, but could have similar traits to popular species
NEW.SPP = subset(COMBO.APNI, AREA_OCCUPANCY > 4000 & LGA.AGG > 66 &
                   Annual_mean_temp_range > 18 & 
                   Number.of.growers < 25)[ c("searchTaxon", "COMBO.count", "AREA_OCCUPANCY", "LGA.AGG",
                                              "Top_200", "Origin", "Annual_mean_temp_range")]


## Reorder DF by species
NEW.SPP = NEW.SPP[with(NEW.SPP, order(searchTaxon)), ] 
dim(NEW.SPP)




#########################################################################################################################
###################################################  TBC ################################################################ 
#########################################################################################################################