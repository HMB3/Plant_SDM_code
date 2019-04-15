#########################################################################################################################
######################################## CREATE SPECIES LISTS FOR HIA MODELLING ######################################### 
#########################################################################################################################


## This code takes the raw list of plants supplied by Matt Plumber and cleaned by Anthony Manea, and then cleans the list 
## as best as possible in R to use the species binomial as the unit for downloading data and analysing niches.
## see the publication in the science for the total environment ::



#########################################################################################################################
## 1). READ IN SPECIES LISTS
#########################################################################################################################


#########################################################################################################################
## Horticultural lists ::
## This evergreen list (CLEAN.list and HIA.lits) derives from all species and varieties sold anywhere in Australia in the last 5 years. 
## Anthony Manea cleaned up the data and cross-linked to growth form and exotic/native status and derived a list of ~1000 
## species that are the Most commonly sold, covering the right ratio of growth forms, regional representation and native/exotic


## Here is the description of the evergreen list from our publication ::


# To analyse horticulturally significant Australian species, we obtained a list of plant species currently grown by the nursery 
# industry across Australia, including a count of nurseries known to be growing each species within each state or territory 
# (M Plumber, pers comm). Duplicate records or records that contained topiary (i.e. shaping) descriptions for a species were 
# condensed into a single record for that species. The remaining species records were designated an origin (native or exotic) 
# using various state and territory government and botanic garden websites. The exotic species records were then removed from 
# the list.


## Taxa in these lists are really "things", that is varities cultivars and some other weird stuff too.
## Ultimately we juset need binomials................................................................


## HIA   = 1k  taxa with <25 growers
## RAW   = 20k uncleaned taxa  
## CLEAN = 13k taxa cleaned by Anthony Manea
HIA.list            = read.csv("./data/base/HIA_LIST/HIA/GREEN_CITIES_DRAFT_LIST_2709_2017.csv", stringsAsFactors = FALSE)
HIA.RAW             = read.csv("./data/base/HIA_LIST/HIA/HIA_ORIGINAL_RAW.csv",                  stringsAsFactors = FALSE)
CLEAN.list          = read.csv("./data/base/HIA_LIST/HIA/HIA.CLEAN.csv",                         stringsAsFactors = FALSE)
TREE.list           = read.csv("./data/base/HIA_LIST/global_tree_search_trees_1_2.csv",          stringsAsFactors = FALSE)
APNI                = readRDS("./data/base/HIA_LIST/ALA/APNI_LIST.rds")


#########################################################################################################################
## After analysing different species lists, we settled on modelling as many of the taxa in the CLEAN.list as possible.
## This is subject to the availability of occurrence data, model performance, etc.


#########################################################################################################################
## This is a list of hollow bearing trees that David Coleman and Alessandro created
## This list of species can be used for testing model setttings, etc.
HOLLOW.SPP          = read.csv("./data/base/HIA_LIST/Tree_hollows/Hollow_species.csv",           stringsAsFactors = FALSE)
HOLLOW.SPP          = HOLLOW.SPP$Hollow_species
HOLLOW.SPP          = as_utf8(HOLLOW.SPP, normalize = TRUE)


#########################################################################################################################
## This is a list of all the species with maxent models rated by Linda so far (i.e. for the Stoten publication). 
## Just provide here as a guide
MAXENT.RATING.LAT   = read.csv("./output/maxent/MAXENT_CHECK_RATING_FEB2019.csv",                stringsAsFactors = FALSE)
MXT.CHECK           = read.csv("./output/maxent/MAXNET_ORIGIN_RESULTS.csv",                      stringsAsFactors = FALSE)


#########################################################################################################################
## This is the list of 176 species analysed for the STOTEN publication
native.good.models  = subset(MXT.CHECK, Origin == "Native" & Check.map <=2)$searchTaxon



#########################################################################################################################
## READ IN ALA/GBIF PLANT DUMP
#########################################################################################################################


#########################################################################################################################
## Now read in tables of all plant records from ALA and GBIF
#ALA.PLANTS  = read.csv("./data/base/HIA_LIST/ALA_plants.csv", stringsAsFactors = FALSE)
#saveRDS(ALA.PLANTS, './data/base/HIA_LIST/ALA_plants.rds')
# redgum.test = ALA.PLANTS[ALA.PLANTS$scientificName %like% "Eucalyptus camaldulensis", ]
# ALA.RG      = ALA.RG[ALA.RG$scientificName %like% "Eucalyptus camaldulensis", ]
# dim(ALA.RG)



#GBIF.PLANTS = read.csv("./data/base/HIA_LIST/GBIF_plants.csv", stringsAsFactors = FALSE)





#########################################################################################################################
## 2). READ IN URBAN TREE INVENTORIES
#########################################################################################################################


#########################################################################################################################
## Now read in the list of species Alessandro created by generating global tree inventories
TI.XY   = readRDS("./data/base/Global_tree_inventories/Global_Inventory_CLEAN_21.03.2019.rds")
TI.XY   = dplyr::select(TI.XY, New_binomial, Lat, Long, Country, City)
TI.CITY = read.csv("./data/base/Global_tree_inventories/GTI_470cities_pop15.csv", stringsAsFactors = FALSE)


## Rename data
names(TI.XY)[names(TI.XY) == 'New_binomial'] <- 'searchTaxon'
names(TI.XY)[names(TI.XY) == 'City'] <-  'INVENTORY'
names(TI.XY)[names(TI.XY) == 'Lat']  <-  'lat'
names(TI.XY)[names(TI.XY) == 'Long'] <-  'lon'
TI.XY$SOURCE = 'INVENTORY'



## Check Ale's data :: why are some councils producing NA lat/lon?
names(TI.XY)
head(TI.XY)
summary(TI.XY)
TI.XY.NA <- TI.XY[rowSums(is.na(TI.XY)) > 0,]
length(unique(TI.XY.NA$INVENTORY))
length(unique(TI.XY.NA$searchTaxon))

TI.XY = na.omit(TI.XY)
length(unique(TI.XY$INVENTORY))
length(unique(TI.XY$searchTaxon))


#########################################################################################################################
## What are the most commonly planted urban trees in Australia, or the world?
## Create a lookup table for the species
TI.LUT = as.data.frame(table(TI.XY$searchTaxon))
names(TI.LUT) = c("searchTaxon", "FREQUENCY")
TI.LUT = TI.LUT[with(TI.LUT, rev(order(FREQUENCY))), ] 
head(TI.LUT);dim(TI.LUT)


## What are the unique binomials?
TI.LIST     = unique(TI.XY$searchTaxon)
TI.SPP      = TI.XY[!duplicated(TI.XY[,c("searchTaxon")]),][c("searchTaxon")]
TI.SPP$Plant_type = "Tree"




#########################################################################################################################
## 3). CLEAN THE EVERGREEN LIST
#########################################################################################################################


#########################################################################################################################
## The problem with the HIA taxa list is that it contains lots of weird varieties and cultivars.
## To create realised niches and SDMs, we need to convert these to binomials.


## Some taxa will always be lost in the process. However, we need to use the GBIF, TPL and ALA taxonomies to get the 
## occurrence data. So the varieties, etc, must be cleaned out. Note that searching GBIF and ALA for the binomials 
## will return lots of varities anyway.


#########################################################################################################################
## First, use sub to convert the HIA taxa names to bionmials
## Also, remove the taxa labelled "spp". We can't do anything with these
HIA.list$Binomial   <- sub('(^\\S+ \\S+).*', '\\1',   HIA.list$Species)   # \\s = white space; \\S = not white space
CLEAN.list$Binomial <- sub('(^\\S+ \\S+).*', '\\1',   CLEAN.list$Species) # \\s = white space; \\S = not white space
CLEAN.list          <- CLEAN.list[!grepl("spp.", CLEAN.list$Binomial),]


## Then remove weird characters - this would be better with a multiple grepl
CLEAN.list$Binomial = gsub(" x",     "",  CLEAN.list$Binomial, perl = TRUE)
CLEAN.list$Binomial = gsub("NA",     "",  CLEAN.list$Binomial, perl = TRUE)
CLEAN.list$Binomial = gsub("  ",     " ", CLEAN.list$Binomial, perl = TRUE)
CLEAN.list$Binomial = gsub(" $",     "",  CLEAN.list$Binomial, perl = TRUE)
CLEAN.list$Binomial = gsub("    $",  "",  CLEAN.list$Binomial, perl = TRUE)


## Then remove those rows with only 1 word, after the above cleaning
CLEAN.list = CLEAN.list[sapply(strsplit(as.character(CLEAN.list$Binomial), " "), length)>1,]


## So there about 4,250 unique binomials in the Evergreen list, if we exclude all the varities, etc.
## This list still has heaps of taxonomic anomalies, etc. But we can use this as the basis for the SDMs
length(unique(CLEAN.list$Binomial))


#########################################################################################################################
## Then sum the number of growers and states across the different varieties and cultivars
number.grow                        <- tapply(CLEAN.list$Number.of.growers, CLEAN.list$Binomial, sum, na.rm = TRUE)
number.state                       <- tapply(CLEAN.list$Number.of.States,  CLEAN.list$Binomial, sum, na.rm = TRUE)
CLEAN.list$Number.of.growers.total <- number.grow[CLEAN.list$Binomial]
CLEAN.list$Number.of.states.total  <- number.state[CLEAN.list$Binomial]


## Create a table of the contextual HIA columns, using the unique binomials
## The "searchTaxon" column is used here, because we will be using this
## to search ALA and GBIF. That variable will be used throughout the analysis
CLEAN.GROW = dplyr::select(CLEAN.list, Binomial, Plant.type, Origin, Number.of.growers.total, Number.of.states.total)
names(CLEAN.GROW) = c("Evergreen_taxon", "Plant_type", "Origin", "Total_growers", "Number_states")
CLEAN.GROW        = CLEAN.GROW [!duplicated(CLEAN.GROW [,c('Evergreen_taxon')]),]
View(CLEAN.GROW)


#########################################################################################################################
## What are the proportions of species from native/exotic and different life forms?
## 64% of this table don't have native/exotic status attributed. Is there a way to look this up?
round(with(CLEAN.GROW, table(Plant_type)/sum(table(Plant_type))*100), 1)
round(with(CLEAN.GROW, table(Origin)/sum(table(Origin))*100), 1)


## Save the Evergreen list to file to check against the GBIF taxonomy
write.csv(CLEAN.GROW, "./data/base/HIA_LIST/HIA/EVERGREEN_LIST_MARCH2019.csv", row.names = FALSE)





#########################################################################################################################
## 4). CHECK TAXONOMY FOR THE HORTICULTURAL SPECIES
#########################################################################################################################


## This species list for the Stoten publication is simply the number of native species with > 50 plantings and 1 grower, 
## with good modelsfrom the previous analyis :: 'native.good.models'
length(intersect(subset(TI.LUT, FREQUENCY > 50)$searchTaxon, CLEAN.GROW$Evergreen_taxon))


#########################################################################################################################
## We will model two list ::


## A). All trees that are planted in Australia and sold (so the intersection of Ale's inventories and CLEAN.GROW) 
##     Native AND exotic. This set can be modlled with ALA, GBIF and Tree inventory data

## B). All the non-trees on the Evergreen list (i.e. CLEAN.GROW)
##     Native AND exotic. This set can be modelled with ALA and GBIF data.



## We need to harmonise the taxonomy between the HIA and Inventory lists 

## 1). Create the inital list by combing the planted trees with evergreen list
##     TREE.HIA.SPP = intersect(subset(TI.LIST, Plantings > 50)$searchTaxon, CLEAN.GROW$searchTaxon) 
##     This is about 400 species
##     Origin
##     NA     Exotic Native 
##     20.5   25.8   53.8 

## 2). Clean this list using the GBIF backbone taxonomy :: use the "species" column in from the GBIF "species lookup" tool
##     https://www.gbif.org/tools/species-lookup

## 3). Run the GBIF "species" list through the TPL taxonomy. Take "New" Species and Genus as the "searchTaxon"
##     Make this column the "searchTaxon" in CLEAN.GROW.

## 4). Searching for native/exotic, etc, should come after the taxonomic check   


#########################################################################################################################
## Read in the list of species that have been checked against the GBIF taxonomy::
CLEAN.GBIF     = read.csv("./data/base/HIA_LIST/HIA/EVERGREEN_LIST_GBIF_RESULT.csv", stringsAsFactors = FALSE)
round(with(CLEAN.GBIF , table(matchType)/sum(table(matchType))*100), 1)
CLEAN.GBIF.SPP = unique(CLEAN.GBIF$species)
GBIF.NA        = CLEAN.GBIF[(CLEAN.GBIF$species == ""), ]$verbatimScientificName  ## Get taxa which returned blank from GBIF taxonomy


## Now join the Evergreen taxonomy to the GBIF taxonomy
CLEAN.GBIF = dplyr::select(CLEAN.GBIF, verbatimScientificName, species)
CLEAN.GROW <- CLEAN.GROW %>%
  dplyr::left_join(., CLEAN.GBIF, by = c("Evergreen_taxon" = "verbatimScientificName")) %>%
  dplyr::select(., Evergreen_taxon, species, Origin, Plant_type, Total_growers, Number_states)
names(CLEAN.GROW)[names(CLEAN.GROW) == 'species'] <- 'GBIF_taxon'


#########################################################################################################################
## Run the TPL function on this list ::
message('Running TPL taxonomy for ', length(CLEAN.GBIF.SPP), ' species in the HIA list')
# HIA.TAXO <- Taxonstand::TPL(CLEAN.GBIF.SPP, infra = TRUE,
#                   corr = TRUE, repeats = 100)  ## to stop it timing out...
# HIA.TAXO$New_binomial = paste(HIA.TAXO$New.Genus, HIA.TAXO$New.Species, sep = " ")
# sum(is.na(HIA.TAXO$New_binomial))
# sum(is.na(HIA.TAXO$New.Taxonomic.status))
# 
# 
# TPL.NA                = HIA.TAXO[(HIA.TAXO$New.Taxonomic.status == ""), ]$Taxon
# names(HIA.TAXO)
# saveRDS(HIA.TAXO, paste0(ALA_path, 'HIA_TPL_TAXO.rds'))
ALA_path = "./data/base/HIA_LIST/ALA/" 
HIA.TAXO = readRDS(paste0(ALA_path, 'HIA_TPL_TAXO.rds'))
TPL.NA   = HIA.TAXO[(HIA.TAXO$New.Taxonomic.status == ""), ]$Taxon


## None of the evergreen or HA taxa are NA
sum(is.na(CLEAN.GROW$Evergreen_taxon))
sum(is.na(HIA.TAXO$GBIF_taxon))



#########################################################################################################################
## Create a New_binomial field from the TPL check. This becomes the "SearchTaxon" column for ALA and GBIF
## There is another list of plants that don't match either GBIF or TPL taxonomy. 
## Some are varieites, some are probably Australian-only taxonomy
CLEAN.TAXO = dplyr::select(HIA.TAXO, Taxon, New_binomial, New.Taxonomic.status)


## Only 220 species should be NA by this join
CLEAN.GROW <- CLEAN.GROW %>%
  dplyr::left_join(., CLEAN.TAXO, by = c("GBIF_taxon" = "Taxon"))

CLEAN.GROW <- dplyr::select(CLEAN.GROW, Evergreen_taxon, GBIF_taxon, New_binomial, 
                            New.Taxonomic.status, Origin, Plant_type, Total_growers, Number_states)

names(CLEAN.GROW)[names(CLEAN.GROW) == 'New_binomial'] <- 'searchTaxon'
names(CLEAN.GROW)[names(CLEAN.GROW) == 'New.Taxonomic.status'] <- 'Taxo_status'
View(CLEAN.GROW)


## Check how the species names match between HIA, GBIF and TPL
CLEAN.NA = CLEAN.GROW[(CLEAN.GROW$searchTaxon == "NA NA"), ]


#########################################################################################################################
## Now create a taxonomic lookup table using the new binomial
CLEAN.CLASS <- unique(CLEAN.GROW$searchTaxon) %>%  
  lookup_table(., by_species = TRUE) 
CLEAN.CLASS  <- setDT(CLEAN.CLASS , keep.rownames = TRUE)[]
colnames(CLEAN.CLASS)[1] <- "searchTaxon"


## Use table later to join ::
GBIF.GROW = join(CLEAN.GROW, CLEAN.CLASS, type = "left")
GBIF.GROW =  dplyr::select(GBIF.GROW, searchTaxon, Evergreen_taxon, GBIF_taxon, Taxo_status, genus, family, 
                           order, group, Origin, Plant_type, Total_growers, Number_states)
View(GBIF.GROW)


#########################################################################################################################
## Save the evergreen table out with botanical columns and context
#write.csv(GBIF.GROW, "./data/base/HIA_LIST/HIA/EVERGREEN_TPL_LIST_MARCH2019.csv", row.names = FALSE)


## Try filling in the trees using Alessandro's data
GBIF.GROW$Plant_type[GBIF.GROW$Plant_type==""]  <- NA 
GBIF.GROW = FillIn(GBIF.GROW, TI.SPP, "Plant_type", "Plant_type", 
                   KeyVar = c("searchTaxon"), allow.cartesian = FALSE, KeepD2Vars = FALSE)
GBIF.GROW.NA <- GBIF.GROW[rowSums(is.na(GBIF.GROW)) > 0,]


## Save the list to file
GBIF.GROW$searchTaxon = gsub("NA NA", "NA", GBIF.GROW$searchTaxon, perl = TRUE)
GBIF.GROW.VALID       = GBIF.GROW[!(GBIF.GROW$searchTaxon == "NA"),]
dim(GBIF.GROW.VALID)


#########################################################################################################################
## All the evergreen plants lack life form attribution, and 64% lack native/exotic data 
round(with(GBIF.GROW.VALID, table(Plant_type)/sum(table(Plant_type))*100), 1)
round(with(GBIF.GROW.VALID, table(Origin)/sum(table(Origin))*100), 1)


## How much overlap is there between Alessandro's inventory data, and the HIA list?
length(subset(GBIF.GROW.VALID , Plant_type == "")$searchTaxon)              ## No    taxa without life form
length(subset(GBIF.GROW.VALID , Origin == "")$searchTaxon)                  ## 2582  taxa without origin
length(intersect(GBIF.GROW.VALID$searchTaxon, unique(TI.XY$searchTaxon)))   ## 1,111 taxa on HIA and Inventory list
TI.HIA = intersect(GBIF.GROW.VALID$searchTaxon, unique(TI.XY$searchTaxon))
summary(GBIF.GROW.VALID)





#########################################################################################################################
## 5). CREATE FINAL SPECIES LISTS FOR ANALYSIS
#########################################################################################################################


#########################################################################################################################
## These are the lists of trees and non-trees to model
WPW.spp         = sort(GBIF.GROW.VALID$searchTaxon)
WPW.tree        = subset(GBIF.GROW.VALID, Plant_type == "Tree")$searchTaxon
WPW.non.tree    = as.character(GBIF.GROW.VALID$searchTaxon[!GBIF.GROW.VALID$searchTaxon %in% WPW.tree ])
TI.spp          = unique(TI.XY$searchTaxon)
TI.HIA          = setdiff(TI.spp, WPW.spp)

WPW.spp         = WPW.spp[lapply(WPW.spp,length)>0]
WPW.NA          = unique(c(TPL.NA, GBIF.NA))          ## These taxa did not match either GBIF or ALA


## This is a list of species we can use to test the models. These spp had poor models without MESS maps
hollow.test.spp = c("Corymbia citriodora",      "Eucalyptus baileyana",  "Eucalyptus baxteri",  "Eucalyptus decorticans", "Eucalyptus dumosa",
                    "Eucalyptus major",         "Eucalyptus megacarpa",  "Eucalyptus moluccana", "Eucalyptus ovata",      "Eucalyptus patens",
                    "Eucalyptus salmonophloia", "Eucalyptus tenuipes")


popular.test    = c("Acer palmatum", "Syzygium smithii", "Magnolia grandiflora", "Callistemon viminalis", 
                    "Gardenia jasminoides", "Pyrus calleryana", "Murraya paniculata", "Ficus microcarpa",
                    "Jacaranda mimosifolia", "Hardenbergia violacea")


popular.test    = intersect(popular.test, TI.spp)
inv.test        = as.character(head(TI.LUT$searchTaxon, 20))
inv.test        = unique(c(inv.test, popular.test))





#########################################################################################################################
## Now we need a list of all the trees, and non trees. Native info not as important.
write.csv(GBIF.GROW.VALID, "./data/base/HIA_LIST/HIA/WPW_MODELLING_LIST_MARCH2019.csv", row.names = FALSE)


#########################################################################################################################
## Check the proportions of species with different maxdent model ratings from the previous run
## 1 = good, 2 = fair, 3 = poor. Current success rate is ~50-60%
table(MAXENT.RATING.LAT$CHECK_MAP)
round(with(MAXENT.RATING.LAT, table(CHECK_MAP)/sum(table(CHECK_MAP))*100), 1)





#########################################################################################################################
## 6). COMBINE THE NICHE RUNS TOGETHER
#########################################################################################################################


# ## Create a list of all dataframes with the extension from this run
# COMBO.NICHE.list  = list.files(DATA_path, pattern = 'COMBO_NICHE_CONTEXT_EVERGREEN',  full.names = TRUE, recursive = TRUE)
SDM.TABLE.list = list.files(DATA_path, pattern = 'SDM_SPAT_OCC_BG_', full.names = TRUE, recursive = TRUE)
  

# #########################################################################################################################
# ## Now combine the niche tables for each species into one table
# COMBO.NICHE.ALL <- COMBO.NICHE.list[1:8] %>%
# 
#   ## pipe the list into lapply
#   lapply(function(x) {
# 
#     ## create the character string
#     f <- paste0(x)
# 
#     ## load each .csv file
#     d <- readRDS(f)
#     d
# 
#   }) %>%
# 
#   ## finally, bind all the rows together
#   bind_rows
# 
# 
# ## Update this
# str(COMBO.NICHE.ALL)
# dim(COMBO.NICHE.ALL)
# 
# 
# ## Make sure the Species are unique
# COMBO.NICHE.ALL = COMBO.NICHE.ALL[!duplicated(COMBO.NICHE.ALL[,c('searchTaxon')]),]
# dim(COMBO.NICHE.ALL)
# length(unique(COMBO.NICHE.ALL$searchTaxon))
# 
# 
# #########################################################################################################################
# ## Now combine the raster tables for each species into one table
# SDM.TABLE.ALL <- SDM.TABLE.list %>%
# 
#   ## pipe the list into lapply
#   lapply(function(x) {
# 
#     ## create the character string
#     f <- paste0(x)
# 
#     ## load each .csv file
#     d <- readRDS(f)
#     d
# 
#   }) %>%
# 
#   ## finally, bind all the rows together
#   bind_rows
# 
# 
# ## This is a summary of maxent output for current conditions
# dim(SDM.TABLE.ALL)
# names(SDM.TABLE.ALL)[1:10]
# 
# 
# ##
# length(unique(SDM.TABLE.ALL$searchTaxon))
# 
# 
# 
# #########################################################################################################################
# ## Save the niche and raster data
# saveRDS(COMBO.NICHE.ALL,  paste0(DATA_path, 'COMBO_NICHE_ALL_',  save_run, '.rds'))
# saveRDS(SDM.TABLE.ALL, paste0(DATA_path, 'SDM_SPAT_OCC_BG_', save_run, '.rds'))


#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
## Once the modelling is finished, these are the columns we can use to rate each species.
## The contextual data should come from the CLEAN taxa list above
results.columns = c("searchTaxon",       ## From the ALA/ GBIF download code 
                    "Origin",            ## native/extoic : from Anthony Manea's spreadsheet, affected by taxonomy....
                    "Plant_type",        ## From Anthony Manea's spreadsheet, will be affected by taxonomy....
                    "Total_growers",     ## From Anthony Manea's spreadsheet.....    
                    "Number_states",     ## From Anthony Manea's spreadsheet.....
                    
                    "Plantings",         ## No. urban plantings :: from urban tree inventory
                    "GLOBAL_RECORDS",    ## No. global records  :: from the R workflow
                    "AUS_RECORDS",       ## No. AUS records     :: from the R workflow
                    "Maxent_records",    ## No. records used in the SDM
                    "SUA_COUNT",         ## No. SUAs each species occurs in :: From the R workflow
                    
                    "Number_var",        ## No. maxent variables :: from Maxent code
                    "Var_pcont",         ## Maxent Variable with highest permutation importance    
                    "Per_cont",          ## The permutaiton importance of that variable
                    "Var_pimp",          ## Maxent Variable with highest permutation importance    
                    "Perm_imp",          ## The permutaiton importance of that variable 
                    "Iterations",               ## No. iterations                                                                    
                    "Training_AUC",             ## training AUC
                    "max_tss",                  ## Maximium True skill statistic
                    "Number_background_points", ## No. background points
                    "Logistic_threshold"        ## Maxent threshold)
)




#########################################################################################################################
## OUTSTANDING LIST TASKS:
#########################################################################################################################


## Find a way to code the infilling of data on origing, life form, etc..................................................






#########################################################################################################################
############################################  END OF HIA LIST CODE ###################################################### 
#########################################################################################################################