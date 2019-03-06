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


## HIA   = 1k taxa with <25 growers
## RAW   = X uncleaned taxa  
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
## 2). READ IN URBAN TREE INVENTORIES
#########################################################################################################################


#########################################################################################################################
## Now read in the list of species Alessandro created by generating global tree inventories
TI.LIST             = read.csv("./data/base/HIA_LIST/COMBO/ALE_TREE_SPP_LIST.csv",               stringsAsFactors = FALSE)
TI.XY               = read.csv("./data/base/HIA_LIST/COMBO/ALE_TREE_SPP_XY.csv",                 stringsAsFactors = FALSE)


## Check Ale's data
TI.LIST = na.omit(TI.LIST)
TI.LIST = TI.LIST[with(TI.LIST, rev(order(Plantings))), ]
head(TI.LIST$searchTaxon, 10)


## Rename tree inventory data
TI.XY = TI.XY[c("searchTaxon", "POINT_X", "POINT_Y", "FILENAME")]
names(TI.XY)[names(TI.XY) == 'FILENAME'] <- 'INVENTORY'
names(TI.XY)[names(TI.XY) == 'POINT_X']  <- 'lon'
names(TI.XY)[names(TI.XY) == 'POINT_Y']  <- 'lat'


## Remove the gunk
TI.XY$INVENTORY = gsub("_trees.shp",      "", TI.XY$INVENTORY)
TI.XY$INVENTORY = gsub("_trees_.shp",     "", TI.XY$INVENTORY)
TI.XY$INVENTORY = gsub("_trees_sign.shp", "", TI.XY$INVENTORY)
unique(TI.XY$INVENTORY)
TI.XY$SOURCE = 'INVENTORY'


## Filter the XY tree inventory data to just the cleaned species
TI.XY  = TI.XY[TI.XY$searchTaxon %in% unique(TI.LIST$searchTaxon), ]
length(unique(TI.XY$searchTaxon))
TI.XY = na.omit(TI.XY)


## Replace the weird shapefile name for hobart
TI.XY$INVENTORY = gsub("_trees_sign_.shp", "", TI.XY$INVENTORY)
head(TI.XY)
unique(TI.XY$INVENTORY)


#########################################################################################################################
## What are the most commonly planted urban trees in Australia, or the world?
## Create a lookup table for the species
TI.LUT = as.data.frame(table(TI.XY$searchTaxon))
names(TI.LUT) = c("searchTaxon", "FREQUENCY")
TI.LUT = TI.LUT[with(TI.LUT, rev(order(FREQUENCY))), ] 
head(TI.LUT);dim(TI.LUT)
length(intersect(CLEAN.SPP$Binomial, TI.LUT $searchTaxon))





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
CLEAN.GROW = select(CLEAN.list, Binomial, Plant.type, Origin, Number.of.growers.total, Number.of.states.total)
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
## 6). CHECK TAXONOMY FOR THE HORTICULTURAL SPECIES
#########################################################################################################################


## This species list for the Stoten publication is simply the number of native species with > 50 plantings and 1 grower, 
## with good modelsfrom the previous analyis :: 'native.good.models'
length(intersect(subset(TI.LIST, Plantings > 50)$searchTaxon, CLEAN.GROW$Evergreen_taxon))


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
CLEAN.GBIF = select(CLEAN.GBIF, verbatimScientificName, species)
CLEAN.GROW <- CLEAN.GROW %>%
  left_join(., CLEAN.GBIF, by = c("Evergreen_taxon" = "verbatimScientificName")) %>%
  select(., Evergreen_taxon, species, Origin, Plant_type, Total_growers, Number_states)
names(CLEAN.GROW)[names(CLEAN.GROW) == 'species'] <- 'GBIF_taxon'


## None of the evergreen or HA taxa are NA
sum(is.na(CLEAN.GROW$Evergreen_taxon))
sum(is.na(HIA.TAXO$GBIF_taxon))


#########################################################################################################################
## Run the TPL function on this list ::
message('Running TPL taxonomy for ', length(CLEAN.GBIF.SPP), ' species in the HIA list')
# HIA.TAXO <- TPL(CLEAN.GBIF.SPP, infra = TRUE,
#                   corr = TRUE, repeats = 100)  ## to stop it timing out...
# HIA.TAXO$New_binomial = paste(HIA.TAXO$New.Genus, HIA.TAXO$New.Species, sep = " ")
# sum(is.na(HIA.TAXO$New_binomial))
# sum(is.na(HIA.TAXO$New.Taxonomic.status))
# 
# 
# TPL.NA                = HIA.TAXO[(HIA.TAXO$New.Taxonomic.status == ""), ]$Taxon
# names(HIA.TAXO)
# saveRDS(HIA.TAXO, paste0(ALA_path, 'HIA_TPL_TAXO.rds'))
HIA.TAXO = readRDS(paste0(ALA_path, 'HIA_TPL_TAXO.rds'))
TPL.NA                = HIA.TAXO[(HIA.TAXO$New.Taxonomic.status == ""), ]$Taxon


#########################################################################################################################
## Create a New_binomial field from the TPL check. This becomes the "SearchTaxon" column for ALA and GBIF
## There is another list of plants that don't match either GBIF or TPL taxonomy. 
## Some are varieites, some are probably Australian-only taxonomy
CLEAN.TAXO = select(HIA.TAXO, Taxon, New_binomial, New.Taxonomic.status)


## Only 220 species should be NA by this join
CLEAN.GROW <- CLEAN.GROW %>%
  left_join(., CLEAN.TAXO, by = c("GBIF_taxon" = "Taxon")) %>%
  select(., Evergreen_taxon, GBIF_taxon, New_binomial, New.Taxonomic.status, Origin, Plant_type, Total_growers, Number_states)
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
GBIF.GROW = join(CLEAN.GROW, CLEAN.CLASS) %>%
  select(., searchTaxon, Evergreen_taxon, GBIF_taxon, Taxo_status, genus, family, 
         order, group, Origin, Plant_type, Total_growers, Number_states)
View(GBIF.GROW)


#########################################################################################################################
## Save the evergreen table out with botanical columns and context
write.csv(GBIF.GROW, "./data/base/HIA_LIST/HIA/EVERGREEN_TPL_LIST_MARCH2019.csv", row.names = FALSE)


#########################################################################################################################
## 14% of all evergreen plants lack life form attribution, and 64% lack native/exotic data 
round(with(GBIF.GROW, table(Plant_type)/sum(table(Plant_type))*100), 1)
round(with(GBIF.GROW, table(Origin)/sum(table(Origin))*100), 1)


## Not much overlap between the global tree list and the species without life form
no.type   = subset(GBIF.GROW, Plant_type == "")$searchTaxon  ## 600 taxa without life form
no.origin = subset(GBIF.GROW, Origin == "")$searchTaxon      ## 2740 taxa without origin
length(intersect(no.type, TREE.list$Tree_taxa))


#########################################################################################################################
## The result of this manual search is listed below
WPW.spp             = sort(GBIF.GROW$searchTaxon)
WPW.tree            = subset(GBIF.GROW, Plant_type == "Tree")$searchTaxon
WPW.non.tree        = subset(GBIF.GROW, Plant_type != "Tree")$searchTaxon

WPW.spp             = WPW.spp[lapply(WPW.spp,length)>0]
WPW.NA              = unique(c(TPL.NA, GBIF.NA))                  ## These taxa did not match either GBIF or ALA


#########################################################################################################################
## Now we need a list of all the trees, and non trees. Native info not as important.L 


## Here is a sample of 10 species on the "Hollow_spp" list that have bad models.
## These can be used as a test bed for different analysis settings
hollow.test.spp = c("Corymbia citriodora",      "Eucalyptus baileyana",  "Eucalyptus baxteri",  "Eucalyptus decorticans", "Eucalyptus dumosa",
                    "Eucalyptus major",         "Eucalyptus megacarpa",  "Eucalyptus moluccana", "Eucalyptus ovata",      "Eucalyptus patens",
                    "Eucalyptus salmonophloia", "Eucalyptus tenuipes")


#########################################################################################################################
## Check the proportions of species with different maxdent model ratings from the previous run
## 1 = good, 2 = fair, 3 = poor.
## Current success rater is ~50-60%
table(MAXENT.RATING.LAT$CHECK_MAP)
round(with(MAXENT.RATING.LAT, table(CHECK_MAP)/sum(table(CHECK_MAP))*100), 1)





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
                    "Number_background_points", ## NO. background points
                    "Logistic_threshold"        ## Maxent threshold)
)




#########################################################################################################################
## OUTSTANDING LIST TASKS:
#########################################################################################################################


## Find a way to code the infilling of data on origing, life form, etc..................................................






#########################################################################################################################
############################################  END OF HIA LIST CODE ###################################################### 
#########################################################################################################################