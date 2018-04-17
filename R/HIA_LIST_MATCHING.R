#########################################################################################################################
########################################  CREATE SPECIES LISTS ########################################################## 
#########################################################################################################################


## This code takes the raw list of plants with 25 or more growers supplied by Matt Plumber and Anthony Manea, and
## then cleans the list as best as possible in R to use the species binomial as the unit for downloading and analysis.

## All the cleaining methods will throw up some anomalies, which need to be tracked, and checked with the team for how
## each case is treated (see outstanding tasks at the bottom)


#########################################################################################################################
## Setup for project 
library(gtools)
library(GISTools)
library(devtools)
library(Rcpp)
library(raster)
library(rgdal)
library(plyr)
library(dplyr)
library(sfsmisc)
#library(spatstat)
library(data.table)
library(rvest)
#library(vegan)

library(SDMTools)
library(rmaxent)
library(dismo)
#library(BiodiversityR)
library(biomod2)
#library(AdaptR)
library(red)
library(ConR)

library(ff)
library(rgeos)
library(sp)
library(raster)
library(rJava)
library(things)

library(ALA4R)
library(rgbif)
library(scrubr)
library(RCurl)
library(httr)

library(taxonlookup)
library(Taxonstand)
#library(speciesgeocodeR)
library(CoordinateCleaner)
library(spatialEco)
library(raster)
library(rnaturalearth)
library(gdalUtils)

#library(knitr)
#library(htmltools)
library(yaml)
library(caTools)
library(bitops)
library(rmarkdown)
library(gsubfn)
library(functional)
library(splitstackshape)

library(tidyverse)
library(stringr)
library(maptools)
#library(ggmap)
library(rgeos)
library(magrittr)
library(datastorr)
#library(baad.data)
#library(rapportools)

library(Cairo)
library(lattice)
library(latticeExtra)
library(PerformanceAnalytics)
library(timetk)


##
p <- c('ff',    'things',         'raster',    'dismo',        'sp',           'latticeExtra', 'data.table', 
       'rgdal', 'rgeos',          'gdalUtils', 'rmaxent',      'readr',        'dplyr',        'tidyr',
       'readr', 'rnaturalearth',  'rasterVis', 'RColorBrewer', 'latticeExtra', 'parallel')


## Require packages
sapply(p, require, character.only = TRUE)


## Source functions
source('./R/HIA_CLEAN_MATCHING.R')
source('./R/GREEN_CITIES_FUNCTIONS.R')
source('./R/MAXENT_FUNCTIONS.R')
source('./R/MAPPING_FUNCTIONS.R')





#########################################################################################################################
## 1). READ IN DRAFT HIA LIST AND CLEAN
#########################################################################################################################


#########################################################################################################################
## Read in the niche data
COMBO.NICHE.CONTEXT = readRDS("./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT_APRIL_2018.rds")
COMBO.NICHE.OLD     = readRDS("./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT_1304_2018.rds")
str(unique(COMBO.NICHE.CONTEXT$searchTaxon))  ## long enough


#########################################################################################################################
## Horticultural lists ::
## This evergreen list (HIA.list) derives from all species and varieties sold anywhere in Australia in the last 5 years. Anthony Manea cleaned 
## up the data and cross-linked to growth form and exotic/native status and derived a list of ~1000 species that are the 
## Most commonly sold, covering the right ratio of growth forms, regional representation and native/exotic
HIA.list            = read.csv("./data/base/HIA_LIST/HIA/GREEN_CITIES_DRAFT_LIST_2709_2017.csv", stringsAsFactors = FALSE)
HIA.RAW             = read.csv("./data/base/HIA_LIST/HIA/HIA_ORIGINAL_RAW.csv",                  stringsAsFactors = FALSE)
CLEAN.list          = read.csv("./data/base/HIA_LIST/HIA/HIA.CLEAN.csv",                         stringsAsFactors = FALSE)
GROWING             = read.csv("./data/base/HIA_LIST/HIA/database_aus_sp_growing.csv",           stringsAsFactors = FALSE)
MISSING             = read.csv("./data/base/HIA_LIST/HIA/MISSING_SPECIES.csv",                   stringsAsFactors = FALSE)
MOD_2               = read.csv("./data/base/HIA_LIST/HIA/MOD2_LIST.csv",                         stringsAsFactors = FALSE)
COMBO.ALL           = read.csv("./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT_ALL_SPP.csv",     stringsAsFactors = FALSE)
KOP.TEST            = read.csv("./data/base/HIA_LIST/COMBO/KOPPEN_TEST_SPP.csv",                 stringsAsFactors = FALSE)
RISK.LIST           = read.csv("./data/base/HIA_LIST/HIA/RISK_LIST.csv",                         stringsAsFactors = FALSE)
RISK.BINOMIAL.CLEAN = read.csv("./data/base/HIA_LIST/HIA/RISK_BINOMIAL_DF.csv",                  stringsAsFactors = FALSE)


#########################################################################################################################
## Experimental lists :: the next round of experiments will be... 
top.200              = read.csv("./data/base/HIA_LIST/HIA/HIA_TOP_200_1309_2017.csv",       stringsAsFactors = FALSE)
renee.full           = read.csv("./data/base/HIA_LIST/HIA/RENEE_FULL_LIST.csv",             stringsAsFactors = FALSE) 
renee.taxa           = read.csv("./data/base/HIA_LIST/HIA/RENEE_TAXA.csv",                  stringsAsFactors = FALSE)
renee.50             = read.csv("./data/base/HIA_LIST/HIA/RENEE_TOP_50.csv",                stringsAsFactors = FALSE)
MQ.glasshouse        = read.csv("./data/base/HIA_LIST/HIA/MQ_glasshouse.csv",               stringsAsFactors = FALSE)
Manuel.experimental  = read.csv("./data/base/HIA_LIST/HIA/Manuel_experimental_species.csv", stringsAsFactors = FALSE)


#########################################################################################################################
## Modelling lists :: which species have bias, and which have enough records in Australia? 
SPP.BIAS             = read.csv("./data/base/HIA_LIST/COMBO/SPP_BIAS_LIST.csv",             stringsAsFactors = FALSE)
SPP.RANDOM           = subset(SPP.BIAS, AUS_BOUND_BIAS == "TRUE")
SPP.TARGET           = subset(SPP.BIAS, AUS_BOUND_BIAS == "FALSE")


## Get a list of species with > 20 records
summary(COMBO.NICHE.CONTEXT$AUS_RECORDS)
SPP.AUS = subset(COMBO.NICHE.CONTEXT, AUS_RECORDS >= 20)
summary(SPP.AUS$COMBO.count);summary(SPP.AUS$AUS_RECORDS)


##
spp.rand             = intersect(SPP.RANDOM$searchTaxon, SPP.AUS$searchTaxon)
spp.target           = intersect(SPP.TARGET$searchTaxon, SPP.AUS$searchTaxon)
str(spp.rand);str(spp.target)


## Test these values, seems very convenient that all these species have > 200 Australian records?
RAND = COMBO.NICHE.CONTEXT[COMBO.NICHE.CONTEXT$searchTaxon %in% spp.rand, ]
TARG = COMBO.NICHE.CONTEXT[COMBO.NICHE.CONTEXT$searchTaxon %in% spp.target, ]


## How do the distributions compare?
summary(COMBO.NICHE.CONTEXT$AUS_RECORDS)
summary(TARG$AUS_RECORDS)
summary(RAND$AUS_RECORDS)


##
head(renee.taxa)
head(MQ.glasshouse)


## have a look
dim(HIA.list)
dim(CLEAN.list)
dim(renee.taxa)
str(HIA.list)
head(HIA.list)


#########################################################################################################################
## intersect(CLEAN.list$Species, GROWING$scientific_name)



#########################################################################################################################
## Create a list of the raw HIA list, but removing the weird characters...
## Just use "TRIM(CELL)" in excel
## trim <- function (x) gsub("^\\s+|\\s+$", "", x) ##  HIA.list$Species <- trim(HIA.list$Species)
RAW.HIA.SPP = gsub("  ",     " ", HIA.list$Species)
RAW.HIA.SPP = gsub(" $",     "",  HIA.list$Species, perl = TRUE)
RAW.HIA.SPP = gsub("    $",  "",  HIA.list$Species, perl = TRUE)
length(RAW.HIA.SPP)


## also, add the "Top 200" species in here
spp.200          = top.200[c("Species")]
spp.200$Species  <- sub('(^\\S+ \\S+).*', '\\1', spp.200$Species) # \\s = white space; \\S = not white space

spp.200$Species  = gsub("  ",     " ", spp.200$Species)
spp.200$Species  = gsub(" $",     "",  spp.200$Species, perl = TRUE)
spp.200$Species  = gsub("    $",  "",  spp.200$Species, perl = TRUE)
spp.200$Top_200  = "TRUE"
spp.200          = dplyr::rename(spp.200, Binomial = Species)


## Just get the species renee selected that are not on the top 1000 or 200
renee.list       = renee.taxa[c("Species", "Growth_Form")]



#########################################################################################################################
## Merge the ~1000 with the top 200
## This merge won't get the ones that match to binomial
HIA.list$Binomial <- sub('(^\\S+ \\S+).*', '\\1', HIA.list$Species) # \\s = white space; \\S = not white space


## Check this reduces the number of Top 200 missing from this list
HIA.list = merge(HIA.list, spp.200, by = "Binomial", all.x = TRUE) 
HIA.list$Top_200[is.na(HIA.list$Top_200)] <- "FALSE"
HIA.list$Origin <- gsub(" ",  "", HIA.list$Origin)


## check
str(HIA.list)
head(HIA.list)
unique(HIA.list$Top_200)
unique(HIA.list$Origin)
length(unique(HIA.list$Binomial)) ## 660 unique binomials
names(HIA.list)


#########################################################################################################################
## WHAT ARE SOME USEFUL MEASURES OF THE DATASET?
#########################################################################################################################


## ~ HALF the total species are native, ~47% of the top 200
dim(subset(HIA.list, Origin == "Native"))[1]/dim(HIA.list)[1]*100
dim(subset(top.200,  Origin == "Native"))[1]/dim(top.200)[1]*100





#########################################################################################################################
## 2). DRAFT LIST PREP: THIS NEEDS TO CHANGE IN CONSULTATION WITH RACH, PAUL, ETC
#########################################################################################################################


## First, taxa with "spp". Return to these later...
DRAFT.HIA.TAXA         = HIA.list
DRAFT.HIA.TAXA         = DRAFT.HIA.TAXA[DRAFT.HIA.TAXA$Species %in% DRAFT.HIA.TAXA$Species[!grepl("spp.", DRAFT.HIA.TAXA$Species)], ]
dim(DRAFT.HIA.TAXA)    ## 948 species after we cut out the "spp."


## Remove weird characters...
DRAFT.HIA.TAXA$Species = gsub(" x",     "",  DRAFT.HIA.TAXA$Species)
DRAFT.HIA.TAXA$Species = gsub("NA",     "",  DRAFT.HIA.TAXA$Species)
DRAFT.HIA.TAXA$Species = gsub("  ",     " ", DRAFT.HIA.TAXA$Species)
DRAFT.HIA.TAXA$Species = gsub(" $",     "",  DRAFT.HIA.TAXA$Species, perl = TRUE)
DRAFT.HIA.TAXA$Species = gsub("    $",  "",  DRAFT.HIA.TAXA$Species, perl = TRUE)


#########################################################################################################################
## Now create a table of how many varieties each species has
# length(unique(HIA.list$Binomial)) 
# length(unique(sub('(^\\S+ \\S+).*', '\\1', DRAFT.HIA.TAXA$Species)))


## Create another binomial column
DRAFT.HIA.TAXA$Binomial <- sub('(^\\S+ \\S+).*', '\\1', DRAFT.HIA.TAXA$Species) # \\s = white space; \\S = not white space


## And count how many varieties each taxa has?
## Before, this needs to change so that every variety that matches a binomial (e.g. magnolia grandiflora) is added to the
## Number of varieties. Also, can we take the variety with the highest number of growers? There are 8 different varieties
## of magnolia, currently I'm not getting the most popular ones.
HIA.VARIETY <- 
  DRAFT.HIA.TAXA$Binomial[DRAFT.HIA.TAXA$Binomial != DRAFT.HIA.TAXA$Species] %>% 
  table %>% 
  as.data.frame %>% 
  setNames(c('Binomial', 'No.of.Varieties')) %>% 
  full_join(DRAFT.HIA.TAXA) %>% 
  dplyr::select(Species, Binomial, No.of.Varieties, Plant.type:WA, Number.of.growers, Number.of.States, Origin, Top_200)

HIA.VARIETY %>% 
  filter(Binomial==Species)


#######################################################################################################################
## Which taxa have the second word capitalised, and are not captured by eliminating the "spp"?
grep('^\\S+ [A-Z]', HIA.VARIETY$Species, val = TRUE)


## Now for GBIF, just get the unique species...
HIA.SPP = HIA.VARIETY[!duplicated(HIA.VARIETY["Binomial"]),]
HIA.SPP = dplyr::rename(HIA.SPP, HIA.Taxa = Species)


## Reorder by species
HIA.SPP = HIA.SPP[with(HIA.SPP, order(Binomial)), ] 
#View(HIA.SPP)


#######################################################################################################################
## Now create list of HIA species. 610 unique species, mius the corrections, etc. 
spp            = unique(as.character(HIA.SPP$Binomial))
spp.renee      = unique(as.character(renee.list$Species)) ## 
length(spp)   


########################################################################################################################
## Try using taxonlookup to check the taxonomy
HIA.SPP.LOOKUP = lookup_table(HIA.SPP[["Binomial"]], by_species = TRUE) ## convert rows to column and merge
HIA.SPP.LOOKUP = setDT(HIA.SPP.LOOKUP , keep.rownames = TRUE)[]
HIA.SPP.LOOKUP = dplyr::rename(HIA.SPP.LOOKUP, Binomial = rn)
head(HIA.SPP.LOOKUP) ## Can merge on the bilogical data here...
## write.csv(HIA.SPP.LOOKUP, "./data/base/HIA_LIST/HIA/HIA_BINOMIAL_LOOKUP.csv", row.names = FALSE)





#########################################################################################################################
## 3). CREATE TEST SPECIES LIST FOR MODEL RUNS 
#########################################################################################################################


## All species on the growers list
spp.all  <- unique(HIA.SPP$Binomial)
str(spp.all)
all.reverse = sort(spp.all, decreasing = TRUE)


##
setdiff(MQ.glasshouse$Species, renee.full$Species)
setdiff(Manuel.experimental$Species, renee.full$Species)
setdiff(Manuel.experimental$Species, spp.all)


## The trial species
glasshouse.spp = trimws(sort(unique(c(renee.full$Species, MQ.glasshouse$Species))))
glasshouse.spp = as.data.frame(glasshouse.spp)

                        
test.spp = trimws(sort(unique(c(renee.full$Species, 
                                Manuel.experimental$Species,
                                MQ.glasshouse$Species,
                                "Betula pendula", "Fraxinus excelsior", "Quercus robur", "Fagus sylvatica"))))

spp.all = unique(sort(c(test.spp, spp.all)))


## Combine with 45 from the main list
HIA.SAMPLE = head(COMBO.NICHE.CONTEXT, 53)[, c("searchTaxon")]
test.spp   = sort(unique(c(test.spp, HIA.SAMPLE)))
str(test.spp)
test.reverse = sort(test.spp, decreasing = TRUE)


##
kop.spp = sort(unique(KOP.TEST$searchTaxon))

##
exp.spp  = c('Swainsona formosa', 'Templetonia retusa', 'Dodonaea baueri', 'Platanus hispanica', 'Kennedia beckxiana')  
exp.rev  = sort(exp.spp, decreasing = TRUE)
miss.spp = c('Corymbia tessellaris', 'Metrosideros excelsa')


## Create lists for the mapping code
all_spp      = gsub(" ", "_", spp.all)
all_reverse  = sort(all_spp, decreasing = TRUE)

test_spp     = gsub(" ", "_", test.spp)
test_reverse = sort(test_spp, decreasing = TRUE)

exp_spp      = gsub(" ", "_", exp.spp)
miss_spp     = gsub(" ", "_", miss.spp)
kop_spp      = gsub(" ", "_", kop.spp)


## Which species are only on the test list?
setdiff(test.spp, spp.all)
save(test.spp, file = paste("./data/base/HIA_LIST/GBIF/test.spp.RData", sep = ""))





#########################################################################################################################
## 4). MATCH EVERGREEN LIST WITH THE PLANT RISK LIST
#########################################################################################################################


## Check if there is difference between the Evergreen Connect list and the Plant risk list. Michelle ::

## Just a cross-check of two lists of most widely sold species, from different sources. There may be some interesting 
## spp on the Nursery Industry list from Anthony Kachenko that aren't on the original list we got from Evergreen, and 
## it's good to keep a look out for these things and make sure we have as good coverage as possible.


## Check for the raw match 
length(intersect(RAW.HIA.SPP, RISK.LIST$Plant.Name))                       ## 232 binomials match
length(intersect(HIA.list$Binomial,  RISK.LIST$Plant.Name))                ## 244 binomials match
length(intersect(COMBO.NICHE.CONTEXT$searchTaxon, RISK.LIST$Plant.Name))   ## 250 binomials match


## Check for the raw difference
length(setdiff(RISK.LIST$Plant.Name, RAW.HIA.SPP))                       ## 713 binomials differ
length(setdiff(RISK.LIST$Plant.Name, HIA.list$Binomial))                 ## 700 binomials differ
length(setdiff(RISK.LIST$Plant.Name, COMBO.NICHE.CONTEXT$searchTaxon))   ## 695 binomials differ


#########################################################################################################################
## Now create a new species list to manipulate
RISK.CLEAN            = RISK.LIST
RISK.CLEAN$Plant.Name = trimws(RISK.CLEAN$Plant.Name, which = c("both"))   ## remove any trailing or leading white space


## Use gsub to find and replace the weirdos :: could use multiple commands, but still messy
RISK.CLEAN$Plant.Name = gsub(" x ",   " ",  RISK.CLEAN$Plant.Name,  perl = TRUE)
RISK.CLEAN$Plant.Name = gsub(" spp ", "",   RISK.CLEAN$Plant.Name,  perl = TRUE)
RISK.CLEAN$Plant.Name = gsub(" spp. ", "",  RISK.CLEAN$Plant.Name,  perl = TRUE)
RISK.CLEAN$Plant.Name = gsub(" sp. ", "",   RISK.LIST$Plant.Name,   perl = TRUE)
RISK.CLEAN$Plant.Name = gsub(" var. ", "",  RISK.LIST$Plant.Name,   perl = TRUE)
RISK.CLEAN$Plant.Name = gsub(" var ", "",   RISK.LIST$Plant.Name,   perl = TRUE)


## Too many weird things to remove
RISK.CLEAN$Plant.Name = gsub(" species and cultivars ", "",  RISK.LIST$Plant.Name,   perl = TRUE)
RISK.CLEAN$Plant.Name = gsub(" hybrids and cultivars ", "",  RISK.LIST$Plant.Name,   perl = TRUE)
RISK.CLEAN$Plant.Name = gsub(" and cultivars ", "",  RISK.LIST$Plant.Name,   perl = TRUE)
RISK.CLEAN$Plant.Name = gsub(" (cultivars) ", "",  RISK.LIST$Plant.Name,   perl = TRUE)


## Again this could be done with some kind of regular expression
RISK.CLEAN$Plant.Name = gsub(" (hybrids) ", "",  RISK.LIST$Plant.Name,   perl = TRUE)
RISK.CLEAN$Plant.Name = gsub(" (Hybrids) ", "",  RISK.LIST$Plant.Name,   perl = TRUE)
RISK.CLEAN$Plant.Name = gsub(" hybrids ", "",  RISK.LIST$Plant.Name,   perl = TRUE)
RISK.CLEAN$Plant.Name = gsub(" Hybrids ", "",  RISK.LIST$Plant.Name,   perl = TRUE)


## Remove words in parenthesis from string
RISK.CLEAN$Plant.Name = trimws(RISK.CLEAN$Plant.Name, which = c("both"))
RISK.CLEAN$Plant.Name = gsub(" x ",   " ", RISK.CLEAN$Plant.Name,  perl = TRUE)
RISK.CLEAN$Plant.Name = gsub("\\s*\\([^\\)]+\\)", "", as.character(RISK.CLEAN$Plant.Name))


## Get the binomials :: causes problems with X, but these are common anyway?
RISK.CLEAN$Binomial <- sub('(^\\S+ \\S+).*', '\\1', RISK.CLEAN$Plant.Name)


## Now get the unique binomials
RISK.BINOMIAL       = unique(RISK.CLEAN$Binomial)
RISK.BINOMIAL.DF    = as.data.frame(RISK.BINOMIAL) 
colnames(RISK.BINOMIAL.DF)[1] = "Plant_name"
names(RISK.BINOMIAL.DF)
#write.csv(RISK.BINOMIAL.DF, "./data/base/HIA_LIST/HIA/RISK_BINOMIAL_DF.csv", row.names = FALSE) 
#RISK.BINOMIAL.CLEAN = read.csv("./data/base/HIA_LIST/HIA/RISK_BINOMIAL_DF.csv", stringsAsFactors = FALSE)


#########################################################################################################################
## Now check the intersection and difference with the cleaned data
## Check for the cleaned match between evergreen and plant risk lists
length(intersect(RAW.HIA.SPP, RISK.BINOMIAL))                     ## 269 binomials match
length(intersect(HIA.list$Binomial,  RISK.BINOMIAL))              ## 311 binomials match
length(intersect(COMBO.NICHE.CONTEXT$searchTaxon, RISK.BINOMIAL)) ## 319 binomials match


## 311/654 41% overlap between the risk list and the evergreen list
(length(intersect(unique(HIA.list$Binomial), RISK.BINOMIAL))/length(unique(HIA.list$Binomial)))*100 


## Check for the cleaned difference  between evergreen and plant risk lists
setdiff(RISK.BINOMIAL, RAW.HIA.SPP)                       ## 605 plants are new
setdiff(RISK.BINOMIAL, unique(HIA.list$Binomial))         ## 564 plants are new
setdiff(RISK.BINOMIAL, COMBO.NICHE.CONTEXT$searchTaxon)   ## 554 plants are new


## 860/1131 (76%) difference between the risk list and the evergreen list
(length(setdiff(RAW.HIA.SPP, RISK.BINOMIAL))/length(RAW.HIA.SPP))*100 


## 555/610 (90%) difference between the risk list and the popular list
(length(setdiff(RISK.BINOMIAL, COMBO.NICHE.CONTEXT$searchTaxon))/length(COMBO.NICHE.CONTEXT$searchTaxon))*100


#########################################################################################################################
## Return the grower info for the species that overlap
## Instead of merge or join, just use %in%
RISK.LOW = RISK.LIST[, c("Plant.Name",
                         "Low.Risk")]
RISK.LOW =  dplyr::rename(RISK.LOW, 
                          searchTaxon = Plant.Name)

## The %in% operator is useful
EVERGREEN.RISK = COMBO.NICHE.CONTEXT[COMBO.NICHE.CONTEXT$searchTaxon %in% RISK.BINOMIAL, ][, c("searchTaxon",
                                                                                               "Origin",
                                                                                               "Top_200",
                                                                                               "Plant.type",
                                                                                               "Number.of.growers", 
                                                                                               "Number.of.States",
                                                                                               "No.of.Varieties")]

EVERGREEN.RISK = merge(EVERGREEN.RISK, RISK.LOW, by = "searchTaxon", all = FALSE)

## Check and save
dim(EVERGREEN.RISK)
#View(EVERGREEN.RISK)
write.csv(EVERGREEN.RISK, "./data/base/HIA_LIST/HIA/EVERGREEN_RISK_MATCH.csv", row.names = FALSE)  





#########################################################################################################################
## 5). MATCH EVERGREEN LIST WITH THE PLANT RISK LIST
#########################################################################################################################


## Get the intersection of all grown spp (that's the CLEAN list), the risky spp and the innovative spp (not sure what this is)
length(intersect(CLEAN.SPP$Binomial, RISK.BINOMIAL)) ## About 500 species in the extra list
EXTRA.SPP = CLEAN.SPP[CLEAN.SPP$Binomial %in% intersect(CLEAN.SPP$Binomial, RISK.BINOMIAL), ]
EXTRA.SPP = EXTRA.SPP[!EXTRA.SPP$Binomial %in% spp.all, ]
write.csv(EXTRA.SPP, "./data/base/HIA_LIST/COMBO/EXTRA_SPP.csv", row.names = FALSE)

## Check and just get the most popular of these?
dim(EXTRA.SPP)
View(EXTRA.SPP)
summary(EXTRA.SPP$Number.of.growers)


##
extra.spp = unique(EXTRA.SPP$Binomial)
length(extra.spp)

## Can we get another random selection of plants from the risky list to plug the gap?



#########################################################################################################################
## LIST EXCEPTIONS:
#########################################################################################################################


## Record each list: Raw top 25 (1135), Varieties (948), Binomials (610) 
## Check exceptions with Paul, Linda and Rach
length(unique(HIA.list$Species))     ## Raw top 25 (1135)
length(unique(HIA.VARIETY$Species))  ## Varieties  (948), excluding "spp.", eg Philodendron spp. Congo, Nandina domestica Moon Bay
length(unique(HIA.SPP$Binomial))     ## Binomials (610), keep Michelia yunnanensis Scented Pearl, exclude Spathiphyllum spp. Assorted


## record the "spp." weirdos
EXCLUDED.SPP         = setdiff(unique(RAW.HIA.SPP), unique(HIA.VARIETY$Species))
EXCLUDED.VARIETIES   = setdiff(unique(HIA.VARIETY$Species), unique(HIA.SPP$HIA.Taxa))   ## Here is the list that spots the exceptions!!!!!!!


## Which species that can't be modelled are on the test list or the top 200 list?
test.glasshouse = unique(c(MQ.glasshouse$Species, renee.full$Species))

intersect(MISSING$searchTaxon, test.spp)
intersect(MISSING$searchTaxon, test.glasshouse)
intersect(MISSING$searchTaxon, top.200$Species)


## What about the mod_2 species
setdiff(MOD_2$Species, HIA.RAW$Species)
setdiff(MOD_2$Species, HIA.list$Binomial)
setdiff(MOD_2$Species, COMBO.NICHE.CONTEXT$searchTaxon)
intersect(MOD_2$Species, test.spp)


## Now restrict the niche dataset to just the MOD2 species
# View(head(COMBO.ALL[COMBO.ALL$searchTaxon %in% MOD_2$Species, ], 13)[, c("searchTaxon",
#                                                                          "Origin",
#                                                                          "Top_200",
#                                                                          "Plant.type",
#                                                                          "Number.of.growers", 
#                                                                          "Number.of.States")])

MOD2.SPP = head(COMBO.ALL[COMBO.ALL$searchTaxon %in% MOD_2$Species, ], 15)[, c("searchTaxon",
                                                                               "Origin",
                                                                               "Top_200",
                                                                               "Plant.type",
                                                                               "Number.of.growers", 
                                                                               "Number.of.States")]


## Remaining anomalies:

## EG: Rhaphiolepis indica has growers for the spp and each variety, should we add them together?
## Magnolia grandiflora has 8 varieties which are being missed by the current code...

## Also, a key point for the future is how to treat the extra varieties, etc.






#########################################################################################################################
############################################  END OF HIA LIST CODE ###################################################### 
#########################################################################################################################