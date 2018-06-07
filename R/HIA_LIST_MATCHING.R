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
rasterOptions(tmpdir = file.path('H:/green_cities_sdm/RTEMP')) 


#########################################################################################################################
## Read in spatial data
aus         = readRDS("F:/green_cities_sdm/data/base/CONTEXTUAL/aus_states.rds")
LAND        = readRDS("F:/green_cities_sdm/data/base/CONTEXTUAL/LAND_world.rds")
areal_unit  = readRDS("F:/green_cities_sdm/data/base/CONTEXTUAL/SUA.rds")
areal_unit  = areal_unit[order(areal_unit$SUA_NAME11),]
Koppen      = readRDS('data/base/CONTEXTUAL/Koppen_1975.rds')

   
# Kopp.future = readOGR("F:/green_cities_sdm/data/base/CONTEXTUAL/CM10_Kop_Shp_V1.2/CM10_Kop_V1.2.shp", 
#                       layer = "CM10_Kop_V1.2")


## Set definitions :: best to minimise the number of projection used in this project
CRS.MOL      <- CRS('+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs')
CRS.WGS.84   <- CRS("+init=epsg:4326")
CRS.AUS.ALB  <- CRS("+init=EPSG:3577")
ALB.CONICAL  <- CRS('+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')




#########################################################################################################################
## 1). READ IN DRAFT HIA LIST AND CLEAN
#########################################################################################################################


#########################################################################################################################
## Read in the niche data
COMBO.NICHE.CONTEXT = readRDS("./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT_APRIL_2018_STANDARD_CLEAN.rds")
CLEAN.NICHE.CONTEXT = readRDS("./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT_APRIL_2018_COORD_CLEAN.rds")
MAXENT.CHECK        = read.csv("./output/maxent/MAXENT_CHECK_RATING.csv",       stringsAsFactors = FALSE)
OVERALL.LOSS        = read.csv("./output/tables/OVERALL_LOSS.csv",              stringsAsFactors = FALSE)
APPENDIX            = read.csv("./data/base/HIA_LIST/COMBO/Appendix_table.csv", stringsAsFactors = FALSE)
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
MOD.2.3             = read.csv("./data/base/HIA_LIST/HIA/MODULE_2_3.csv",                        stringsAsFactors = FALSE)
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
Manuel.group         = read.csv("./MANUEL/SUA_by_SPP.csv",                                  stringsAsFactors = FALSE)
TREE.NETWORK         = read.csv("./MANUEL/COMBO_APNI.csv",                                  stringsAsFactors = FALSE)
NURSE.MATCH          = read.csv("./MANUEL/nurseries.csv",                                   stringsAsFactors = FALSE)
campbelltown         = read.csv("./data/base/HIA_LIST/HIA/campbelltown_species.csv",        stringsAsFactors = FALSE)


##
intersect(campbelltown$Species, CLEAN.NICHE.CONTEXT$searchTaxon)
intersect(NURSE.MATCH$species, CLEAN.NICHE.CONTEXT$searchTaxon)
new.spp  = trimws(unique((sort(c(NURSE.MATCH$species, campbelltown$Species)))))
camp.spp = trimws(campbelltown$Species)


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
spp.model            = 
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


## Now sum the number of growers across multiple varieties
## Just get the count for each unique binomial
n <- tapply(HIA.list$Number.of.growers, HIA.list$Binomial, sum, na.rm = TRUE)
HIA.list$Number.of.growers.total <- n[HIA.list$Binomial]
TOT.GROW = HIA.list[c("Binomial",
                        "Number.of.growers.total")]
names(TOT.GROW) = c("searchTaxon", "Total.growers")
TOT.GROW        = TOT.GROW[!duplicated(TOT.GROW[,c('searchTaxon')]),] 


## Join the total growers to the NICHE data
COMBO.NICHE.CONTEXT = join(COMBO.NICHE.CONTEXT, TOT.GROW)
COMBO.NICHE.CONTEXT =  COMBO.NICHE.CONTEXT[, c(1:14, 199, 16:198)] 
names(COMBO.NICHE.CONTEXT[1:15])


CLEAN.NICHE.CONTEXT = join(CLEAN.NICHE.CONTEXT, TOT.GROW)
CLEAN.NICHE.CONTEXT =  CLEAN.NICHE.CONTEXT[, c(1:14, 199, 16:198)] 
names(CLEAN.NICHE.CONTEXT[1:15])


## Important, join on the total growers here
names(COMBO.NICHE.CONTEXT)


#########################################################################################################################
## Update Manuel's data
TREE.NETWORK = TREE.NETWORK[,c("searchTaxon",      "Plant.type",      "Origin",    "APNI",    "PLANTED_GROWING",  "Top_200",             
                               "ACT",              "NSW",             "NT",        "QLD",     "SA",    "VIC",    "WA",              
                               "Number.of.States", "No.of.Varieties", "Number.of.growers",    "AREA_OCCUPANCY",  "LGA.AGG")]
names(TREE.NETWORK)

## Merge Manules data with the new climate data
CLEAN.CLIMATE = COMBO.NICHE.CONTEXT[, c(1, 15, 19:198)] 
CLEAN.CLIMATE = join(TREE.NETWORK, CLEAN.CLIMATE)
names(CLEAN.CLIMATE)
write.csv(CLEAN.CLIMATE, "./data/base/HIA_LIST/COMBO/TREE_NETWORK_NICHES.csv", row.names = FALSE)


## Check this reduces the number of Top 200 missing from this list
HIA.list = merge(HIA.list, spp.200, by = "Binomial", all.x = TRUE) 
HIA.list$Top_200[is.na(HIA.list$Top_200)] <- "FALSE"
HIA.list$Origin <- gsub(" ",  "", HIA.list$Origin)


## check
str(HIA.list)
head(HIA.list)
unique(HIA.list$Top_200)
unique(HIA.list$Origin)
length(unique(HIA.list$Binomial)) ## 654 unique binomials
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
## of magnolia, currently I'm not getting the most popular variety.


## Some spp/varieties are not working..............................................................................................
HIA.VARIETY <- 
  DRAFT.HIA.TAXA$Binomial[DRAFT.HIA.TAXA$Binomial != DRAFT.HIA.TAXA$Species] %>% 
  table %>% 
  as.data.frame %>% 
  setNames(c('Binomial', 'No.of.Varieties')) %>% 
  full_join(DRAFT.HIA.TAXA) %>% 
  dplyr::select(Species, Binomial, No.of.Varieties, Plant.type:WA, Number.of.growers.total, Number.of.States, Origin, Top_200)

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
#glasshouse.spp = as.data.frame(glasshouse.spp)

                        
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
                                                                                               "Total.growers", 
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
#View(EXTRA.SPP)
summary(EXTRA.SPP$Number.of.growers)


##
spp.extra = unique(EXTRA.SPP$Binomial)
length(spp.extra)




#########################################################################################################################
## 6). MILESTONE JUNE 2018 LIST
#########################################################################################################################


## We want all the current experimental spp (so the 50 as of April 2018). 
## The 16 species Manuel is monitoring
## Plus another 70 of the right sort
MOD.2 = COMBO.NICHE.CONTEXT[COMBO.NICHE.CONTEXT$searchTaxon %in% MOD_2$Species, ][,  c("searchTaxon",
                                                                                       "COMBO.count",
                                                                                       "AUS_RECORDS",
                                                                                       "Total.growers",
                                                                                       "Top_200")]


## So find 80 species which have the most growers, have the most Aus records and with at least 50 of them trees
dim(subset(COMBO.NICHE.CONTEXT, Total.growers > 25 & AUS_RECORDS > 20 & COMBO.count > 100 & Top_200 == "TRUE"))

MILE.1         = subset(COMBO.NICHE.CONTEXT, Total.growers > 25 & AUS_RECORDS > 20 & COMBO.count > 100 & Top_200 == "TRUE")
MILE.CLEAN     = subset(CLEAN.NICHE.CONTEXT, Total.growers > 25 & AUS_RECORDS > 20 & COMBO.count > 100 & Top_200 == "TRUE")
MILE.1.EXTRA   = subset(CLEAN.NICHE.CONTEXT, Total.growers > 25 & AUS_RECORDS > 20 & COMBO.count > 100)
MILE.1.EXTRA   = MILE.1.EXTRA[with(MILE.1.EXTRA, order(-Total.growers)), ]
MILE.1.EXTRA   = MILE.1.EXTRA [!MILE.1.EXTRA$searchTaxon %in% MILE.CLEAN$searchTaxon, ]
summary(MILE.1.EXTRA$AUS_RECORDS)
summary(MILE.1.EXTRA$COMBO.count)

spp.mile.extra     = head(MILE.1.EXTRA$searchTaxon, 84)
spp.mile.1         = unique(sort(c(MOD.2.3$Species, MILE.1$searchTaxon, spp.mile.extra)))
length(spp.mile.1)


## Join on the column which shows the kind of data bias
MILE.1.SPP     =  COMBO.NICHE.CONTEXT[COMBO.NICHE.CONTEXT$searchTaxon %in% spp.mile.1, ]
MILE.1.SPP     =  join(MILE.1.SPP, SPP.BIAS, type = "left")

MILE.CLEAN.SPP =  CLEAN.NICHE.CONTEXT[CLEAN.NICHE.CONTEXT$searchTaxon %in% spp.mile.1, ]
MILE.CLEAN.SPP =  join(MILE.CLEAN.SPP, SPP.BIAS, type = "left")

MILE.1.SPP = MILE.1.SPP[,  c(1, 15, 16, 18, 7, 17, 3, 4, 5, 199, 19:198)]
MILE.1.SPP = MILE.1.SPP[with(MILE.1.SPP, order(-Total.growers)), ]

MILE.CLEAN.SPP = MILE.CLEAN.SPP[,  c(1, 15, 16, 18, 7, 17, 3, 4, 5, 199, 19:198)]
MILE.CLEAN.SPP = MILE.CLEAN.SPP[with(MILE.CLEAN.SPP, order(-Total.growers)), ]

#View(MILE.1.SPP)
#View(MILE.CLEAN.SPP)


## What is the distribution?
dim(MILE.1.SPP)
summary(MILE.1.SPP$Total.growers)
summary(MILE.1.SPP$AUS_RECORDS)
summary(MILE.1.SPP$COMBO.count)


## Also how should these be model
with(MILE.1.SPP, table(MILE.1.SPP$AUS_BOUND_BIAS))
spp.mile.targ = subset(MILE.1.SPP, AUS_BOUND_BIAS == "FALSE")$searchTaxon
spp.mile.rand = subset(MILE.1.SPP, AUS_BOUND_BIAS == "TRUE")$searchTaxon
spp.mile      = sort(unique(MILE.1.SPP$searchTaxon))
spp.mile.rev  = sort(spp.mile, decreasing = TRUE)
   
spp_mile      = gsub(" ", "_", spp.mile)
spp_mile_rev  = sort(spp_mile, decreasing = TRUE)

spp_new      = gsub(" ", "_", new.spp)
spp_new_rev  = sort(spp_new, decreasing = TRUE)


spp_camp     = gsub(" ", "_", camp.spp)


## Add plant type data for the missing species not sure here
unique(MILE.CLEAN.SPP$Plant.type)

## Rename the species with missing functional types
MILE.CLEAN.SPP[193, "Plant.type"] = "Tree"  ; MILE.CLEAN.SPP[193, "Origin"] = "Native"
MILE.CLEAN.SPP[194, "Plant.type"] = "Tree"  ; MILE.CLEAN.SPP[194, "Origin"] = "Native"
MILE.CLEAN.SPP[195, "Plant.type"] = "Shrub" ; MILE.CLEAN.SPP[195, "Origin"] = "Native"
MILE.CLEAN.SPP[196, "Plant.type"] = "Shrub" ; MILE.CLEAN.SPP[196, "Origin"] = "Native"
MILE.CLEAN.SPP[197, "Plant.type"] = "Tree"  ; MILE.CLEAN.SPP[197, "Origin"] = "Native"
MILE.CLEAN.SPP[198, "Plant.type"] = "Tree"  ; MILE.CLEAN.SPP[198, "Origin"] = "Native"


## And create a table of the functional types
with(MILE.1, table(Plant.type,     useNA = "always"))
with(MILE.1.SPP, table(Plant.type, useNA = "always"))


## Origin 
MILE.CLEAN.SPP[180, "Orign"] = "Native"
MILE.CLEAN.SPP[182, "Orign"] = "Native"
MILE.CLEAN.SPP[183, "Orign"] = "Native"
MILE.CLEAN.SPP[184, "Orign"] = "Native"
MILE.CLEAN.SPP[188, "Orign"] = "Exotic"
MILE.CLEAN.SPP[192, "Orign"] = "Native"
View(MILE.CLEAN.SPP)


## Check the join was ok
length(unique(MILE.CLEAN.SPP$searchTaxon))


########################################################################################################################
## Proportions
round(with(MILE.CLEAN.SPP, table(Plant.type)/sum(table(Plant.type))*100), 1)


## 
spp.combo      = sort(unique(c(spp.all, spp.extra, spp.mile, MILE.CLEAN.SPP$searchTaxon)))
combo.rev      = sort(spp.combo, decreasing = TRUE)
combo_spp      = gsub(" ", "_", spp.combo)
combo_reverse  = sort(combo_spp, decreasing = TRUE)
length(spp.combo)


## Add the module to the species output
renee.full$Module = 3
MOD_2$Module      = 2
spp.modules = join(renee.full, MOD_2, type = "full")
names(spp.modules) = c("searchTaxon", "Module")


## Join modules to milestones
MILE.CLEAN.MODULE = join(MILE.CLEAN.SPP, spp.modules, type = "left")
MILE.CLEAN.MODULE$Module[is.na(MILE.CLEAN.MODULE$Module)] <- 1
MILE.CLEAN.MODULE = MILE.CLEAN.MODULE[MILE.CLEAN.MODULE$searchTaxon %in% MILE.CLEAN.SPP$searchTaxon, ]
MILE.CLEAN.MODULE = MILE.CLEAN.MODULE[!duplicated(MILE.CLEAN.MODULE$searchTaxon),]
 

## Bar plot 
module.counts <- table(MILE.CLEAN.MODULE$Module)
barplot(module.counts, main = "Species per module",
        xlab = "Module", ylab = "number of species", col = c("lightblue", "pink", "orange"))
        #legend = rownames(counts))


## Pie chart
pie(module.counts, col = c("lightblue", "pink", "orange"))


########################################################################################################################
# Check the overlap between the models that have been checked, and the most popular
setdiff(MILE.CLEAN.MODULE$searchTaxon, MAXENT.CHECK$searchTaxon)
top.losers = head(unique(OVERALL.LOSS$SPECIES), 20)
top.losers = gsub("_",  " ", top.losers)


##
MILE.LOSE  = MILE.CLEAN.MODULE[MILE.CLEAN.MODULE$searchTaxon %in% top.losers, ]
MILE.LOSE  = MILE.LOSE[rev(order(MILE.LOSE$Total.growers)),]
View(MILE.LOSE)


## Intersect existing list with the number of growers
APP.GROW = COMBO.NICHE.CONTEXT[,  c("searchTaxon",
                                    "Total.growers")]
APPENDIX = join(APPENDIX, APP.GROW)


## Proportions for this appendix
with(APPENDIX , table(Plant.type))
round(with(APPENDIX , table(Plant.type)/sum(table(Plant.type))*100), 1)


## Save ::
## setdiff(MOD.2.3$Species, MILE.1.SPP$searchTaxon)
write.csv(MILE.1.SPP,          "./data/base/HIA_LIST/COMBO/MILESTONE_TAXA_APRIL_2018_STANDARD_CLEAN.csv",      row.names = FALSE)
write.csv(MILE.CLEAN.MODULE,   "./data/base/HIA_LIST/COMBO/MILESTONE_TAXA_APRIL_2018_COORD_CLEAN_MODULES.csv", row.names = FALSE)
write.csv(MOD.2,               "./data/base/HIA_LIST/COMBO/MOD_2_SPP_RECORDS.csv",                             row.names = FALSE)
write.csv(COMBO.NICHE.CONTEXT, "./data/base/HIA_LIST/COMBO/COMBO_TOT_GROWERS.csv",                             row.names = FALSE)
write.csv(APPENDIX,            "./data/base/HIA_LIST/COMBO/APPENDIX_TOT_GROWERS.csv",                          row.names = FALSE)



########################################################################################################################
## Final counts
dim(COMBO.NICHE.CONTEXT) 
HIA.COUNT   = COMBO.NICHE.CONTEXT[COMBO.NICHE.CONTEXT$searchTaxon %in% HIA.SPP$Binomial, ]
LOW.RECORDS = dim(subset (HIA.COUNT, COMBO.count < 20))[1] + (dim(HIA.SPP)[1] - dim(HIA.COUNT)[1]) 
LOW.NO      = dim(HIA.SPP)[1] - 455


#########################################################################################################################
## LIST EXCEPTIONS:
#########################################################################################################################


## Record each list: Raw top 25 (1135), Varieties (948), Binomials (610) 
## Check exceptions with Paul, Linda and Rach
length(unique(HIA.list$Species))     ## Raw top 25 (1135)
length(unique(HIA.VARIETY$Species))  ## Varieties  (948), excluding "spp.", eg Philodendron spp. Congo, Nandina domestica Moon Bay
length(unique(HIA.SPP$Binomial))     ## Binomials  (610), keep Michelia yunnanensis Scented Pearl, exclude Spathiphyllum spp. Assorted


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
# View(head(COMBO.NICHE.CONTEXT[COMBO.NICHE.CONTEXT$searchTaxon %in% MOD_2$Species, ], 13)[, c("searchTaxon",
#                                                                          "Origin",
#                                                                          "Top_200",
#                                                                          "Plant.type",
#                                                                          "Total.growers", 
#                                                                          "Number.of.States")])

MOD2.SPP = COMBO.NICHE.CONTEXT[COMBO.NICHE.CONTEXT$searchTaxon %in% MOD_2$Species, ]
MOD2.SPP = MOD2.SPP[rev(order(MOD2.SPP$AUS_RECORDS)),]
#View(MOD2.SPP)


## Remaining anomalies:

## EG: Rhaphiolepis indica has growers for the spp and each variety, should we add them together?
## Magnolia grandiflora has 8 varieties which are being missed by the current code...

## Also, a key point for the future is how to treat the extra varieties, etc.






#########################################################################################################################
############################################  END OF HIA LIST CODE ###################################################### 
#########################################################################################################################