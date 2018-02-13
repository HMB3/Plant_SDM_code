#########################################################################################################################
########################################  CREATE SPECIES LIST ########################################################### 
#########################################################################################################################


#########################################################################################################################
# ## Setup for project 
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

#library("biglm")
#library("bigmemory")
#library("biganalytics")
#library("bigtabulate")


##
p <- c('ff',    'things',         'raster',    'dismo',        'sp',           'latticeExtra', 'data.table', 
       'rgdal', 'rgeos',          'gdalUtils', 'rmaxent',      'readr',        'dplyr',        'tidyr',
       'readr', 'rnaturalearth',  'rasterVis', 'RColorBrewer', 'latticeExtra', 'parallel')


## Require packages
sapply(p, require, character.only = TRUE)


## source functions
source('./R/GREEN_CITIES_FUNCTIONS.R')
source('./R/SDM_FUNCTIONS.R')
source('./R/MAPPING_FUNCTIONS.R')





#########################################################################################################################
## 1). READ IN DRAFT HIA LIST AND CLEAN
#########################################################################################################################

## This code takes the raw list of plants with 25 or more growers supplied by Matt Plumber and Anthony Manea, and
## then cleans the list as best as possible in R to use the species binomial as the unit for downloading and analysis.

## All the cleaining methods will throw up some anomalies, which need to be tracked, and checked with the team for how
## each case is treated (see outstanding tasks at the bottom)


#########################################################################################################################
## This list derives from all species and varieties sold anywhere in Australia in the last 5 years. Anthony Manea cleaned 
## Up the data and cross-linked to growth form and exotic/native status and derived a list of ~1000 species that are the 
## Most commonly sold, covering the right ratio of growth forms, regional representation and native/exotic
load("./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT_1601_2018.RData")
HIA.list   = read.csv("./data/base/HIA_LIST/HIA/GREEN_CITIES_DRAFT_LIST_2709_2017.csv", stringsAsFactors = FALSE)
CLEAN.list = read.csv("./data/base/HIA_LIST/HIA/HIA.CLEAN.csv",                         stringsAsFactors = FALSE)
GROWING    = read.csv("./data/base/HIA_LIST/HIA/database_aus_sp_growing.csv",           stringsAsFactors = FALSE)
MISSING    = read.csv("./data/base/HIA_LIST/HIA/MISSING_SPECIES.csv",                   stringsAsFactors = FALSE)


top.200              = read.csv("./data/base/HIA_LIST/HIA/HIA_TOP_200_1309_2017.csv",       stringsAsFactors = FALSE)
renee.full           = read.csv("./data/base/HIA_LIST/HIA/RENEE_FULL_LIST.csv",             stringsAsFactors = FALSE) 
renee.taxa           = read.csv("./data/base/HIA_LIST/HIA/RENEE_TAXA.csv",                  stringsAsFactors = FALSE)
renee.50             = read.csv("./data/base/HIA_LIST/HIA/RENEE_TOP_50.csv",                stringsAsFactors = FALSE)
MQ.glasshouse        = read.csv("./data/base/HIA_LIST/HIA/MQ_glasshouse.csv",               stringsAsFactors = FALSE)
Manuel.experimental  = read.csv("./data/base/HIA_LIST/HIA/Manuel_experimental_species.csv", stringsAsFactors = FALSE)


##
head(renee.taxa)
head(MQ.glasshouse)


## have a look
dim(HIA.list)
dim(CLEAN.list)
dim(renee.taxa)
str(HIA.list)
head(HIA.list)

##
#intersect(CLEAN.list$Species, GROWING$scientific_name)

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
test.spp = trimws(sort(unique(c(renee.full$Species, 
                                Manuel.experimental$Species,
                                MQ.glasshouse$Species,
                                "Betula pendula", "Fraxinus excelsior", "Quercus robur", "Fagus sylvatica"))))


## Combine with 45 from the main list
HIA.SAMPLE = head(COMBO.NICHE.CONTEXT, 53)[, c("searchTaxon")]
test.spp   = sort(unique(c(test.spp, HIA.SAMPLE)))
str(test.spp)
test.reverse = sort(test.spp, decreasing = TRUE)

##
exp.spp  = c('Swainsona formosa', 'Templetonia retusa', 'Dodonaea baueri', 'Platanus hispanica', 'Kennedia beckxiana')  
exp.rev  = sort(exp.spp, decreasing = TRUE)
miss.spp = c('Corymbia tessellaris', 'Metrosideros excelsa')


## Create lists for the mapping code
all_spp     = gsub(" ", "_", spp.all)
all_reverse = sort(all_spp, decreasing = TRUE)

test_spp     = gsub(" ", "_", test.spp)
test_reverse = sort(test_spp, decreasing = TRUE)

exp_spp      = gsub(" ", "_", exp.spp)
miss_spp     = gsub(" ", "_", miss.spp)


## Which species are only on the test list?
setdiff(test.spp, spp.all)

  


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


## Remaining anomalies:

## EG: Rhaphiolepis indica has growers for the spp and each variety, should we add them together?
## Magnolia grandiflora has 8 varieties which are being missed by the current code...

## Also, a key point for the future is how to treat the extra varieties, etc.






#########################################################################################################################
############################################  END OF HIA LIST CODE ###################################################### 
#########################################################################################################################