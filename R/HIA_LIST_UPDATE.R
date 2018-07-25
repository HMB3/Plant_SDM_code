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
library(spThin)


##
p <- c('ff',    'things',         'raster',    'dismo',        'sp',           'latticeExtra', 'data.table', 
       'rgdal', 'rgeos',          'gdalUtils', 'rmaxent',      'readr',        'dplyr',        'tidyr',
       'readr', 'rnaturalearth',  'rasterVis', 'RColorBrewer', 'latticeExtra', 'parallel')


## Require packages
sapply(p, require, character.only = TRUE)


## Source functions
source('./R/GREEN_CITIES_FUNCTIONS.R')
source('./R/MAXENT_FUNCTIONS.R')
source('./R/MAPPING_FUNCTIONS.R')
source('./R/HIA_CLEAN_MATCHING.R')
rasterOptions(tmpdir = file.path('H:/green_cities_sdm/RTEMP')) 


#########################################################################################################################
## Read in spatial data once, rather than in each script
aus           = readRDS("F:/green_cities_sdm/data/base/CONTEXTUAL/aus_states.rds")
LAND          = readRDS("F:/green_cities_sdm/data/base/CONTEXTUAL/LAND_world.rds")
areal_unit    = readRDS("F:/green_cities_sdm/data/base/CONTEXTUAL/SUA.rds")
areal_unit    = areal_unit[order(areal_unit$SUA_NAME11),]
Koppen        = readRDS('data/base/CONTEXTUAL/Koppen_1975.rds')
Koppen_zones  = unique(readOGR('data/base/CONTEXTUAL/WC05_1975H_Koppen_Shapefile/WC05_1975H_Koppen_Kriticos_2012.shp')@data[, 1:2])
Koppen_1975   = raster('data/Koppen_1000m_Mollweide54009.tif')
AUS_RAIN      = readRDS('data/base/CONTEXTUAL/BOM/BOM_RAIN_AGG.rds')

SUA           = readRDS("F:/green_cities_sdm/data/base/CONTEXTUAL/IN_SUA_AUS.rds")
LGA           = readRDS("F:/green_cities_sdm/data/base/CONTEXTUAL/LGA.rds")
AUS           = readRDS("F:/green_cities_sdm/data/base/CONTEXTUAL/aus_states.rds")
ALL.SUA.POP   = read.csv("./data/base/CONTEXTUAL/ABS_SUA_POP.csv", stringsAsFactors = FALSE)


##
template.raster = raster("./data/template_hasData.tif")
template.cells  = readRDS("./data/hasData_cells.rds")
load("./data/base/CONTEXTUAL/urbanareas.rda")


# AUS_RAIN = readOGR("H:/green_cities_sdm/data/base/CONTEXTUAL/BOM/BOM_RAIN_AGG.shp",
#                    layer = "BOM_RAIN_AGG")
# saveRDS(AUS_RAIN, file = paste("./data/base/CONTEXTUAL/BOM/BOM_RAIN_AGG.rds"))


## Set definitions :: best to minimise the number of projection used in this project
CRS.MOL      <- CRS('+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs')
CRS.MOL.SDM  <- CRS('+init=ESRI:54009 +proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs84=0,0,0')
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
MAXENT.CHECK        = read.csv("./output/maxent/MAXENT_CHECK_1707_2018.csv",    stringsAsFactors = FALSE)
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
MAXENT.RATING       = read.csv("./output/maxent/MAXENT_RATING_26_2018.csv",                      stringsAsFactors = FALSE)
TI.LIST             = read.csv("./data/base/HIA_LIST/COMBO/ALE_TREE_SPP_LIST.csv",               stringsAsFactors = FALSE)
TI.XY               = read.csv("./data/base/HIA_LIST/COMBO/ALE_TREE_SPP_XY.csv",                 stringsAsFactors = FALSE)


## Check Ale's data
names(TI.LIST)
names(TI.XY)


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
head(TI.XY)


## What is the intersection betwen the evergreen list and the Tree inventories?
TREE.HIA = intersect(head(TI.LIST, 426)$searchTaxon, CLEAN.SPP$Binomial)
TREE_HIA = gsub(" ", "_", TREE.HIA)


EVERGREEN = CLEAN.SPP[c("Binomial", "Number.of.growers", "Number.of.States", "Origin")]
names(EVERGREEN )[names(EVERGREEN) == 'Binomial'] <- 'searchTaxon'
TREE.EVERGREEN = merge(EVERGREEN, TI.LIST)
TREE.EVERGREEN = TREE.EVERGREEN[TREE.EVERGREEN$searchTaxon %in% TREE.HIA, ]
TREE.EVERGREEN = TREE.EVERGREEN [with(TREE.EVERGREEN, rev(order(Plantings))), ]


## What does the dataset look like?
round(with(TREE.EVERGREEN, table(Origin)/sum(table(Origin))*100), 1)


## Test the new urban data on a subset of species
test.exotics = c("Platanus acerifolia", "Pyrus calleryana",  "Jacaranda mimosifolia")
test_exotics = gsub(" ", "_", test.exotics)


## The list of species with checked maxent maps
SPP.BIAS      = read.csv("./output/maxent/SPP_BOUNDARY_BIAS.csv", stringsAsFactors = FALSE)
SPP.BIAS      = merge(SPP.BIAS, COMBO.NICHE.CONTEXT[c("searchTaxon", "Origin")])
SPP.BIAS      = subset(SPP.BIAS, AUS_BOUND_BIAS == "TRUE" & Origin == "Native")
SPP.BIAS      = subset(SPP.BIAS)$searchTaxon
SPP_BIAS      = gsub(" ", "_", SPP.BIAS)


#########################################################################################################################
## European garden nursery lists
EURO.SPP            = read.csv("./data/base/HIA_LIST/URBAN/Euro_garden_flora_spp.csv",            stringsAsFactors = FALSE)
EURO.NURSE          = read.csv("./data/base/HIA_LIST/URBAN/Euro_garden_flora_spp_nurse.csv",      stringsAsFactors = FALSE)
EURO.NURSE.LOC      = read.csv("./data/base/HIA_LIST/URBAN/Euro_garden_flora_nurse_location.csv", stringsAsFactors = FALSE)


#########################################################################################################################
## Experimental trait lists ::
HEAT.RISK  = read.csv("./data/base/HIA_LIST/RENEE/MOD3_HEAT_RANKS_072018.csv",                   stringsAsFactors = FALSE)
TRAIT.SPP  = read.csv("./data/base/HIA_LIST/RENEE/MOD3_TRAIT_SPP.csv",                           stringsAsFactors = FALSE)


## Now find the match between the trait species and the trait species... 
colnames(HEAT.RISK)[colnames(HEAT.RISK)=="Species"] <- "searchTaxon"
colnames(TRAIT.SPP)[colnames(TRAIT.SPP)=="Species"] <- "searchTaxon"
TRAIT.SPP = merge(TRAIT.SPP, HEAT.RISK, by = "searchTaxon")
trait.spp = unique(TRAIT.SPP$searchTaxon)
trait_spp = gsub(" ", "_", trait.spp)


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


## The species for the SUA analysis
SUA.spp = read.csv("./output/maxent/MAXENT_SUA_SPP.csv", stringsAsFactors = FALSE)
SUA.spp = SUA.spp[order(SUA.spp$searchTaxon),] 
SUA.spp = unique(SUA.spp$searchTaxon)
SUA_spp = gsub(" ", "_", SUA.spp)


##
intersect(campbelltown$Species, CLEAN.NICHE.CONTEXT$searchTaxon)
intersect(NURSE.MATCH$species, CLEAN.NICHE.CONTEXT$searchTaxon)
new.spp  = trimws(unique((sort(c(NURSE.MATCH$species, campbelltown$Species)))))
camp.spp = trimws(campbelltown$Species)


## Get a list of species with > 20 records
summary(COMBO.NICHE.CONTEXT$AUS_RECORDS)
SPP.AUS = subset(COMBO.NICHE.CONTEXT, AUS_RECORDS >= 20)
summary(SPP.AUS$COMBO.count);summary(SPP.AUS$AUS_RECORDS)


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


#########################################################################################################################
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
#write.csv(CLEAN.CLIMATE, "./data/base/HIA_LIST/COMBO/TREE_NETWORK_NICHES.csv", row.names = FALSE)


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
## 2). DRAFT LIST PREP
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
#save(test.spp, file = paste("./data/base/HIA_LIST/GBIF/test.spp.RData", sep = ""))





#########################################################################################################################
## 4). WHAT IS THE MATCH BETWEEN THE EVERGREEN LIST AND THE EUROPEAN NURSERY LIST
#########################################################################################################################


#########################################################################################################################
## What is the unique species list of the European nurseries?
EURO.SPP = unique(sort(EURO.SPP$Updated_species_name))


## Check for the raw match 
length(intersect(RAW.HIA.SPP, EURO.SPP))                       ## 258  binomials match
length(intersect(HIA.list$Binomial,  EURO.SPP))                ## 321  binomials match
length(intersect(top.200$Species,    EURO.SPP))                ## 127  binomials match
length(intersect(CLEAN.NICHE.CONTEXT$searchTaxon, EURO.SPP))   ## 1575 binomials match


## Check for the raw difference
length(setdiff(RAW.HIA.SPP, EURO.SPP))                         
length(intersect(HIA.list$Binomial,  EURO.SPP))                
length(intersect(top.200$Species,    EURO.SPP))                
length(intersect(CLEAN.NICHE.CONTEXT$searchTaxon, EURO.SPP))   


#########################################################################################################################
## Now create a table of species * lat/long
names(EURO.NURSE)
names(EURO.NURSE.LOC)
head(EURO.NURSE);head(EURO.NURSE.LOC)
names(EURO.NURSE.LOC)[names(EURO.NURSE.LOC) == 'Nursery_Location'] <- 'Nursery'
identical(length(unique(EURO.NURSE$Nursery)), length(unique(EURO.NURSE.LOC$Nursery)))


## Join the lat/long onto the species names
EURO.NURSE    = EURO.NURSE [, c("Species", "Nursery")]
SPP.NURSE.LOC = join(EURO.NURSE, EURO.NURSE.LOC, type = "left")


## Check
dim(SPP.NURSE.LOC)
dim(EURO.NURSE)
names(SPP.NURSE.LOC)


## How many species are singletons?
NURSE.count   = as.data.frame(table(SPP.NURSE.LOC$Species))
names(NURSE.count) = c("Species", "Count")
NURSE.count = NURSE.count[with(NURSE.count, rev(order(Count))), ]
summary(NURSE.count)

## Then, check the distribution across Europe!





#########################################################################################################################
## 6). MILESTONE JUNE 2018 LIST
#########################################################################################################################


## We want all the current experimental spp (so the 50 as of April 2018).
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


#########################################################################################################################
## Now just get the species with good models and enough data
MODEL.CHECK = merge(CLEAN.NICHE.CONTEXT, MAXENT.CHECK)
MODEL.CHECK = MODEL.CHECK[, c("searchTaxon",
                              "COMBO.count",
                              "AUS_RECORDS",
                              "Total.growers",
                              "Origin",
                              "Plant.type",
                              "Top_200",                                                                        
                              "CHECK_MAP")]


## Add a column for module two and 3
non.HIA = intersect(setdiff(MOD.2.3$Species, top.200$Species), SUA.spp)
MODEL.CHECK$MODULE = ifelse(MODEL.CHECK$searchTaxon %in% MOD.2.3$Species, "Two_three", "One")
MODEL.CHECK = MODEL.CHECK[with(MODEL.CHECK, order(-Total.growers)), ]
View(MODEL.CHECK)


#########################################################################################################################
## Now how many trees do we have? 68 modelled, another 86 with ok data.
## This will increase with Alessandro's data. Re-create niche file, and re-run the criteria.
## Let's get another 40 exotic trees with good models
MODEL.CHECK[1, "Plant.type"]
checked.trees = subset(MODEL.CHECK, Total.growers >= 25 & CHECK_MAP <= 2 & Plant.type == "Tree")$searchTaxon ## 67
subset(MODEL.CHECK, Total.growers >= 25 & CHECK_MAP <= 2 & Plant.type == "Tree" & Origin == "Exotic")$searchTaxon
exotic.trees = subset(COMBO.NICHE.CONTEXT, Total.growers >= 25 & COMBO.count > 300 & Plant.type == "Tree" & Origin == "Exotic")$searchTaxon
exotic_trees = gsub(" ", "_", exotic.trees)


#########################################################################################################################
## Get the match between ALE's tree data and HIA data
extra.trees   = subset(MILE.1.EXTRA, Plant.type == "Tree")$searchTaxon
MS.trees      = unique(c(checked.trees, extra.trees))
new.trees     = setdiff(extra.trees, checked.trees)
new_trees     = gsub(" ", "_", new.trees)

# ALE.TREE = merge(TI.LIST, CLEAN.NICHE.CONTEXT[c("searchTaxon",  "Plant.type", "Origin", "Total.growers", "COMBO.count")])
# ALE.SPP  = subset(ALE.TREE, Total.growers >= 25 & Frequency > 300 & Plant.type == "Tree")$searchTaxon
# write.csv(ALE.TREE,  "./data/base/HIA_LIST/COMBO/ALE_TREE_MATCH.csv",    row.names = FALSE)

## Now create a table for the supplementary material
## Spp, type, origin, growers, Count
# MS.spp = sort(unique(c(checked.trees, exotic.trees, ALE.SPP)))
# MS.COL = CLEAN.NICHE.CONTEXT[c("searchTaxon",  "Plant.type", "Origin", "Total.growers", "COMBO.count")]
# MS.SPP = MS.COL[MS.COL$searchTaxon %in% MS.spp, ]
# MS.SPP = MS.SPP[with(MS.SPP, rev(order(Total.growers))), ] 
# dim(MS.SPP)
# 
# 
# ## Now create a table of the exotics v natives
# with(MS.SPP, table(Origin))
# round(with(MS.SPP, table(Origin)/sum(table(Origin))*100), 1)
# write.csv(MS.SPP,  "./data/base/HIA_LIST/COMBO/MS_SPP_TABLE.csv", row.names = FALSE)





#########################################################################################################################
############################################  END OF HIA LIST CODE ###################################################### 
#########################################################################################################################