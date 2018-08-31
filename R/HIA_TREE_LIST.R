#########################################################################################################################
########################################  CREATE SPECIES LISTS ########################################################## 
#########################################################################################################################


## This code takes the raw list of plants with 25 or more growers supplied by Matt Plumber and Anthony Manea, and
## then cleans the list as best as possible in R to use the species binomial as the unit for downloading and analysis.

## All the cleaining methods will throw up some anomalies, which need to be tracked, and checked with the team for how
## each case is treated (see outstanding tasks at the bottom)


#########################################################################################################################
## Setup for project 
# library(gtools)
# library(GISTools)
# library(devtools)
# library(Rcpp)
# library(raster)
# library(rgdal)
# library(plyr)
# library(dplyr)
# library(sfsmisc)
# #library(spatstat)
# library(data.table)
# library(rvest)
# #library(vegan)
# 
# library(SDMTools)
# library(rmaxent)
# library(dismo)
# #library(BiodiversityR)
# library(biomod2)
# #library(AdaptR)
# library(red)
# library(ConR)
# library(dat)
# 
# library(ff)
# library(rgeos)
# library(sp)
# library(raster)
# library(rJava)
# library(things)
# library(digest)
# 
# library(ALA4R)
# library(rgbif)
# library(scrubr)
# library(RCurl)
# library(httr)
# 
# library(taxonlookup)
# library(Taxonstand)
# #library(speciesgeocodeR)
# library(CoordinateCleaner)
# library(spatialEco)
# library(raster)
# library(rnaturalearth)
# library(gdalUtils)
# 
# #library(knitr)
# #library(htmltools)
# library(yaml)
# library(caTools)
# library(bitops)
# library(rmarkdown)
# library(gsubfn)
# library(functional)
# library(splitstackshape)
# 
# library(tidyverse)
# library(stringr)
# library(maptools)
# #library(ggmap)
# library(rgeos)
# library(magrittr)
# library(datastorr)
# #library(baad.data)
# #library(rapportools)
# 
# library(Cairo)
# library(lattice)
# library(latticeExtra)
# library(PerformanceAnalytics)
# library(timetk)
# library(spThin)


## Load only the packages needed for the analysis
p <- c('ff',    'things', 'raster',        'dismo',        'sp',           'latticeExtra', 'data.table', 
       'rgdal', 'rgeos',  'gdalUtils',     'rmaxent',      'readr',        'plyr',         'dplyr',        
       'tidyr', 'readr',  'rnaturalearth', 'rasterVis',    'RColorBrewer', 'latticeExtra', 'parallel',     
       'taxonlookup',     'ALA4R',         'stringr',      'Taxonstand',   'CoordinateCleaner')


## Require packages
sapply(p, require, character.only = TRUE)


## Source functions
source('./R/GREEN_CITIES_FUNCTIONS.R')
source('./R/MAXENT_FUNCTIONS.R')
source('./R/MAPPING_FUNCTIONS.R')
source('./R/HIA_CLEAN_MATCHING.R')
rasterOptions(tmpdir = file.path('H:/green_cities_sdm/RTEMP')) 





#########################################################################################################################
## 1). READ IN SPATIAL DATA
#########################################################################################################################


#########################################################################################################################
## Read in spatial data once, rather than in each script
aus           = readRDS("./data/base/CONTEXTUAL/aus_states.rds")
LAND          = readRDS("./data/base/CONTEXTUAL/LAND_world.rds")
areal_unit    = readRDS("./data/base/CONTEXTUAL/SUA.rds")
areal_unit    = areal_unit[order(areal_unit$SUA_NAME11),]
Koppen        = readRDS('data/base/CONTEXTUAL/Koppen_1975.rds')
Koppen_aus    = readRDS('data/base/CONTEXTUAL/KOPPEN_AUS.rds')
Koppen_zones  = unique(readOGR('data/base/CONTEXTUAL/WC05_1975H_Koppen_Shapefile/WC05_1975H_Koppen_Kriticos_2012.shp')@data[, 1:2])
Koppen_1975   = raster('data/Koppen_1000m_Mollweide54009.tif')
AUS_RAIN      = readRDS('data/base/CONTEXTUAL/BOM/BOM_RAIN_AGG.rds')

SUA           = readRDS("./data/base/CONTEXTUAL/IN_SUA_AUS.rds")
LGA           = readRDS("./data/base/CONTEXTUAL/LGA.rds")
AUS           = readRDS("./data/base/CONTEXTUAL/aus_states.rds")
ALL.SUA.POP   = read.csv("./data/base/CONTEXTUAL/ABS_SUA_POP.csv", stringsAsFactors = FALSE)


##
template.raster = raster("./data/template_hasData.tif")
template.cells  = readRDS("./data/hasData_cells.rds")
load("./data/base/CONTEXTUAL/urbanareas.rda")


# Koppen_aus = readOGR("H:/green_cities_sdm/data/base/CONTEXTUAL/KOPPEN_AUS.shp",
#                    layer = "KOPPEN_AUS")
# saveRDS(Koppen_aus, file = paste("./data/base/CONTEXTUAL/KOPPEN_AUS.rds"))


## Set definitions :: best to minimise the number of projection used in this project
CRS.MOL      <- CRS('+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs')
CRS.MOL.SDM  <- CRS('+init=ESRI:54009 +proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs84=0,0,0')
CRS.WGS.84   <- CRS("+init=epsg:4326")
CRS.AUS.ALB  <- CRS("+init=EPSG:3577")
ALB.CONICAL  <- CRS('+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')





#########################################################################################################################
## 2). READ IN SPECIES LISTS
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
MOD_2               = read.csv("./data/base/HIA_LIST/HIA/MOD2_LIST.csv",                         stringsAsFactors = FALSE)
MOD.2.3             = read.csv("./data/base/HIA_LIST/HIA/MODULE_2_3.csv",                        stringsAsFactors = FALSE)
RISK.LIST           = read.csv("./data/base/HIA_LIST/HIA/RISK_LIST.csv",                         stringsAsFactors = FALSE)
RISK.BINOMIAL.CLEAN = read.csv("./data/base/HIA_LIST/HIA/RISK_BINOMIAL_DF.csv",                  stringsAsFactors = FALSE)
MAXENT.RATING       = read.csv("./output/maxent/MAXENT_RATING_26_2018.csv",                      stringsAsFactors = FALSE)
TI.LIST             = read.csv("./data/base/HIA_LIST/COMBO/ALE_TREE_SPP_LIST.csv",               stringsAsFactors = FALSE)
TI.XY               = read.csv("./data/base/HIA_LIST/COMBO/ALE_TREE_SPP_XY.csv",                 stringsAsFactors = FALSE)
SPP.BIAS            = read.csv("./output/maxent/SPP_BOUNDARY_BIAS.csv", stringsAsFactors = FALSE)


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


## The list of species with checked maxent maps
SPP.BIAS      = merge(SPP.BIAS, COMBO.NICHE.CONTEXT[c("searchTaxon", "Origin")])
SPP.BIAS      = subset(SPP.BIAS, AUS_BOUND_BIAS == "TRUE" & Origin == "Native")
SPP.BIAS      = subset(SPP.BIAS)$searchTaxon



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
camp_spp = gsub(" ", "_", camp.spp) 




#########################################################################################################################
## 3). CLEAN THE HIA LIST
#########################################################################################################################


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
## 4). CREATE LIST OF BINOMIALS : might not need this now
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
length(spp)   


########################################################################################################################
## Try using taxonlookup to check the taxonomy
HIA.SPP.LOOKUP = lookup_table(HIA.SPP[["Binomial"]], by_species = TRUE) ## convert rows to column and merge
HIA.SPP.LOOKUP = setDT(HIA.SPP.LOOKUP , keep.rownames = TRUE)[]
HIA.SPP.LOOKUP = dplyr::rename(HIA.SPP.LOOKUP, Binomial = rn)
head(HIA.SPP.LOOKUP) ## Can merge on the bilogical data here...
## write.csv(HIA.SPP.LOOKUP, "./data/base/HIA_LIST/HIA/HIA_BINOMIAL_LOOKUP.csv", row.names = FALSE)





#########################################################################################################################
## 7). MANUSCRIPT LIST
#########################################################################################################################


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


## Now how many trees do we have? 68 modelled, another 86 with ok data.
## This will increase with Alessandro's data. Re-create niche file, and re-run the criteria.
## Let's get another 40 exotic trees with good models
checked.trees = subset(MODEL.CHECK, Total.growers >= 25 & CHECK_MAP <= 2 & Plant.type == "Tree")$searchTaxon ## 67
subset(MODEL.CHECK, Total.growers >= 25 & CHECK_MAP <= 2 & Plant.type == "Tree" & Origin == "Exotic")$searchTaxon


#########################################################################################################################
## What is the intersection betwen the evergreen list and the Tree inventories?
TREE.HIA.SPP = intersect(subset(TI.LIST, Plantings > 50)$searchTaxon, CLEAN.SPP$Binomial) 
TREE.200.SPP = intersect(subset(TI.LIST, Plantings > 500)$searchTaxon, CLEAN.SPP$Binomial) 
##intersect(head(TI.LIST, 600)$searchTaxon, CLEAN.SPP$Binomial)
summary(TI.LIST[TI.LIST$searchTaxon %in% TREE.HIA.SPP, ])


## Create a list of synonyms using the TPL search
synonyms     = c("Fraxinus angustifolia", "Ficus microcarpa", "Callistemon viminalis", "Schinus molle", 
                 "Alnus acuminata",       "Schinus areira",   "Ficus microcarpa",      "Alnus jorullensis", 
                 "Quercus coccinea",      "Libidibia ferrea", "Fraxinus oxycarpa",     "Prunus blireiana",
                 "Cupressus leylandii",   "Elaeocarpus angustifolius",                 "Eucalyptus globulus",
                 "Corymbia ficifolia",    "Eucalyptus bicolor",                        "Ulmus minor")


## Get the intersection of the synonyms
TREE.HIA.SPP = sort(unique(c(TREE.HIA.SPP, synonyms)))
TREE_HIA_SPP = gsub(" ", "_", TREE.HIA.SPP)
TREE_HIA_SPP = gsub("_trees_sign_.shp", "", TREE.HIA.SPP)


## Replace the weird shapefile name for hobart
TI.XY$INVENTORY = gsub("_trees_sign_.shp", "", TI.XY$INVENTORY)
head(TI.XY)
unique(TI.XY$INVENTORY)


## Only process biased species from the target list, and remove those which are not really biased...........
SPP.BIAS    = intersect(SPP.BIAS, TREE.HIA.SPP)
SPP.BIAS    = setdiff(SPP.BIAS, c("Banksia integrifolia", "Brachychiton acerifolius", "Callistemon viminalis",
                                  "Corymbia maculata", "Elaeocarpus grandis", "Eucalyptus camaldulensis",
                                  "Flindersia schottiana", "Leptospermum petersonii", "Lophostemon suaveolens",
                                  "Stenocarpus sinuatus", "Syzygium luehmannii", "Acacia pendula"))
SPP_BIAS      = gsub(" ", "_", SPP.BIAS)


## Reverse the model list
TREE.HIA.REV = sort(TREE.HIA.SPP, decreasing = TRUE) 
TREE.HIA.REV = sort(TREE.HIA.SPP, decreasing = TRUE) 


## Create a table of the species modelled - join this on in step 8
EVERGREEN = CLEAN.SPP[c("Binomial", "Number.of.growers", "Number.of.States", "Origin")]
names(EVERGREEN )[names(EVERGREEN) == 'Binomial'] <- 'searchTaxon'
TREE.EVERGREEN = merge(EVERGREEN, TI.LIST)
TREE.EVERGREEN = TREE.EVERGREEN[TREE.EVERGREEN$searchTaxon %in% TREE.HIA.SPP, ]
TREE.EVERGREEN = TREE.EVERGREEN [with(TREE.EVERGREEN, rev(order(Plantings))), ]
TREE.EVERGREEN = TREE.EVERGREEN[c("searchTaxon", "Origin", "Plantings", "Number.of.States", "Number.of.growers")]
TREE.EVERGREEN = TREE.EVERGREEN[!duplicated(TREE.EVERGREEN[,c('searchTaxon')]),]


## What does the dataset look like?
TREE.EVERGREEN$Origin = trimws(TREE.EVERGREEN$Origin)
dim(TREE.EVERGREEN)
head(TREE.EVERGREEN)
round(with(TREE.EVERGREEN, table(Origin)/sum(table(Origin))*100), 1)


## How many trees have already been modelled well
length(intersect(TREE.EVERGREEN$searchTaxon, checked.trees))
write.csv(TREE.EVERGREEN, "./data/base/HIA_LIST/COMBO/TREE_EVERGREEN.csv", row.names = FALSE)







#########################################################################################################################
## 8). CHECK TAXONOMY AND ORIGIN OF THE CHOSEN SPECIES
#########################################################################################################################


# ## Run taxonstand for the manuscript tree species
# TREE.TAXO <- TPL(unique(TREE.EVERGREEN$searchTaxon), infra = TRUE,
#                  corr = TRUE, repeats = 100)  ## to stop it timing out...
# 
# 
# ## Most of the species are resolved 
# sort(names(TREE.TAXO))
# table(TREE.TAXO$New.Taxonomic.status)


## Unique(TREE.EVERGREEN$New.Taxonomic.status)
#saveRDS(TREE.EVERGREEN, file = paste("./data/base/HIA_LIST/COMBO/TREE_LIST.rds"))
#SUA.SPP             = read.csv("./data/base/HIA_LIST/COMBO/SUA_SPP_TPL.csv",                     stringsAsFactors = FALSE)
SUA.SPP             = read.csv("./data/base/HIA_LIST/COMBO/GBIF_SPP.csv",                     stringsAsFactors = FALSE)
SUA.SPP             = sort(unique(SUA.SPP$GBIF_SPP))
SUA_SPP             = gsub(" ", "_", SUA.SPP)
setdiff(TREE.HIA.SPP, SUA.SPP)


# TPL.SUA <- TPL(unique(SUA.SPP), infra = TRUE, corr = TRUE, repeats = 100)
# View(TPL.SUA[c("Taxon", "Taxonomic.status", "New.Taxonomic.status", "New.Genus", "New.Species")])
# saveRDS(TPL.SUA, file = paste("./data/base/HIA_LIST/GBIF/TPL_SUA.rds"))
#write.csv(TPL.SUA, "./data/base/HIA_LIST/GBIF/TPL_SUA.csv", row.names = FALSE)


## Note that TPL is not always right for Australian species :: Eucalyptus largiflorens is the accepted Australian name
## bicolor is not recognised by the ALA
TPL.SUA = readRDS("./data/base/HIA_LIST/GBIF/TPL_SUA.rds")
TPL.SPP = paste(TPL.SUA$New.Genus, TPL.SUA$New.Species, sep = " ")
TPL_SPP = gsub(" ", "_", TPL.SPP)
setdiff(TREE.HIA.SPP, SUA.SPP)


# test.gbif = c("Elaeocarpus grandis", "Melaleuca viminalis",       "Melaleuca pallida", 
#               "Platanus hybrida",    "Triadica sebifera",         "Eugenia biflora", 
#               "Cupressus leylandii", "Elaeocarpus angustifolius", "Eucalyptus globulus", 
#               "Corymbia ficifolia",  "Fraxinus angustifolia",     "Callistemon viminalis", 
#               "Prunus blireiana",    "Talipariti tiliaceum",      "Quercus coccinea", 
#               "Quercus palustris")



#########################################################################################################################
############################################  END OF HIA LIST CODE ###################################################### 
#########################################################################################################################