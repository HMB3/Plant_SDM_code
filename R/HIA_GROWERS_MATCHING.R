#########################################################################################################################
########################################  CREATE SPECIES LIST ########################################################### 
#########################################################################################################################


#########################################################################################################################
## 1). READ IN DRAFT HIA LIST AND CLEAN
#########################################################################################################################

## The aim here is to take the raw list of plants with 25 or more growers supplied by Matt Plumber and Anthony Manea, and
## then clean the list as best as possible in R to use the species binomial as the unit of downloading and analysis.

## All the cleaining methods will throw up some anomalies, which need to be tracked, and checked with the team for how
## each case is treated (see outstanding tasks at the bottom)


#########################################################################################################################
## This list derives from all species and varieties sold anywhere in Australia in the last 5 years. Anthony Manea cleaned 
## Up the data and cross-linked to growth form and exotic/native status and derived a list of ~1000 species that are the 
## Most commonly sold, covering the right ratio of growth forms, regional representation and native/exotic
CLEAN.list = read.csv("./data/base/HIA_LIST/HIA/HIA.CLEAN.csv",                         stringsAsFactors = FALSE)
GROW.list  = read.csv("./data/base/HIA_LIST/HIA/database_aus_sp_growing.csv",           stringsAsFactors = FALSE) 
top.200    = read.csv("./data/base/HIA_LIST/HIA/HIA_TOP_200_1309_2017.csv",             stringsAsFactors = FALSE)


#########################################################################################################################
## Dim
dim(GROW.list)
GROW.list = dplyr::rename(GROW.list, Species = scientific_name)


## Create a list of the raw HIA list, but removing the weird characters...
## Just use "TRIM(CELL)" in excel
## trim <- function (x) gsub("^\\s+|\\s+$", "", x) ##  GROW.list$Species <- trim(GROW.list$Species)
RAW.HIA.SPP = gsub("  ",     " ", GROW.list$Species)
RAW.HIA.SPP = gsub(" $",     "",  GROW.list$Species, perl = TRUE)
RAW.HIA.SPP = gsub("    $",  "",  GROW.list$Species, perl = TRUE)
length(RAW.HIA.SPP)


#########################################################################################################################
## Taxon check on the original list
#intersect(GROW.list$Species, GROWING$scientific_name)
CLEAN.UNIQUE = as.character(unique(RAW.HIA.SPP))
FIRST.LOOKUP = lookup_table(CLEAN.UNIQUE, by_species = TRUE, missing_action = "NA") ## convert rows to column and merge
FIRST.LOOKUP = setDT(FIRST.LOOKUP, keep.rownames = TRUE)[]
FIRST.LOOKUP = dplyr::rename(FIRST.LOOKUP, Binomial = rn)
head(FIRST.LOOKUP)
View(FIRST.LOOKUP)
setdiff(FIRST.LOOKUP$Binomial, CLEAN.UNIQUE)


## also, add the "Top 200" species in here
spp.200          = top.200[c("Species", "t200_MATCH_25")]
spp.200$Species  <- sub('(^\\S+ \\S+).*', '\\1', spp.200$Species) # \\s = white space; \\S = not white space

spp.200$Species  = gsub("  ",     " ", spp.200$Species)
spp.200$Species  = gsub(" $",     "",  spp.200$Species, perl = TRUE)
spp.200$Species  = gsub("    $",  "",  spp.200$Species, perl = TRUE)
spp.200$Top_200  = "TRUE"
spp.200          = dplyr::rename(spp.200, Binomial = Species)


#########################################################################################################################
## Merge the ~1000 with the top 200
## This merge won't get the ones that match to binomial
GROW.list$Binomial <- sub('(^\\S+ \\S+).*', '\\1', GROW.list$Species) # \\s = white space; \\S = not white space


## Check this reduces the number of Top 200 missing from this list
GROW.list = merge(GROW.list, spp.200, by = "Binomial", all.x = TRUE) 
GROW.list$Top_200[is.na(GROW.list$Top_200)] <- "FALSE"
GROW.list$t200_MATCH_25[is.na(GROW.list$t200_MATCH_25)] <- "TRUE"
#GROW.list$Origin <- gsub(" ",  "", GROW.list$Origin)


## check
str(GROW.list)
head(GROW.list)
unique(GROW.list$Top_200)
unique(GROW.list$t200_MATCH_25)
unique(GROW.list$Origin)
length(unique(GROW.list$Binomial)) ## 1077 unique binomials
names(GROW.list)





#########################################################################################################################
## 2). DRAFT LIST PREP: THIS NEEDS TO CHANGE IN CONSULTATION WITH RACH, PAUL, ETC
#########################################################################################################################


## First, taxa with "spp". Return to these later...
DRAFT.HIA.TAXA         = GROW.list
DRAFT.HIA.TAXA         = DRAFT.HIA.TAXA[DRAFT.HIA.TAXA$Species %in% DRAFT.HIA.TAXA$Species[!grepl("spp.", DRAFT.HIA.TAXA$Species)], ]
DRAFT.HIA.TAXA         = DRAFT.HIA.TAXA[DRAFT.HIA.TAXA$Species %in% DRAFT.HIA.TAXA$Species[!grepl("sp.", DRAFT.HIA.TAXA$Species)], ]
dim(DRAFT.HIA.TAXA)    ## 948 species after we cut out the "spp."


## Remove weird characters...
DRAFT.HIA.TAXA$Species = gsub(" x",     "",  DRAFT.HIA.TAXA$Species)
DRAFT.HIA.TAXA$Species = gsub("NA",     "",  DRAFT.HIA.TAXA$Species)
DRAFT.HIA.TAXA$Species = gsub("  ",     " ",  DRAFT.HIA.TAXA$Species)
DRAFT.HIA.TAXA$Species = gsub(" $",     "",  DRAFT.HIA.TAXA$Species, perl = TRUE)
DRAFT.HIA.TAXA$Species = gsub("    $",  "",  DRAFT.HIA.TAXA$Species, perl = TRUE)


#########################################################################################################################
## Now create a table of how many varieties each species has
# length(unique(GROW.list$Binomial)) 
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
  select(Binomial, state, SUA, LGA, Species, common_name, ref, Top_200)

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
View(HIA.SPP)


#######################################################################################################################
## Now create list of HIA species. 610 unique species, mius the corrections, etc. 
spp.grow = unique(as.character(HIA.SPP$Binomial))
length(spp.grow)   


########################################################################################################################
## Try using taxonlookup to check the taxonomy
HIA.SPP.LOOKUP = lookup_table(HIA.SPP[["Binomial"]], by_species = TRUE, missing_action = "NA") ## convert rows to column and merge
HIA.SPP.LOOKUP = setDT(HIA.SPP.LOOKUP , keep.rownames = TRUE)[]
HIA.SPP.LOOKUP = dplyr::rename(HIA.SPP.LOOKUP, Binomial = rn)
head(HIA.SPP.LOOKUP) ## Can merge on the bilogical data here...
View(HIA.SPP.LOOKUP)
## setdiff(spp.grow, HIA.SPP.LOOKUP$Binomial)
## write.csv(HIA.SPP.LOOKUP, "./data/base/HIA_LIST/HIA/HIA_BINOMIAL_LOOKUP.csv", row.names = FALSE)





#########################################################################################################################
## LIST EXCEPTIONS:
#########################################################################################################################


## Record each list: Raw top 25 (1135), Varieties (948), Binomials (610) 
## Check exceptions with Paul, Linda and Rach
length(unique(GROW.list$Species))   ## Raw top 25 (13736)
length(unique(HIA.VARIETY$Species))  ## Varieties  (8822), excluding "spp.", eg Philodendron spp. Congo, Nandina domestica Moon Bay
length(unique(HIA.SPP$Binomial))     ## Binomials (4306), keep Michelia yunnanensis Scented Pearl, exclude Spathiphyllum spp. Assorted


## record the "spp." weirdos
EXCLUDED.SPP         = setdiff(unique(RAW.HIA.SPP), unique(HIA.VARIETY$Species))
EXCLUDED.VARIETIES   = setdiff(unique(HIA.VARIETY$Species), unique(HIA.SPP$HIA.Taxa))   ## Here is the list that spots the exceptions!!!!!!!


## Remaining anomalies:

## EG: Rhaphiolepis indica has growers for the spp and each variety, should we add them together?
## Magnolia grandiflora has 8 varieties which are being missed by the current code...








#########################################################################################################################
############################################  END OF HIA LIST CODE ###################################################### 
#########################################################################################################################