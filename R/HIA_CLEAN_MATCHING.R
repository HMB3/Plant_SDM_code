#########################################################################################################################
########################################  CREATE SPECIES LIST ########################################################### 
#########################################################################################################################


#########################################################################################################################
## 1). READ IN DRAFT HIA LIST
#########################################################################################################################

## The aim here is to take the raw list of plants grown in Aus supplied by Matt Plumber and Anthony Manea, and
## then clean the list as best as possible in R to use the species binomial as the unit of downloading and analysis.

## All the cleaining methods will throw up some anomalies, which need to be tracked, and checked with the team for how
## each case is treated (see outstanding tasks at the bottom)


#########################################################################################################################
## This list derives from all species and varieties sold anywhere in Australia in the last 5 years. Anthony Manea cleaned 
## up the data and cross-linked to growth form and exotic/native status and derived a list of ~1300 species that are the 
## Most commonly sold, covering the right ratio of growth forms, regional representation and native/exotic
CLEAN.list = read.csv("./data/base/HIA_LIST/HIA/HIA.CLEAN.csv",                         stringsAsFactors = FALSE)
GROW.list  = read.csv("./data/base/HIA_LIST/HIA/HIA.CLEAN.csv",                         stringsAsFactors = FALSE)
top.200    = read.csv("./data/base/HIA_LIST/HIA/HIA_TOP_200_1309_2017.csv",             stringsAsFactors = FALSE)


#########################################################################################################################
## Dim
dim(CLEAN.list)


## Create a list of the raw HIA list, but removing the weird characters...
## Can use "TRIM(CELL)" in excel, OR:
## trim <- function (x) gsub("^\\s+|\\s+$", "", x) ##  CLEAN.list$Species <- trim(CLEAN.list$Species)
RAW.HIA.SPP = gsub("  ",     " ", CLEAN.list$Species)
RAW.HIA.SPP = gsub(" $",     "",  CLEAN.list$Species, perl = TRUE)
RAW.HIA.SPP = gsub("    $",  "",  CLEAN.list$Species, perl = TRUE)
length(RAW.HIA.SPP)





#########################################################################################################################
## 2). MERGE ON THE TOP 200 MOST POPULAR SPECIES
#########################################################################################################################


#########################################################################################################################
## Clean the top 200 for character errors
spp.200          = top.200[c("Species")]
spp.200$Species  <- sub('(^\\S+ \\S+).*', '\\1', spp.200$Species) # \\s = white space; \\S = not white space

spp.200$Species  = gsub("  ",     " ", spp.200$Species)
spp.200$Species  = gsub(" $",     "",  spp.200$Species, perl = TRUE)
spp.200$Species  = gsub("    $",  "",  spp.200$Species, perl = TRUE)
spp.200$Top_200  = "TRUE"
spp.200          = dplyr::rename(spp.200, Binomial = Species)


## Merge the ~13000 with the top 200
CLEAN.list$Binomial <- sub('(^\\S+ \\S+).*', '\\1', CLEAN.list$Species) # \\s = white space; \\S = not white space


## Merge the 13,000 witth the top 200
CLEAN.list = merge(CLEAN.list, spp.200, by = "Binomial", all.x = TRUE) 
CLEAN.list$Top_200[is.na(CLEAN.list$Top_200)] <- "FALSE"
CLEAN.list$t200_MATCH_25[is.na(CLEAN.list$t200_MATCH_25)] <- "TRUE"
CLEAN.list$Origin <- gsub(" ",  "", CLEAN.list$Origin)





#########################################################################################################################
## 3). DRAFT LIST PREP: REMOVE SUBSPECIES, VARIETIES, ETC
#########################################################################################################################


## First, remove taxa with "spp", etc
DRAFT.HIA.TAXA         = CLEAN.list
DRAFT.HIA.TAXA         = DRAFT.HIA.TAXA[DRAFT.HIA.TAXA$Species %in% DRAFT.HIA.TAXA$Species[!grepl("spp.", DRAFT.HIA.TAXA$Species)], ]
dim(DRAFT.HIA.TAXA)    ## 8836 species after we cut out the "spp."


## Remove weird characters...
DRAFT.HIA.TAXA$Species = gsub(" x",     "",  DRAFT.HIA.TAXA$Species)
DRAFT.HIA.TAXA$Species = gsub("NA",     "",  DRAFT.HIA.TAXA$Species)
DRAFT.HIA.TAXA$Species = gsub("  ",     " ",  DRAFT.HIA.TAXA$Species)
DRAFT.HIA.TAXA$Species = gsub(" $",     "",  DRAFT.HIA.TAXA$Species, perl = TRUE)
DRAFT.HIA.TAXA$Species = gsub("    $",  "",  DRAFT.HIA.TAXA$Species, perl = TRUE)


#########################################################################################################################
## Now check the taxonomy
CLEAN.UNIQUE = as.character(unique(DRAFT.HIA.TAXA$Species))
FIRST.LOOKUP = lookup_table(CLEAN.UNIQUE, by_species = TRUE, missing_action = "NA") ## convert rows to column and merge
FIRST.LOOKUP = setDT(FIRST.LOOKUP, keep.rownames = TRUE)[]
FIRST.LOOKUP = dplyr::rename(FIRST.LOOKUP, Binomial = rn)
head(FIRST.LOOKUP)

## Just get the NA rows
CLEAN_TAXO_ERRORS <- FIRST.LOOKUP[rowSums(is.na(FIRST.LOOKUP)) > 0,]
View(CLEAN_TAXO_ERRORS)  ## lots of taxonomy problems here - but I've fixed most of the problems in the popular species


## Create another binomial column
DRAFT.HIA.TAXA$Binomial <- sub('(^\\S+ \\S+).*', '\\1', DRAFT.HIA.TAXA$Species) # \\s = white space; \\S = not white space


#########################################################################################################################
## And count how many varieties each taxa has...
## This needs to change so that every variety that matches a binomial is added to the Number of varieties. Also, can we 
## take the variety with the highest number of growers? There are 8 different varieties of magnolia, currently I'm not 
## getting the most popular ones (grandiflora).
HIA.VARIETY <- 
  DRAFT.HIA.TAXA$Binomial[DRAFT.HIA.TAXA$Binomial != DRAFT.HIA.TAXA$Species] %>% 
  table %>% 
  as.data.frame %>% 
  setNames(c('Binomial', 'No.of.Varieties')) %>% 
  full_join(DRAFT.HIA.TAXA) %>% 
  select(Species, Binomial, No.of.Varieties, Plant.type:WA, Number.of.growers, Number.of.States, Origin, Top_200)

HIA.VARIETY %>% 
  filter(Binomial==Species)




#########################################################################################################################
## 4). CREATE FINAL UNIQUE LIST FOR GBIF SEARCHING
#########################################################################################################################


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
spp.clean            = unique(as.character(HIA.SPP$Binomial))
length(spp.clean)   


########################################################################################################################
## Try using taxonlookup to check the taxonomy
HIA.SPP.LOOKUP = lookup_table(HIA.SPP[["Binomial"]], by_species = TRUE, missing_action = "NA") ## convert rows to column and merge
HIA.SPP.LOOKUP = setDT(HIA.SPP.LOOKUP , keep.rownames = TRUE)[]
HIA.SPP.LOOKUP = dplyr::rename(HIA.SPP.LOOKUP, Binomial = rn)
head(HIA.SPP.LOOKUP) ## Can merge on the bilogical data here...
View(HIA.SPP.LOOKUP)
## setdiff(spp, HIA.SPP.LOOKUP$Binomial)
## write.csv(HIA.SPP.LOOKUP, "./data/base/HIA_LIST/HIA/HIA_BINOMIAL_LOOKUP.csv", row.names = FALSE)





#########################################################################################################################
## LIST EXCEPTIONS:
#########################################################################################################################


## Record each list: Raw top 25 (1135), Varieties (948), Binomials (610) 
## Check exceptions with Paul, Linda and Rach
length(unique(CLEAN.list$Species))   ## All species being grown (13736)
length(unique(HIA.VARIETY$Species))  ## Varieties  (8822), excluding "spp.", eg Philodendron spp. Congo, Nandina domestica Moon Bay
length(unique(HIA.SPP$Binomial))     ## Binomials  (4306), keep Michelia yunnanensis Scented Pearl, exclude Spathiphyllum spp. Assorted


## record the "spp." weirdos
EXCLUDED.SPP         = setdiff(unique(RAW.HIA.SPP), unique(HIA.VARIETY$Species))
EXCLUDED.VARIETIES   = setdiff(unique(HIA.VARIETY$Species), unique(HIA.SPP$HIA.Taxa))   


## These list includes exceptions


## Remaining anomalies:

## EG: Rhaphiolepis indica has growers for both the species and each variety, should we add them together?
## Magnolia grandiflora has 8 varieties, which are being missed by the current code...will ask john to help fix it








#########################################################################################################################
############################################  END OF HIA LIST CODE ###################################################### 
#########################################################################################################################