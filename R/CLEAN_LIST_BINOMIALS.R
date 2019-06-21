#########################################################################################################################
########################################  'CLEAN' SPECIES LIST ########################################################## 
#########################################################################################################################


#########################################################################################################################
## 1). READ IN DRAFT HIA LIST
#########################################################################################################################

## The aim here is to take the raw list of plants grown in Aus supplied by Matt Plumber and Anthony Manea (the 'clean tab'), and
## then clean the list as best as possible in R to use the species binomial as the unit of downloading and analysis.

## All the cleaining methods will throw up some anomalies, which need to be tracked, and checked with the team for how
## each case is treated (see outstanding tasks at the bottom)

## PAUL:
## Merge the species native in the SUA with the clean records from ALA.


## Ultimately, we need a table of all the species modelled, with the data for 


#########################################################################################################################
## This list derives from all species and varieties sold anywhere in Australia in the last 5 years. Anthony Manea cleaned 
## up the data and cross-linked to growth form and exotic/native status and derived a list of ~1300 species that are the 
## Most commonly sold, covering the right ratio of growth forms, regional representation and native/exotic
CLEAN.list = read.csv("./data/base/HIA_LIST/HIA/HIA.CLEAN.csv",                               stringsAsFactors = FALSE)
HIA.list   = read.csv("./data/base/HIA_LIST/HIA/GREEN_CITIES_DRAFT_LIST_2709_2017.csv",       stringsAsFactors = FALSE)



#########################################################################################################################
## Dim
dim(CLEAN.list);dim(HIA.list)


## Create a list of the raw HIA list, but removing the weird characters...
## Can use "TRIM(CELL)" in excel, OR:
## trim <- function (x) gsub("^\\s+|\\s+$", "", x) ##  CLEAN.list$Species <- trim(CLEAN.list$Species)
CLEAN.SPECIES = gsub("  ",     " ", CLEAN.list$Species)
CLEAN.SPECIES = gsub(" $",     "",  CLEAN.list$Species, perl = TRUE)
CLEAN.SPECIES = gsub("    $",  "",  CLEAN.list$Species, perl = TRUE)
length(CLEAN.SPECIES)





#########################################################################################################################
## 2). MERGE ON THE TOP 200 MOST POPULAR SPECIES
#########################################################################################################################


#########################################################################################################################
## Merge the ~13000 with the top 200
CLEAN.list$Binomial <- sub('(^\\S+ \\S+).*', '\\1', CLEAN.list$Species) # \\s = white space; \\S = not white space





#########################################################################################################################
## 3). DRAFT 'CLEAN' LIST PREP: REMOVE SUBSPECIES, VARIETIES, ETC
#########################################################################################################################


## First, remove taxa with "spp", etc: see setdiff at the end for the ones which are removed
DRAFT.HIA.TAXA         = CLEAN.list
DRAFT.HIA.TAXA         = DRAFT.HIA.TAXA[DRAFT.HIA.TAXA$Species %in% DRAFT.HIA.TAXA$Species[!grepl("spp.", DRAFT.HIA.TAXA$Species)], ]
dim(DRAFT.HIA.TAXA)    ## 8836 species after we cut out the "spp."


## Remove weird characters...
DRAFT.HIA.TAXA$Species = gsub(" x",     "",  DRAFT.HIA.TAXA$Species)
DRAFT.HIA.TAXA$Species = gsub("NA",     "",  DRAFT.HIA.TAXA$Species)
DRAFT.HIA.TAXA$Species = gsub("  ",     " ", DRAFT.HIA.TAXA$Species)
DRAFT.HIA.TAXA$Species = gsub(" $",     "",  DRAFT.HIA.TAXA$Species, perl = TRUE)
DRAFT.HIA.TAXA$Species = gsub("    $",  "",  DRAFT.HIA.TAXA$Species, perl = TRUE)


#########################################################################################################################
## Now check the taxonomy
CLEAN.UNIQUE = as.character(unique(DRAFT.HIA.TAXA$Species))
CLEAN.LOOKUP = lookup_table(CLEAN.UNIQUE, by_species = TRUE, missing_action = "NA") ## convert rows to column and merge
CLEAN.LOOKUP = setDT(CLEAN.LOOKUP, keep.rownames = TRUE)[]
CLEAN.LOOKUP = dplyr::rename(CLEAN.LOOKUP, Binomial = rn)
head(CLEAN.LOOKUP);dim(CLEAN.LOOKUP)


## Just get the NA rows
CLEAN_TAXO_ERRORS <- CLEAN.LOOKUP[rowSums(is.na(CLEAN.LOOKUP)) > 0,]
#View(CLEAN_TAXO_ERRORS)  ## lots of taxonomy problems here - but I've fixed most of the problems in the popular species


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
  select(Species, Binomial, No.of.Varieties, Plant.type:WA, Number.of.growers, Number.of.States, Origin)

HIA.VARIETY %>% 
  filter(Binomial==Species)

dim(HIA.VARIETY);head(HIA.VARIETY)




#########################################################################################################################
## 4). CREATE FINAL UNIQUE LIST FOR GBIF SEARCHING
#########################################################################################################################


#######################################################################################################################
## Which taxa have the second word capitalised, and are not captured by eliminating the "spp"?
grep('^\\S+ [A-Z]', HIA.VARIETY$Species, val = TRUE)


## Now for GBIF, just get the unique species...
CLEAN.SPP = HIA.VARIETY[!duplicated(HIA.VARIETY["Binomial"]),]
CLEAN.SPP = dplyr::rename(CLEAN.SPP, HIA.Taxa = Species)
dim(CLEAN.SPP);head(CLEAN.SPP)

#########################################################################################################################
## Now sum the number of growers across multiple varieties
## Just get the count for each unique binomial
n.clean <- tapply(CLEAN.SPP$Number.of.growers, CLEAN.SPP$Binomial, sum, na.rm = TRUE)

CLEAN.SPP$Number.of.growers.total <- n.clean[CLEAN.SPP$Binomial]
TOT.GROW = CLEAN.SPP[c("Binomial", "Number.of.growers.total")]
names(TOT.GROW) = c("searchTaxon", "Total.growers")
TOT.GROW        = TOT.GROW[!duplicated(TOT.GROW[,c('searchTaxon')]),]
CLEAN.GROW      = join(CLEAN.SPP, TOT.GROW)



## Reorder by species
CLEAN.SPP = CLEAN.SPP[with(CLEAN.SPP, order(HIA.Taxa)), ] 
#write.csv(CLEAN.SPP , "./data/base/HIA_LIST/HIA/CLEAN_HIA_BIONOMIAL.csv", row.names = FALSE) 
#View(CLEAN.SPP)


#######################################################################################################################
## Now create list of HIA species. 610 unique species, mius the corrections, etc. 
spp.clean            = unique(as.character(CLEAN.SPP$Binomial))
length(spp.clean)   


########################################################################################################################
## Try using taxonlookup to check the taxonomy
CLEAN.LOOKUP = lookup_table(CLEAN.SPP[["Binomial"]], by_species = TRUE, missing_action = "NA") ## convert rows to column and merge
CLEAN.LOOKUP = setDT(CLEAN.LOOKUP , keep.rownames = TRUE)[]
CLEAN.LOOKUP = dplyr::rename(CLEAN.LOOKUP, Binomial = rn)
head(CLEAN.LOOKUP) ## Can merge on the bilogical data here...
#View(CLEAN.LOOKUP)
## setdiff(spp, CLEAN.SPP.LOOKUP$Binomial)
## write.csv(CLEAN.SPP.LOOKUP, "./data/base/HIA_LIST/HIA/HIA_BINOMIAL_LOOKUP.csv", row.names = FALSE)





#########################################################################################################################
## LIST EXCEPTIONS:
#########################################################################################################################


## Record each list: Raw top 25 (1135), Varieties (948), Binomials (610) 
## Check exceptions with team
length(unique(CLEAN.list$Species))    ## All 'things' being grown (13736)
length(unique(HIA.VARIETY$Species))   ## Varieties  (8822), excluding "spp.", eg Philodendron spp. Congo, Nandina domestica Moon Bay
length(unique(CLEAN.SPP$Binomial))    ## Binomials  (4303), keep Michelia yunnanensis Scented Pearl, exclude Spathiphyllum spp. Assorted


## Record the "spp." weirdos
EXCLUDED.SPP         = setdiff(unique(CLEAN.SPECIES), unique(HIA.VARIETY$Species))
EXCLUDED.VARIETIES   = setdiff(unique(HIA.VARIETY$Species), unique(CLEAN.SPP$HIA.Taxa))


#########################################################################################################################
## What is the difference between the growers list and the 'clean'HIA list, without excluding the varieties, etc........?
source('./R/HIA_GROWERS_MATCHING.R')
length(setdiff(CLEAN.SPECIES, GROW.SPECIES))
length(intersect(CLEAN.SPECIES, GROW.SPECIES))



#########################################################################################################################
## What is the difference between the growers list and the 'clean'HIA list, excluding the varieties, etc.?
length(setdiff(spp.grow, spp.clean))
length(intersect(spp.grow, spp.clean))



## How many species need the data to be filled in?
round(with(CLEAN.SPP,   table(Origin)/sum(table(Origin))*100), 1)
round(with(CLEAN.SPP,   table(Plant.type)/sum(table(Plant.type))*100), 1)




#########################################################################################################################
############################################  END OF HIA LIST CODE ###################################################### 
#########################################################################################################################