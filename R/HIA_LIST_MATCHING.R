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
HIA.list   = read.csv("./data/base/HIA_LIST/HIA/GREEN_CITIES_DRAFT_LIST_0809_2017.csv", stringsAsFactors = FALSE)
top.200    = read.csv("./data/base/HIA_LIST/HIA/HIA_TOP_200.csv",                       stringsAsFactors = FALSE)
renee.taxa = read.csv("./data/base/HIA_LIST/HIA/RENEE_TAXA.csv",                        stringsAsFactors = FALSE)


## have a look
dim(HIA.list)
dim(renee.taxa)
str(HIA.list)
head(HIA.list)


## also, add the "Top 200" species in here
spp.200          = top.200[c("Species", "t200_MATCH_25")]
spp.200$Species  <- sub('(^\\S+ \\S+).*', '\\1', spp.200$Species) # \\s = white space; \\S = not white space

spp.200$Species  = gsub("  ",     " ", spp.200$Species)
spp.200$Species  = gsub(" $",     "",  spp.200$Species, perl = TRUE)
spp.200$Species  = gsub("    $",  "",  spp.200$Species, perl = TRUE)
spp.200$Top_200  = "TRUE"
spp.200          = rename(spp.200, Binomial = Species)


## Just get the species renee selected that are not on the top 1000 or 200
renee.list       = renee.taxa[c("Species", "Growth_Form")]


#########################################################################################################################
## Merge the ~1000 with the top 200
## This merge won't get the ones that match to binomial
HIA.list$Binomial <- sub('(^\\S+ \\S+).*', '\\1', HIA.list$Species) # \\s = white space; \\S = not white space


## Check this reduces the number of Top 200 missing from this list
HIA.list = merge(HIA.list, spp.200, by = "Binomial", all.x = TRUE) 
HIA.list$Top_200[is.na(HIA.list$Top_200)] <- "FALSE"
HIA.list$Top_200[is.na(HIA.list$t200_MATCH_25)] <- "TRUE"
HIA.list$Origin <- gsub(" ",  "", HIA.list$Origin)


## check
str(HIA.list)
head(HIA.list)
unique(HIA.list$Top_200)
unique(HIA.list$Origin)
length(unique(HIA.list$Binomial)) ## 660 unique binomials


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
DRAFT.HIA.TAXA$Species = gsub("  ",     " ",  DRAFT.HIA.TAXA$Species)
DRAFT.HIA.TAXA$Species = gsub(" $",     "",  DRAFT.HIA.TAXA$Species, perl = TRUE)
DRAFT.HIA.TAXA$Species = gsub("    $",  "",  DRAFT.HIA.TAXA$Species, perl = TRUE)


#########################################################################################################################
## Now create a table of how many varieties each species has
# length(unique(HIA.list$Binomial)) 
# length(unique(sub('(^\\S+ \\S+).*', '\\1', DRAFT.HIA.TAXA$Species)))
# setdiff(unique(HIA.list$Binomial), unique(sub('(^\\S+ \\S+).*', '\\1', DRAFT.HIA.TAXA$Species)))


## Create another binomial column
DRAFT.HIA.TAXA$Binomial <- sub('(^\\S+ \\S+).*', '\\1', DRAFT.HIA.TAXA$Species) # \\s = white space; \\S = not white space


## And count how many varieties each taxa has? 
HIA.VARIETY <- 
  DRAFT.HIA.TAXA$Binomial[DRAFT.HIA.TAXA$Binomial != DRAFT.HIA.TAXA$Species] %>% 
  table %>% 
  as.data.frame %>% 
  setNames(c('Binomial', 'No.of.Varieties')) %>% 
  full_join(DRAFT.HIA.TAXA) %>% 
  select(Species, Binomial, No.of.Varieties, Plant.type:WA, Number.of.growers, Number.of.States, Origin, Top_200)

HIA.VARIETY %>% 
  filter(Binomial==Species)


#######################################################################################################################
## Which taxa have the second word capitalised, and are not captured by eliminating the "spp"?
grep('^\\S+ [A-Z]', HIA.VARIETY$Species, val = TRUE)


## Now for GBIF, just get the unique species...
HIA.SPP = HIA.VARIETY[!duplicated(HIA.VARIETY["Binomial"]),]
HIA.SPP = rename(HIA.SPP, HIA.Taxa = Species)


## Reorder by species
HIA.SPP = HIA.SPP[with(HIA.SPP, order(Species)), ] 
View(HIA.SPP)


## This still leaves single genera?
# vapply(lapply(strsplit(HIA.SPP$binomial, " "),
#               unique), paste, character(1L), collapse = " ")


#######################################################################################################################
## Now create list of HIA species. 496 unique species, mius the corrections, etc. 
spp            = unique(as.character(HIA.SPP$Binomial))
spp.renee      = unique(as.character(renee.list$Species)) ## 
length(spp)   ## why 660? later check on what happens with the different queries


########################################################################################################################
## Try using taxonlookup to check the taxonomy
HIA.SPP.LOOKUP = lookup_table(HIA.SPP[["Binomial"]], by_species = TRUE) ## convert rows to column and merge
HIA.SPP.LOOKUP = setDT(HIA.SPP.LOOKUP , keep.rownames = TRUE)[]
HIA.SPP.LOOKUP = rename(HIA.SPP.LOOKUP, Binomial = rn)
head(HIA.SPP.LOOKUP) ## Can merge on the bilogical data here...





#########################################################################################################################
## 3). NOW MATCH THE LISTS
#########################################################################################################################


#########################################################################################################################
## Five spp have more than 200K records...
load("./data/base/HIA_LIST/GBIF/GBIF_NICHE_CONTEXT.RData")
load("./data/base/HIA_LIST/GBIF/skipped_species.RData")
skipped.species.df[ which(skipped.species.df$Reason_skipped == "Number of records > 200,000"), ]
View(skipped.species.df)
View(GBIF.NICHE.CONTEXT)


## Get the difference between the original list and the processed list
missed.HIA.processed = setdiff(HIA.SPP$Binomial, GBIF.NICHE.CONTEXT$searchTaxon)        ## return elements beloning to HIA only
missed.processed.HIA = setdiff(GBIF.NICHE.CONTEXT$searchTaxon, HIA.SPP$Binomial)        ## beloning to processed only


## Plus the difference between the top 200 and the processed list
missed.t200.processed = setdiff(spp.200$Binomial, subset(HIA.SPP, Top_200 == "TRUE")[["Binomial"]])
missed.t200.processed = setdiff(spp.200$Binomial, GBIF.NICHE.CONTEXT$searchTaxon) 


## Need 660 rows in the processed data with contextual data
length(missed.HIA.processed) + length(missed.t200.processed)                            ## 133 missing taxa... 660 -559
missing.taxa = unique(c(missed.HIA.processed, missed.t200.processed))
missing.taxa = gsub("    $",  "",  missing.taxa, perl = TRUE)
missing.taxa


#########################################################################################################################
## Get the unique list from the final data
length(GBIF.NICHE.CONTEXT[ which(GBIF.NICHE.CONTEXT$Number.of.growers >= 25), ][["searchTaxon"]])
GBIF.NICHE.CONTEXT.UNIQUE = unique(GBIF.NICHE.CONTEXT[ which(GBIF.NICHE.CONTEXT$Number.of.growers >= 25), ][["searchTaxon"]])
setdiff(GBIF.NICHE.CONTEXT.UNIQUE, missing.taxa)


#########################################################################################################################
## List tasks:
## Record each list: Raw top 25 (1135), Varieties (948), Binomials (610) 
length(unique(HIA.list$Species))     ## Raw top 25 (1135)
length(unique(HIA.VARIETY$Species))  ## Varieties  (948), excluding "spp.", eg Philodendron spp. Congo, Nandina domestica Moon Bay
 length(unique(HIA.SPP$Binomial))     ## Binomials (610), keep Michelia yunnanensis Scented Pearl, exclude Spathiphyllum spp. Assorted



## All the unique species that are still on the list, especially the top 200
## Find a way to store the number of varieties for each species
## Keep a list of the exceptions

#########################################################################################################################
############################################  END OF HIA LIST CODE ###################################################### 
#########################################################################################################################