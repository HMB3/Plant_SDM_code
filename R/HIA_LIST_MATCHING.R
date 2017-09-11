#########################################################################################################################
########################################  CREATE SPECIES LIST ########################################################### 
#########################################################################################################################


#########################################################################################################################
## 1). READ IN DRAFT HIA LIST AND CLEAN
#########################################################################################################################


## This list derives from all species and varieties sold anywhere in Australia in the last 5 years. Anthony Manea cleaned 
## Up the data and cross-linked to growth form and exotic/native status and derived a list of ~1000 species that are the 
## Most commonly sold, covering the right ratio of growth forms, regional representation and native/exotic
HIA.list   = read.csv("./data/base/HIA_LIST/HIA/GREEN_CITIES_DRAFT_LIST_0809_2017.csv", stringsAsFactors = FALSE)
top.200    = read.csv("./data/base/HIA_LIST/HIA/HIA_TOP_200.csv",                       stringsAsFactors = FALSE)
renee.taxa = read.csv("./data/base/HIA_LIST/RENEE/RENEE_TAXA.csv",                      stringsAsFactors = FALSE)


## have a look
dim(HIA.list)
dim(renee.taxa)
str(HIA.list)
head(HIA.list)


## also, add the "Top 200" species in here
spp.200          = top.200[c("Species")]
spp.200$Species  = gsub(" $", "", spp.200$Species, perl = TRUE)
spp.200$Top_200  = "TRUE"


## Just get the species renee selected that are not on the top 1000 or 200
renee.list       = renee.taxa[c("Species", "Growth_Form")]


## Merge the ~1000 with the top 200
HIA.list = merge(HIA.list, spp.200, by = "Species", all.x = TRUE) 
HIA.list$Top_200[is.na(HIA.list$Top_200)] <- "FALSE"
HIA.list$Origin <- gsub(" ",  "", HIA.list$Origin)


## check
str(HIA.list)
head(HIA.list)
unique(HIA.list$Top_200)
unique(HIA.list$Origin)





#########################################################################################################################
## WHAT ARE SOME USEFUL MEASURES OF THE DATASET?
#########################################################################################################################


## ~ HALF the total species are native, ~47% of the top 200
dim(subset(HIA.list, Origin == "Native"))[1]/dim(HIA.list)[1]*100
dim(subset(top.200,  Origin == "Native"))[1]/dim(top.200)[1]*100


#########################################################################################################################
## 2). DRAFT LIST PREP: THIS NEEDS TO CHANGE IN CONSULTATION WITH RACH, PAUL, ETC
#########################################################################################################################


## First, get rid of all the lines with "spp".
DRAFT.HIA.TAXA         = HIA.list
EXCLUDE.SPP            = DRAFT.HIA.TAXA$Species[!grepl("spp.", DRAFT.HIA.TAXA$Species)]
DRAFT.HIA.TAXA         = DRAFT.HIA.TAXA[DRAFT.HIA.TAXA$Species %in% EXCLUDE.SPP, ]
dim(DRAFT.HIA.TAXA)


## For now, remove all the subspecies, varieties etc.
## But check if some of them work on GBIF, e.g. a separate list of varities
DRAFT.HIA.TAXA$Species = gsub(" x",   "",  DRAFT.HIA.TAXA$Species)
DRAFT.HIA.TAXA$Species = gsub("NA",   "",  DRAFT.HIA.TAXA$Species)
DRAFT.HIA.TAXA$Species = gsub(" $",   "",  DRAFT.HIA.TAXA$Species, perl = TRUE)


##

DRAFT.HIA.TAXA$binomial <- sub('(^\\S+ \\S+).*', '\\1', DRAFT.HIA.TAXA$Species) # \\s = white space; \\S = not white space
DRAFT.HIA.TAXA <- 
  DRAFT.HIA.TAXA$binomial[DRAFT.HIA.TAXA$binomial != DRAFT.HIA.TAXA$Species] %>% 
  table %>% 
  as.data.frame %>% 
  setNames(c('binomial', 'n_infraspecific')) %>% 
  full_join(DRAFT.HIA.TAXA) %>% 
  select(Species, binomial, n_infraspecific, Plant.type:WA)

DRAFT.HIA.TAXA %>% 
  filter(binomial==Species)

grep('^\\S+ [A-Z]', DRAFT.HIA.TAXA$Species, val=T)


## Then remove the varieties and subspecies
DRAFT.HIA.TAXA$Species = vapply(lapply(strsplit(DRAFT.HIA.TAXA$Species, " "),
                                 unique), paste, character(1L), collapse = " ")


## Then get just the first two words (again cleaning up the subspecies, and single genera)
DRAFT.HIA.TAXA$Species = vapply(lapply(strsplit(DRAFT.HIA.TAXA$Species, " "), 
                                 string_fun_first_two_words), paste, character(1L), collapse = " ")


## Reorder by species
DRAFT.HIA.TAXA = DRAFT.HIA.TAXA[with(DRAFT.HIA.TAXA, order(Species)), ] 


## check
str(DRAFT.HIA.TAXA)
View(DRAFT.HIA.TAXA)


## in here, try to count how many varieties each species has... 


## Now create list of HIA taxa. 768 unique species, mius the corrections, etc. 
DRAFT.HIA.TAXA = DRAFT.HIA.TAXA[with(DRAFT.HIA.TAXA, order(Species)), ]
spp            = unique(as.character(DRAFT.HIA.TAXA$Species))
spp.renee      = unique(as.character(renee.list$Species)) ## 
str(spp)   ## why 660? later check on what happens with the different queries
head(spp, 50)
tail(spp, 50)


# ## also make a genus list?
# genera = unique(vapply(lapply(strsplit(HIA.list$Species, " "), 
#                               string_fun_first_word), paste, character(1L), collapse = " "))
# 
# 
# ## check
# str(genera)
# head(genera, 50)
# tail(genera, 50)


## also, for Stuarts code EG, I need a df not a list. Get just the rows of HIA.list which have
## unique species names. 
DRAFT.HIA.TAXA = DRAFT.HIA.TAXA[!duplicated(DRAFT.HIA.TAXA$Species), ]
dim(DRAFT.HIA.TAXA)
head(DRAFT.HIA.TAXA)
View(DRAFT.HIA.TAXA)


## get rid of the gunk again
HIA.list$Species = gsub("spp.", "",  HIA.list$Species)
HIA.list$Species = gsub(" x",  "",   HIA.list$Species)
HIA.list$Species = gsub("NA",  "",   HIA.list$Species)
HIA.list$Species = gsub(" $","",     HIA.list$Species, perl = TRUE)


## then remove the varieties and subspecies
DRAFT.HIA.TAXA$Species = vapply(lapply(strsplit(DRAFT.HIA.TAXA$Species, " "),
                                 unique), paste, character(1L), collapse = " ")

## then get just the first two words (again cleaning up the subspecies, and single genera)
DRAFT.HIA.TAXA$Species = vapply(lapply(strsplit(DRAFT.HIA.TAXA$Species, " "), 
                             string_fun_first_two_words), paste, character(1L), collapse = " ")

## remove NA
unique.HIA.Species = DRAFT.HIA.TAXA$Species[!grepl(paste0("NA", collapse = "|"), DRAFT.HIA.TAXA$Species)]
length(unique.HIA.Species)


########################################################################################################################
## Try using taxonlookup to check the taxonomy
DRAFT.TAXA.LOOKUP = lookup_table(DRAFT.HIA.TAXA[["Species"]], by_species = TRUE) ## convert rows to column and merge
DRAFT.TAXA.LOOKUP = setDT(DRAFT.TAXA.LOOKUP , keep.rownames = TRUE)[]
DRAFT.TAXA.LOOKUP = rename(DRAFT.TAXA.LOOKUP, searchTaxon = rn)
head(DRAFT.TAXA.LOOKUP)





#########################################################################################################################
## So some of the character replacements are not working as intended. 
## But not many of those species are on the top 200. 


## Get the difference between the original list and the processed list
missed.HIA.processed = setdiff(HIA.list$Species, GBIF.NICHE.CONTEXT$searchTaxon)        ## return elements beloning to HIA only
missed.processed.HIA = setdiff(GBIF.NICHE.CONTEXT$searchTaxon, unique.HIA.Species)        ## return elements beloning to processed only
missed.spp.processed = setdiff(unique.HIA.Species, GBIF.NICHE.CONTEXT$searchTaxon)    ## return elements beloning to spp only


## Plus the difference between the top 200 and the processed list
missed.t200.processed = setdiff(spp.200$Species, GBIF.NICHE.CONTEXT$searchTaxon)        ## return elements beloning to 200 only
missed.t200.processed = setdiff(spp.200$Species, missed.spp.processed)                  ## return elements beloning to 200 only


## Plus the difference between the setdiff, and the top 200
setdiff.t200          = setdiff(missed.spp.processed, spp.200$Species)                  ## return elements beloning to setdiff only
t200.setdiff          = setdiff(missed.spp.processed, spp.200$Species)                  ## return elements beloning to setdiff only


## Need 660 rows in the processed data with contextual data
length(missed.spp.processed) + length(missed.t200.processed)                            ## 133 missing taxa... 660 -559
missing.taxa = unique(c(missed.spp.processed, missed.t200.processed))


# ## get rid of the gunk again
# missing.taxa = gsub("spp.", "",  missing.taxa)
# missing.taxa = gsub(" x",  "",   missing.taxa)
# missing.taxa = gsub("NA",  "",   missing.taxa)
# missing.taxa = gsub(" $","",     missing.taxa, perl = TRUE)
# 
# 
# ## then get just the first two words (again cleaning up the subspecies, and single genera)
# missing.taxa = vapply(lapply(strsplit(missing.taxa, " "), 
#                              string_fun_first_two_words), paste, character(1L), collapse = " ")
# 
# 
# ## remove NA
# test = missing.taxa[!grepl(paste0("NA", collapse = "|"), missing.taxa)]


## Get the unique list from the final data
length(GBIF.NICHE.CONTEXT[ which(GBIF.NICHE.CONTEXT$Number.of.growers >= 25), ][["searchTaxon"]])
GBIF.NICHE.CONTEXT.UNIQUE = unique(GBIF.NICHE.CONTEXT[ which(GBIF.NICHE.CONTEXT$Number.of.growers >= 25), ][["searchTaxon"]])
setdiff(GBIF.NICHE.CONTEXT.UNIQUE, missing.taxa)





#########################################################################################################################
## So some of the character replacements are not working as intended. 
## But not many of those species are on the top 200. 

## The lists are getting messy! Outstanding tasks:
## Get all the unique species that are still on the list, especially the top 200
## Find a way to store the number of varieties for each species




