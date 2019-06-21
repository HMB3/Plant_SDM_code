#########################################################################################################################
########################################  CLEAN URBAN LIST ############################################################## 
#########################################################################################################################


## This code cleans the urban tree inventory data


## Read in Alessandro's data
ALE.LIST            = read.csv("./data/base/HIA_LIST/COMBO/AleTreeInventory summary_20_7_18.csv", stringsAsFactors = FALSE)
ALE.DF              = read.csv("./data/base/HIA_LIST/COMBO/ALL_trees_WGS84_24.07.18.csv",         stringsAsFactors = FALSE)
CLEAN.NICHE.CONTEXT = readRDS("./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT_APRIL_2018_COORD_CLEAN.rds")


## What are the common exceptions
setdiff(ALE.LIST$searchTaxon, CLEAN.NICHE.CONTEXT$searchTaxon)



#########################################################################################################################
## 1). REMOVE EXCEPTIONS
#########################################################################################################################


#########################################################################################################################
## Trim white space
ALE.DF$Species <- trim(ALE.DF$SPECIES)


#########################################################################################################################
## Replace multiple words at once with a space
to.replace = list("spp" = " ", "spp." = " ", 
                  "Spp" = " ", "Spp." = " ", 
                  "Sp" = " ",  "Sp." = " ", 
                  "sp" = " ",  "sp." = " ",
                  " spp " = " ",
                  "species" = " ",
                  "Dead species" = " ",
                  "Dead tree" = " ",
                  
                  
                  ## Vague
                  "Specified" = " ",
                  "Not Specified" = " ",
                  "UNKNOWN" = " ",
                  "unknown" = " ",
                  "Known" = " ",
                  "Unknown" = " ",
                  "Unidentified" = " ",
                  "UNIDENTIFIED" = "",
                  "check species" = "",
                  "Vacant Planting" = " ",
                  "Unknown unknown" = " ",
                  "Not Known"       = " ",
                  "natives"         = " ",
                  
                  
                  ## Replace
                  "Euc." = "Eucalyptus", 
        
                  # ## varieties, etc
                  " x "      = " ",
                  " X "      = " ",
                  " * "      = " ",
                  " v "      = " ",
                  "XL"       = " ",
                  "Natives"  = " ",
                  "Special" = " ",
                  "special" = " ",
                  "M." = " ",
                  "Kings Park" = " ",

                  " was "    = " ",
                   "Yes"     = " ",
                  "test"     = " ",
                  
                  "subsp."   = " ",
                  " p. "     = " ",
                  " P. "     = " ",
                  " f. "     = " ",
                  "var."     = " ",
                  "Vacant2"  = " ",
                  "Vacant"   = " ",
                  "Not Specified" = " ",
                  "Not Specified Not Specified" = " ",
                  "To Define" = " ",
                  "PINK"    = " ",

                  "ssp." = " ",
                  "Planting" = " ",
                  "cv" = " ",
                  "stump" = " ", 
                  "Stump" = " ",
                  "no tree" = " ",
                  "Recorded" = " ",
                  "NOT_TREE" = " ",
                  "No tree" = " ",
                  "Shrub" = " ",
                  "Mixed" = " ",
                  "mixed" = " ",
                  "not applicable" = " ",
                  "Not Applicable" = " ",
                  "Not applicable Not applicable" = " ")


## Now use gsubfn to replace multiple terms
ALE.DF$Species = gsubfn(paste(names(to.replace),collapse = "|"), to.replace, ALE.DF$Species)
"spp" %in% ALE.DF$Species;"Spp" %in% ALE.DF$Species;"spp." %in% ALE.DF$Species


#########################################################################################################################
##  Use a regular expression to create a binomial column
ALE.DF$Species<- trim(ALE.DF$Species)
ALE.DF$Species = gsub("  ",     " ", ALE.DF$Species)
ALE.DF$Species = gsub("   ",     " ", ALE.DF$Species)
ALE.DF$Species = gsub("     ",     " ", ALE.DF$Species)

ALE.DF$Species = gsub("-  -",     "", ALE.DF$Species)
ALE.DF$Species = gsub("  -",      "", ALE.DF$Species)
setdiff(ALE.DF$Species, CLEAN.NICHE.CONTEXT$searchTaxon)


ALE.DF$Binomial <- sub('(^\\S+ \\S+).*', '\\1', ALE.DF$Species) # \\s = white space; \\S = not white space
head(ALE.DF)


## Make another column for just the genus
ALE.DF$Binomial<- trim(ALE.DF$Binomial)
ALE.DF$Spp.length <- sapply(strsplit(ALE.DF$Binomial, " "), length)
table(ALE.DF$Spp.length)
head((subset(ALE.DF, Spp.length > 2)))
ALE.SPP = subset(ALE.DF, Spp.length == 2)
length(unique(ALE.SPP$Binomial))


#########################################################################################################################
## Subset Ale's dataframe to just the species on this binomial list
ALE.DF.SUBSET        = ALE.SPP[ALE.SPP$Binomial %in% unique(ALE.SPP$Binomial), ]
dim(DRAFT.HIA.TAXA)    ## 948 species after we cut out the "spp."


#######################################################################################################################
## Which taxa have the second word capitalised, and are not captured by eliminating the "spp"?
grep('^\\S+ [A-Z]', ALE.DF$Binomial, val = TRUE)




#########################################################################################################################
#######################################################  TBC ############################################################ 
#########################################################################################################################
