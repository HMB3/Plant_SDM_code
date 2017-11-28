#########################################################################################################################
########################################  'CLEAN' canonicalName LIST ########################################################## 
#########################################################################################################################


#########################################################################################################################
## 1). READ IN DRAFT HIA LIST
#########################################################################################################################

## The aim here is to take the raw list of plants grown in Aus supplied by Matt Plumber and Anthony Manea (the 'clean tab'), and
## then clean the list as best as possible in R to use the canonicalName binomial as the unit of downloading and analysis.


#########################################################################################################################
## This list derives from
APC.list = read.csv("./data/base/TRAITS/apc_taxon_distribution.csv", stringsAsFactors = FALSE)
dim(APC.list)
APC.TAXA        = APC.list
str(unique(APC.TAXA$canonicalName))


## Remove weird characters, subspecies and varieties
APC.TAXA$canonicalName = gsub(" x",     "",  APC.TAXA$canonicalName)
APC.TAXA$canonicalName = gsub("NA",     "",  APC.TAXA$canonicalName)
APC.TAXA$canonicalName = gsub("  ",     " ", APC.TAXA$canonicalName)
APC.TAXA$canonicalName = gsub(" $",     "",  APC.TAXA$canonicalName, perl = TRUE)
APC.TAXA$canonicalName = gsub("    $",  "",  APC.TAXA$canonicalName, perl = TRUE)
APC.TAXA$canonicalName = gsub("subsp.", "",  APC.TAXA$canonicalName, perl = TRUE)
APC.TAXA$canonicalName = gsub("var.",   "",  APC.TAXA$canonicalName, perl = TRUE)


## Create another binomial column: this removes taxa with > 2 words
APC.TAXA$Binomial <- sub('(^\\S+ \\S+).*', '\\1', APC.TAXA$canonicalName) # \\s = white space; \\S = not white space


## Which taxa have the second word capitalised, and are not captured by eliminating the "spp"?
grep('^\\S+ [A-Z]', APC.TAXA$Binomial, val = TRUE)
length(unique(APC.TAXA$Binomial));length(unique(APC.TAXA$canonicalName))
tail(unique(APC.TAXA$Binomial));tail(unique(APC.TAXA$canonicalName))



#########################################################################################################################
## 2). CREATE FINAL UNIQUE LIST FOR GBIF SEARCHING
#########################################################################################################################


#######################################################################################################################
## Now for GBIF, just get the unique canonicalName...
APC.SPP = APC.TAXA[!duplicated(APC.TAXA["Binomial"]),]


## Reorder by canonicalName
APC.SPP = APC.SPP[with(APC.SPP, order(Binomial)), ] 
SPP.APC = unique(as.character(APC.SPP$Binomial))
length(SPP.APC)   





#########################################################################################################################
############################################  END OF APC LIST CODE ###################################################### 
#########################################################################################################################