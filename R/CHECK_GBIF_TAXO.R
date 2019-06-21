#########################################################################################################################
## CHECK TAXONOMY RETURNED BY GBIF USING TAXONSTAND
######################################################################################################################### 


## The problems is the mismatch between what we searched, and what GBIF returned. We can check this by taking the 
## scientificName, and run that through TPL. Then, of the names in this list which are accepted, but which don't match 
## our list, get rid of them.



#########################################################################################################################
## Use "Taxonstand" to check the taxonomy.
GBIF.TAXO <- TPL(unique(GBIF.TRIM$scientificName), infra = TRUE,
                 corr = TRUE, repeats = 100)  ## to stop it timing out...
sort(names(GBIF.TAXO))



#########################################################################################################################
## The procedure used for taxonomic standardization is based on function TPLck. A progress bar
## indicates the proportion of taxon names processed so far. 


## Check the taxonomy by taking scientificName, and run that through TPL. Then, of the names in this list which are 
## accepted, but which don't match our list Get rid of them.


## Then join the GBIF data to the taxonomic check, using "scientificName" as the join field...
GBIF.TRIM.TAXO <- GBIF.TRIM %>%
  left_join(., GBIF.TAXO, by = c("scientificName" = "Taxon"))
names(GBIF.TRIM.TAXO)


## However, the scientificName string and the searchTaxon string are not the same. So to test which SN are accepted but
## not on the ST list, we need string matching using regular expressions, or the like.
## currently using str_detect
Match.SN = GBIF.TRIM.TAXO  %>%
  mutate(Match.SN.ST = 
           str_detect(scientificName, searchTaxon)) %>%
  
  select(one_of(c("scientificName",
                  "searchTaxon",
                  "Taxonomic.status",
                  "New.Taxonomic.status",
                  "New.Genus",
                  "New.Species",
                  "country",
                  "Match.SN.ST")))


## How many records don't match? 44k or 1.5%
unique(Match.SN$Taxonomic.status)
with(Match.SN, table(Match.SN.ST))
with(Match.SN, table(Taxonomic.status))


## Make a list of the species names which didn't match
accepted.SN   = unique(subset(Match.SN, Taxonomic.status == "Accepted")$scientificName)
synonym.SN    = unique(subset(Match.SN, Taxonomic.status == "Synonym")$scientificName)
unresolved.SN = unique(subset(Match.SN, Taxonomic.status == "Unresolved")$scientificName)
misapp.SN     = unique(subset(Match.SN, Taxonomic.status == "Misapplied")$scientificName)
blank.SN      = unique(subset(Match.SN, Taxonomic.status == "")$scientificName)



## What are the proportions
round(with(Match.SN, table(Match.SN.ST)/sum(table(Match.SN.ST))*100), 2)
round(with(Match.SN, table(Taxonomic.status)/sum(table(Taxonomic.status))*100), 2)
View(Match.SN)


#########################################################################################################################
## Get the subset of species which are accpeted, but not on our list
acc.no.match = subset(Match.SN, Taxonomic.status == "Accepted" & Match.SN.ST == "FALSE")
acc.mis.sn  = unique(acc.no.match$scientificName)    ## the unique scientific names we didn't ask for
acc.mis.st  = unique(acc.no.match$searchTaxon)       ## The unique species returning dodgy scientific names
View(acc.no.match)



#########################################################################################################################
###################################################### TBC ############################################################## 
#########################################################################################################################