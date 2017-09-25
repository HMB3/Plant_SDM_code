#########################################################################################################################
##########################################  GBFI API SPECIES  ########################################################### 
#########################################################################################################################



#########################################################################################################################
## 1). READ IN API SPECIES
#########################################################################################################################


## Read in each file as text?
Magnolia.g = gbif('Magnolia grandiflora', download = TRUE)
Magnolia.g.names = sort(names(Magnolia.g))
unique(Magnolia.g$taxonRank)


##
Magnolia.gbif = gbif(genus = "Magnolia", species = "Magnolia grandiflora", geo = TRUE,
                     removeZeros = TRUE, download = TRUE)


##
test.occ  = occ_search(scientificName = 'Magnolia grandiflora')$data
occ.names = sort(names(test.occ))
unique(test.occ$issues)


##
Fagus.s = read.csv("./data/base/HIA_LIST/GBIF/SPECIES/Fagus_sylvatica.csv", stringsAsFactors = FALSE)
Fagus.s.names = sort(names(Fagus.s))


## How different are the manually downloaded species?
setdiff(Magnolia.g.names, occ.names)
setdiff(occ.names, Magnolia.g.names)

setdiff(Magnolia.g.names, Fagus.s.names)
setdiff(Fagus.s.names, Magnolia.g.names)
setdiff(gbif.keep, Fagus.s.names)


##
unique(Magnolia.g$country)
unique(Betula.p$countryCode)
unique(Betula.p$coordinateUncertaintyInMeters)
unique(Betula.p$continent)


## Save all files as .Rdata, then try concatenating them
Betula.pendula     = read.csv("./data/base/HIA_LIST/GBIF/SPECIES/Betula_pendula.csv",     stringsAsFactors = FALSE)
Fagus.sylvatica    = read.csv("./data/base/HIA_LIST/GBIF/SPECIES/Fagus_sylvatica.csv",    stringsAsFactors = FALSE)
Fraxinus.excelsior = read.csv("./data/base/HIA_LIST/GBIF/SPECIES/Fraxinus_excelsior.csv", stringsAsFactors = FALSE)
Hedera.helix       = read.csv("./data/base/HIA_LIST/GBIF/SPECIES/Fraxinus_excelsior.csv", stringsAsFactors = FALSE)
Quercus.robur      = read.csv("./data/base/HIA_LIST/GBIF/SPECIES/Quercus_robur.csv",      stringsAsFactors = FALSE)


## Change column names
Betula.pendula     = rename(Betula.pendula,     lat  = decimalLatitude, lon  = decimalLongitude)
Fagus.sylvatica    = rename(Fagus.sylvatica,    lon  = long)
Fraxinus.excelsior = rename(Fraxinus.excelsior, lon  = long)
Hedera.helix       = rename(Hedera.helix,       lon  = long)
Quercus.robur      = rename(Quercus.robur,      lon  = long)


## Save
save(Betula.pendula,     file = paste("./data/base/HIA_LIST/GBIF/SPECIES/Betula pendula_GBIF_records.RData",     sep = ""))
save(Fagus.sylvatica,    file = paste("./data/base/HIA_LIST/GBIF/SPECIES/Fagus sylvatica_GBIF_records.RData",    sep = ""))
save(Fraxinus.excelsior, file = paste("./data/base/HIA_LIST/GBIF/SPECIES/Fraxinus excelsior_GBIF_records.RData", sep = ""))
save(Hedera.helix,       file = paste("./data/base/HIA_LIST/GBIF/SPECIES/Hedera helix_GBIF_records.RData",       sep = ""))
save(Quercus.robur,      file = paste("./data/base/HIA_LIST/GBIF/SPECIES/Quercus robur_GBIF_records.RData",      sep = ""))


#########################################################################################################################
## 2). TRY RENAMING AND CLEANING
#########################################################################################################################


## Remove impossible coordiates
Fagus.coord = Fagus.sylvatica

dframe(Fagus.coord)  %>%
  coord_impossible() %>%
  coord_incomplete() %>%
  coord_unlikely()

dim(Fagus.coord[1]) - dim(Fagus.sylvatica[1])


## remove duplicates
Fagus.1000 = Fagus.sylvatica [1:1000,]

dp <- dframe(Fagus.1000) %>% dedup()
dim(dp[1]) - dim(Fagus.1000[1])



## GBIF ISSUES
(key <- name_suggest(q = 'Helianthus annuus', rank = 'species')$key[1])
(res <- occ_search(taxonKey=key, limit=100))
head(gbif_issues())



## How would you test if this has worked? Try Geoclean
Fagus.coord.rename = rename(Fagus.coord, identifier = species, 
                           XCOOR = decimalLatitude, YCOOR = decimalLongitude, country = countryCode)

test = GeoClean(Fagus.coord.rename)



## Try getting spatial duplicates
library(spatstat)
x <- GBIF.RASTER.CONTEXT$lon ; y<-GBIF.RASTER.CONTEXT$lat
w <- ripras(x, y)
wp <- ppp(x,y, window = w)
dupv<-duplicated.ppp(wp)

x2<-x[which(dupv==FALSE)] 
y2<-y[which(dupv==FALSE)]


##
# function gbifapi { 
#   
#   curl -i –user popple_1500:Popple1500 -H "Content-Type: 
#   application/json" -H "Accept: application/json" -X POST -d "{\"creator\":
#   \”yourGbifUserName\", \"notification_address\": [\”hugh.burley@mq.edu.au\"], 
#   \"predicate\": {\"type\":\"and\", \"predicates\": [{\"type\":\"equals\",\"key\":
#   \"HAS_COORDINATE\",\"value\":\"true\"}, {\"type\":\"equals\", \"key\":
#   \"TAX O N _KEY\", \"value\":\"$1\"}] }}" 
#   http://api.gbif.org/v1/occurrence/download/request  >> log.txt 
#   echo -e "\r\n$1 $2\r\n\r\n----------------\r\n\r\n" >> log.txt
#   
# }  
# 
# $ gbifapi 4140730 "Betula pendula" 
# $ gbifapi 2882316 "Fagus sylvatica" 
# $ gbifapi 3172358 "Fraxinus excelsior" 
# $ gbifapi 8351737 "Hedera helix" 
# $ gbifapi 2878688 "Quercus robur"


