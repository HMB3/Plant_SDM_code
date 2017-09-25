#########################################################################################################################
##########################################  GBFI API SPECIES  ########################################################### 
#########################################################################################################################



#########################################################################################################################
## 1). READ IN API SPECIES
#########################################################################################################################


## Read in each file as text?
Magnolia.g = gbif('Magnolia grandiflora', download = TRUE)
Magnolia.g.names = sort(names(Magnolia.g))


Betula.p = read.csv("./data/base/HIA_LIST/GBIF/SPECIES/Betula_pendula.csv", stringsAsFactors = FALSE)
Betula.p.names = sort(names(Betula.p))


## How different are the manually downloaded species?
setdiff(Magnolia.g.names, Betula.p.names)
setdiff(Betula.p.names, Magnolia.g.names)
setdiff(gbif.keep, Betula.p.names)


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





#########################################################################################################################
## 2). TRY RENAMING AND CLEANING
#########################################################################################################################


## Remove impossible coordiates
Betula.coord = Betula.pendula

dframe(Betula.coord) %>%
  coord_impossible() %>%
  coord_incomplete() %>%
  coord_unlikely()

dim(Betula.coord[1]) - dim(Betula.pendula[1])


## remove duplicates
dp <- dframe(Betula.pendula) %>% dedup()
dim(dp[1]) - dim(Betula.pendula[1])


## How would you test if this has worked?


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


