#########################################################################################################################
##############################  DOWNLOAD GLOBAL SPECIES RECORDS FOR HORT AUS LIST ####################################### 
#########################################################################################################################


## This code downloads all occurrence records for species on the Horticulture Australia list, from the GBIF database. 
## The trick here will be substituing the functions from the ALA4R package for the rgbif functions...
## Matching the core fields could be tricky for some taxa. Also, 


## setup 
library(gtools)
library(devtools)
#devtools::install_github("ropensci/rgbif")
library(data.table)
library(RCurl)
library(Rcpp)
library(raster)
library(rgdal)
library(plyr)

library(SDMTools)
library(rmaxent)
library(dismo)
library(ALA4R)
library(rgbif)
library(speciesgeocodeR)
library(raster)

library(knitr)
library(htmltools)
library(yaml)
library(caTools)
library(bitops)
library(rmarkdown)
library(gsubfn)
library(functional)
library(splitstackshape)


library(tidyverse)
library(stringr)
library(maptools)
library(rgeos)
library(magrittr)


## source functions
source('./R/GREEN_CITIES_FUNCTIONS.R')





#########################################################################################################################
## 1). READ IN DRAFT HIA LIST AND CLEAN
#########################################################################################################################


## this list derives from all species and varieties sold anywhere in Australia in the last 5 years. Anthony Maneahas cleaned 
## up the data and cross-linked to growth form and exotic/native status and derived a list of ~1000 species that are the most 
## commonly sold, covering the right ratio of growth forms, regional representation and native/exotic
spp.list = read.csv("./data/base/HIA_LIST/HIA/GREEN_CITIES_DRAFT_LIST.csv", stringsAsFactors = FALSE)
top.200  = read.csv("./data/base/HIA_LIST/HIA/HIA_TOP_200.csv", stringsAsFactors = FALSE)


## have a look
dim(spp.list)
str(spp.list)
head(spp.list)


## also, add the "Top 200 species in here"
spp.200          = top.200[c("Species")]
spp.200$Top_200  = "TRUE"

spp.list = merge(spp.list, spp.200, by = "Species", all.x = TRUE) 
spp.list$Top_200[is.na(spp.list$Top_200)] <- "FALSE"
spp.list$Origin <- gsub(" ",  "", spp.list$Origin)


## check
str(spp.list)
head(spp.list)
unique(spp.list$Top_200)
unique(spp.list$Origin)


#########################################################################################################################
## WHAT ARE SOME USEFUL MEASURES OF THE DATASET?
#########################################################################################################################


## ~ HALF the total species are native, ~47% of the top 200
dim(subset(spp.list, Origin == "Native"))[1]/dim(spp.list)[1]*100
dim(subset(top.200, Origin == "Native"))[1]/dim(top.200)[1]*100


#########################################################################################################################
## DRAFT CLEAN: THIS NEEDS TO CHANGE IN CONSULTATION WITH RACH, PAUL, ETC
#########################################################################################################################


## for now, remove all the subspecies, varieties etc.
## but check if some of them work on GBIF, e.g. a separate list of varities
spp.list$Species <- gsub("spp.", "", spp.list$Species)
spp.list$Species <- gsub(" x ",  " ", spp.list$Species)
spp.list$Species <- gsub("  ",  " ", spp.list$Species)


## then remove the varieties and subspecies
spp.list$Species = vapply(lapply(strsplit(spp.list$Species, " "),
                                                 unique), paste, character(1L), collapse = " ")


## then get just the first two words (again cleaning up the subspecies, and single genera)
spp.list$Species = vapply(lapply(strsplit(spp.list$Species, " "), 
                                 string_fun_first_two_words), paste, character(1L), collapse = " ")


## if exclude NA, spp.list = spp.list[!grepl("NA", spp.list$Species),]
spp.list = spp.list[with(spp.list, order(Species)), ] 


## check
str(spp.list)
View(spp.list)


## Now create list of HIA taxa. 768 unique species, mius the corrections, etc. 
spp.list = spp.list[with(spp.list, order(Species)), ]
spp = unique(as.character(spp.list$Species))
str(spp)   ## why 680? later check on what happens with the different queries
head(spp, 50)
tail(spp, 50)


## also make a genus list?
genera = unique(vapply(lapply(strsplit(spp.list$Species, " "), 
                              string_fun_first_word), paste, character(1L), collapse = " "))


## check
str(genera)
head(genera, 50)
tail(genera, 50)





#########################################################################################################################
## 2). DOWNLOAD RECORDS FROM GBIF USING HIA LIST
#########################################################################################################################


## now run lOOPs to dowload species in the "spp" list from GBIF
## not including any data quality checks here, just downloading everything ##
## 


## set a few global variables to be used inside the functions...
#GBIF.download.limit = 200000
skip.spp.list       = list()
skip.gen.list       = list()


## run the download function on the species and genera lists
skipped.species = download_GBIF_all_species(spp)    ## saves each spp as .Rdata file, returning list of skipped spp 
skipped.genera  = download_GBIF_all_genera(genera)  ## saves each gen as .Rdata file, returning list of skipped genera 


## now convert the lists of skipped species and genera into a dataframe
str(skipped.species)
str(skipped.genera)


## converting the lists of skipped species and genera into a dataframe
skipped.species <- data.frame(matrix(unlist(skipped.species), nrow = length(skipped.species), byrow = TRUE))
skipped.genera  <- data.frame(matrix(unlist(skipped.genera),  nrow = length(skipped.spp), byrow = TRUE))


## split the reason and the species into separate columns
skipped.species <- cSplit(skipped.species, 1:ncol(skipped.species), sep = "|", stripWhite = TRUE, type.convert = FALSE)


## update names
colnames(skipped.species)[1] <- "Reason_skipped"
colnames(skipped.species)[2] <- "Taxa"


## check
str(skipped.species)
head(skipped.species)
tail(skipped.species)
View(skipped.species)


## get subset for each type
max.records  <- skipped.species[ which(skipped.species$Reason_skipped == "Number of records > 200,000"), ]
name.records <- skipped.species[ which(skipped.species$Reason_skipped == "Possible incorrect nomenclature"), ]
no.records   <- skipped.species[ which(skipped.species$Reason_skipped == "No GBIF records"), ]


## create lists for each category
max.records.spp  = unique(as.character(max.records$Taxa))
name.records.spp = unique(as.character(name.records$Taxa))
no.records.spp   = unique(as.character(no.records$Taxa))





#########################################################################################################################
## 2). LOAD GBIF RECORDS FOR EACH SPECIES INTO ONE DATAFRAME
#########################################################################################################################


## could also do the cleaing for each species individually?
load("./data/base/HIA_LIST/GBIF/Viburnum suspensum_GBIF_records.RData")
str(GBIF)


## Load a 1km grid of Australia (should be 250m grid?)
map = raster("./data/base/world_clim/wc2.0_30s_prec_01.tif")
plot(map)


## Simple table of columns needed from GBIF to check data quality
## problem here is that not all species have the same columns!
## can we get all the columns at once, even if they have a different number?
## look at all the GBIF fields for an example
## do these GBIF names match the ALA names?


## Stu: As to binding multiple data frames that may not all have the same columns: dplyr::bind_rows() will 
## do this. I do this once all the individual files are downloaded for each species. This can get very 
## slow for a large number of data frames if you do them one by one, I found a quicker way was to read 
## all data frames to a list, then use a single bind_rows() to do them all in one hit.


## check GBIF field names for a key species...
Magnolia.grandiflora = gbif('Magnolia grandiflora', download = TRUE)
GBIF.names = sort(names(Magnolia.grandiflora))
GBIF.names


## can only do this once/if fields are decided?
GBIF.HIA.RECORDS = data.frame()
  
  ## taxonomic resolution fields
  # "occurrenceID"     = character(),
  # "gbifID"           = numeric(),
  # #"basisOfRecord"    = character(),
  # #"taxonRank"        = character(),
  # #"taxonRankID"      = character(),
  # "type"             = character(),
  # "institutionCode"  = character(),
  # "year"             = numeric(),
  # "eventDate"        = numeric(),
  # "identifiedBy"     = character(),
  # 
  # ## taxonomic fields
  # "scientificName"   = character(),
  # "order"            = character(),
  # "family"           = character(),
  # "genus"            = character(),
  # "species"          = character(),
  # 
  # ## spatial fields
  # "lon"              = numeric(),
  # "lat"              = numeric(),
  # #"geodeticDatum"    = numeric(),
  # "coordinateUncertaintyInMeters" = numeric(),
  # #"elevation"        = numeric(),
  # "country"          = character(),
  # "locality"         = character()
  # 
)



## spp doesn't match the downloaded species
## now create a list from the downloaded files
spp.download = list.files("./data/base/HIA_LIST/GBIF/SPECIES/", pattern = ".RData")
spp.download = gsub("_GBIF_records.RData", "", spp.download)
str(spp.download)



## for all the species in the HIA list
for(sp.n in spp.download) {
  
  # Load GBIF records for each species as a .CSV file
  filename = paste0("./data/base/HIA_LIST/GBIF/SPECIES/", sp.n, "_GBIF_records.RData")
  load(filename)
  
  ## the restrict the dataframe to only those in common...can olny do this once fields are decided....
  GBIF = GBIF
  # GBIF = GBIF[,c(
  # 
  #   ## taxonomic resolution fields
  #   "occurrenceID",
  #   "gbifID",
  #   #"basisOfRecord",
  #   #"taxonRank",
  #   #"taxonRankID",
  #   "type",
  #   "institutionCode",
  #   "year",
  #   "eventDate",
  #   "identifiedBy",
  #   
  #   ## taxonomic fields
  #   "scientificName",
  #   "order",
  #   "family",
  #   "genus",
  #   "species",
  #   
  #   ## spatial fields
  #   "lon",
  #   "lat",
  #   #"geodeticDatum",
  #   "coordinateUncertaintyInMeters",
  #   #"elevation",
  #   "country",
  #   "locality" 
  #   
  #   )]

  
  ## then bind the different species dataframes together
  #GBIF.HIA.RECORDS  = rbind(GBIF.HIA.RECORDS, GBIF)
  GBIF.HIA.RECORDS  = bind_rows(GBIF.HIA.RECORDS, GBIF)
  
}


## now check the combined data frame
dim(GBIF.HIA.RECORDS)
str(GBIF.HIA.RECORDS)
View(GBIF.HIA.RECORDS)
head(GBIF.HIA.RECORDS, 9)
tail(GBIF.HIA.RECORDS, 10)


## get numeric spatial error for later processing, somehow it was converted to character
summary(GBIF.HIA.RECORDS$coordinateUncertaintyInMeters)


## now plot records, just to see what they look like
GBIF.points       = GBIF.HIA.RECORDS[,c("lon", "lat")]
GBIF.points       = data.matrix(GBIF.points)

plot(map)
points(GBIF.points, cex = 0.1, col = "blue")


## Save tables
save(GBIF.HIA.RECORDS,   file = paste("./data/base/HIA_LIST/GBIF/GLOBAL_GBIF_HIA_records_RAW_DATE.RData",       sep = ""))
write.csv(GBIF.HIA.RECORDS,       "./data/base/HIA_LIST/GBIF/GLOBAL_GBIF_HIA_records_RAW_DATE.csv",       row.names = FALSE)





#########################################################################################################################
## 3). CLEAN GBIF RECORDS
#########################################################################################################################

## GBIF records could also be cleaned one at a time...



## which fields do 

## read the GBIF data back in
GBIF.HIA.RECORDS = read.csv("./data/GLOBAL_GBIF_HIA_records_RAW_DATE.csv", stringsAsFactors = FALSE)


## have a look, 1951283 - 1933051 = 18,232
dim(GBIF.HIA.RECORDS)
head(GBIF.HIA.RECORDS)
View(GBIF.HIA.RECORDS)
summary(GBIF.HIA.RECORDS)





#########################################################################################################################
## now how do we eliminate bad records?
## TYPE
unique(GBIF.HIA.RECORDS$type)


## remvove species with certain observation type?
# records <- subset(records, 
#                   recordType == 'PhysicalObject'    | 
#                   recordType == 'specimen'          | 
#                   recordType == 'PreservedSpecimen' ) 


## most records don't have a type? 
dim(subset(GBIF.HIA.RECORDS, type != "NA"))[1] - dim(GBIF.HIA.RECORDS)[1]


#########################################################################################################################
## DATE
summary(GBIF.HIA.RECORDS$year)


## histogram
hist(GBIF.HIA.RECORDS$year, 
     xlab = 'Collection year', 
     ylab = 'Number of records',
     main = 'Collection year (all years)',
     border = NA, col = "grey", breaks = 50)


## X recrods recorded before 1950
dim(GBIF.HIA.RECORDS[GBIF.HIA.RECORDS$year >= 1950, ])[1] - dim(GBIF.HIA.RECORDS)[1]


#########################################################################################################################
## COORDINATES: X records with no lat/long
nrow(GBIF.HIA.RECORDS[-which(is.na(GBIF.HIA.RECORDS$lon) | is.na(GBIF.HIA.RECORDS$lat)), ]) -
  nrow(GBIF.HIA.RECORDS)


## now plot records, just to see what they look like
GBIF.points       = GBIF.HIA.RECORDS[,c("lon", "lat")]
GBIF.points       = data.matrix(GBIF.points)

plot(map)
points(GBIF.points, cex = 0.1, col = "blue")


## X rows with 0,0
nrow(GBIF.points[-which(GBIF.HIA.RECORDS$lon == 0 & GBIF.HIA.RECORDS$lat == 0), ])


## now look at coordinate uncertainty
hist(GBIF.HIA.RECORDS$coordinateUncertaintyInMeters, 
     breaks = 50,
     border = NA,
     xlab   = 'Coordinate uncertainty (m)', 
     ylab   = 'Number of records', main = '', 
     col    = 'red')


## where are the records with errors >5000 m, or which have no stated uncertainty?
GBIF.UNCERT <- GBIF.HIA.RECORDS[GBIF.HIA.RECORDS$coordinateUncertaintyInMeters > 5000 | 
                                  is.na(GBIF.HIA.RECORDS$coordinateUncertaintyInMeters), ]


## most of the records are uncertain?
nrow(GBIF.UNCERT)


## plot all the records, then the uncertain ones...
plot(map)
points(GBIF.points, pch = 21, bg = 'mediumseagreen') 
points(GBIF.UNCERT$lon, GBIF.UNCERT$lat, pch = 2, col = 'red') 


## add a legend?
legend('bottomleft', 
       legend = c('OK', 'Uncertain'), 
       pch    = c(16, 2), 
       col    = c('mediumseagreen', 'red'), 
       bg     = 'white')



# key = (key <- name_suggest(q = 'Magnolia grandiflora', rank = 'species')$key[1])
# res = (res <- occ_search(taxonKey=key, limit=100))







#########################################################################################################################
############################################  END OF GBIF DOWLOAD CODE ################################################## 
#########################################################################################################################