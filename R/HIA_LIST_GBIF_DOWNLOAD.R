#########################################################################################################################
##############################  DOWNLOAD GLOBAL SPECIES RECORDS FOR HORT AUS LIST ####################################### 
#########################################################################################################################


# This code downloads all occurrence records for species on Horticulture Australia list from the GBIF database. First we week
# can get all the records for just one list, then once that works, we can create a hieararchy as Rachel outlined.
# The trick here will be substituing the functions from the ALA4R package for the rgbif functions...


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


## source
source('./R/GREEN_CITIES_FUNCTIONS.R')


#########################################################################################################################
## 1). DOWNLOAD RECORDS FROM GBIF USING HIA LIST
#########################################################################################################################


## have a look at GBIF: e.g. return the total no. of records on the database?
occ_count(basisOfRecord = 'OBSERVATION')
occ_count(georeferenced = TRUE)


## search for a genus
head(name_lookup(query = 'Magnolia', rank = "genus", return = "data"), 20)


## how bout the no. 1 species on the HIA list?
dim(name_lookup(query  = 'Magnolia grandiflora', rank = "species", return = 'data'))
head(name_lookup(query = 'Magnolia grandiflora', rank = "species", return = 'data'), 20)


## now read in the HIA list: this could be one list with different criteria to create the hierarchy, 
## or several lists
spp.list = read.csv("./data/base/HIA_LIST/HIA/HIA_DATE.csv", stringsAsFactors = FALSE)


## also consider how to download, what is the taxonomic resolution? If the varieties are troublesome, are we best
## downloading all the records for a genus, rather than a species? Discuss with Rach once I've got more work done...
dim(spp.list)
str(spp.list)
head(spp.list)


## now remove all the subspecies, etc,
# test = gsubfn(".", list("subsp." = "", 
#                         "var."   = ""), paste(spp.list$Matched_Scientific_name))
# spp.list$Matched_Scientific_name = gsub("subsp.", "",      paste(spp.list$Matched_Scientific_name))
# spp.list$Matched_Scientific_name = gsub("sect.", "",       paste(spp.list$Matched_Scientific_name))
# spp.list$Matched_Scientific_name = gsub("var.", "",        paste(spp.list$Matched_Scientific_name))
# spp.list$Matched_Scientific_name = gsub("Plantago a", "",  paste(spp.list$Matched_Scientific_name))
# spp.list$Matched_Scientific_name = gsub("  ", " ",         paste(spp.list$Matched_Scientific_name))


# test = gsubfn(".", list("subsp.", "", 
#                         "sect.", "",
#                         "var.", "",
#                         "Plantago a", "",
#                         "  ", " "), spp.list$Matched_Scientific_name)


## now remove all the subspecies, etc,
## pipe through?
spp.list$Matched_Scientific_name %<>%
  
  ## list all the character strings and what to replace them with...
  gsub("subsp.", "", .) %>%
  gsub("sect.", "", .) %>%
  gsub("var.", "", .) %>%
  gsub("Plantago a", "", .) %>%
  # gsub("#Agrostis capillaris", "", .) %>%
  # gsub("#Agrostis capillaris", "", .) %>%
  gsub("  ", " ", .)


## check the list 
str(spp.list$Matched_Scientific_name)


## then remove the varieties and subspecies
spp.list$Matched_Scientific_name = vapply(lapply(strsplit(spp.list$Matched_Scientific_name, " "), 
                                                 unique), paste, character(1L), collapse = " ")


## then get just the first two words (again cleaning up the subspecies, and single genera)
spp.list$Matched_Scientific_name = vapply(lapply(strsplit(spp.list$Matched_Scientific_name, " "), 
                                                 string_fun), paste, character(1L), collapse = " ")


## this creates "NA" species where only a genus was given, so just exclude these rows
spp.list = spp.list[!grepl("NA", spp.list$Matched_Scientific_name),]
spp.list = spp.list[!(spp.list$Matched_Scientific_name == ""), ]
spp.list = spp.list[with(spp.list, order(Matched_Scientific_name)), ] 

str(spp.list)
View(spp.list)


## now what is the simplest way to do the download? EG: download all fields for a species, then do the culling later?
## n = 5 * ~ 200 species, so records for ~1000 taxa. So not too many to download at once! 
Magnolia.grandiflora = gbif('Magnolia grandiflora', download = TRUE)


## look at all the GBIF fields for an example
## do these GBIF names match the ALA names?
GBIF.names = sort(names(Magnolia.grandiflora))
GBIF.names


## Might want to have this as "refresh" but trying off first. 
## The cache should clear every time an R session is started anyway.
## R is unstable with this package and frequently crashes
#ala_config(caching = "off")


## Now create list of HIA taxa
spp = unique(as.character(spp.list$Matched_Scientific_name))
str(spp)   ## why 680? later check on what happens with the different queries
head(spp, 20)
tail(spp, 20)





#########################################################################################################################
## now run a lOOP to dowload species in the "spp" list from GBIF
## not including any data quality checks here, just downloading everything
GBIF.download.limit = 200000


## for every species in the list
for(sp.n in spp){
  
  
  ## 1). First, check if the f*&%$*# file exists
  file = paste0("./data/base/HIA_LIST/GBIF/", sp.n, "_GBIF_records.RData")
  
  
  ## need this for longer list, could be up to thousands of species in total
  ## also GBIF, ALA, AVH will all fail intermittently 
  if (file.exists (file)) {
    
    print (paste ("file exists for species", sp.n, "skipping"))
    next
    
  }
  
  ## 2). Check the spelling...incorrect nomenclature will return NULL result
  if (is.null(occ_search(scientificName = sp.n, limit = 1)$meta$count) == TRUE) {
    
    print (paste ("Incorrect nomenclature", sp.n, "skipping"))
    next
    
  }
  
  ## 3). Skip species with no records
  if (occ_search(scientificName = sp.n)$meta$count == 0) {
    
    print (paste ("No GBIF records for", sp.n, "skipping"))
    next
    
  }
  
  ## 4). Check how many records there are, and skip if there are over 200k
  if (occ_search(scientificName = sp.n, limit = 1)$meta$count > GBIF.download.limit) {
    
    print (paste ("Number of records > max for GBIF download via R (200,000)", sp.n, "skipping"))
    next
    
  }
  
  ## 5). Download ALL records from GBIF
  GBIF = gbif(sp.n, download = TRUE)   ## could use more arguments here, download_reason_id = 7, etc.
  
  ## 6). save records to .Rdata file, note that using .csv files seemed to cause problems...
  save(GBIF, file = paste("./data/base/HIA_LIST/GBIF/", sp.n, "_GBIF_records.RData", sep = ""))
  
} 


## Errors...
## The number of records is larger than the maximum for download via this service (200,000)





#########################################################################################################################
## 2). LOAD GBIF RECORDS FOR EACH SPECIES INTO ONE DATAFRAME
#########################################################################################################################


## could also do the cleaing for each species individually?
load("./data/base/HIA_LIST/GBIF/Agonis flexuosa_GBIF_records.RData")
str(GBIF)

test            = read.csv("./data/base/HIA_LIST/GBIF/Corymbia maculata_GBIF_records.csv")
Opuntia.stricta = read.csv("./data/base/HIA_LIST/GBIF/Opuntia stricta_GBIF_records.csv")
str(Opuntia.stricta)


## Load a 1km grid of Australia (should be 250m grid?)
map = raster("./data/base/world_clim/wc2.0_30s_prec_01.tif")
plot(map)


## Simple table of columns needed from GBIF to check data quality
## problem here is that not all species have the same columns!
## can we get all the columns at once, even if they have a different number?
GBIF.names


## can only do this once/if fields are decided?
GBIF.HIA.RECORDS = data.frame(
  
  ## taxonomic resolution fields
  "occurrenceID"     = character(),
  "gbifID"           = numeric(),
  #"basisOfRecord"    = character(),
  #"taxonRank"        = character(),
  #"taxonRankID"      = character(),
  "type"             = character(),
  "institutionCode"  = character(),
  "year"             = numeric(),
  "eventDate"        = numeric(),
  "identifiedBy"     = character(),
  
  ## taxonomic fields
  "scientificName"   = character(),
  "order"            = character(),
  "family"           = character(),
  "genus"            = character(),
  "species"          = character(),
  
  ## spatial fields
  "lon"              = numeric(),
  "lat"              = numeric(),
  #"geodeticDatum"    = numeric(),
  "coordinateUncertaintyInMeters" = numeric(),
  #"elevation"        = numeric(),
  "country"          = character(),
  "locality"         = character()
  
)


## for all the species in the HIA list
for(sp.n in spp){
  

  # Load GBIF records for each species as a .CSV file
  filename = paste0("./data/base/HIA_LIST/GBIF/", sp.n, "_GBIF_records.csv")
  load(filename)
  
  ## the restrict the dataframe to only those in common...can olny do this once fields are decided....
  GBIF = GBIF[,c(

    ## taxonomic resolution fields
    "occurrenceID",
    "gbifID",
    #"basisOfRecord",
    #"taxonRank",
    #"taxonRankID",
    "type",
    "institutionCode",
    "year",
    "eventDate",
    "identifiedBy",
    
    ## taxonomic fields
    "scientificName",
    "order",
    "family",
    "genus",
    "species",
    
    ## spatial fields
    "lon",
    "lat",
    #"geodeticDatum",
    "coordinateUncertaintyInMeters",
    #"elevation",
    "country",
    "locality" 
    
    )]

  
  ## then bind the different species dataframes together
  GBIF.HIA.RECORDS  = rbind(GBIF.HIA.RECORDS, GBIF)
  
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
## Maxent models will be run one at a time...

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