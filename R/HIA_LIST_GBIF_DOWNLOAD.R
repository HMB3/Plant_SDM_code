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
library(dismo)
library(ALA4R)
library(rgbif)
library(raster)





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


## now read in the HIA list: this could be one list with different criteria to create the hierarchy, or several lists
#spp.list = read.csv("./data/base/HIA_LIST_date.csv", stringsAsFactors = FALSE)


## also consider how to download, what is the taxonomic resolution? If the varieties are troublesome, are we best
## downloading all the records for a genus, rather than a species? Discuss with Rach once I've got more work done...


## just use a dummy list...
spp.list = c("Magnolia grandiflora", 
             "Macadamia integrifolia", 
             "Macadamia jansenii",
             "Macadamia ternifolia",
             "Macadamia tetraphylla")


## have a look
dim(spp.list)
str(spp.list)
head(spp.list)

 
## now what is the simplest way to do the download? EG: download all fields for a species, then do the culling later?
## n = 5 * ~ 200 species, so records for ~1000 taxa. So not too many to download at once! 
magnolia.test           = gbif('Magnolia',               download = TRUE)
Magnolia.grandiflora    = gbif('Magnolia grandiflora',   download = TRUE)
Macadamia.integrifolia  = gbif('Macadamia integrifolia', download = TRUE)
Macadamia.jansenii      = gbif('Macadamia jansenii',     download = TRUE)
Macadamia.ternifolia    = gbif('Macadamia ternifolia',   download = TRUE)



## look at all the GBIF fields...
GBIF.names = sort(names(Magnolia.grandiflora))
GBIF.names


## Might want to have this as "refresh" but trying off first. The cache should clear every time an R session is started anyway.
## R is unstable with this package and frequently crashes
ala_config(caching = "off")


## Now create list of HIA taxa
#spp = unique(as.character(spp.list$binomial))
spp = spp.list
str(spp)
head(spp, 20)
tail(spp, 20)





#########################################################################################################################
## now run a lOOP to dowload species in the "spp" list from the GBIF
for(sp.n in spp){


  ## First, check if the f*&%$*# file exists
  file = paste0("./data/base/HIA_LIST/GBIF/", sp.n, "_GBIF_records.RData")
  
  
  if (file.exists (file)) {
    
    print (paste ("file exists for species", sp.n, "skipping"))
    next
    
  }

  ## 1). download all records from GBIF
  #GBIF = occurrences(taxon = sp.n, download_reason_id = 7)
  GBIF = gbif(sp.n, download = TRUE)
    

  # 8. save.
  save(GBIF, file = paste("./data/base/HIA_LIST/GBIF/", sp.n, "_GBIF_records.csv", sep = ""))
  
} 





#########################################################################################################################
## 2). LOAD GBIF RECORDS FOR EACH SPECIES INTO ONE DATAFRAME
#########################################################################################################################


## Load a 1km grid of Australia (should be 250m grid?)
map = raster("./data/base/world_clim/wc2.0_30s_prec_01.tif")
plot(map)


## Simple table of details needed from GBIF
## problem here is that not all species have the same columns...
## can only do this once/if fields are decided

GBIF.HIA.RECORDS = data.frame(
  
  ## taxonomic resolution fields
  "occurrenceID"     = character(),
  #"basisOfRecord"    = character(),
  #"taxonRank"        = character(),
  #"taxonRankID"      = character(),
  "institutionCode"  = character(),
  "year"             = numeric(),
  
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
  "country"          = character()
  
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
    #"basisOfRecord",
    #"taxonRank",
    #"taxonRankID",
    "institutionCode",
    "year",
    
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
    "country"          
    
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
GBIF.points   = GBIF.HIA.RECORDS[,c("lon", "lat")]
GBIF.points   = data.matrix(GBIF.points)
plot(map)
points(GBIF.points, cex = 0.1, col = "blue")


## Save tables
save(GBIF.HIA.RECORDS,   file = paste("./data/base/HIA_LIST/GBIF/GLOBAL_GBIF_HIA_records_DATE.RData",       sep = ""))
write.csv(GBIF.HIA.RECORDS,       "./data/base/HIA_LIST/GBIF/GLOBAL_GBIF_HIA_records_DATE.csv",       row.names = FALSE)





#########################################################################################################################
## 3). CLEAN GBIF RECORDS
#########################################################################################################################


## read the GBIF data back in
GBIF.HIA.RECORDS = read.csv("./data/GLOBAL_GBIF_HIA_records_DATE.csv", stringsAsFactors = FALSE)


## have a look, 1951283 - 1933051 = 18,232
dim(GBIF.HIA.RECORDS)
head(GBIF.HIA.RECORDS)
View(GBIF.HIA.RECORDS)
summary(GBIF.HIA.RECORDS)





#########################################################################################################################
## now how do we eliminate bad records?
Magnolia.grandiflora
key = (key <- name_suggest(q = 'Magnolia grandiflora', rank = 'species')$key[1])
res = (res <- occ_search(taxonKey=key, limit=100))





#########################################################################################################################
############################################  END OF GBIF DOWLOAD CODE ################################################## 
#########################################################################################################################