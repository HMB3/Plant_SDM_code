#########################################################################################################################
#############################################  FUNCTIONS FOR HORT AUS LIST ############################################## 
#########################################################################################################################


#########################################################################################################################
## Package functions
#########################################################################################################################


ipak <- function(pkg){
  
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  
  sapply(pkg, require, character.only = TRUE)
  
}


#########################################################################################################################
## sorting functions
#########################################################################################################################

string_fun_first_two_words <- function(x) {
  
  ul = unlist(strsplit(x, split = "\\s+"))[1:2]
  paste(ul, collapse = " ") 
  
}


string_fun_first_word <- function(x) {
  
  ul = unlist(strsplit(x, split = "\\s+"))[1]
  paste(ul, collapse = " ") 
  
}





#########################################################################################################################
## downloading functions
#########################################################################################################################


#########################################################################################################################
## SPECIES
#########################################################################################################################


download_GBIF_all_species = function (list) {
  
  ## create variables
  skip.spp.list       = list()
  GBIF.download.limit = 200000
  
  ## for every species in the list
  for(sp.n in spp){
    
    ## 1). First, check if the f*&%$*# file exists
    file = paste0("./data/base/HIA_LIST/GBIF/", sp.n, "_GBIF_records.RData")
    
    ## If it's already downloaded, skip
    if (file.exists (file)) {
      
      print (paste ("file exists for species", sp.n, "skipping"))
      next
      
    }
    
    ## 2). Then check the spelling...incorrect nomenclature will return NULL result
    if (is.null(occ_search(scientificName = sp.n, limit = 1)$meta$count) == TRUE) {
      
      ## now append the species which had incorrect nomenclature to the skipped list
      ## this is slow, but it works for now
      print (paste ("Possible incorrect nomenclature", sp.n, "skipping"))
      nomenclature = paste ("Possible incorrect nomenclature |", sp.n)
      skip.spp.list <- c(skip.spp.list, nomenclature)
      next
      
    }
    
    ## 3). Skip species with no records
    if (occ_search(scientificName = sp.n)$meta$count == 0) {
      
      ## now append the species which had no records to the skipped list
      print (paste ("No GBIF records for", sp.n, "skipping"))
      records = paste ("No GBIF records |", sp.n)
      skip.spp.list <- c(skip.spp.list, records)
      next
      
    }
    
    ## 4). Check how many records there are, and skip if there are over 200k
    if (occ_search(scientificName = sp.n, limit = 1)$meta$count > GBIF.download.limit) {
      
      ## now append the species which had > 200k records to the skipped list
      print (paste ("Number of records > max for GBIF download via R (200,000)", sp.n, "skipping"))
      max =  paste ("Number of records > 200,000 |", sp.n)
      skip.spp.list <- c(skip.spp.list, max)
      next
      
    }
    
    ## 5). Download ALL records from GBIF
    ## ala = occurrences(taxon = sp.n, download_reason_id = 7)
    GBIF = gbif(sp.n, download = TRUE)   ## could use more arguments here, download_reason_id = 7, etc.
    
    ## 6). save records to .Rdata file, note that using .csv files seemed to cause problems...
    save(GBIF, file = paste("./data/base/HIA_LIST/GBIF/", sp.n, "_GBIF_records.RData", sep = ""))
    return(skip.spp.list)
    
  }
  
}





#########################################################################################################################
## GENERA
#########################################################################################################################


## Function to download all genera
download_GBIF_all_genera = function (list) {
  
  ## create variables
  skip.gen.list       = list()
  GBIF.download.limit = 200000
  
  ## for every unique genus in the list
  for(gen.n in genera){
    
    
    ## 1). First, check if the f*&%$*# file exists
    file = paste0("./data/base/HIA_LIST/GBIF/GENERA/", gen.n, "_GBIF_records.RData")
    
    
    ## If it's already downloaded, skip
    if (file.exists (file)) {
      
      print (paste ("file exists for genera", gen.n, "skipping"))
      next
      
    }
    
    ## 2). Then check the spelling...incorrect nomenclature will return NULL result
    if (is.null(occ_search(scientificName = gen.n, limit = 1)$meta$count) == TRUE) {
      
      print (paste ("Possible incorrect nomenclature", gen.n, "skipping"))
      nomenclature = paste ("Possible incorrect nomenclature |", gen.n)
      skip.gen.list <- c(skip.gen.list, nomenclature)
      next
      
    }
    
    ## 3). Skip species with no records
    if (occ_search(scientificName = gen.n)$meta$count == 0) {
      
      print (paste ("No GBIF records for", gen.n, "skipping"))
      records = paste ("No GBIF records for |", gen.n)
      skip.gen.list <- c(skip.gen.list, records)
      next
      
    }
    
    ## 4). Check how many records there are, and skip if there are over 200k
    if (occ_search(scientificName = gen.n, limit = 1)$meta$count > GBIF.download.limit) {
      
      print (paste ("Number of records > max for GBIF download via R (200,000)", gen.n, "skipping"))
      max =  paste ("Number of records > 200,000 |", gen.n)
      skip.gen.list <- c(skip.gen.list, max)
      next
      
    }
    
    ## 5). Download ALL records from GBIF
    ## ala = occurrences(taxon = sp.n, download_reason_id = 7)
    GBIF.GEN = gbif(gen.n, download = TRUE)   ## could use more arguments here, download_reason_id = 7, etc.
    
    ## 6). save records to .Rdata file, note that using .csv files seemed to cause problems...
    save(GBIF.GEN, file = paste("./data/base/HIA_LIST/GBIF/GENERA/", gen.n, "_GBIF_records.RData", sep = ""))
    return(skip.gen.list)
    
  }
  
}





#########################################################################################################################
## GBIF fields we don't need
#########################################################################################################################


gbifColsToDrop <- c("cloc",
                    "crawlId",
                    "disposition",
                    "dynamicProperties",
                    "elevationAccuracy",
                    "endDayOfYear",
                    "startDayOfYear",
                    "eventRemarks",
                    "eventTime",
                    "eventID",
                    "month",
                    
                    "class",
                    "kingdom",
                    "phylum",
                    
                    "familyKey",
                    "genusKey",
                    "kingdomKey",
                    "phylumKey",
                    "publishingOrgKey",
                    "taxonomicStatus",
                    "lowestBiostratigraphicZone",
                    "reproductiveCondition",
                    "parentNameUsage",
                    "preparations",
                    "specificEpithet",
                    
                    "coordinatePrecision",
                    "georeferenceVerificationStatus",
                    "eventID",
                    "lastCrawled",
                    "higherGeography",
                    "municipality",
                    "verbatimCoordinateSystem",
                    "verbatimLocality",
                    "georeferenceRemarks",
                    "georeferencedBy",
                    "georeferencedDate",
                    "georeferenceProtocol",
                    "georeferenceSources",
                    "higherGeographyID",
                    "http://unknown.org/occurrenceDetails",
                    "lowestBiostratigraphicZone",
                    "highestBiostratigraphicZone",
                    "waterBody",
                    "verbatimSRS",
                    "depth",
                    "depthAccuracy",
                    "disposition",
                    
                    "informationWithheld",
                    "language",
                    "protocol",
                    "rightsHolder",
                    "publishingCountry",
                    "accessRights",
                    "associatedSequences",
                    "dynamicProperties",
                    "earliestEpochOrLowestSeries",
                    "earliestPeriodOrLowestSystem",
                    
                    "http...unknown.org.classs", 
                    "identificationID", 
                    "identificationQualifier", 
                    "identificationRemarks",                   
                    "identifiedBy", 
                    "identifier", 
                    "individualCount",
                    
                    "namePublishedInYear",  
                    "http...unknown.org.organismQuantity", 
                    "http...unknown.org.organismQuantityType", 
                    "http...unknown.org.sampleSizeUnit",       
                    "http...unknown.org.sampleSizeValue", 
                    "pointRadiusSpatialFit", 
                    "samplingEffort", 
                    "taxonConceptID",  
                    "nameAccordingToID", 
                    "parentEventID", 
                    "earliestEraOrLowestErathem", 
                    "formation",                               
                    "http...unknown.org.layer", 
                    "latestEpochOrHighestSeries", 
                    "latestPeriodOrHighestSystem", 
                    "sampleSizeUnit",
                    
                    "sampleSizeValue",                          
                    "sex",                                      
                    "footprintSpatialFit",                      
                    "group",                                   
                    "latestEraOrHighestErathem",                
                    "associatedOccurrences",                    
                    "identificationReferences",                 
                    "http...unknown.org.organismID",  
                    
                    "earliestEonOrLowestEonothem",              
                    "latestEonOrHighestEonothem",               
                    "verbatimLabel",                            
                    "verbatimDepth",                           
                    "furtherInformationURL",                    
                    "materialSampleID",                         
                    "parentNameUsageID",                        
                    "earliestAgeOrLowestStage", 
                    
                    "latestAgeOrHighestStage",                  
                    "member",                                   
                    "acceptedNameUsageID",                      
                    "bed",                                     
                    "lithostratigraphicTerms",                  
                    "behavior",                                 
                    "http...unknown.org.furtherInformationURL", 
                    "minimumDistanceAboveSurfaceInMeters",
                    "organismScope")





#########################################################################################################################
## NICHE BREADTH CALCULATIONS
#########################################################################################################################


## 
nicheBreadth <- function(data, species, r, verbose = TRUE) {

  # handle Inf/-Inf better! Output something more meaningful
  # remove NA before stats, so don't need na.rm everywhere
  # - if species not specified then look for taxa list in data
  # - do not presume column names
  
  ## transform raster simply to get it into memory
  r <- t(t(r))
  
  ## get data columns of interest
  data <- data[, c("species", "longitude_raw", "latitude_raw")]
  
  ## use data.table package for improved speed
  myDT <- data.table(data)
  setkey(myDT, species)
  
  ## setup output data frame ()
  out <- data.frame(matrix(ncol = 11, nrow = length(species)))
  
  colnames(out) <- c("species", "n", "min", "max", "median",
                     "perc02", "perc05", "perc95", "perc98", 
                     "breadth", "breadth95", "breadth98")
  
  ## for all species
  for (i in seq_len(length(species))) {
    
    speciesName <- species[i]
    
    if (verbose) {
      
      message(paste(i, ':', speciesName))
      
    }
    
    ## get unique occurrence data for species
    speciesData <- data.frame(myDT[speciesName, ])
    pts         <- unique(speciesData[, 2:3])
    
    ## extract raster values at each point
    x <- extract(r, pts)
    
    ## get niche stats for species
    out[i, "species"] <- speciesName
    out[i, "n"]       <- length(x)
    out[i, "min"]     <- min(x, na.rm = TRUE)
    out[i, "max"]     <- max(x, na.rm = TRUE)
    out[i, "median"]  <- median(x, na.rm = TRUE)
    out[i, "perc02"]  <- quantile(x, 0.02, names = FALSE, na.rm = TRUE)
    out[i, "perc05"]  <- quantile(x, 0.05, names = FALSE, na.rm = TRUE)
    out[i, "perc95"]  <- quantile(x, 0.95, names = FALSE, na.rm = TRUE)
    out[i, "perc98"]  <- quantile(x, 0.98, names = FALSE, na.rm = TRUE)
    
  }
  
  ## compute niche breadths
  out$breadth   <- out$max - out$min
  out$breadth95 <- out$perc95 - out$perc05
  out$breadth98 <- out$perc98 - out$perc02
  
  ## this returns "out"?
  return(out)
  
}





#########################################################################################################################
## GBIF data
## so here, x would be a "GBIF" object that comes from the previous loop
# load("./data/base/HIA_LIST/GBIF/Agonis flexuosa_GBIF_records.RData")
# x = GBIF
## run the function as 
## occ_data_clean = CLEAN_ALA(GBIF)



CLEAN_GBIF_MACQU <- function(x) {
  
  
  ## for this function to work, every record from GBIF needs to have the same column names... 
  ## need to consider how to deal with this problem for each taxa
  
  ## remove records that are older than 1950
  yr        <- ifelse(is.na(x$year) | x$year < 1950, 'year', '')           ## unique(yr) length(yr)
  coords.na <- is.na(x$lon) | is.na(x$lat) | x$lon == '' | x$lat == ''     ## unique(coords.na) length(coords.na) 
  # summary(x$lon)
  # length(coords.na[coords.na==TRUE])
  
  
  ## remove records with no coordinates, with inverted coordinates, and with a third condition
  ## check 'Latitude' is 'lat', etc
  coords.0.lat  <- if('zeroLatitude' %in% names(x)) x$zeroLatitude else FALSE                ## unique(coords.0.lat)  length(coords.0.lat)
  coords.0.lon  <- if('zeroLongitude' %in% names(x)) x$zeroLongitude else FALSE              ## unique(coords.0.lon)
  coords.inv    <- if('invertedCoordinates' %in% names(x)) x$invertedCoordinates else FALSE  ## unique(coords.inv)
  coords        <- ifelse(coords.na|coords.0.lat|coords.0.lon|coords.inv, 'coords', '')      ## unique(coords)  length(coords)
  # length(coords[coords == ''])
  
  
  ## remove records with spatial accuracy < 1km ## summary(x$coordinateUncertaintyInMeters)
  coorduncert <- ifelse(!((x$coordinateUncertaintyInMeters <= 1000) | is.na(x$coordinateUncertaintyInMeters)), 'coorduncert', '')
  ## length(coorduncert[coorduncert == ''])
  
  
  ## remove records from zoos (only needed for animals...)
  #loc <- ifelse(grepl('\\bzoo\\b', x$locality, ignore.case = TRUE), 'locality', '')
  
  ## remove records with taxonomic problems, or that are fossils
  #taxonid <- ifelse(x$taxonIdentificationIssue == 'questionSpecies', 'taxonid', '')
  basis   <- ifelse(x$basisOfRecord == 'FOSSIL_SPECIMEN',  'basis', '')   ## length(x$basisOfRecord[x$basisOfRecord == "UNKNOWN"])
  
  
  
  ## remove records that don't meet an AUS collection code criteria...can't use this on GBIF?
  # bionet <- ifelse(x$collectionCode %in% c(
  #   'BioNet Atlas of NSW Wildlife', 
  #   'NSW Office of Environment and Heritage BioNet Atlas of NSW Wildlife'), 
  #   'bionet', '')
  
  
  ## remove cultivated records: again, we want these?
  # cultiv_esc <- if('occCultivatedEscapee' %in% names(x)) {
  #   ifelse(x$occCultivatedEscapee, 'cultiv_esc', '') 
  #   
  # } else ''
  
  
  ## the list of reasons for excluding records...
  reason <- gsub('^,+|(?<=,),+|,+$',
                 '',
                 paste(yr, 
                       coords, 
                       coorduncert, 
                       #loc, 
                       #taxonid, 
                       basis, 
                       #bionet, 
                       #cultiv_esc, 
                       sep = ','),
                 
                 perl = TRUE)
  
  
  ## now exclude records that don't meet each condition as NA
  x <- x %>% 
    
    mutate(exclude = ifelse(reason !='', TRUE, FALSE),
           reason  = ifelse(exclude, reason, NA))
  
  ## return (x) at the end 
  x
  #return(x)
  
}





#########################################################################################################################
#############################################  FUNCTIONS FOR HORT AUS LIST ############################################## 
#########################################################################################################################