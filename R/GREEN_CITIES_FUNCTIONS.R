#########################################################################################################################
#############################################  FUNCTIONS FOR HORT AUS LIST ############################################## 
#########################################################################################################################


#########################################################################################################################
## sorting functions
#########################################################################################################################

string_fun <- function(x) {
  
  ul = unlist(strsplit(x, split = "\\s+"))[1:2]
  paste(ul, collapse = " ") 
  
}





#########################################################################################################################
## data cleaning functions
#########################################################################################################################


#########################################################################################################################
## GBIF data
## so here, x would be a "GBIF" object that comes from the previous loop
# load("./data/base/HIA_LIST/GBIF/Agonis flexuosa_GBIF_records.RData")
# x = GBIF


CLEAN_GBIF <- function(x) {
  
  
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
  basis   <- ifelse(x$basisOfRecord == 'FOSSIL_SPECIMEN',  'basis', '')                      ## length(x$basisOfRecord[x$basisOfRecord == "UNKNOWN"])
  

  
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
                       loc, 
                       #taxonid, 
                       basis, 
                       #bionet, 
                       #cultiv_esc, 
                       sep = ','),
                 
                 perl = TRUE)
  
  
  ## now exclude records that don't meet each condition as NA
  x <- x %>% 
    
    mutate(exclude = ifelse(reason = '', TRUE, FALSE),
           reason  = ifelse(exclude, reason, NA))
  
  ## return (x) at the end 
  return(x)
  
}





#########################################################################################################################
## ALA data
CLEAN_ALA <- function(x) {
  
  
  ## remove records that are older than 1950
  yr      <- ifelse(is.na(x$year) | x$year < 1950, 'year', '')
  coords1 <- is.na(x$longitude) | is.na(x$latitude) | x$longitude == '' | x$latitude == ''
  
  
  ## remove records with no coordinates, with inverted coordinates, and with a third condition
  coords2 <- if('zeroLatitude' %in% names(x)) x$zeroLatitude else FALSE
  coords3 <- if('zeroLongitude' %in% names(x)) x$zeroLongitude else FALSE
  coords4 <- if('invertedCoordinates' %in% names(x)) x$invertedCoordinates else FALSE # could just flip them back though...
  coords  <- ifelse(coords1|coords2|coords3|coords4, 'coords', '')
  
  
  ## remove records with spatial accuracy < 1km
  coorduncert <- ifelse(!((x$coordinateUncertaintyInMetres <= 1000) | is.na(x$coordinateUncertaintyInMetres)), 'coorduncert', '')
  
  
  ## remove records from zoos (only needed for animals...)
  #loc <- ifelse(grepl('\\bzoo\\b', x$locality, ignore.case = TRUE), 'locality', '')
  
  
  ## remove records with taxonomic problems, or that are fossils
  taxonid <- ifelse(x$taxonIdentificationIssue == 'questionSpecies', 'taxonid', '')
  basis   <- ifelse(x$basisOfRecord == 'FossilSpecimen', 'basis', '')
  
  
  ## remove records that don't meet an AUS collection code criteria (can't use this on GBIF)
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
                       loc, 
                       taxonid, 
                       basis, 
                       bionet, 
                       cultiv_esc, 
                       sep=','),
                 
                 perl = TRUE)
  
  
  ## now exclude records that don't meet each condition as NA
  x <- x %>% 
    mutate(exclude = ifelse(reason = '', TRUE, FALSE),
           reason  = ifelse(exclude, reason, NA))
  
  ## return x at the end 
  x
  
}





#########################################################################################################################
#############################################  FUNCTIONS FOR HORT AUS LIST ############################################## 
#########################################################################################################################