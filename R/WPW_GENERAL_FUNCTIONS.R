#########################################################################################################################
#############################################  FUNCTIONS FOR HORT AUS LIST ############################################## 
#########################################################################################################################


#########################################################################################################################
## DOWNLOADING FUNCTIONS
#########################################################################################################################


#########################################################################################################################
## SPECIES
#########################################################################################################################


#########################################################################################################################
## GBIF PARSE


# GBIF has recently made a bunch of handy tools available via their revamped API. These tools include a species name parser, 
# which seems very useful for cleaning long lists of taxon names.
# 
# Here’s a simple R function that takes a vector of taxon names and parses them using GBIF’s API, extracting, among other 
# details, the genus, species, infraspecific rank and epithet, nothorank (i.e., indicating the taxonomic rank of hybridisation), 
# and authorship


# gbif_parse <- function(x) {
#   # x is a vector of species names
#   library(RJSONIO)
#   library(RCurl)
#   library(plyr)
#   u <- "http://api.gbif.org/v1/parser/name"
#   res <- fromJSON(
#     postForm(u,
#              .opts = list(postfields = RJSONIO::toJSON(x),
#                           httpheader = c('Content-Type' = 'application/json')))
#   )
#   do.call(rbind.fill, lapply(res, as.data.frame))  
# }





#########################################################################################################################
## GBIF
download_GBIF_all_species = function (species_list, path) {
  
  ## create variables
  skip.spp.list       = list()
  GBIF.download.limit = 200000
  
  ## for every species in the list
  for(sp.n in species_list){
    
    ## 1). First, check if the f*&%$*# file exists
    ## data\base\HIA_LIST\GBIF\SPECIES
    file_name = paste0(path, sp.n, "_GBIF_records.RData")
    
    ## If it's already downloaded, skip
    if (file.exists (file_name)) {
      
      print(paste ("file exists for species", sp.n, "skipping"))
      next
      
    }
    #  create a dummy file
    dummy = data.frame()
    save (dummy, file = file_name)
    
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
    if (occ_search(scientificName = sp.n)$meta$count <= 2) {
      
      ## now append the species which had no records to the skipped list
      print (paste ("No GBIF records for", sp.n, "skipping"))
      records = paste ("No GBIF records |", sp.n)
      skip.spp.list <- c(skip.spp.list, records)
      next
      
    }
    
    ## 4). Check how many records there are, and skip if there are over 200k
    if (occ_search(scientificName = sp.n, limit = 1)$meta$count > GBIF.download.limit) {
      
      ## now append the species which had > 200k records to the skipped list
      print (paste ("Number of records > max for GBIF download via R (200,000)", sp.n))
      max =  paste ("Number of records > 200,000 |", sp.n)
      
      ## and send a request to GBIF for download
      # message("Sending request to GBIF to download ", sp.n, " using rgbif :: occ_download")
      # key  <- name_backbone(name = sp.n, rank = 'species')$usageKey
      # GBIF = occ_download(paste('taxonKey = ', key),  user = "popple_1500", pwd = "Popple1500", email = "hugh.burley@mq.edu.au")
      # save(GBIF, file = paste(path, sp.n, "_GBIF_request.RData", sep = ""))
      # skip.spp.list <- c(skip.spp.list, max)
      
    } else {
      
      ## 5). Download ALL records from GBIF
      message("Downloading GBIF records for ", sp.n, " using rgbif :: occ_data")
      key <- name_backbone(name = sp.n, rank = 'species')$usageKey
      # x = name_lookup(sp.n)
      # keys = x$data$key
      # 
      # GBIF <- keys %>%         
      #   
      #   ## pipe the list into lapply
      #   lapply(function(x) {
      #     
      #     ## Create the character string
      #     f <- occ_data(taxonKey = x, limit = GBIF.download.limit)
      #     f =  as.data.frame(f$data)
      #     ## Load each .RData file
      #     
      #   }) %>%
      #   
      #   ## Finally, bind all the rows together
      #   dplyr::bind_rows
      
      GBIF <- occ_data(taxonKey = key, limit = GBIF.download.limit)
      GBIF <- as.data.frame(GBIF$data)
      
      cat("Synonyms returned for :: ",  sp.n, unique(GBIF$scientificName), sep="\n")
      cat("Names returned for :: ", sp.n, unique(GBIF$name),               sep="\n")
      cat("Takonkeys returned for :: ", sp.n, unique(GBIF$taxonKey),       sep="\n")
      
      ## Could also only use the key searched, but that could knock out a lot of species
      #GBIF = GBIF[GBIF$taxonKey %in% key, ]
      #View(GBIF[c("name", "scientificName", "taxonKey")])
      
      message(dim(GBIF[1]), " Records returned for ", sp.n)
      
      ## 6). save records to .Rdata file, note that using .csv files seemed to cause problems...
      #save(GBIF, file = paste(path, sp.n, "_GBIF_records.RData", sep = ""))
      save(GBIF, file = file_name)
      
    }
    
  }
  
  return(skip.spp.list)
  
}


#########################################################################################################################
## ALA
## sp.n = "Polyalthia longifolia"
## path    = ALA_path
download_ALA_all_species = function (species_list, path) {
  
  ## create variables
  skip.spp.list       = list()
  #ALA.download.limit  = 200000
  
  ## for every species in the list
  for(sp.n in species_list){
    
    ## Get the ID?
    lsid <- ALA4R::specieslist(sp.n)$taxonConceptLsid
    
    ## 1). First, check if the f*&%$*# file exists
    file_name = paste0(path, sp.n, "_ALA_records.RData")
    
    ## If it's already downloaded, skip
    if (file.exists (file_name)) {
      
      print (paste ("file exists for species", sp.n, "skipping"))
      next
      
    }
    #  create a dummy file
    dummy = data.frame()
    save (dummy, file = file_name)
    
    ## 2). Then check the spelling...incorrect nomenclature will return NULL result
    ## Also, add a condtion to only get species with not that many synonyms?
    ## length(unique(ALA4R::occurrences(taxon = paste('taxon_name:\"', sp.n, '\"',sep=""), 
    ##                    download_reason_id = 7)$data$scientificName)) > 30
    if (is.null(ALA4R::occurrences(taxon = paste('taxon_name:\"', sp.n, '\"',sep=""), download_reason_id = 7)$data) == TRUE) {
      
      ## Now, append the species which had incorrect nomenclature to the skipped list
      ## this is slow, but it works for now.........................................
      print (paste ("Possible incorrect nomenclature", sp.n, "skipping"))
      nomenclature = paste ("Possible incorrect nomenclature |", sp.n)
      skip.spp.list <- c(skip.spp.list, nomenclature)
      next
      
    }
    
    ## 3). Skip species with no records
    if (nrow(ALA4R::occurrences(taxon = paste('taxon_name:\"', sp.n, '\"',sep=""), download_reason_id = 7)$data) <= 2) {
      
      ## now append the species which had no records to the skipped list
      print (paste ("No ALA records for", sp.n, "skipping"))
      records = paste ("No ALA records |", sp.n)
      skip.spp.list <- c(skip.spp.list, records)
      next
      
    }
    
    ## 4). Download ALL records from ALA :: 
    ## Try testing the searchtaxon for a few species 
    
    
    message("Downloading ALA records for ", sp.n, " using ALA4R :: occurrences")
    # taxon="taxon_name:\"Alaba vibex\""
    # paste('modelCheck(var"',i,'_d.bug")',sep="")
    # paste('taxon_name:\', sp.n, '\"', sep="")
    ALA = ALA4R::occurrences(taxon = paste('taxon_name:\"', sp.n, '\"',sep=""), download_reason_id = 7)   
    ALA = ALA[["data"]]
    
    cat("Synonyms returned for :: ", sp.n, unique(ALA$scientificName), sep="\n")
    message(dim(ALA[1]), " Records returned for ", sp.n)
    
    ## 5). save records to .Rdata file, note that using .csv files seemed to cause problems...
    save(ALA, file = file_name)
    
  }
  
}


#########################################################################################################################
## GENERA
#########################################################################################################################


## Function to download all genera
# download_GBIF_all_genera = function (list) {
#   
#   ## create variables
#   skip.gen.list       = list()
#   GBIF.download.limit = 200000
#   
#   ## for every unique genus in the list
#   for(gen.n in genera){
#     
#     
#     ## 1). First, check if the f*&%$*# file exists
#     file = paste0("./data/base/HIA_LIST/GBIF/GENERA/", gen.n, "_GBIF_records.RData")
#     
#     
#     ## If it's already downloaded, skip
#     if (file.exists (file)) {
#       
#       print (paste ("file exists for genera", gen.n, "skipping"))
#       next
#       
#     }
#     
#     ## 2). Then check the spelling...incorrect nomenclature will return NULL result
#     if (is.null(occ_search(scientificName = gen.n, limit = 1)$meta$count) == TRUE) {
#       
#       print (paste ("Possible incorrect nomenclature", gen.n, "skipping"))
#       nomenclature = paste ("Possible incorrect nomenclature |", gen.n)
#       skip.gen.list <- c(skip.gen.list, nomenclature)
#       next
#       
#     }
#     
#     ## 3). Skip species with no records
#     if (occ_search(scientificName = gen.n)$meta$count == 0) {
#       
#       print (paste ("No GBIF records for", gen.n, "skipping"))
#       records = paste ("No GBIF records for |", gen.n)
#       skip.gen.list <- c(skip.gen.list, records)
#       next
#       
#     }
#     
#     ## 4). Check how many records there are, and skip if there are over 200k
#     if (occ_search(scientificName = gen.n, limit = 1)$meta$count > GBIF.download.limit) {
#       
#       print (paste ("Number of records > max for GBIF download via R (200,000)", gen.n, "skipping"))
#       max =  paste ("Number of records > 200,000 |", gen.n)
#       skip.gen.list <- c(skip.gen.list, max)
#       next
#       
#     }
#     
#     ## 5). Download ALL records from GBIF
#     ## ala = occurrences(taxon = sp.n, download_reason_id = 7)
#     GBIF.GEN = gbif(gen.n, download = TRUE)   ## could use more arguments here, download_reason_id = 7, etc.
#     
#     ## 6). save records to .Rdata file, note that using .csv files seemed to cause problems...
#     save(GBIF.GEN, file = paste("./data/base/HIA_LIST/GBIF/GENERA/", gen.n, "_GBIF_records.RData", sep = ""))
#     return(skip.gen.list)
#     
#   }
#   
# }





#########################################################################################################################
## Package functions
#########################################################################################################################


# ipak <- function(pkg){
#   
#   new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
#   
#   if (length(new.pkg)) 
#     install.packages(new.pkg, dependencies = TRUE)
#   
#   sapply(pkg, require, character.only = TRUE)
#   
# }


#########################################################################################################################
## Sorting functions
#########################################################################################################################


## Trim trailin and leading white space
# trim <- function (x) gsub("^\\s+|\\s+$", "", x)
# 
# 
# ## Reorder column names
# movetolast <- function(data, move) {
#   
#   data[c(setdiff(names(data), move), move)]
#   
# }



## Get the first two words of a string
# string_fun_first_two_words <- function(x) {
#   
#   ul = unlist(strsplit(x, split = "\\s+"))[1:2]
#   paste(ul, collapse = " ") 
#   
# }


## Get the first word of a string
# string_fun_first_word <- function(x) {
#   
#   ul = unlist(strsplit(x, split = "\\s+"))[1]
#   paste(ul, collapse = " ") 
#   
# }


##
# make.true.NA <- function(x) if(is.character(x)||is.factor(x)){
#   is.na(x) <- x=="NA"; x} else {
#     x}


## Get a complete df
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
  
}


# round_df <- function(x, digits) {
#   # round all numeric variables
#   # x: data frame 
#   # digits: number of digits to round
#   numeric_columns <- sapply(x, mode) == 'numeric'
#   x[numeric_columns] <-  round(x[numeric_columns], digits)
#   x
# }


# ## Make the first letter of a string lower case
# firs_lett_lower <- function(x) {
#   substr(x, 1, 1) <- tolower(substr(x, 1, 1))
#   x
# }


# convert_ALA_columns <- function(ds) {
#   
#   
#   ## Make the first letter lowercase
#   firs_lett_lower <- function(x) {
#     substr(x, 1, 1) <- tolower(substr(x, 1, 1))
#     x
#   }
#   
#   ## 
#   names(ds) <- firs_lett_lower(names(ds))
#   
#   ## Convert any number of consecutive dots to a single space.
#   names(ds) <- gsub(x = names(ds),
#                     pattern = "(\\.)+",
#                     replacement = "")
#   
#   ## Drop the trailing spaces.
#   names(ds) <- gsub(x = names(ds),
#                     pattern = "( )+$",
#                     replacement = "")
#   
#   ## Remove dash
#   names(ds) <- gsub(x = names(ds),
#                     pattern = " - ",
#                     replacement = "")
#   
#   ## Remove parsed 
#   names(ds) <- gsub(x = names(ds),
#                     pattern = " - parsed",
#                     replacement = "")
#   
#   ## Remove white space.
#   names(ds) <- gsub(x = names(ds),
#                     pattern = " ",
#                     replacement = "")
#   
#   # ## for loop solution
#   # for(i in myobject){
#   #   mycolumns[grepl(i, mycolumns)] <- i
#   # }
#   
#   ds
# }





#########################################################################################################################
## CORRELATION PLOTTING FUNCTIONS
#########################################################################################################################


## if(missing(cex.cor)) cex = 0.8/strwidth(txt)
## text(0.5, 0.5, txt, cex = 5) 
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
  
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = 2)
}


panel.loess = function(x, y, digits=2, prefix="", cex.cor) {
  
  usr = par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  
  r = (loess(x, y))
  
  txt = format(c(r, 0.123456789), digits=digits)[1]
  txt = paste(prefix, txt, sep = "")
  
  if(missing(cex.cor)) cex = 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = 5)
  
}


panel.sig <- function(x, y, digits = 2, cex.cor)
  
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  test <- cor.test(x,y)
  Signif <- ifelse(round(test$p.value,3)<0.001,"p<0.001",paste("p=",round(test$p.value,3)))  
  text(0.5, 0.25, paste("r=",txt))
  text(.5, .75, Signif)
  
}


## tweak these arguments: lwd = 8, cex, etc.
panel.smooth <- function (x, y, col = "blue", bg = NA, pch = 19, 
                          cex = 1.2, col.smooth = "red",
                          cex.axis = 3,
                          span = 2/3, iter = 3, ...) 
  
{
  points(x, y, pch = pch, col = col, bg = bg, cex = cex, cex.axis = cex.axis)
  
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
          col = "orange", lwd = 4, cex.axis = 3, ...)
  
}


## tweak these arguments: border = NA, colour = grey, etc.
## cannot figure out how to get the axes of these plots to enlarge...
panel.hist <- function(x, ...)
  
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, #col = "grey", 
       border = "NA", ...)
  
}


# Mode <- function(x) {
#   
#   ux <- unique(x)
#   ux[which.max(tabulate(match(x, ux)))]
#   
# }





#########################################################################################################################
## COMBINE FUNCTIONS
#########################################################################################################################


# run_function_concatenate = function (list, DF, exp) {
#   
#   
#   COMBO.TABLE <- env.variables[c(1:length(list))] %>% 
#     
#     ## Pipe the list into lapply
#     lapply(function(x) {
#       
#       ## Now use the niche width function on each colname (so 8 environmental variables)
#       ## Also, need to figure out how to make the aggregating column generic (species, genus, etc.)
#       ## currently it only works hard-wired
#       func (exp)
#       
#       ## would be good to remove the duplicate columns here
#       
#     }) %>% 
#     
#     ## finally, create one dataframe for all niches
#     as.data.frame
#   
# }


## Read in a list of tables, bind them together and output one dataframe
# read_bind_tables = function (table.list, path) {
#   
#   READ.BIND.TABLE <- table.list[c(1:length(table.list))] %>% 
#     
#     ## pipe the list into lapply
#     lapply(function(x) {
#       
#       ## create the character string
#       f <- paste0(path, x)
#       
#       ## read each .csv file
#       dat <- read.csv(f, stringsAsFactors = FALSE)
#       
#       ## now drop the columns which we don't need
#       if(!is.character(dat$Catalog_Nu)) {
#         
#         dat$Catalog_Nu <- as.character(dat$Catalog_Nu)
#         
#       }
#       
#       dat
#       
#     }) %>%
#     
#     ## Finally, bind all the rows together
#     dplyr::bind_rows
#   
# }
# 
# 
# ## Loop over a list of subfolders
# read_bind_maxent = function (table.list, path) {
#   
#   READ.BIND.TABLE <- table.list[c(1:length(table.list))] %>% 
#     
#     ## pipe the list into lapply
#     lapply(function(x) {
#       
#       ## create the character string
#       f <- paste0(path, x, "/full/maxentResults.csv")
#       
#       ## read each .csv file
#       d <- read.csv(f)
#       
#       ## now add a model column
#       cbind(GBIF_Taxon = x,
#             Model_run  = path, 
#             d)
#       
#       ## Remove path gunk, and species
#       d$GBIF_Taxon = gsub("_", " ", d$GBIF_Taxon)
#       d$Model_run  = gsub("./output/maxent/", "", d$Model_run)
#       d$Model_run  = gsub("/", "", d$Model_run)
#       
#       
#     }) %>% 
#     
#     ## finally, bind all the rows together
#     dplyr::bind_rows
#   
# }



#########################################################################################################################
## GBIF FIELDS WE DON'T NEED
#########################################################################################################################


## Keep
gbifColsToDrop <- c(
  "crawlId",
  "disposition",
  "dynamicProperties",
  "elevationAccuracy",
  "endDayOfYear",
  "startDayOfYear",
  
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

## 
gbif.keep <- c(## TAXONOMY
  "searchTaxon",
  "species",
  "scientificName",
  "taxonRank",
  "taxonKey",
  "genus",
  "family",
  
  ## CULTIVATION
  "cloc",
  "basisOfRecord",
  "locality",
  "establishmentMeans",
  "institutionCode",
  "datasetName",
  "habitat",
  "eventRemarks",
  
  ## RECORD ID
  "recordedBy",
  "identifiedBy",
  "gbifID",
  "catalogNumber",
  
  ## PLACE/TIME
  "lat",
  "lon",
  "decimalLatitude",
  "decimalLongitude",
  "country",
  "coordinateUncertaintyInMeters",
  "geodeticDatum",
  "year",
  "month",
  "day",
  "eventDate",
  "eventID")


ALA.keep <- c(## TAXONOMY
  "searchTaxon",
  "scientificName",
  "scientificNameOriginal",
  "species",
  "taxonRank",
  "rank",
  "genus",
  "family",
  
  ## CULTIVATION
  "occCultivatedEscapee",
  "basisOfRecord",
  "locality",
  "establishmentMeans",
  "institutionCode",
  "datasetName",
  "habitat",
  "eventRemarks",
  "taxonomicQuality",
  
  ## RECORD ID
  "recordedBy",
  "id",
  "catalogNumber",
  "identificationID",
  "identifiedBy",
  "occurrenceID",
  "basisOfRecord",
  "institutionCode",
  
  ## PLACE/TIME
  "latitude",
  "longitude",
  "lat",
  "lon",
  "coordinateUncertaintyInMetres",
  "coordinateUncertaintyInMeters",
  "zeroCoordinates",
  "country",
  "state",
  "IBRA7Regions",
  "IBRA.7.Subregions",
  "localGovernmentAreas",
  "locality",
  "geodeticDatum",
  "year",
  "month",
  "day",
  "eventDate",
  "eventID")


TPL.keep <- c(## GBIF TAXONOMY
  "searchTaxon",
  "scientificName",
  "taxonRank",
  "genus",
  "family",
  
  ## TPL fields
  "Taxonomic.status",
  "Infraspecific.rank",
  "New.Taxonomic.status",
  "New.ID",
  "TPL_binomial",
  "taxo_agree",
  
  ## CULTIVATION
  "cloc",
  "basisOfRecord",
  "locality",
  "establishmentMeans",
  "institutionCode",
  "datasetName",
  "habitat",
  "eventRemarks",
  
  ## RECORD ID
  "recordedBy",
  "identifiedBy",
  "gbifID",
  "catalogNumber",
  
  ## PLACE/TIME
  "lat",
  "lon",
  "country",
  "coordinateUncertaintyInMeters",
  "geodeticDatum",
  "year",
  "month",
  "day",
  "eventDate",
  "eventID")





#########################################################################################################################
## NICHE BREADTH CALCULATIONS
#########################################################################################################################





#########################################################################################################################
## NICHE BREADTH CALCULATIONS
#########################################################################################################################


########################################################################################################################
## SHAWN'S NICHE
########################################################################################################################


## This is much simpler than Stuart's EG
niche_estimate = function (DF, 
                           colname) {
  
  ## Use ddply inside a function to create niche widths and medians for each species
  ## This syntax is tricky, maybe ask John and Stu what they think
  
  ## Also, need to figure out how to make the aggregating column generic (species, genus, etc.)
  summary = ddply(DF, 
                  .(searchTaxon),           ## currently grouping column only works hard-wired
                  .fun = function (xx, col) {
                    
                    ## All the different columns
                    min      = min(xx[[col]])
                    max      = max(xx[[col]])
                    
                    q02      = quantile(xx[[col]], .02)
                    q05      = quantile(xx[[col]], .05)
                    q95      = quantile(xx[[col]], .95)
                    q98      = quantile(xx[[col]], .98)
                    
                    median   = median(xx[[col]])
                    mean     = mean(xx[[col]])
                    mode     = Mode(xx[[col]])
                    range    = max - min
                    q95_q05  = (q95 - q05)
                    q98_q02  = (q98 - q02)
                    
                    ## Then crunch them together
                    c(min, max, median, mode, mean, range, q05, q95,  q95_q05, q98_q02)
                    
                  },
                  
                  colname
                  
  )
  
  ## Concatenate output
  ## Also, need to figure out how to make the aggregating column generic (species, genus, etc.)
  ## currently it only works hard-wired
  colnames(summary) = c("searchTaxon", 
                        paste0(colname,  "_min"),
                        paste0(colname,  "_max"),
                        paste0(colname,  "_median"),
                        paste0(colname,  "_mode"),
                        paste0(colname,  "_mean"),
                        paste0(colname,  "_range"),
                        paste0(colname,  "_q05"),
                        paste0(colname,  "_q95"),
                        paste0(colname,  "_q95_q05"),
                        paste0(colname,  "_q98_q02"))
  
  ## return the summary of niche width and median
  return (summary)
  
}






#########################################################################################################################
#############################################  FUNCTIONS FOR HORT AUS LIST ############################################## 
#########################################################################################################################