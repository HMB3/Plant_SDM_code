download_GBIF_all_species = function (species_list, path) {
  
  ## create variables
  skip.spp.list       = list()
  GBIF.download.limit = 200000
  
  ## for every species in the list
  for(sp.n in species_list){
    
    ## 1). First, check if the f*&%$*# file exists
    file = paste0(path, sp.n, "_GBIF_records.RData")
    
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
    
    ## 4). Check how many records there are, and try to request them from GBIF if over 200k
    if (occ_search(scientificName = sp.n, limit = 1)$meta$count > GBIF.download.limit) {
      
      ## now append the species which had > 200k records to the skipped list
      print (paste ("Number of records > max for GBIF download via R (200,000)", sp.n))
      max =  paste ("Number of records > 200,000 |", sp.n)
      
      ## and send a request to GBIF for download
      message("Sending request to GBIF to download ", sp.n, " using rgbif :: occ_download")
      key  <- name_backbone(name = sp.n, rank = 'species')$usageKey
      GBIF = occ_download(paste('taxonKey = ', key),  user = "popple_1500", pwd = "Popple1500", email = "hugh.burley@mq.edu.au")
      save(GBIF, file = paste(path, sp.n, "_GBIF_request.RData", sep = ""))
      skip.spp.list <- c(skip.spp.list, max)
      
    } else {
      
      ## 5). Download ALL records from GBIF
      message("Downloading GBIF records for ", sp.n, " using rgbif :: occ_data")
      key  <- name_backbone(name = sp.n, rank = 'species')$usageKey

      GBIF <- occ_data(taxonKey = key, limit = GBIF.download.limit)
      GBIF <- as.data.frame(GBIF$data)
      cat("Synonyms returned for :: ", sp.n, unique(GBIF$scientificName), sep="\n")
      cat("Takonkeys returned for :: ", sp.n, unique(GBIF$taxonKey), sep="\n")
      message(dim(GBIF[1]), " Records returned for ", sp.n)
      
      
      ## 6). save records to .Rdata file
      save(GBIF, file = paste(path, sp.n, "_GBIF_records.RData", sep = ""))
      
    }
    
  }
  
}