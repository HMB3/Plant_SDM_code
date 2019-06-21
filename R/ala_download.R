## path    = ALA_path
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

