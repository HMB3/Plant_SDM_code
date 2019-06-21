#########################################################################################################################
#######################################  COMBINE ALL SPECIES DATA FRAMES INTO ONE ####################################### 
#########################################################################################################################


#########################################################################################################################
## This code combines all downloaded records for each species into a single database.


#########################################################################################################################
## 1). CREATE LIST OF TAXA FROM DOWNLOADED FILES
#########################################################################################################################


#########################################################################################################################
## Read in all data to run the SDM code :: species lists, shapefile, rasters & tables
#source('./R/HIA_LIST_MATCHING.R')
spp.download = list.files("./data/base/HIA_LIST/ALA/TREE_SPECIES/", pattern = ".RData")



#########################################################################################################################
## 2). COMBINE SPECIES DATAFRAMES INTO ONE
#########################################################################################################################


## memory is a problem. So we need more RAM. Can these sorts of operations be run in parallel?
memory.limit()
gc()

 
#########################################################################################################################
## Combine all the taxa into a single dataframe at once
ALA.TREES <- spp.download %>%   ## spp.download[c(1:length(spp.download))] 

  ## Pipe the list into lapply
  lapply(function(x) {

    ## Create a character string of each .RData file
    f <- sprintf("./data/base/HIA_LIST/ALA/TREE_SPECIES/%s", x)

    ## Load each file
    d <- get(load(f))

    ## Now drop the columns which we don't need
    dat <- data.frame(searchTaxon = x, d[, colnames(d) %in% ALA.keep],
                      stringsAsFactors = FALSE)
    
    if(!is.character(dat$id)) {

      dat$id <- as.character(dat$id)

    }
    
    ## Need to print the object within the loop
    dat$coordinateUncertaintyInMetres = as.numeric(dat$coordinateUncertaintyInMetres)
    dat$year        = as.numeric(dat$year)
    dat$month       = as.numeric(dat$month)
    dat$searchTaxon = gsub("_ALA_records.RData", "", dat$searchTaxon)
    dat
    
  }) %>%

  ## Finally, bind all the rows together
  bind_rows


## Check the output :: how does this compare to John's 
sort(names(ALA.TREES))
dim(ALA.TREES)
length(unique(ALA.TREES$scientificName))
length(unique(ALA.TREES$searchTaxon))
sort(unique(ALA.TREES$scientificName))


## How can this be cleaned?
## How do the searched and returned items compare?
head(ALA.TREES, 100)[, c("scientificName",
                         "searchTaxon")]

tail(ALA.TREES, 100)[, c("scientificName",
                         "searchTaxon")]





#########################################################################################################################
## 2). CHECK TAXONOMY RETURNED BY ALA USING TAXONSTAND
######################################################################################################################### 


#########################################################################################################################
## Use "Taxonstand" to check the taxonomy. However, this also assumes that the ALA data is clean
ALA.TAXO <- TPL(unique(ALA.TREES$scientificName), infra = TRUE,
                corr = TRUE, repeats = 100)  ## to stop it timing out...
sort(names(ALA.TAXO))





#########################################################################################################################
#####################################################  END ############################################################## 
#########################################################################################################################