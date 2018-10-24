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
## Sorting functions
#########################################################################################################################


## Trim trailin and leading white space
trim <- function (x) gsub("^\\s+|\\s+$", "", x)


## Reorder column names
movetolast <- function(data, move) {
  
  data[c(setdiff(names(data), move), move)]
  
}



## Get the first two words of a string
string_fun_first_two_words <- function(x) {
  
  ul = unlist(strsplit(x, split = "\\s+"))[1:2]
  paste(ul, collapse = " ") 
  
}


## Get the first word of a string
string_fun_first_word <- function(x) {
  
  ul = unlist(strsplit(x, split = "\\s+"))[1]
  paste(ul, collapse = " ") 
  
}


completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}



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


Mode <- function(x) {
  
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
  
}





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


gbif_parse <- function(x) {
  # x is a vector of species names
  library(RJSONIO)
  library(RCurl)
  library(plyr)
  u <- "http://api.gbif.org/v1/parser/name"
  res <- fromJSON(
    postForm(u,
             .opts = list(postfields = RJSONIO::toJSON(x),
                          httpheader = c('Content-Type' = 'application/json')))
  )
  do.call(rbind.fill, lapply(res, as.data.frame))  
}


#########################################################################################################################
## GBIF
# download_GBIF_all_species = function (species_list, path) {
#   
#   ## create variables
#   skip.spp.list       = list()
#   GBIF.download.limit = 200000
#   
#   ## for every species in the list
#   for(sp.n in species_list){
#     
#     ## 1). First, check if the f*&%$*# file exists
#     ## data\base\HIA_LIST\GBIF\SPECIES
#     file = paste0(path, sp.n, "_GBIF_records.RData")
#     
#     ## If it's already downloaded, skip
#     if (file.exists (file)) {
#       
#       print (paste ("file exists for species", sp.n, "skipping"))
#       next
#       
#     }
#     
#     ## 2). Then check the spelling...incorrect nomenclature will return NULL result
#     if (is.null(occ_search(scientificName = sp.n, limit = 1)$meta$count) == TRUE) {
#       
#       ## now append the species which had incorrect nomenclature to the skipped list
#       ## this is slow, but it works for now
#       print (paste ("Possible incorrect nomenclature", sp.n, "skipping"))
#       nomenclature = paste ("Possible incorrect nomenclature |", sp.n)
#       skip.spp.list <- c(skip.spp.list, nomenclature)
#       next
#       
#     }
#     
#     ## 3). Skip species with no records
#     if (occ_search(scientificName = sp.n)$meta$count == 0) {
#       
#       ## now append the species which had no records to the skipped list
#       print (paste ("No GBIF records for", sp.n, "skipping"))
#       records = paste ("No GBIF records |", sp.n)
#       skip.spp.list <- c(skip.spp.list, records)
#       next
#       
#     }
#     
#     ## 4). Check how many records there are, and skip if there are over 200k
#     if (occ_search(scientificName = sp.n, limit = 1)$meta$count > GBIF.download.limit) {
#       
#       ## now append the species which had > 200k records to the skipped list
#       print (paste ("Number of records > max for GBIF download via R (200,000)", sp.n, "skipping"))
#       max =  paste ("Number of records > 200,000 |", sp.n)
#       skip.spp.list <- c(skip.spp.list, max)
#       next
#       
#     }
#     
#     ## 5). Download ALL records from GBIF
#     ## ala = occurrences(taxon = sp.n, download_reason_id = 7)
#     print (paste (sp.n))
#     GBIF = gbif(sp.n, download = TRUE)   ## could use more arguments here, download_reason_id = 7, etc.
#     
#     ## 6). save records to .Rdata file, note that using .csv files seemed to cause problems...
#     save(GBIF, file = paste(path, sp.n, "_GBIF_records.RData", sep = ""))
#     #return(skip.spp.list)
#     
#   }
#   
# }


download_GBIF_all_species = function (species_list, path) {
  
  ## create variables
  skip.spp.list       = list()
  GBIF.download.limit = 200000
  
  ## for every species in the list
  for(sp.n in species_list){
    
    ## 1). First, check if the f*&%$*# file exists
    ## data\base\HIA_LIST\GBIF\SPECIES
    file = paste0(path, sp.n, "_GBIF_records.RData")
    
    ## If it's already downloaded, skip
    if (file.exists (file)) {
      
      print(paste ("file exists for species", sp.n, "skipping"))
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
      #   bind_rows
      
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
      save(GBIF, file = paste(path, sp.n, "_GBIF_records.RData", sep = ""))
      #return(skip.spp.list)
      
    }
    
  }
  
}


#########################################################################################################################
## ALA
download_ALA_all_species = function (species_list, path) {
  
  ## create variables
  skip.spp.list       = list()
  #ALA.download.limit  = 200000
  
  ## for every species in the list
  for(sp.n in species_list){
    
    ## 1). First, check if the f*&%$*# file exists
    file = paste0(path, sp.n, "_ALA_records.RData")
    
    ## If it's already downloaded, skip
    if (file.exists (file)) {
      
      print (paste ("file exists for species", sp.n, "skipping"))
      next
      
    }
    
    ## 2). Then check the spelling...incorrect nomenclature will return NULL result
    if (is.null(occurrences(taxon = sp.n, download_reason_id = 7)$data) == TRUE) {
      
      ## now append the species which had incorrect nomenclature to the skipped list
      ## this is slow, but it works for now
      print (paste ("Possible incorrect nomenclature", sp.n, "skipping"))
      nomenclature = paste ("Possible incorrect nomenclature |", sp.n)
      skip.spp.list <- c(skip.spp.list, nomenclature)
      next
      
    }
    
    ## 3). Download ALL records from ALA :: 
    message("Downloading ALA records for ", sp.n, " using ALA4R :: occurrences")
    ALA = occurrences(taxon = sp.n, download_reason_id = 7)   ## could use more arguments here, download_reason_id = 7, etc.
    ALA = ALA[["data"]]
    cat("Synonyms returned for :: ", sp.n, unique(ALA$scientificName), sep="\n")
    message(dim(ALA[1]), " Records returned for ", sp.n)
    
    ## 4). save records to .Rdata file, note that using .csv files seemed to cause problems...
    save(ALA, file = paste(path, sp.n, "_ALA_records.RData", sep = ""))
    #return(skip.spp.list)
    
  }
  
}





# ## remaining species
# download_GBIF_setdiff_species = function (list) {
#   
#   ## create variables
#   skip.spp.list       = list()
#   GBIF.download.limit = 200000
#   
#   ## for every species in the list
#   for(sp.n in list){
#     
#     ## 1). First, check if the f*&%$*# file exists
#     file = paste0("./data/base/HIA_LIST/GBIF/SETDIFF/", sp.n, "_GBIF_records.RData")
#     
#     ## If it's already downloaded, skip
#     if (file.exists (file)) {
#       
#       print (paste ("file exists for species", sp.n, "skipping"))
#       next
#       
#     }
#     
#     ## 2). Then check the spelling...incorrect nomenclature will return NULL result
#     if (is.null(occ_search(scientificName = sp.n, limit = 1)$meta$count) == TRUE) {
#       
#       ## now append the species which had incorrect nomenclature to the skipped list
#       ## this is slow, but it works for now
#       print (paste ("Possible incorrect nomenclature", sp.n, "skipping"))
#       nomenclature = paste ("Possible incorrect nomenclature |", sp.n)
#       skip.spp.list <- c(skip.spp.list, nomenclature)
#       next
#       
#     }
#     
#     ## 3). Skip species with no records
#     if (occ_search(scientificName = sp.n)$meta$count == 0) {
#       
#       ## now append the species which had no records to the skipped list
#       print (paste ("No GBIF records for", sp.n, "skipping"))
#       records = paste ("No GBIF records |", sp.n)
#       skip.spp.list <- c(skip.spp.list, records)
#       next
#       
#     }
#     
#     ## 4). Check how many records there are, and skip if there are over 200k
#     if (occ_search(scientificName = sp.n, limit = 1)$meta$count > GBIF.download.limit) {
#       
#       ## now append the species which had > 200k records to the skipped list
#       print (paste ("Number of records > max for GBIF download via R (200,000)", sp.n, "skipping"))
#       max =  paste ("Number of records > 200,000 |", sp.n)
#       skip.spp.list <- c(skip.spp.list, max)
#       next
#       
#     }
#     
#     ## 5). Download ALL records from GBIF
#     ## ala = occurrences(taxon = sp.n, download_reason_id = 7)
#     print (paste (sp.n))
#     GBIF = gbif(sp.n, download = TRUE)   ## could use more arguments here, download_reason_id = 7, etc.
#     
#     ## 6). save records to .Rdata file, note that using .csv files seemed to cause problems...
#     save(GBIF, file = paste("./data/base/HIA_LIST/GBIF/SETDIFF/", sp.n, "_GBIF_records.RData", sep = ""))
#     #return(skip.spp.list)
#     
#   }
#   
# }
# 
# 
# 
# 
# 
# ## remaining species
# download_GBIF_Renee_species = function (list) {
#   
#   ## create variables
#   skip.spp.list       = list()
#   GBIF.download.limit = 200000
#   
#   ## for every species in the list
#   for(sp.n in list){
#     
#     ## 1). First, check if the f*&%$*# file exists
#     file = paste0("./data/base/HIA_LIST/GBIF/RENEE/", sp.n, "_GBIF_records.RData")
#     
#     ## If it's already downloaded, skip
#     if (file.exists (file)) {
#       
#       print (paste ("file exists for species", sp.n, "skipping"))
#       next
#       
#     }
#     
#     ## 2). Then check the spelling...incorrect nomenclature will return NULL result
#     if (is.null(occ_search(scientificName = sp.n, limit = 1)$meta$count) == TRUE) {
#       
#       ## now append the species which had incorrect nomenclature to the skipped list
#       ## this is slow, but it works for now
#       print (paste ("Possible incorrect nomenclature", sp.n, "skipping"))
#       nomenclature = paste ("Possible incorrect nomenclature |", sp.n)
#       skip.spp.list <- c(skip.spp.list, nomenclature)
#       next
#       
#     }
#     
#     ## 3). Skip species with no records
#     if (occ_search(scientificName = sp.n)$meta$count == 0) {
#       
#       ## now append the species which had no records to the skipped list
#       print (paste ("No GBIF records for", sp.n, "skipping"))
#       records = paste ("No GBIF records |", sp.n)
#       skip.spp.list <- c(skip.spp.list, records)
#       next
#       
#     }
#     
#     ## 4). Check how many records there are, and skip if there are over 200k
#     if (occ_search(scientificName = sp.n, limit = 1)$meta$count > GBIF.download.limit) {
#       
#       ## now append the species which had > 200k records to the skipped list
#       print (paste ("Number of records > max for GBIF download via R (200,000)", sp.n, "skipping"))
#       max =  paste ("Number of records > 200,000 |", sp.n)
#       skip.spp.list <- c(skip.spp.list, max)
#       next
#       
#     }
#     
#     ## 5). Download ALL records from GBIF
#     ## ala = occurrences(taxon = sp.n, download_reason_id = 7)
#     print (paste (sp.n))
#     GBIF = gbif(sp.n, download = TRUE)   ## could use more arguments here, download_reason_id = 7, etc.
#     
#     ## 6). save records to .Rdata file, note that using .csv files seemed to cause problems...
#     save(GBIF, file = paste("./data/base/HIA_LIST/GBIF/RENEE/", sp.n, "_GBIF_records.RData", sep = ""))
#     #return(skip.spp.list)
#     
#   }
#   
# }




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
## COMBINE FUNCTIONS
#########################################################################################################################


run_function_concatenate = function (list, DF, exp) {


COMBO.TABLE <- env.variables[c(1:length(list))] %>% 
  
  ## Pipe the list into lapply
  lapply(function(x) {
    
    ## Now use the niche width function on each colname (so 8 environmental variables)
    ## Also, need to figure out how to make the aggregating column generic (species, genus, etc.)
    ## currently it only works hard-wired
    func (exp)
    
    ## would be good to remove the duplicate columns here
    
  }) %>% 
  
  ## finally, create one dataframe for all niches
  as.data.frame

}


## Read in a list of tables, bind them together and output one dataframe
read_bind_tables = function (table.list, path) {
  
  READ.BIND.TABLE <- table.list[c(1:length(table.list))] %>% 
    
    ## pipe the list into lapply
    lapply(function(x) {
      
      ## create the character string
      f <- paste0(path, x)
      
      ## read each .csv file
      dat <- read.csv(f, stringsAsFactors = FALSE)
      
      ## now drop the columns which we don't need
      if(!is.character(dat$Catalog_Nu)) {
        
        dat$Catalog_Nu <- as.character(dat$Catalog_Nu)
        
      }
      
      dat
      
    }) %>%
    
    ## Finally, bind all the rows together
    bind_rows
  
}


## Loop over a list of subfolders
read_bind_maxent = function (table.list, path) {
  
  READ.BIND.TABLE <- table.list[c(1:length(table.list))] %>% 
    
    ## pipe the list into lapply
    lapply(function(x) {
      
      ## create the character string
      f <- paste0(path, x, "/full/maxentResults.csv")
      
      ## read each .csv file
      d <- read.csv(f)
      
      ## now add a model column
      cbind(GBIF_Taxon = x,
            Model_run  = path, 
            d)
      
      ## Remove path gunk, and species
      d$GBIF_Taxon = gsub("_", " ", d$GBIF_Taxon)
      d$Model_run  = gsub("./output/maxent/", "", d$Model_run)
      d$Model_run  = gsub("/", "", d$Model_run)
      
      
    }) %>% 
    
    ## finally, bind all the rows together
    bind_rows
  
}


# dfl <- list(df1,df2)
# 
# # This would generate your error
# # df <- bind_rows(dfl)
# 
# # This would solve it while keeping the factors
# df <- dfl %>%
#   lapply(function(x) mutate_each(x, funs('as.character'))) %>%
#   bind_rows() %>% 
#   mutate_each(funs('as.factor'))

# for(i in LGA.list){
#   f <- paste0("./data/base/HIA_LIST/LGA/", i)
#   assign(i, read.csv(f, stringsAsFactors = FALSE))
#   
# }


# list = function (table.list, path) {
#   
#   for(i in LGA.list){
#     f <- paste0("./data/base/HIA_LIST/LGA/", i)
#     assign(i, read.csv(f, stringsAsFactors = FALSE))
#     
#   }
#   
#   return(list)
# }  



# ## Read in a list of tables, bind them together and output one dataframe
# read_LGA_list = function (table.list, path) {
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
#       d <- read.csv(f, stringsAsFactors = FALSE)
#       
#       ## now drop the columns which we don't need
#       assign(i, read.delim(filepath,
#                            colClasses=c("character","factor",rep("numeric",4)),
#                            sep = "\t"))
#       
#     })
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
  "identificationID",
  "identifiedBy",
  "occurrenceID",
  "catalogNumber",
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


# cultivated.synonyms <- "garden|afgrøde|bebauen|bebouwen|beskärs|çiftçilik|crescătorie|cultiv|Ernte|frø|gefokt|gepropageerd|însămânțat|kırpılmış|kylvetään|odlad|oogsten|opdrættede|propagate|propagé|propagiert|récolt|seribaşı|viljellyn|viljelykasvi|yayılır|yönetilen|выращиваемых|размножают|allevat|ausgesät|avl|avlet|beskjæres|bewerkt|bred|coltiv|conduce|cortad|crescut|criad|crop|cultiva|cultivé|decupată|dirigid|dyrk|ekili|élevé|ensemencé|farm|fikk til|former|forplantet|fortplantas|förvaltade|frø|gestit|gezaaid|gezüchtet|gospodarować|hodowlany|kasvatetaan|kultiv|manag|manage|manejad|oppdretts|o'stiriladi|propaga|seed|sembrad|semead|seminat|siewny|uppfödda|upraw|uprawiać ziemię|viljel|wychowany|ympades|засевали|культура|культурный|разводятся|удалось|לִזרוֹעַ|מְתוּרבָּת|מופץ|מעובדים|תְבוּאָה|المزارع|م|حصول|مزروع|مزروع|نشر|ولدت|교양 있는|씨를 뿌린|양식장|자란|전파 된|伝播|传播|养殖|孕育|栽培された|種まき|繁殖した"



#########################################################################################################################
## NICHE BREADTH CALCULATIONS
#########################################################################################################################


## Loop over a list of subfolders
read_bind_tables = function (table.list, path) {
  
  READ.BIND.TABLE <- table.list[c(1:length(table.list))] %>% 
    
    ## pipe the list into lapply
    lapply(function(x) {
      
      ## create the character string
      f <- paste0(path, x)
      
      ## read each .csv file
      d <- read.csv(f)
      
      ## now add the new columns
      d <- d[, columns]
      
    }) %>% 
    
    ## finally, bind all the rows together
    bind_rows
  
}




#########################################################################################################################
## NICHE BREADTH CALCULATIONS
#########################################################################################################################


########################################################################################################################
## SHAWN'S NICHE
########################################################################################################################


## This is much simpler than Stuart's EG, but we can add in more stuff later
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
  
  ## concatenate output
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
## STU'S CALC
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
## CLEANING FUNCTIONS
#########################################################################################################################


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





########################################################################################################################
## PLOTTING FUNCTIONS
########################################################################################################################


########################################################################################################################
## MAPPING FUNCTIONS
########################################################################################################################


## create simple maps of all GBIF records for selected taxa, in Australia and overseas
map_GBIF_records = function (taxa.list, DF) {
  
  ###############################
  ## for all the taxa in the list
  for (taxa.n in taxa.list) {
    
    ################################################################
    ## If the dim = 0 for that taxa subset to Australia, skip to next
    if (dim(DF[ which(DF$searchTaxon == taxa.n 
                      & DF$country == "Australia"), ][, c("lon", "lat")])[1] == 0) {
      
      print (paste ("Possible no Australian records for ", taxa.n, "skipping"))
      
      ## First, check if the file exists
      file   = paste("./output/Figures/niche_summary/maps/", taxa.n, "_WORLD_GBIF_map.png", sep = "")
      
      ## If it's already downloaded, skip
      if (file.exists (file)) {
        
        print (paste ("World map exists for", taxa.n, "skipping"))
        next
        
      }
      
      
      ########################################
      ## Plot global occurences for taxa.n
      plot(LAND, col = 'grey', bg = 'sky blue')
      title(paste0("Global occurrences for ", taxa.n))
      
      ## add points
      points(DF[ which(DF$searchTaxon == taxa.n), ][, c("lon", "lat")], 
             pch = ".", col = "red", cex = 3, asp = 1)
      
      
      ###############################
      ## Save global maps to file
      ## start Cairo device
      CairoPNG(width  = 16180, height = 10000, 
               file   = paste("./output/Figures/niche_summary/maps/", taxa.n, "_WORLD_GBIF_map.png", sep = ""), 
               canvas = "white", bg = "white", units = "px", dpi = 600)
      
      par(mgp      = c(10, 4, 0), 
          oma      = c(1.5, 1.5, 1.5, 1.5),
          font.lab = 2) 
      
      ## Add land
      plot(LAND, #add = TRUE, 
           lwd = 1.8, asp = 1, col = 'grey', bg = 'sky blue')
      
      ## add points
      points(DF[ which(DF$searchTaxon == taxa.n), ][, c("lon", "lat")], 
             pch = ".", cex = 7, col = "red", cex.lab = 3, cex.main = 4, cex.axis = 2, 
             main = paste0("Global occurrences for ", taxa.n), 
             xlab = "", ylab = "", asp = 1)
      
      ## title 
      title(paste0("Global occurrences for ", taxa.n),
            cex.main = 4,   font.main = 4, col.main = "blue")
      
      ## finsh the device
      dev.off()
      
    }
    
    else {
      
      ## First, check if the file exists
      file   = paste("./output/Figures/niche_summary/maps/", taxa.n, "_WORLD_GBIF_map.png", sep = "")
      
      ## If it's already downloaded, skip
      if (file.exists (file)) {
        
        print (paste ("Global map exists for", taxa.n, "skipping"))
        next
        
      }
      
      
      ##############################################
      ## Plot global occurences for taxa.n to screen
      plot(LAND, col = 'grey', bg = 'sky blue')
      title(paste0("Global occurrences for ", taxa.n))
      
      ## add points
      points(DF[ which(DF$searchTaxon == taxa.n), ][, c("lon", "lat")], 
             pch = ".", col = "red", cex = 3, asp = 1)
      
      ## PLot Australian occurrences to screen
      plot(DF[ which(DF$searchTaxon == taxa.n 
                                      & DF$country == "Australia"), ][, c("lon", "lat")], 
           pch = ".", cex = 5, col = "red", asp = 1)
      
      ## add title
      title(paste0("Australian occurrences for ", taxa.n))
      plot(LAND, add = TRUE, asp = 1)
      
      ##################################
      ## Save global record maps to file
      ## start Cairo device
      CairoPNG(width  = 16180, height = 10000, 
               file = paste("./output/Figures/niche_summary/maps/", taxa.n, "_WORLD_GBIF_map.png", sep = ""), 
               canvas = "white", bg = "white", units = "px", dpi = 600)
      
      ## set dimensions
      par(mgp      = c(10, 4, 0), 
          oma      = c(1.5, 1.5, 1.5, 1.5),
          font.lab = 2)
      
      ## Add land
      plot(LAND, #add = TRUE, 
           lwd = 1.8, asp = 1, col = 'grey', bg = 'sky blue')
      
      ## add points
      points(DF[ which(DF$searchTaxon == taxa.n), ][, c("lon", "lat")], 
             pch = ".", cex = 7, col = "red", cex.lab = 3, cex.main = 4, cex.axis = 2, 
             main = paste0("Global occurrences for ", taxa.n), 
             xlab = "", ylab = "", asp = 1)
      
      ## title 
      title(paste0("Global occurrences for ", taxa.n),
            cex.main = 4,   font.main = 4, col.main = "blue")
      
      ## finsh the device
      dev.off()
      
      ################################
      ## Save Aus records maps to file
      ## start Cairo device
      CairoPNG(width  = 16180, height = 10000, 
               file   = paste("./output/Figures/niche_summary/maps/", taxa.n, "_AUS_GBIF_map.png", sep = ""), 
               canvas = "white", bg = "white", units = "px", dpi = 600)
      
      ## set par
      par(mgp      = c(10, 3, 0), 
          oma      = c(1.5, 1.5, 1.5, 1.5),
          font.lab = 2)
      
      ## plot
      plot(DF[ which(DF$searchTaxon == taxa.n 
                                      & DF$country == "Australia"), ][, c("lon", "lat")], 
           pch = ".", cex = 10, col = "red", cex.lab = 3, cex.main = 4, cex.axis = 2.5,
           font.main = 4, col.main = "blue",
           main = paste0("Australian occurrences for ", taxa.n), 
           xlab = "", ylab = "", asp = 1)
      
      ## Add land
      plot(LAND, add = TRUE, lwd = 1.8, asp = 1) # col = 'grey', bg = 'sky blue')
      
      ## finsh the device
      dev.off()
      
    }
    
  }
  
}



########################################################################################################################
## Print simple maps of all GBIF records for selected taxa, in Australia and overseas
print_occurrence_records = function (taxa.list, DF) {
  
  ###############################
  ## for all the taxa in the list
  for (taxa.n in taxa.list) {
    
    ################################################################
    ## If the dim = 0 for that taxa subset to Australia, skip to next
    if (dim(DF[ which(DF$searchTaxon == taxa.n 
                                       & DF$country == "Australia"), ][, c("lon", "lat")])[1] == 0) {
      
      print (paste ("Possible no Australian records for ", taxa.n, "skipping"))
      
      ########################################
      ## Plot global occurences for taxa.n
      plot(LAND, col = 'grey', bg = 'sky blue')
      title(paste0("Global occurrences for ", taxa.n))
      
      ## add points
      points(DF[ which(DF$searchTaxon == taxa.n), ][, c("lon", "lat")], 
             pch = ".", col = "red", cex = 3, asp = 1)
      
    }
    
    else {
      
      ##############################################
      ## Plot global occurences for taxa.n to screen
      plot(LAND, col = 'grey', bg = 'sky blue')
      title(paste0("Global occurrences for ", taxa.n))
      
      ## add points
      points(DF[ which(DF$searchTaxon == taxa.n), ][, c("lon", "lat")], 
             pch = ".", col = "red", cex = 3, asp = 1)
      
      ## PLot Australian occurrences to screen
      plot(DF[ which(DF$searchTaxon == taxa.n 
                                      & DF$country == "Australia"), ][, c("lon", "lat")], 
           pch = ".", cex = 5, col = "red", asp = 1)
      
      ## add title
      title(paste0("Australian occurrences for ", taxa.n))
      plot(LAND, add = TRUE, asp = 1)
      
      ##############################################
      ## Plot Histograms to the screen
      ## Env.1: create the plot dimensions
      # nf <- layout(mat = matrix(c(1,2),2,1, byrow = TRUE),  height = c(1,3))
      # par(mar = c(3.1, 3.1, 1.1, 2.1),
      #     oma = c(1.5, 1.5, 1.5, 1.5)) 
      # 
      # ## create vector for the environmental dimenson, and set min and max for plotting
      # env.1 = DF[ which(DF$searchTaxon == taxa.n), ][[env.var.1]]
      # min.1 = min(env.1)
      # max.1 = max(env.1)
      # 
      # ## now create the boxplot
      # boxplot(env.1, horizontal = TRUE,  outline = TRUE, ylim  = c(min.1, max.1), frame = FALSE, col = env.col.1, axes = FALSE)
      # 
      # ## and the histogram
      # hist(env.1, xlim = c(min.1, max.1),
      #      breaks = 50, border = NA, col = env.col.1, main = taxa.n,
      #      xlab = paste0("Worldclim ", env.var.1, " ", env.units.1, sep = " "))
      
    }
    
  }
  
}




########################################################################################################################
## HISTOGRAM FUNCTIONS
########################################################################################################################


## create simple maps of all GBIF records for selected taxa, in Australia and overseas
Print_global_histogram = function (taxa.list, DF, env.var.1, env.col.1, env.units.1,
                                   env.var.2, env.col.2, env.units.2) {
  
  ###############################
  ## for all the taxa in the list
  for (taxa.n in taxa.list) {
    
    #####################################################
    ## 1). Plot histograms for global occurences of taxa.n
    
    #####################################################
    ## Env.1: create the plot dimensions
    nf <- layout(mat = matrix(c(1,2),2,1, byrow = TRUE),  height = c(1,3))
    par(mar = c(5, 1.6, 1.6, 1.6),
        oma = c(0.5, 0.5, 0.5, 0.5)) 
    
    ## Create vector for the environmental dimenson, and set min and max for plotting
    env.1 = DF[ which(DF$searchTaxon == taxa.n), ][[env.var.1]]
    min.1 = min(env.1, na.rm = TRUE)
    max.1 = max(env.1, na.rm = TRUE)
    
    ## now create the boxplot
    boxplot(env.1, horizontal = TRUE,  outline = TRUE, 
            ylim  = c(min.1, max.1), frame = FALSE, 
            col = env.col.1, axes = FALSE)
    
    ## and the histogram
    hist(env.1, xlim = c(min.1, max.1),
         breaks = 50, border = NA, col = env.col.1, main = taxa.n,
         xlab = paste0("Worldclim ", env.var.1, " ", env.units.1, sep = " "))
    
    #####################################################
    ## Env.2: create the plot dimensions
    nf <- layout(mat = matrix(c(1,2),2,1, byrow = TRUE),  height = c(1,3))
    par(mar = c(5, 1.6, 1.6, 1.6),
        oma = c(0.5, 0.5, 0.5, 0.5)) 
    
    ## Set min and max
    env.2 = DF[ which(DF$searchTaxon == taxa.n), ][[env.var.2]]
    min.2  = min(env.2, na.rm = TRUE)
    max.2  = max(env.2, na.rm = TRUE)
    
    ## Now create the boxplot
    boxplot(env.2, horizontal = TRUE,  outline = TRUE, 
            #ylim  = c(min.2, max.2), 
            frame = FALSE, 
            col = env.col.2, axes = FALSE)
    
    ## And the histogram
    hist(env.2, 
         #xlim = c(min.2, max.2),
         breaks = 50, border = NA, col = env.col.2, main = taxa.n,
         xlab = paste0("AWAP ", env.var.2, " ", env.units.2, sep = " "))
    
  }
  
}



## create simple maps of all GBIF records for selected taxa, in Australia and overseas
histogram_GBIF_records = function (DF, taxa.list, env.var.1, env.col.1, env.units.1,
                                   env.var.2, env.col.2, env.units.2) {
  
  ###############################
  ## for all the taxa in the list
  for (taxa.n in taxa.list) {
    
    ## First, check if the file exists
    file  = paste("./output/figures/Traits_v_glasshouse/", 
                  taxa.n, "_", env.var.1, "_world_GBIF_histo.png", sep = "")
    
    ## If it's already downloaded, skip
    if (file.exists (file)) {
      
      print (paste ("Global histogram exists for", taxa.n, "skipping")) 
      next
      
    }
    
    #####################################################
    ## 1). Plot histograms for global occurences of taxa.n
    
    #####################################################
    ## Env.1: create the plot dimensions
    nf <- layout(mat = matrix(c(1,2),2,1, byrow = TRUE),  height = c(1,3))
    par(mar = c(3.1, 3.1, 1.1, 2.1),
        oma = c(1.5, 1.5, 1.5, 1.5)) 
    
    ## create vector for the environmental dimenson, and set min and max for plotting
    env.1 = DF[ which(DF$searchTaxon == taxa.n), ][[env.var.1]]
    min.1 = min(env.1, na.rm = TRUE)
    max.1 = max(env.1, na.rm = TRUE)
    
    ## now create the boxplot
    boxplot(env.1, horizontal = TRUE,  outline = TRUE, 
            ylim  = c(min.1, max.1), 
            frame = FALSE, col = env.col.1, axes = FALSE)
    
    ## and the histogram
    hist(env.1, xlim = c(min.1, max.1),
         breaks = 50, border = NA, col = env.col.1, main = taxa.n,
         xlab = paste0("Worldclim ", env.var.1, " ", env.units.1, sep = " "))
    
    #####################################################
    ## Env.2: create the plot dimensions
    nf <- layout(mat = matrix(c(1,2),2,1, byrow = TRUE),  height = c(1,3))
    par(mar = c(3.1, 3.1, 1.1, 2.1),
        oma = c(1.5, 1.5, 1.5, 1.5)) 
    
    ## set min and max
    env.2 = DF[ which(DF$searchTaxon == taxa.n), ][[env.var.2]]
    min.2  = min(env.2, na.rm = TRUE)
    max.2  = max(env.2, na.rm = TRUE)
    
    ## now create the boxplot
    boxplot(env.2, horizontal = TRUE,  outline = TRUE, 
            #ylim  = c(min.2, max.2), 
            frame = FALSE, col = env.col.2, axes = FALSE)
    
    ## and the histogram
    hist(env.2, 
         #xlim = c(min.2, max.2),
         breaks = 50, border = NA, col = env.col.2, main = taxa.n,
         xlab = paste0("Worldclim ", env.var.2, " ", env.units.2, sep = " "))
    
    #################################
    ## 2). Save env.1 histogram to file
    message("Saving histogram for ", taxa.n)
    CairoPNG(width  = 16180, height = 12000,
             file   = paste("./output/figures/Traits_v_glasshouse/", 
                            taxa.n, "_", env.var.1, "_world_GBIF_histo.png", sep = ""),
             canvas = "white", bg = "white", units = "px", dpi = 600)
    
    ## create layout
    nf <- layout(mat = matrix(c(1,2), 2,1, byrow = TRUE), height = c(1,3))
    par(mar = c(7.5, 5.5, 4, 2),
        oma = c(4, 4, 4, 4),
        mgp = c(6, 3, 0),
        font.lab = 2, lwd = 2)
    
    ## print the boxplot
    boxplot(env.1, horizontal = TRUE,  outline = TRUE, 
            ylim  = c(min.1, max.1), 
            frame = FALSE, col = env.col.1, axes = FALSE)
    
    ## and print the
    hist(env.1, xlim = c(min.1, max.1),
         breaks = 50, border = NA, col = env.col.1, main = "",
         xlab = paste0(taxa.n, " ", env.var.1, " ", "(", env.units.1, ")", sep = ""), ylab = "",
         cex.lab = 5, cex.axis = 4)
    
    ## finsh the device
    dev.off()
    
    #################################
    ## 3). Save env.2 histogram to file
    CairoPNG(width  = 16180, height = 12000,
             file   = paste("./output/figures/Traits_v_glasshouse/", 
                            taxa.n, "_", env.var.2, "_world_GBIF_histo.png", sep = ""),
             canvas = "white", bg = "white", units = "px", dpi = 600)
    
    ## create layout
    nf <- layout(mat = matrix(c(1,2), 2,1, byrow = TRUE), height = c(1,3))
    par(mar = c(7.5, 5.5, 4, 2),
        oma = c(4, 4, 4, 4),
        mgp = c(6, 3, 0),
        font.lab = 2, lwd = 2)
    
    ## print the boxplot
    boxplot(env.2, horizontal = TRUE,  outline = TRUE, 
            #ylim  = c(min.2, max.2), 
            frame = FALSE, col = env.col.2, axes = FALSE)
    
    ## and print the
    hist(env.2, 
         #xlim = c(min.2, max.2),
         breaks = 50, border = NA, 
         col = env.col.2, main = "",
         xlab = paste0(taxa.n, " ", env.var.2, " ", "(", env.units.2, ")", sep = ""), 
         ylab = "", cex.lab = 5, cex.axis = 4)
    
    ## finsh the device
    dev.off()
    
  }
  
}





########################################################################################################################
## TABLE FUNCTIONS
########################################################################################################################


########################################################################################################################
## create simple summaries for selected taxa
COMBO_check_records = function (taxa.list, 
                                #columns, 
                                DF) {
  
  ###############################
  ## for all the taxa in the list
  for (taxa in taxa.list) {
    
    ## slice the table 
    #summary.table <- DF[, columns][ which(DF[["searchTaxon"]] == taxa ), ]
    summary.table <- subset(DF, searchTaxon == taxa)
    
    ## print to screen
    View(summary.table)
    
  }
  
}


##
GBIF_summary_slice = function (taxa.list, env.cols, GBIF) {
  
  ###############################
  ## for all the taxa in the list
  for (taxa.n in taxa.list) {
    
    ## slice the table 
    summary.table <- GBIF[, env.cols][ which(GBIF[["searchTaxon"]] == taxa.n ), ]

    ## print to screen
    print(kable(summary.table, row.names = FALSE))
    
  }
  
}


########################################################################################################################
## Create lists of species in each LGA
LGA_lists = function (taxa.list, env.cols, GBIF) {
  
  ###############################
  ## for all the taxa in the list
  for (taxa.n in taxa.list) {
    
    ## slice the table 
    summary.table <- GBIF[, env.cols][ which(GBIF[["searchTaxon"]] == taxa.n ), ]
    
    ## print to screen
    print(kable(summary.table, row.names = FALSE))
    
  }
  
}





########################################################################################################################
## AREA FUNCTIONS
########################################################################################################################


########################################################################################################################
## create simple summaries for selected taxa
GBIF_AOO = function (taxa.list) {
  
  ###############################
  ## for all the taxa in the list
  for (taxa.n in taxa.list) {
    
    ## slice the table 
    test      = subset(data, searchTaxon == sp.n)[, c("lon", "lat")]
    
    ## print to screen
    aoo(test)
    
  }
  
}


#########################################################################################################################
#############################################  FUNCTIONS FOR HORT AUS LIST ############################################## 
#########################################################################################################################