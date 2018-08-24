#########################################################################################################################
############################################# UPDATE ALA  ############################################################### 
#########################################################################################################################


#########################################################################################################################
## This code updates the ALA data
ala.download = list.files(ALA_path, pattern = ".rds")


#########################################################################################################################
## 2). COMBINE ALA SPECIES INTO ONE DATASET
#########################################################################################################################


#########################################################################################################################
## Combine all the taxa into a single dataframe at once
ALA.ALL <- ala.download %>%   ## spp.download[c(1:length(spp.download))] 
  
  ## Pipe the list into lapply
  lapply(function(x) {
    
    ## Create a character string of each .RData file
    #f <- sprintf("./data/base/HIA_LIST/ALA/TREE_SPECIES/%s", x)
    f <- sprintf(paste0(ALA_path, "%s"), x)
    
    ## Load each file
    d <- get(load(f))
    
    ## Now drop the columns which we don't need
    dat <- data.frame(searchTaxon = x,
                      d[, colnames(d) %in% ALA.keep],
                      stringsAsFactors = FALSE)
    
    if(!is.character(dat$id)) {
      
      dat$id <- as.character(dat$id)
      
    }
    
    ## This is a list of columns in different files which have weird characters
    dat$coordinateUncertaintyInMetres = as.numeric(dat$coordinateUncertaintyInMetres)
    dat$year             = as.numeric(dat$year)
    dat$month            = as.numeric(dat$month)
    #dat$latitudeOriginal = as.numeric(dat$latitudeOriginal)
    #dat$individualCount  = as.numeric(dat$individualCount)
    
    dat$searchTaxon = gsub("_ALA_records.RData", "", dat$searchTaxon)
    names(dat)[names(dat) == 'latitude']  <- 'lat'
    names(dat)[names(dat) == 'longitude'] <- 'lon'
    return(dat)
    
  }) %>%
  
  ## Finally, bind all the rows together
  bind_rows


#########################################################################################################################
## Now check the ALA species names : which has the least holes?
(sum(is.na(ALA.ALL$scientificName))          + dim(subset(ALA.ALL, scientificName == ""))[1])/dim(ALA.ALL)[1]*100
(sum(is.na(ALA.ALL$scientificNameOriginal))  + dim(subset(ALA.ALL, scientificNameOriginal == ""))[1])/dim(ALA.ALL)[1]*100
(sum(is.na(ALA.ALL$species))                 + dim(subset(ALA.ALL, species == ""))[1]/dim(ALA.ALL))[1]*100


## What is the match between 'species' and 'scientificNameOriginal'?
## What is the match between scientificNameOriginal, and the searchTaxon?
Match.ALA = ALA.ALL  %>%
  mutate(Match.SN.ST = 
           str_detect(scientificNameOriginal, searchTaxon)) %>%
  
  mutate(Match.SP.ST = 
           str_detect(species, searchTaxon)) %>%
  
  mutate(Match.SN.SP = 
           str_detect(scientificNameOriginal, species)) %>%
  
  select(one_of(c("scientificName",
                  "scientificNameOriginal",
                  "species",
                  "searchTaxon",
                  "Match.SN.ST",
                  "Match.SP.ST",
                  "Match.SN.SP")))

View(Match.ALA)

## So for 15-12% of the records, neither the scientificNameOriginal or the species match the search taxon.
dim(subset(Match.ALA,  Match.SN.ST == "FALSE"))[1]/dim(Match.ALA)[1]*100
dim(subset(Match.ALA,  Match.SP.ST == "FALSE"))[1]/dim(Match.ALA)[1]*100



## So rename 'scientificNameOriginal' to 'scientificName', to match GBIF
ALA.RENAME = ALA.ALL
ALA.RENAME$scientificName = NULL 
ALA.RENAME$scientificName = ALA.RENAME$scientificNameOriginal 
ALA.RENAME$scientificNameOriginal = NULL 
(sum(is.na(ALA.RENAME$scientificName)) + dim(subset(ALA.RENAME, scientificName == ""))[1])/dim(ALA.ALL)[1]*100


########################################################################################################################
## Check names
## What names get returned?
sort(names(ALA.RENAME))
ALA.TRIM <- ALA.RENAME %>% 
  select(one_of(ALA.keep))

dim(ALA.TRIM)
names(ALA.TRIM)


#########################################################################################################################
## Just get the newly downloaded species................................................................................
ALA.TRIM = ALA.TRIM[ALA.TRIM$searchTaxon %in% GBIF.spp, ]
dim(ALA.TRIM)


## What are the unique species?
length(unique(ALA.TRIM$searchTaxon))
length(unique(ALA.TRIM$scientificName)) 
length(unique(ALA.TRIM$species))
(sum(is.na(ALA.TRIM$scientificName)) + dim(subset(ALA.TRIM, scientificName == ""))[1])/dim(ALA.TRIM)[1]*100




#########################################################################################################################
## 2). CHECK TAXONOMY RETURNED BY ALA USING TAXONSTAND
######################################################################################################################### 

## The problems is the mismatch between what we searched, and what GBIF returned. 

## 1). Create the inital list by combing the planted trees with evergreen list
##     TREE.HIA.SPP = intersect(subset(TI.LIST, Plantings > 50)$searchTaxon, CLEAN.SPP$Binomial)

## 2). Clean this list using the GBIF backbone taxonomy :: use the "species" column

## 3). Run the GBIF "species" list through the TPL taxonomy. Take "New" Species and Genus as the "searchTaxon"

## 4). Use rgbif and ALA4R to download occurence data. ALA is ok, because the taxonomy is resolved.
##     For GBIF, we use
##     key  <- name_backbone(name = sp.n, rank = 'species')$usageKey
##     GBIF <- occ_data(taxonKey = key, limit = GBIF.download.limit)

##     This returns multiple keys and synonyms, but there is no simple way to skip these


## 5). Join the TPL taxonomy to the "scientificName" field. We can't use "name" (the equivalent of "species", it seems),
##     because name is always the same as the searchTaxon and not reliable (i.e. they will always match, and we know that
##     no one has gone through and checked each one).

##     Exclude records where the "scientificName" both doesn't match the "searchTaxon", and, also is not a synonym according to TPL
##     This is the Same as taking the SNs which are accepted, but which don't match ST.

##     Then we model these records as before. Of 390 species we downloaded, we will pick the 200 with acceptable maps.
##     One line in the MS : we matched the GBIF backbone taxo against the TPL taxo, and searched the currently accepted names
##     in GBIF and ALA, excluding incorreclty matching records (probably don't say this).



#########################################################################################################################
## Use "Taxonstand" to check the taxonomy :: which field to use?
ALA.TREES.TAXO <- TPL(unique(ALA.TRIM$scientificName), infra = TRUE,
                 corr = TRUE, repeats = 100)  ## to stop it timing out...
sort(names(ALA.TREES.TAXO))
saveRDS(ALA.TREES.TAXO, 'data/base/HIA_LIST/COMBO/ALA_TAXO_400.rds')





#########################################################################################################################
## 3). CREATE TABLE OF PRE-CLEAN FLAGS AND FILTER RECORDS
#########################################################################################################################


#########################################################################################################################
## Create a table which counts the number of records meeting each criteria:
## Note that TRUE indicates there is a problem (e.g. if a record has no lat/long, it will = TRUE)
ALA.TREES.PROBLEMS <- with(ALA.TREES.TAXO,
                      
                      table(
                        
                        ## Note this list is incomplete
                        
                        ## No coordinates
                        is.na(lon)|is.na(lat),
                        
                        ## Taxon rank is genus/form?
                        #taxonRank == 'GENUS' & 'FORM',
                        
                        ## Taxonomic status: consider if this is the right filter, it seems too restrictive
                        #New.Taxonomic.status == 'Unresolved',
                        
                        ## Cultivated
                        #CULTIVATED == 'CULTIVATED',
                        
                        ## Establishment means is "MANAGED" is included in the above
                        #establishmentMeans == 'MANAGED' & !is.na(establishmentMeans),
                        
                        ## Collected before 1950 or na.year
                        year < 1950 & !is.na(year),
                        
                        ## No year
                        is.na(year),
                        
                        ## Coordinate uncertainty is > 100 or is not NA
                        coordinateUncertaintyInMetres > 1000 & 
                          !is.na(coordinateUncertaintyInMetres)
                        
                        ## Other checks using coordinateCleaner
                        
                      )
                      
) %>% 
  
  ## Create a data frame and set the names
  as.data.frame %>%  
  
  setNames(c('NO_COORD',     
             #'TAXON_STATUS', 
             #'CULTIVATED', 
             'PRE_1950',     'NO_YEAR', 
             'COORD_UNCERT', 'COUNT'))


#########################################################################################################################
## Now filter the ALA records using conditions which are not too restrictive
ALA.TREES.CLEAN <- ALA.TREES.TAXO %>% 
  
  ## Note that these filters are very forgiving...
  ## Unless we include the NAs, very few records are returned!
  filter(!is.na(lon) & !is.na(lat),
         
         ## CULTIVATED == 'CULTIVATED', 
         ## #establishmentMeans!='MANAGED' | is.na(establishmentMeans),
         ## New.Taxonomic.status == 'Accepted'
         
         year >= 1950 & !is.na(year))


## Check
names(ALA.TREES.CLEAN)





#########################################################################################################################
## 4). REMOVE POINTS OUTSIDE WORLDCLIM LAYERS...
#########################################################################################################################


## Can use WORLDCIM rasters to get only records where wordlclim data is. 
## First, get one of the BIOCLIM variables
world.temp = raster("./data/base/worldclim/world/0.5/bio/current/bio_01")
#plot(world.temp)


## Now get the XY centroids of the unique 1km * 1km WORLDCLIM blocks where ALA records are found
## Get cell number(s) of WORLDCLIM raster from row and/or column numbers. Cell numbers start at 1 in the upper left corner, 
## and increase from left to right, and then from top to bottom. The last cell number equals the number of raster cells 
xy <- cellFromXY(world.temp, ALA.TREES.CLEAN[c("lon", "lat")]) %>% 
  
  ## get the unique raster cells
  unique %>% 
  
  ## Get coordinates of the center of raster cells for a row, column, or cell number of WORLDCLIM raster
  xyFromCell(world.temp, .)


## Take a look at xy: NA's should be removed...
# summary(xy)
# str(xy)
# points(xy, pch = ".", col = "red")


## For some reason, we need to convert the xy coords to a spatial points data frame, in order to avoid this error:
## 'NAs introduced by coercion to integer range'
xy <- SpatialPointsDataFrame(coords = xy, data = as.data.frame(xy),
                             proj4string = CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))


## Now extract the temperature values for the unique 1km centroids which contain ALA data
class(xy)
z   = extract(world.temp, xy)

# Warning message:
#   In .doExtract(x, i, ..., drop = drop) :
#   some indices are invalid (NA returned)

hist(z, border = NA, col = "orange", breaks = 50, main = "", xlab = "Worldclim Annual temp")


## Then track which values of Z are on land or not
onland = z %>% is.na %>%  `!` # %>% xy[.,]  cells on land or not
#summary(onland)


## Finally, filter the cleaned ALA data to only those points on land. 
## This is achieved with the final [onland]
ALA.TREES.LAND = filter(ALA.TREES.CLEAN, cellFromXY(world.temp, ALA.TREES.CLEAN[c("lon", "lat")]) %in% 
                     unique(cellFromXY(world.temp, ALA.TREES.CLEAN[c("lon", "lat")]))[onland])


## how many records were on land?
records.ocean = dim(ALA.TREES.CLEAN)[1] - dim(ALA.TREES.LAND)[1]  ## 91575 records are in the ocean   


## Print the dataframe dimensions to screen
dim(ALA.TREES.LAND)
length(unique(ALA.TREES.LAND$searchTaxon))


## Add a source column
ALA.TREES.LAND$SOURCE = 'ALA'


## Free some memory
gc()



#########################################################################################################################
## save data
saveRDS(ALA.TREES.LAND, file = paste("./data/base/HIA_LIST/GBIF/ALA_TREES_LAND.rds"))




#########################################################################################################################
## OUTSTANDING OCCURRENCE TASKS:
#########################################################################################################################


##


#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################