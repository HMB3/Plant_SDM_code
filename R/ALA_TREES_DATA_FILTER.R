#########################################################################################################################
############################################# UPDATE ALA  ############################################################### 
#########################################################################################################################


#########################################################################################################################
## This code updates the ALA data
spp.download = list.files("./data/base/HIA_LIST/ALA/TREE_SPECIES/", pattern = ".RData")


## Load ALA update :: the coordinates are coming in as characters
AVH.TREES = read.csv("./data/base/HIA_LIST/ALA/AVH_tree_spp.csv", stringsAsFactors = FALSE)
# ALA.UPDATE = read.csv("./data/base/HIA_LIST/ALA/Plants.csv", stringsAsFactors = FALSE)
# saveRDS(ALA.UPDATE, file = paste("./data/base/HIA_LIST/ALA/ALA_PLANTS_UPDATE.rds"))
#ALA.UPDATE = readRDS("./data/base/HIA_LIST/GBIF/ALA_PLANTS_UPDATE.rds")


# dim(ALA.UPDATE)
# dim(AVH.TREES)
# dim(ALA.LAND)
# 
# length(unique(ALA.LAND$scientificname))
# length(unique(ALA.UPDATE$scientificName))  
# length(unique(AVH.TREES$scientificName))  





#########################################################################################################################
## 1). LOAD TAXA, ADD SEARCH TAXA, DROP COLUMNS AND COMBINE SPECIES DATAFRAMES INTO ONE
#########################################################################################################################


## What are the names?
# sort(names(ALA.UPDATE))
# intersect(names(ALA.UPDATE), names(ALA.LAND))
# intersect(names(ALA.UPDATE), ALA.keep)
# 
# 
# ## Assuming decimal lat and lon are the right columns 
# ALA.UPDATE$lat = as.numeric(ALA.UPDATE$decimalLatitude)
# ALA.UPDATE$lon = as.numeric(ALA.UPDATE$decimalLongitude)
# summary(ALA.UPDATE$lat)
# summary(ALA.UPDATE$lon)
# 
# 
# #########################################################################################################################
# ## Now, get just the columns we want to keep. Note gc() frees up RAM
# ALA.TRIM <- ALA.UPDATE %>% 
#   select(one_of(ALA.keep))
# 
# 
# ## Check names
# dim(ALA.TRIM)
# names(ALA.TRIM)
# setdiff(ALA.keep, names(ALA.TRIM))
# 
# 
# ## Just get the species we need - this should be the intersection of Ale's list and the HIA list
# ALA.OLD    = ALA.LAND[ALA.LAND$scientificname %in% GBIF.spp, ]
# 
# 
# ## How big are the datasets?
# dim(ALA.TRIM);dim(ALA.OLD);dim(AVH.UPDATE)
# length(unique(ALA.TRIM$scientificName))    ## not many records per species
# length(unique(AVH.UPDATE$scientificName)) 
# length(unique(ALA.OLD$scientificname))   




#########################################################################################################################
## 2). COMBINE ALA SPECIES INTO ONE DATASET
#########################################################################################################################


#########################################################################################################################
## Combine all the taxa into a single dataframe at once
ALA.TREES.TRIM <- spp.download %>%   ## spp.download[c(1:length(spp.download))] 
  
  ## Pipe the list into lapply
  lapply(function(x) {
    
    ## Create a character string of each .RData file
    f <- sprintf("./data/base/HIA_LIST/ALA/TREE_SPECIES/%s", x)
    
    ## Load each file
    d <- get(load(f))
    
    ## Now drop the columns which we don't need
    dat <- data.frame(searchTaxon = x,
                      d[, colnames(d) %in% ALA.keep],
                      stringsAsFactors = FALSE)
    
    if(!is.character(dat$id)) {
      
      dat$id <- as.character(dat$id)
      
    }
    
    ## Need to print the object within the loop
    dat$coordinateUncertaintyInMetres = as.numeric(dat$coordinateUncertaintyInMetres)
    dat$year        = as.numeric(dat$year)
    dat$month       = as.numeric(dat$month)
    dat$searchTaxon = gsub("_ALA_records.RData", "", dat$searchTaxon)
    names(dat)[names(dat) == 'latitude']  <- 'lat'
    names(dat)[names(dat) == 'longitude'] <- 'lon'
    dat
    
  }) %>%
  
  ## Finally, bind all the rows together
  bind_rows


## Check the output :: how does this compare to John's 
sort(names(ALA.TREES.TRIM))
dim(ALA.TREES.TRIM)
length(unique(ALA.TREES.TRIM$scientificName))
length(unique(ALA.TREES.TRIM$searchTaxon))
sort(unique(ALA.TREES.TRIM$scientificName))


## How can this be cleaned
head(ALA.TREES.TRIM, 10)[, c("scientificName",
                             "searchTaxon")]

tail(ALA.TREES.TRIM, 10)[, c("scientificName",
                             "searchTaxon")]

## Check the lat/lon
class(ALA.TREES.TRIM$lat)
summary(ALA.TREES.TRIM$lat)
summary(ALA.TREES.TRIM$lon)


#########################################################################################################################
## Now, get just the columns we want to keep. Note gc() frees up RAM
## Check names
dim(ALA.TREES.TRIM)
names(ALA.TREES.TRIM)
setdiff(ALA.keep, names(ALA.TREES.TRIM))





#########################################################################################################################
## 2). CHECK TAXONOMY RETURNED BY ALA USING TAXONSTAND
######################################################################################################################### 


#########################################################################################################################
## Use "Taxonstand" to check the taxonomy. However, this also assumes that the ALA data is clean
ALA.TAXO <- TPL(unique(ALA.TREES.TRIM$scientificName), infra = TRUE,
                corr = TRUE, repeats = 100)  ## to stop it timing out...
sort(names(ALA.TAXO))


# #########################################################################################################################
# ## Use "Taxonstand" to check the taxonomy. However, this also assumes that the ALA data is clean
# ALA.TREES.TRIM.TAXO <- TPL(unique(ALA.TRIM$scientificName), infra = TRUE,
#                  corr = TRUE, repeats = 100)  ## to stop it timing out...
# sort(names(ALA.TREES.TRIM.TAXO))


#########################################################################################################################
## The procedure used for taxonomic standardization is based on function TPLck. A progress bar
## indicates the proportion of taxon names processed so far. In case the TPL website cannot be reached
## temporarily, the function returns an error but repeats trying to match the given taxon multiple times
## (see repeats). If standardization is still not successful, the input taxon is returned in field ’Taxon’
## with NA in all other fields


## Then join the ALA data to the taxonomic check, using "scientificName" as the join field...
ALA.TREES.TAXO <- ALA.TREES.TRIM %>%
  left_join(., ALA.TAXO, by = c("scientificName" = "Taxon"))


## So we can filter by the agreement between "scientificName", and "New.Taxonomic.status"?
## A better way could be to create a new filter for agreement between the search taxon genus and species,
## and the new genus and species. Which of these columns to keep?


#########################################################################################################################
## Now create a column for the agreement between the new genus and the old genus
## First, trim the spaces out
ALA.TREES.TAXO$scientificName  = trimws(ALA.TREES.TAXO$scientificName)


## Then combine the genus and species returned by TPL into
ALA.TREES.TAXO$TPL_binomial  = with(ALA.TREES.TAXO, paste(New.Genus, New.Species, sep = " "))


## Now match the searchTaxon with the binomial returned by TPL :: this would be the best field to filter on
ALA.TREES.TAXO$taxo_agree <- ifelse(
  ALA.TREES.TAXO$scientificName == ALA.TREES.TAXO$TPL_binomial, TRUE, FALSE)


## How many species agree?
round(with(ALA.TREES.TAXO, table(taxo_agree)/sum(table(taxo_agree))*100), 1)
round(with(ALA.TREES.TAXO, table(New.Taxonomic.status)/sum(table(New.Taxonomic.status))*100), 1)


## Also keep the unresolved records:
ALA.TAXO.UNRESOLVED <- ALA.TREES.TAXO %>%
  
  ## Note that these filters are very forgiving...
  ## Unless we include the NAs, very few records are returned!
  filter(New.Taxonomic.status == 'Unresolved')


## Also keep the unresolved records:
ALA.TAXO.RESOLVED <- ALA.TREES.TAXO %>%
  
  ## Note that these filters are very forgiving...
  ## Unless we include the NAs, very few records are returned!
  filter(New.Taxonomic.status == 'Accepted')


## Also keep the unresolved records:
ALA.TAXO.DISAGREE <- ALA.TREES.TAXO %>%
  
  ## Note that these filters are very forgiving...
  ## Unless we include the NAs, very few records are returned!
  filter(taxo_agree == 'FALSE')



#########################################################################################################################
## What do examples of resolved and unresolved records look like?
head(ALA.TREES.TAXO, 10)[, c("scientificName", 
                             "TPL_binomial",
                             "taxo_agree",
                             "New.Taxonomic.status")]


head(ALA.TAXO.RESOLVED, 10)[, c("scientificName", 
                                "TPL_binomial",
                                "taxo_agree",
                                "New.Taxonomic.status")]


head(ALA.TAXO.DISAGREE, 10)[, c("scientificName", 
                                "TPL_binomial",
                                "taxo_agree",
                                "New.Taxonomic.status")]


## Also keep the managed records:
unique(ALA.TAXO.UNRESOLVED$New.Taxonomic.status)
dim(ALA.TAXO.UNRESOLVED)   ##





#########################################################################################################################
## 2). CREATE TABLE OF PRE-CLEAN FLAGS AND FILTER RECORDS
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
## 6). REMOVE POINTS OUTSIDE WORLDCLIM LAYERS...
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
## OUTSTANDING NICHE TASKS:
#########################################################################################################################


#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################