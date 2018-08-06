#########################################################################################################################
####################################################  FILTER GBIF DATA ################################################## 
#########################################################################################################################


#########################################################################################################################
## This code combines the GBIF data into one file, and filters to the most reliable recrods. 
## Two main sources of uncertainty:


## Taxonomic
## Spatial 


#########################################################################################################################
## Read in all data to run the SDM code :: species lists, shapefile, rasters & tables
#source('./R/HIA_LIST_MATCHING.R')
rasterTmpFile()


#########################################################################################################################
## 1). LOAD TAXA, ADD SEARCH TAXA, DROP COLUMNS AND COMBINE SPECIES DATAFRAMES INTO ONE
#########################################################################################################################


#########################################################################################################################
## Create a list of species from the downloaded files
spp.download = list.files("./data/base/HIA_LIST/GBIF/SPECIES/", pattern = ".RData")
spp.download = gsub("_GBIF_records.RData", "", spp.download)
spp.download = trimws(spp.download)
bind.spp     = trimws(GBIF.spp)
spp.download = intersect(bind.spp, spp.download)


## memory is a problem. Can these sorts of operations be run in parallel?
memory.limit()
gc()


## Check GBIF column names:
sort(unique(gbifColsToDrop))
sort(unique(gbif.keep))
intersect(unique(gbifColsToDrop), unique(gbif.keep))     ## Check we don't get rid of any of the columns we want to keep


#########################################################################################################################
## Combine all the taxa into a single dataframe at once
GBIF.ALL <- spp.download[c(1:length(spp.download))] %>%   ## spp.download[c(1:length(spp.download))] 
  
  ## Pipe the list into lapply
  lapply(function(x) {
    
    ## Create a character string of each .RData file
    f <- sprintf("./data/base/HIA_LIST/GBIF/SPECIES/%s_GBIF_records.RData", x)
    
    ## Load each file
    d <- get(load(f))
    
    ## Now drop the columns which we don't need
    dat <- data.frame(searchTaxon = x, d[, !colnames(d) %in% gbifColsToDrop],
                      stringsAsFactors = FALSE)
    
    if(!is.character(dat$gbifID)) {
      
      dat$gbifID <- as.character(dat$gbifID)
      
    }
    
    ## Need to print the object within the loop
    dat
    
  }) %>%
  
  ## Finally, bind all the rows together
  bind_rows


#########################################################################################################################
## CHECK SEARCHED AND RETURNED TAXONOMIC NAMES FROM GBIF
## What are the names? Can't get date for GBIF across 
sort(names(GBIF.ALL))
str(GBIF.ALL$recordedBy)
str(GBIF.ALL$dateIdentified)
identical(unique(GBIF.ALL$searchTaxon), GBIF.spp)


## Now get just the columns we want to keep. Note gc() frees up RAM
GBIF.TRIM <- GBIF.ALL %>% 
  select(one_of(gbif.keep))


## Check names
dim(GBIF.TRIM)
names(GBIF.TRIM)
setdiff(names(GBIF.TRIM), gbif.keep)
intersect(names(GBIF.TRIM), gbif.keep)


#########################################################################################################################
## Just get the newly downloaded species................................................................................
GBIF.TRIM = GBIF.TRIM[GBIF.TRIM$searchTaxon %in% GBIF.spp, ]
dim(GBIF.TRIM)


## What are the unique species
length(unique(GBIF.TRIM$scientificName))  
length(unique(GBIF.TRIM$searchTaxon))  


## How do the searched and returned items compare?
head(GBIF.TRIM, 100)[, c("scientificName",
                         "searchTaxon")]

head(GBIF.TRIM, 100)[, c("scientificName",
                         "searchTaxon")]



#########################################################################################################################
## 2). CHECK TAXONOMY RETURNED BY GBIF USING TAXONSTAND
######################################################################################################################### 


#########################################################################################################################
## Use "Taxonstand" to check the taxonomy. However, this also assumes that the ALA data is clean
GBIF.TAXO <- TPL(unique(GBIF.TRIM$scientificName), infra = TRUE,
                 corr = TRUE, repeats = 100)  ## to stop it timing out...
sort(names(GBIF.TAXO))


## Now create table for Alessandro.......................................................................................
GBIF.LUT = as.data.frame(table(GBIF.TRIM$scientificName))
names(GBIF.LUT) = c("scientificName", "FREQUENCY")
GBIF.LUT = GBIF.LUT[with(GBIF.LUT, rev(order(FREQUENCY))), ] 
head(GBIF.LUT);dim(GBIF.LUT)


#########################################################################################################################
## The procedure used for taxonomic standardization is based on function TPLck. A progress bar
## indicates the proportion of taxon names processed so far. In case the TPL website cannot be reached
## temporarily, the function returns an error but repeats trying to match the given taxon multiple times
## (see repeats). If standardization is still not successful, the input taxon is returned in field ’Taxon’
## with NA in all other fields


## Then join the GBIF data to the taxonomic check, using "scientificName" as the join field...
GBIF.TRIM <- GBIF.TRIM %>%
  left_join(., GBIF.TAXO, by = c("scientificName" = "Taxon"))


## So we can filter by the agreement between "scientificName", and "New.Taxonomic.status"?
## A better way could be to create a new filter for agreement between the search taxon genus and species,
## and the new genus and species. Which of these columns to keep?


#########################################################################################################################
## Now create a column for the agreement between the new genus and the old genus
## First, trim the spaces out
GBIF.TRIM$searchTaxon  = trimws(GBIF.TRIM$searchTaxon)


## Then combine the genus and species returned by TPL into
GBIF.TRIM$TPL_binomial  = with(GBIF.TRIM, paste(New.Genus, New.Species, sep = " "))


## Now match the searchTaxon with the binomial returned by TPL :: this would be the best field to filter on
GBIF.TRIM$taxo_agree <- ifelse(
  GBIF.TRIM$searchTaxon == GBIF.TRIM$TPL_binomial, TRUE, FALSE)


## How many species agree?
round(with(GBIF.TRIM, table(taxo_agree)/sum(table(taxo_agree))*100), 1)
round(with(GBIF.TRIM, table(New.Taxonomic.status)/sum(table(New.Taxonomic.status))*100), 1)


## Also keep the unresolved records:
GBIF.UNRESOLVED <- GBIF.TRIM %>%

  ## Note that these filters are very forgiving...
  ## Unless we include the NAs, very few records are returned!
  filter(New.Taxonomic.status == 'Unresolved')


## Also keep the unresolved records:
GBIF.RESOLVED <- GBIF.TRIM %>%
  
  ## Just get the accepted taxa
  filter(New.Taxonomic.status == 'Accepted')


#########################################################################################################################
## What do examples of resolved and unresolved records look like?
head(GBIF.UNRESOLVED, 50)[, c("searchTaxon",
                              "scientificName", 
                              "TPL_binomial",
                              "taxo_agree")]


head(GBIF.RESOLVED, 50)[, c("searchTaxon",
                            "scientificName", 
                            "TPL_binomial",
                            "taxo_agree")]


## Also keep the managed records:
unique(GBIF.UNRESOLVED$New.Taxonomic.status)
dim(GBIF.UNRESOLVED)   ## 1.2 million unresolved records, quite a lot!



## Unique(GBIF.UNRESOLVED$New.Taxonomic.status)
saveRDS(GBIF.UNRESOLVED, file = paste("./data/base/HIA_LIST/GBIF/SUA_TREE_GBIF_UNRESOLVED.rds"))


## Just keep these columns:
## Taxonomic.status, Infraspecific.rank, New.Taxonomic.status, New.ID, New_binomial, taxo_agree
GBIF.TRIM.TAXO <- GBIF.TRIM %>%
  select(one_of(TPL.keep))


## Unique(GBIF.UNRESOLVED$New.Taxonomic.status)
saveRDS(GBIF.TAXO,       file = paste("./data/base/HIA_LIST/GBIF/SUA_TREES_GBIF_TAXO.rds"))
saveRDS(GBIF.TRIM.TAXO,  file = paste("./data/base/HIA_LIST/GBIF/SUA_TREES_GBIF_TRIM_TAXO.rds"))
write.csv(GBIF.TRIM.TAXO,             "./data/base/HIA_LIST/GBIF/SUA_TREES_GBIF_TRIM_TAXO.csv", row.names = FALSE)
write.csv(GBIF.LUT,                   "./data/base/HIA_LIST/GBIF/SUA_TREES_GBIF_LUT.csv",       row.names = FALSE)




# #########################################################################################################################
# ## 3). MARK CULTIVATED RECORDS
# #########################################################################################################################
# 
# 
# ## This only works on species lists with cultivated recrods..............................................................
# 
# 
# ## To create search terms, we can look for keywords in differemt languages: e.g. garden cultivated, etc.
# GBIF.COUNTRY = GBIF.TRIM.TAXO[, c("country")]
# GBIF.COUNTRY = as.data.frame(table(GBIF.COUNTRY))
# names(GBIF.COUNTRY) <-  c('country', 'count')
# GBIF.COUNTRY = GBIF.COUNTRY[order(GBIF.COUNTRY$count, decreasing = TRUE), ]
# 
# 
# ## Look for synonyms in the countries which have the most GBIF records
# #GBIF.COUNTRY
# #write.csv(GBIF.COUNTRY, "./data/base/HIA_LIST/GBIF/SPECIES/GBIF_COUNTRY.csv", row.names = FALSE)
# 
# 
# ## This creates a list of synonyms:
# #cultivated.synonyms
# 
# 
# #########################################################################################################################
# ## Now search for these words in particular columns
# ## Can these terms be searched for across the whole data frame (i.e. any column)?
# ## Also lot's of Australian species don't have these columns, which might make it tricky to run the clean 
# ## unique(GBIF.TRIM.TAXO$country)
# 
# ## Try using the big list of synonyms across all the data
# ## grepl("garden|cultiva",   GBIF.TRIM.TAXO$locality,           ignore.case = TRUE) | 
# GBIF.TRIM.TAXO$CULTIVATED <- ifelse(grepl(cultivated.synonyms,   GBIF.TRIM.TAXO$locality,           ignore.case = TRUE) |
#                                       grepl(cultivated.synonyms, GBIF.TRIM.TAXO$habitat,            ignore.case = TRUE) |
#                                       grepl(cultivated.synonyms, GBIF.TRIM.TAXO$eventRemarks,       ignore.case = TRUE) |
#                                       grepl(cultivated.synonyms, GBIF.TRIM.TAXO$cloc,               ignore.case = TRUE) |
#                                       grepl("managed",           GBIF.TRIM.TAXO$establishmentMeans, ignore.case = TRUE),
# 
#                                     "CULTIVATED", "UNKNOWN")
#  
#  
# # ## How many records are knocked out by using this definition?
# # ## This is probably a bit strict, in that for some of the fields, garden doesn't = cultivated
# GBIF.CULTIVATED = subset(GBIF.TRIM.TAXO, CULTIVATED == "CULTIVATED")
# dim(GBIF.CULTIVATED)[1]
 
 
## Still very few records being returned as "cultivated"?
# dim(GBIF.CULTIVATED)[1]/dim(GBIF.TRIM.TAXO)[1]
#View(GBIF.CULTIVATED)


## Also keep the cultivated records:
# GBIF.CULTIVATED <- GBIF.TRIM.TAXO %>% 
#   
#   ## Note that these filters are very forgiving...
#   filter(CULTIVATED == "CULTIVATED")


## Unique(GBIF.UNRESOLVED$New.Taxonomic.status)
# dim(GBIF.CULTIVATED)
# names(GBIF.CULTIVATED)
# unique(GBIF.CULTIVATED$CULTIVATED)
#saveRDS(GBIF.CULTIVATED, file = paste("./data/base/HIA_LIST/GBIF/GBIF_CULTIVATED.rds"))





#########################################################################################################################
## 4). CREATE TABLE OF PRE-CLEAN FLAGS AND FILTER RECORDS
#########################################################################################################################


#########################################################################################################################
## Create a table which counts the number of records meeting each criteria:
## Note that TRUE indicates there is a problem (e.g. if a record has no lat/long, it will = TRUE)
GBIF.TRIM.TAXO = GBIF.TRIM
GBIF.PROBLEMS <- with(GBIF.TRIM.TAXO,
                      
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
                        coordinateUncertaintyInMeters > 1000 & 
                          !is.na(coordinateUncertaintyInMeters)
                        
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


## Print table to screen: in .Rmd file, no need to save
## Note that TRUE indicates there is a problem (e.g. no lat/long = TRUE)
# str(GBIF.PROBLEMS)
# kable(GBIF.PROBLEMS)


## Quickly check the total record number matches the count of problems
#Total.count = sum(GBIF.PROBLEMS$COUNT)
#identical(dim(GBIF.TRIM.TAXO)[1], Total.count)  ## identical matches two objects


#########################################################################################################################
## Now filter the GBIF records using conditions which are not too restrictive
GBIF.CLEAN <- GBIF.TRIM.TAXO %>% 
  
  ## Note that these filters are very forgiving...
  ## Unless we include the NAs, very few records are returned!
  filter(!is.na(lon) & !is.na(lat),
         
         ## CULTIVATED == 'CULTIVATED', 
         ## #establishmentMeans!='MANAGED' | is.na(establishmentMeans),
         ## New.Taxonomic.status == 'Accepted'
         
         year >= 1950 & !is.na(year))


## The table above gives the details, but worth documenting how many records are knocked out by each filter
## What does this dataframe look like?
names(GBIF.CLEAN)





#########################################################################################################################
## 6). REMOVE POINTS OUTSIDE WORLDCLIM LAYERS...
#########################################################################################################################


## Can use WORLDCIM rasters to get only records where wordlclim data is. 
## At the global scale, there probably is no alterntive to using WORLDCLIM...


## First, get one of the BIOCLIM variables
world.temp = raster("./data/base/worldclim/world/0.5/bio/current/bio_01")
#plot(world.temp)


## Now get the XY centroids of the unique 1km * 1km WORLDCLIM blocks where GBIF records are found
## Get cell number(s) of WORLDCLIM raster from row and/or column numbers. Cell numbers start at 1 in the upper left corner, 
## and increase from left to right, and then from top to bottom. The last cell number equals the number of raster cells 
xy <- cellFromXY(world.temp, GBIF.CLEAN[c("lon", "lat")]) %>% 
  
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


## Now extract the temperature values for the unique 1km centroids which contain GBIF data
class(xy)
z   = extract(world.temp, xy)

# Warning message:
#   In .doExtract(x, i, ..., drop = drop) :
#   some indices are invalid (NA returned)

hist(z, border = NA, col = "orange", breaks = 50, main = "", xlab = "Worldclim Annual temp")


## Then track which values of Z are on land or not
onland = z %>% is.na %>%  `!` # %>% xy[.,]  cells on land or not
#summary(onland)


## Finally, filter the cleaned GBIF data to only those points on land. 
## This is achieved with the final [onland]
GBIF.LAND = filter(GBIF.CLEAN, cellFromXY(world.temp, GBIF.CLEAN[c("lon", "lat")]) %in% 
                     unique(cellFromXY(world.temp, GBIF.CLEAN[c("lon", "lat")]))[onland])


## how many records were on land?
records.ocean = dim(GBIF.CLEAN)[1] - dim(GBIF.LAND)[1]  ## 91575 records are in the ocean   


## Print the dataframe dimensions to screen
dim(GBIF.LAND)
length(unique(GBIF.LAND$searchTaxon))


## Add a source column
GBIF.LAND$SOURCE = 'GBIF'


## Free some memory
gc()


## note that this still leaves lots of points in the Islands. So we need to decide if those are legitimate.
# WORLD <- readOGR("./data/base/CONEXTUAL/TM_WORLD_BORDERS-0.3.shp", layer = "TM_WORLD_BORDERS-0.3")
# LAND  <- readOGR("./data/base/CONTEXTUAL/ne_10m_land.shp", layer = "ne_10m_land")
# plot(WORLD)
# plot(LAND)


#########################################################################################################################
## save data
#saveRDS(GBIF.LAND, file = paste("./data/base/HIA_LIST/GBIF/GBIF_LAND_POINTS.rds"))


## Now save .rds file for the next session
save.image("STEP_3_GBIF_CLEAN.RData")
#load("STEP_3_GBIF_CLEAN.RData")



#########################################################################################################################
## OUTSTANDING CLEANING TASKS:
#########################################################################################################################


## Check distribution maps for each species for spatial outliers  - Keep a spreasheet of all species...

## Estimate native/naturalised ranges as a separate colum         - APC data + Ubran polygons for AUS, USA, EU: only a subset

## GBIF taxonomic errors                                          - Use TPL

## Keep cultivated records as a separate column/file              - Get cultivated column from ALA data...

## Duplicates between GBIF and ALA                                - See email from CSIRO - only a problem for niches...





#########################################################################################################################
###################################################### TBC ############################################################## 
#########################################################################################################################