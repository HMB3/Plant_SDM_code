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
## 1). COMBINE GBIF DATAFRAMES INTO ONE
#########################################################################################################################


#########################################################################################################################
## Create a list of species from the downloaded files
gbif.download = list.files(GBIF_path, pattern = ".RData")


## what is causing the NA lat problem?
## Some of the species have lat/lon, and some have decimal lat/lon. Need to re-download the latest species


#########################################################################################################################
## Combine all the taxa into a single dataframe at once
GBIF.ALL <- gbif.download %>%   ## spp.download[c(1:length(spp.download))] 
  
  ## Pipe the list into lapply
  lapply(function(x) {
    
    ## Create a character string of each .RData file
    #f <- sprintf(paste0(GBIF_path, "%s_GBIF_records.RData"), x)
    f <- sprintf(paste0(GBIF_path, "%s"), x)
    
    ## Load each file
    d <- get(load(f))
    
    ## Now drop the columns which we don't need
    dat <- data.frame(searchTaxon = x, d[, colnames(d) %in% gbif.keep],
                      stringsAsFactors = FALSE)
    
    if(!is.character(dat$gbifID)) {
      
      dat$gbifID <- as.character(dat$gbifID)
      
    }
    
    ## Need to print the object within the loop
    names(dat)[names(dat) == 'decimalLatitude']  <- 'lat'
    names(dat)[names(dat) == 'decimalLongitude'] <- 'lon'
    dat$searchTaxon = gsub("_GBIF_records.RData", "", dat$searchTaxon)
    dat
    
  }) %>%
  
  ## Finally, bind all the rows together
  bind_rows


#########################################################################################################################
## what proportion of the dataset has no lat/lon?
dim(GBIF.ALL)
(sum(is.na(GBIF.ALL$lat))            + dim(subset(GBIF.ALL, lat == ""))[1])/dim(GBIF.ALL)[1]*100


## Almost none of the GBIF data has no scientificName. This is the right field to use for matching taxonomy
(sum(is.na(GBIF.ALL$species))        + dim(subset(GBIF.ALL, species == ""))       [1])/dim(GBIF.ALL)[1]*100
(sum(is.na(GBIF.ALL$name))           + dim(subset(GBIF.ALL, name == ""))          [1])/dim(GBIF.ALL)[1]*100
(sum(is.na(GBIF.ALL$scientificName)) + dim(subset(GBIF.ALL, scientificName == ""))[1])/dim(GBIF.ALL)[1]*100


## Now get just the columns we want to keep.
GBIF.TRIM <- GBIF.ALL %>% 
  select(one_of(gbif.keep))
gc()


## Check names
dim(GBIF.TRIM)
names(GBIF.TRIM)
intersect(names(GBIF.TRIM), gbif.keep)


#########################################################################################################################
## Just get the newly downloaded species
GBIF.TRIM = GBIF.TRIM[GBIF.TRIM$searchTaxon %in% GBIF.spp, ]
dim(GBIF.TRIM)


## What are the unique species?
length(unique(GBIF.TRIM$name))
length(unique(GBIF.TRIM$species))
length(unique(GBIF.TRIM$searchTaxon)) 
length(unique(GBIF.TRIM$scientificName))  
 




#########################################################################################################################
## 2). CHECK TAXONOMY RETURNED BY GBIF USING TAXONSTAND
######################################################################################################################### 


## The problems is the mis-match between what we searched, and what GBIF returned. 

## 1). Create the inital list by combing the planted trees with evergreen list
##     TREE.HIA.SPP = intersect(subset(TI.LIST, Plantings > 50)$searchTaxon, CLEAN.SPP$Binomial)

## 2). Clean this list using the GBIF backbone taxonomy :: use the "species" column

## 3). Run the GBIF "species" list through the TPL taxonomy. Take "New" Species and Genus as the "searchTaxon"

## 4). Use rgbif and ALA4R to download occurence data, using "searchTaxon". ALA is ok, because the taxonomy is resolved.
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
GBIF.TAXO <- TPL(unique(GBIF.TRIM$scientificName), infra = TRUE,
                 corr = TRUE, repeats = 100)  ## to stop it timing out...
sort(names(GBIF.TAXO))
saveRDS(GBIF.TAXO, 'data/base/HIA_LIST/COMBO/GBIF_TAXO_400.rds')
GBIF.TAXO = readRDS('data/base/HIA_LIST/COMBO/GBIF_TAXO_400.rds')


#########################################################################################################################
## The procedure used for taxonomic standardization is based on function TPLck. A progress bar
## indicates the proportion of taxon names processed so far. 


## Check the taxonomy by running scientificName through TPL. Then join the GBIF data to the taxonomic check, using 
## "scientificName" as the join field
GBIF.TRIM.TAXO <- GBIF.TRIM %>%
  left_join(., GBIF.TAXO, by = c("scientificName" = "Taxon"))
names(GBIF.TRIM.TAXO)


#########################################################################################################################
## However, the scientificName string and the searchTaxon string are not the same. 
## currently using 'str_detect'
Match.SN = GBIF.TRIM.TAXO  %>%
  mutate(Match.SN.ST = 
           str_detect(scientificName, searchTaxon)) %>%
  
  select(one_of(c("scientificName",
                  "searchTaxon",
                  "Taxonomic.status",
                  "New.Taxonomic.status",
                  "New.Genus",
                  "New.Species",
                  "country",
                  "Match.SN.ST")))


## How many records don't match?
dim(Match.SN)
unique(Match.SN$Taxonomic.status)
unique(Match.SN$New.Taxonomic.status)


#########################################################################################################################
## Incude records where the "scientificName" and the "searchTaxon" match, and where the taxonomic status is 
## accepted, synonym or unresolved


## Also include records where the "scientificName" and the "searchTaxon" don't match, but status is synonym
## This is the aame as the subset of species which are accpeted, but not on our list
match.true  = unique(subset(Match.SN, Match.SN.ST == "TRUE")$scientificName)
match.false = unique(subset(Match.SN, Match.SN.ST == "FALSE" &
                         Taxonomic.status == "Synonym")$scientificName)  
keep.SN     = unique(c(match.true, match.false))
length(keep.SN)


## Now create the match column in the main table of records
# GBIF.TRIM.TAXO$Match.SN = GBIF.TRIM.TAXO  %>%
#   mutate(Match.SN.ST = 
#            str_detect(scientificName, searchTaxon))
# unique(GBIF.TRIM.TAXO$Match.SN)


#########################################################################################################################
## Now remove these from the GBIF dataset?
GBIF.TRIM.MATCH = GBIF.TRIM.TAXO[GBIF.TRIM.TAXO$scientificName %in% keep.SN, ]


## How many records were removed?
message(dim(GBIF.TRIM.TAXO)[1] - dim(GBIF.TRIM.MATCH)[1], " records removed")
message(round((dim(GBIF.TRIM.MATCH)[1])/dim(GBIF.TRIM.TAXO)[1]*100, 2), 
        " % records retained using TPL mismatch")


## Check the taxonomic status of the updated table
unique(GBIF.TRIM.MATCH$Taxonomic.status)
unique(GBIF.TRIM.MATCH$New.Taxonomic.status)
round(with(GBIF.TRIM.MATCH, table(Taxonomic.status)/sum(table(Taxonomic.status))*100), 2)
round(with(GBIF.TRIM.MATCH, table(New.Taxonomic.status)/sum(table(New.Taxonomic.status))*100), 2)
#View(GBIF.TRIM.MATCH[c("searchTaxon", "scientificName", "Taxonomic.status", "New.Taxonomic.status")])





#########################################################################################################################
## 5). CREATE TABLE OF PRE-CLEAN FLAGS AND FILTER RECORDS
#########################################################################################################################


#########################################################################################################################
## Create a table which counts the number of records meeting each criteria:
## Note that TRUE indicates there is a problem (e.g. if a record has no lat/long, it will = TRUE)
#GBIF.TRIM.TAXO = GBIF.TRIM.TAXO
GBIF.PROBLEMS <- with(GBIF.TRIM.MATCH,
                      
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


#########################################################################################################################
## Now filter the GBIF records using conditions which are not too restrictive
GBIF.CLEAN <- GBIF.TRIM.MATCH %>% 
  
  ## Note that these filters are very forgiving...
  ## Unless we include the NAs, very few records are returned!
  filter(!is.na(lon) & !is.na(lat),
         
         ## CULTIVATED == 'CULTIVATED', 
         ## #establishmentMeans!='MANAGED' | is.na(establishmentMeans),
         ## New.Taxonomic.status == 'Accepted'
         
         year >= 1950 & !is.na(year))


## How many species are there?
names(GBIF.CLEAN)
length(unique(GBIF.TRIM.TAXO$searchTaxon))
length(unique(GBIF.CLEAN$searchTaxon))





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
names(GBIF.LAND)

## Free some memory
gc()


## note that this still leaves lots of points in the Islands. So we need to decide if those are legitimate.
# WORLD <- readOGR("./data/base/CONEXTUAL/TM_WORLD_BORDERS-0.3.shp", layer = "TM_WORLD_BORDERS-0.3")
# LAND  <- readOGR("./data/base/CONTEXTUAL/ne_10m_land.shp", layer = "ne_10m_land")
# plot(WORLD)
# plot(LAND)




#########################################################################################################################
## save data
saveRDS(GBIF.LAND, file = paste("./data/base/HIA_LIST/GBIF/GBIF_TREES_LAND.rds"))


## Now save .rds file for the next session
gc()




#########################################################################################################################
## OUTSTANDING CLEANING TASKS:
#########################################################################################################################


## Check taxonomy 


#########################################################################################################################
###################################################### TBC ############################################################## 
#########################################################################################################################