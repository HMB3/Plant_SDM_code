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


## Print the species run to the screen
message('Combining GBIF occurrence data for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")


#########################################################################################################################
## Create a list of species from the downloaded files
gbif.download = list.files(GBIF_path, pattern = ".RData")
length(gbif.download)


## Now these lists are getting too long for the combine step. Can we restrict them to just the strings that partially 
## match the  species list for each run?
gbif.spp.download <- paste(GBIF.spp, "_ALA_records.RData", sep = "")
gbif.download     = gbif.download[gbif.download %in% gbif.spp.download ] 
message('downloaded species ', length(ala.download), ' analyzed species ', length(GBIF.spp))



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
    message ('Reading GBIF data for ', x)
    dat <- data.frame(searchTaxon = x, d[, colnames(d) %in% gbif.keep],
                      stringsAsFactors = FALSE)
    
    if(!is.character(dat$gbifID)) {
      
      dat$gbifID <- as.character(dat$gbifID)
      
    }
    
    ## Need to print the object within the loop
    names(dat)[names(dat) == 'decimalLatitude']  <- 'lat'
    names(dat)[names(dat) == 'decimalLongitude'] <- 'lon'
    dat$searchTaxon = gsub("_GBIF_records.RData", "", dat$searchTaxon)
    return(dat)
    
  }) %>%
  
  ## Finally, bind all the rows together
  bind_rows


#########################################################################################################################
## What proportion of the dataset has no lat/lon? Need to check this so we know the latest download is working
formatC(dim(GBIF.ALL)[1], format = "e", digits = 2)
(sum(is.na(GBIF.ALL$lat))            + dim(subset(GBIF.ALL, year < 1950))[1])/dim(GBIF.ALL)[1]*100


## Almost none of the GBIF data has no scientificName. This is the right field to use for matching taxonomy
(sum(is.na(GBIF.ALL$species))        + dim(subset(GBIF.ALL, species == ""))       [1])/dim(GBIF.ALL)[1]*100
(sum(is.na(GBIF.ALL$scientificName)) + dim(subset(GBIF.ALL, scientificName == ""))[1])/dim(GBIF.ALL)[1]*100


## Now get just the columns we want to keep.
GBIF.TRIM <- GBIF.ALL %>% 
  select(one_of(gbif.keep))
names(GBIF.TRIM)
gc()


#########################################################################################################################
## Just get the newly downloaded species
GBIF.TRIM = GBIF.TRIM[GBIF.TRIM$searchTaxon %in% GBIF.spp, ]
formatC(dim(GBIF.TRIM)[1], format = "e", digits = 2)


## What are the unique species?
length(unique(GBIF.TRIM$species))
length(unique(GBIF.TRIM$searchTaxon)) 
length(unique(GBIF.TRIM$scientificName))  
 




#########################################################################################################################
## 2). CHECK TAXONOMY RETURNED BY GBIF USING TAXONSTAND
######################################################################################################################### 


## The problems is the mis-match between what we searched, and what GBIF returned. 

## 1). Create the inital list by combing the planted trees with evergreen list
##     TREE.HIA.SPP = intersect(subset(TI.LIST, Plantings > 50)$searchTaxon, CLEAN.SPP$Binomial) 
##     This is about 400 species
##     Origin
##     NA     Exotic Native 
##     20.5   25.8   53.8 

## 2). Clean this list using the GBIF backbone taxonomy :: use the "species" column in from the GBIF "species lookup" tool
##     https://www.gbif.org/tools/species-lookup

## 3). Run the GBIF "species" list through the TPL taxonomy. Take "New" Species and Genus as the "searchTaxon"

## 4). Use rgbif and ALA4R to download occurence data, using "searchTaxon".
##     For GBIF, we use
##     key  <- name_backbone(name = sp.n, rank = 'species')$usageKey
##     GBIF <- occ_data(taxonKey = key, limit = GBIF.download.limit)
##     ALA  <- occurrences(taxon = sp.n, download_reason_id = 7)
    

##     This returns multiple keys and synonyms, but there is no simple way to skip these


## 5). Join the TPL taxonomy to the "scientificName" field. We can't use "name" (the equivalent of "species", it seems),
##     because name is always the same as the searchTaxon and not reliable (i.e. they will always match, and we know that
##     no one has gone through and checked each one.
     
##     Exclude records where the "scientificName" both doesn't match the "searchTaxon", and, also is not a synonym according to TPL
##     The remaining records are either "accepted" "synonym" or "uresolved", with 97% of searched records matching returned records.
     
##     Then we model these records as before. Of ~400 species we downloaded, we will pick the 200 with the best maps.
##     One line in the MS : we matched the GBIF backbone taxo against the TPL taxo, and searched the ALA and GBIF for the currently 
##     accepted names, etc.



#########################################################################################################################
## Use "Taxonstand" to check the taxonomy :: which field to use?
GBIF.TAXO <- TPL(unique(GBIF.TRIM$scientificName), infra = TRUE,
                 corr = TRUE, repeats = 100)  ## to stop it timing out...
sort(names(GBIF.TAXO))

saveRDS(GBIF.TAXO, paste0('data/base/HIA_LIST/COMBO/GBIF_TAXO_', save_run, '.rds'))
#GBIF.TAXO = readRDS('data/base/HIA_LIST/COMBO/GBIF_TAXO_200.rds')


#########################################################################################################################
## The procedure used for taxonomic standardization is based on function TPLck. A progress bar
## indicates the proportion of taxon names processed so far. 


## Check the taxonomy by running scientificName through TPL. Then join the GBIF data to the taxonomic check, using 
## "scientificName" as the join field
GBIF.TRIM.TAXO <- GBIF.TRIM %>%
  left_join(., GBIF.TAXO, by = c("scientificName" = "Taxon"))
names(GBIF.TRIM.TAXO)


## Check NAs again
(sum(is.na(GBIF.TRIM.TAXO$scientificName)) + dim(subset(GBIF.TRIM.TAXO, scientificName == ""))[1])/dim(GBIF.TRIM)[1]*100


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


#########################################################################################################################
## Now remove these from the GBIF dataset?
GBIF.TRIM.MATCH = GBIF.TRIM.TAXO[GBIF.TRIM.TAXO$scientificName %in% keep.SN, ]
Match.record    = Match.SN[Match.SN$scientificName %in% keep.SN, ]


round(with(Match.record, table(Match.SN.ST)/sum(table(Match.SN.ST))*100), 2)
round(with(Match.record, table(Taxonomic.status)/sum(table(Taxonomic.status))*100), 2)
round(with(Match.record, table(New.Taxonomic.status)/sum(table(New.Taxonomic.status))*100), 2)


## How many records were removed?
message(dim(GBIF.TRIM.TAXO)[1] - dim(GBIF.TRIM.MATCH)[1], " records removed")
message(round((dim(GBIF.TRIM.MATCH)[1])/dim(GBIF.TRIM.TAXO)[1]*100, 2), 
        " % records retained using TPL mismatch")


## Check the taxonomic status of the updated table
unique(GBIF.TRIM.MATCH$Taxonomic.status)
unique(GBIF.TRIM.MATCH$New.Taxonomic.status)

round(with(GBIF.TRIM.MATCH, table(Taxonomic.status)/sum(table(Taxonomic.status))*100), 2)
round(with(GBIF.TRIM.MATCH, table(New.Taxonomic.status)/sum(table(New.Taxonomic.status))*100), 2)


## Check NAs again
(sum(is.na(GBIF.TRIM.MATCH$scientificName)) + dim(subset(GBIF.TRIM.MATCH, scientificName == ""))[1])/dim(GBIF.TRIM.MATCH)[1]*100
View(GBIF.TRIM.MATCH[is.na(GBIF.TRIM.MATCH$scientificName),])





#########################################################################################################################
## 3). FILTER RECORDS TO THOSE WITH COORDINATES, AND AFTER 1950
#########################################################################################################################


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
## 4). REMOVE POINTS OUTSIDE WORLDCLIM LAYERS...
#########################################################################################################################


## Can use WORLDCIM rasters to get only records where wordlclim data is. 
## At the global scale, there probably is no alterntive to using WORLDCLIM...


## First, get one of the BIOCLIM variables
world.temp = raster("./data/base/worldclim/world/0.5/bio/current/bio_01")


## Now get the XY centroids of the unique 1km * 1km WORLDCLIM blocks where GBIF records are found
## Get cell number(s) of WORLDCLIM raster from row and/or column numbers. Cell numbers start at 1 in the upper left corner, 
## and increase from left to right, and then from top to bottom. The last cell number equals the number of raster cells 
xy <- cellFromXY(world.temp, GBIF.CLEAN[c("lon", "lat")]) %>% 
  
  ## get the unique raster cells
  unique %>% 
  
  ## Get coordinates of the center of raster cells for a row, column, or cell number of WORLDCLIM raster
  xyFromCell(world.temp, .)


## For some reason, we need to convert the xy coords to a spatial points data frame, in order to avoid this error:
## 'NAs introduced by coercion to integer range'
xy <- SpatialPointsDataFrame(coords = xy, data = as.data.frame(xy),
                             proj4string = CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))


## Now extract the temperature values for the unique 1km centroids which contain GBIF data
message('Removing GBIF points in the ocean for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")
class(xy)
z   = raster::extract(world.temp, xy)
hist(z, border = NA, col = "orange", breaks = 50, main = "", xlab = "Worldclim Annual temp")


## Then track which values of Z are on land or not
onland = z %>% is.na %>%  `!` # %>% xy[.,]  cells on land or not


## Finally, filter the cleaned GBIF data to only those points on land. 
## This is achieved with the final [onland]
GBIF.LAND = filter(GBIF.CLEAN, cellFromXY(world.temp, GBIF.CLEAN[c("lon", "lat")]) %in% 
                     unique(cellFromXY(world.temp,    GBIF.CLEAN[c("lon", "lat")]))[onland])


## how many records were on land?
records.ocean = dim(GBIF.CLEAN)[1] - dim(GBIF.LAND)[1]  ## 91575 records are in the ocean   


## Print the dataframe dimensions to screen
dim(GBIF.LAND)
length(unique(GBIF.LAND$searchTaxon))


## Add a source column
GBIF.LAND$SOURCE = 'GBIF'
unique(GBIF.LAND$SOURCE)
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
saveRDS(GBIF.LAND, paste0('data/base/HIA_LIST/GBIF/GBIF_TREES_LAND_', save_run, '.rds'))


## Now save .rds file for the next session
gc()




#########################################################################################################################
## OUTSTANDING CLEANING TASKS:
#########################################################################################################################


## Increase the taxonomic check to include all species on HIA list


#########################################################################################################################
###################################################### TBC ############################################################## 
#########################################################################################################################