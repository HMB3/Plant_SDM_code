#########################################################################################################################
############################################# UPDATE ALA  ############################################################### 
#########################################################################################################################


#########################################################################################################################
## This code updates the ALA data
ala.download = list.files(ALA_path, pattern = ".RData")


## For now, use the old version of the ALA
## the base url for biocache downloads (used by offline occurrence downloads



#########################################################################################################################
## 1). COMBINE ALA SPECIES INTO ONE DATASET
#########################################################################################################################


## Some species fail because ALA does not recognise their taxonomy. E.G.:
## Tabebuia rosea, Roystonea regia


#########################################################################################################################
## Combine all the taxa into a single dataframe at once
ALA.ALL <- ala.download %>% 
  
  ## Pipe the list into lapply
  lapply(function(x) {
    
    ## Create a character string of each .RData file
    f <- sprintf(paste0(ALA_path, "%s"), x)
    
    ## Load each file - check if some are already dataframes 
    d <- get(load(f))
    if (length(class(d)) > 1) {
      
      d <- d[["data"]]
      
    }
    
    ## Create the searchTaxon column
    d[,"searchTaxon"] = x
    d[,"searchTaxon"] = gsub("_ALA_records.RData", "", d[,"searchTaxon"])
    
    ## Choose onl 
    d = d %>%
      select(one_of(ALA.keep))
    
    if(!is.character(d["id"])) {
      
      d["id"] <- as.character(d["id"])
      
    }
    
    ## This is a list of columns in different ALA files which have weird characters
    d[,"coordinateUncertaintyInMetres"] = as.numeric(unlist(d["coordinateUncertaintyInMetres"]))
    d["year"]  = as.numeric(unlist(d["year"]))
    d["month"] = as.numeric(unlist(d["month"]))
    d["id"]    = as.character(unlist(d["id"]))
    
    names(d)[names(d) == 'latitude']  <- 'lat'
    names(d)[names(d) == 'longitude'] <- 'lon'
    return(d)
    
  }) %>%
  
  ## Finally, bind all the rows together
  bind_rows


#########################################################################################################################
## Just get the newly downloaded species
ALA.ALL = ALA.ALL[ALA.ALL$searchTaxon %in% GBIF.spp, ]
length(unique(ALA.ALL$searchTaxon))


#########################################################################################################################
## What proportion of the ALA dataset has no lat/lon? 
sort(names(ALA.ALL))
dim(ALA.ALL)
(sum(is.na(ALA.ALL$lat))            + dim(subset(ALA.ALL, year < 1950))[1])/dim(ALA.ALL)[1]*100


## How many records have no species label? Need to check this so we know the latest download is working
(sum(is.na(ALA.ALL$scientificNameOriginal))  + dim(subset(ALA.ALL, scientificNameOriginal == ""))[1])/dim(ALA.ALL)[1]*100
(sum(is.na(ALA.ALL$scientificName))          + dim(subset(ALA.ALL, scientificName == ""))[1])/dim(ALA.ALL)[1]*100
(sum(is.na(ALA.ALL$species))                 + dim(subset(ALA.ALL, species == ""))[1])/dim(ALA.ALL)[1]*100


## What is the match between scientificNameOriginal, and the searchTaxon?
# Match.ALA = ALA.ALL  %>%
#   mutate(Match.SNO.ST =
#            str_detect(scientificNameOriginal, searchTaxon)) %>%
# 
#   mutate(Match.SN.ST =
#            str_detect(scientificName, searchTaxon)) %>%
# 
#   mutate(Match.SP.ST =
#            str_detect(species, searchTaxon)) %>%
#   
#   select(one_of(c("searchTaxon",
#                   "scientificName",
#                   "scientificNameOriginal",
#                   "species",
#                   "Match.SNO.ST",
#                   "Match.SN.ST",
#                   "Match.SP.ST")))
# View(Match.ALA)
# 
# 
# ## So for 15-12% of the records, neither the scientificNameOriginal or the species match the search taxon.
# dim(subset(Match.ALA,  Match.SNO.ST == "FALSE"))[1]/dim(Match.ALA)[1]*100
# dim(subset(Match.ALA,  Match.SN.ST  == "FALSE"))[1]/dim(Match.ALA)[1]*100
# dim(subset(Match.ALA,  Match.SN.ST  == "FALSE"))[1]/dim(Match.ALA)[1]*100


## So rename 'scientificNameOriginal' to 'scientificName', to match GBIF
# ALA.RENAME = ALA.ALL
# ALA.RENAME$scientificName = NULL 
# ALA.RENAME$scientificName = ALA.RENAME$scientificNameOriginal
# ALA.RENAME$species = NULL 
# (sum(is.na(ALA.RENAME$scientificName)) + dim(subset(ALA.RENAME, scientificName == ""))[1])/dim(ALA.ALL)[1]*100


########################################################################################################################
## Check names
## What names get returned?
sort(names(ALA.ALL))
ALA.TRIM <- ALA.ALL%>% 
  select(one_of(ALA.keep))

dim(ALA.TRIM)
sort(names(ALA.TRIM))


## What are the unique species?
length(unique(ALA.TRIM$searchTaxon))
length(unique(ALA.TRIM$scientificNameOriginal)) 
length(unique(ALA.TRIM$scientificName)) 
length(unique(ALA.TRIM$species))
(sum(is.na(ALA.TRIM$scientificName)) + dim(subset(ALA.TRIM, scientificName == ""))[1])/dim(ALA.TRIM)[1]*100





#########################################################################################################################
## 2). CHECK TAXONOMY RETURNED BY ALA USING TAXONSTAND
######################################################################################################################### 


## The problems is the mismatch between what we searched, and what GBIF returned. 


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
saveRDS(ALA.TREES.TAXO, paste0('data/base/HIA_LIST/COMBO/ALA_TAXO_', save_run, '.rds'))
ALA.TREES.TAXO = readRDS('data/base/HIA_LIST/COMBO/ALA_TAXO_OLD_ALA.rds')


## Check the taxonomy by running scientificName through TPL. Then join the GBIF data to the taxonomic check, using 
## "scientificName" as the join field
ALA.TRIM.TAXO <- ALA.TRIM %>%
  left_join(., ALA.TREES.TAXO, by = c("scientificName" = "Taxon"))
names(ALA.TRIM.TAXO)


## Check NAs again
(sum(is.na(ALA.TRIM.TAXO$scientificName)) + dim(subset(ALA.TRIM.TAXO, scientificName == ""))[1])/dim(ALA.TRIM)[1]*100


#########################################################################################################################
## However, the scientificName string and the searchTaxon string are not the same. 
## currently using 'str_detect'
Match.SN = ALA.TRIM.TAXO  %>%
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
## Now remove these from the ALA dataset
ALA.TRIM.MATCH  = ALA.TRIM.TAXO[ALA.TRIM.TAXO$scientificName %in% keep.SN, ]
Match.record    = Match.SN[Match.SN$scientificName %in% keep.SN, ]


## Check the taxonomic status
round(with(Match.record, table(Match.SN.ST)/sum(table(Match.SN.ST))*100), 2)
round(with(Match.record, table(Taxonomic.status)/sum(table(Taxonomic.status))*100), 2)
round(with(Match.record, table(New.Taxonomic.status)/sum(table(New.Taxonomic.status))*100), 2)


## How many records were removed by taxonomic filtering?
message(dim(ALA.TRIM.TAXO)[1] - dim(ALA.TRIM.MATCH)[1], " records removed")
message(round((dim(ALA.TRIM.MATCH)[1])/dim(ALA.TRIM.TAXO)[1]*100, 2), 
        " % records retained using TPL mismatch")


## Check the taxonomic status of the updated table
unique(ALA.TRIM.MATCH$Taxonomic.status)
unique(ALA.TRIM.MATCH$New.Taxonomic.status)

round(with(ALA.TRIM.MATCH, table(Taxonomic.status)/sum(table(Taxonomic.status))*100), 2)
round(with(ALA.TRIM.MATCH, table(New.Taxonomic.status)/sum(table(New.Taxonomic.status))*100), 2)


## Check NAs again
(sum(is.na(ALA.TRIM.MATCH$scientificName)) + dim(subset(ALA.TRIM.MATCH, scientificName == ""))[1])/dim(ALA.TRIM.MATCH)[1]*100





#########################################################################################################################
## 3). FILTER RECORDS TO THOSE WITH COORDINATES, AND AFTER 1950
#########################################################################################################################


#########################################################################################################################
## Now filter the ALA records using conditions which are not too restrictive
ALA.TREES.CLEAN <- ALA.TRIM.MATCH %>% 
  
  ## Note that these filters are very forgiving...
  ## Unless we include the NAs, very few records are returned!
  filter(!is.na(lon) & !is.na(lat),
         
         ## CULTIVATED == 'CULTIVATED', 
         ## #establishmentMeans!='MANAGED' | is.na(establishmentMeans),
         ## New.Taxonomic.status == 'Accepted'
         
         year >= 1950 & !is.na(year))


## Check
names(ALA.TREES.CLEAN)


## How many records were removed by filtering?
message(dim(ALA.TRIM.MATCH)[1] - dim(ALA.TREES.CLEAN)[1], " records removed")
message(round((dim(ALA.TREES.CLEAN)[1])/dim(ALA.TRIM.MATCH)[1]*100, 2), 
        " % records retained using spatially valid records")






#########################################################################################################################
## 4). REMOVE POINTS OUTSIDE WORLDCLIM LAYERS...
#########################################################################################################################


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


## For some reason, we need to convert the xy coords to a spatial points data frame, in order to avoid this error:
## 'NAs introduced by coercion to integer range'
xy <- SpatialPointsDataFrame(coords = xy, data = as.data.frame(xy),
                             proj4string = CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))


## Now extract the temperature values for the unique 1km centroids which contain ALA data
class(xy)
z   = raster::extract(world.temp, xy)
hist(z, border = NA, col = "orange", breaks = 50, main = "", xlab = "Worldclim Annual temp")


## Then track which values of Z are on land or not
onland = z %>% is.na %>%  `!` # %>% xy[.,]  cells on land or not


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
message(round((dim(ALA.TREES.LAND)[1])/dim(ALA.TREES.CLEAN)[1]*100, 2), 
        " % records retained using spatially valid records")


## Free some memory
gc()


#########################################################################################################################
## save data
dim(ALA.TREES.LAND)
length(unique(ALA.TREES.LAND$searchTaxon))


## One of the ALA columns is causing trouble. Reduce the ALA dataset to a minimum set
intersect(sort(names(GBIF.LAND)), sort(names(ALA.TREES.LAND)))
ALA.TREES.LAND = ALA.TREES.LAND[c("searchTaxon",      "scientificName", "SOURCE", 
                                  "lon", "lat",       "coordinateUncertaintyInMetres", "geodeticDatum", "year", "locality", "country", 
                                  "basisOfRecord",    "institutionCode", "rank",
                                  "Taxonomic.status", "New.Taxonomic.status", "New.Genus", "New.Species")]
saveRDS(ALA.TREES.LAND, paste0('data/base/HIA_LIST/ALA/ALA_TREES_LAND_', save_run, '.rds'))





#########################################################################################################################
## OUTSTANDING ALA TASKS:
#########################################################################################################################


## Re-run with updated ALA data


#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################