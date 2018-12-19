#########################################################################################################################
############################################# UPDATE ALA  ############################################################### 
#########################################################################################################################


## This code combines all the the ALA downloads for each species into one file
## The issue now becomes how to deal with different species runs that overlap. 
## Can all the species records from multiple runs be combined? Otherwise R is very slow generating the big table of species


#########################################################################################################################
## For now, use the old version of the ALA
## the base url for biocache downloads (used by offline occurrence downloads)
# install.packages("devtools")
# devtools::install_github("AtlasOfLivingAustralia/ALA4R")
# devtools::build(vignettes = FALSE)


#########################################################################################################################
## 1). COMBINE ALA SPECIES INTO ONE DATASET
#########################################################################################################################


## Print the species run to the screen
message('Combining ALA occurrence data for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")


#########################################################################################################################
## Create a list of Rdata files
ala.download = list.files(ALA_path, pattern = ".RData")
length(ala.download)


## Now these lists are getting too long for the combine step. 
## Restrict them to just the strings that partially match the  species list for each run
ala.spp.download <- paste(GBIF.spp, "_ALA_records.RData", sep = "")
ala.download     = ala.download[ala.download %in% ala.spp.download ] 
message('downloaded species ', length(ala.download), ' analyzed species ', length(GBIF.spp))


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
      
    } else {
      
      d = data.frame("searchTaxon" = c())
      
    }
    
    ## Check if the dataframes have data
    if (nrow(d) <= 2) {
      
      ## If the species has < 2 records, escape the loop
      #print (paste ("No ALA records for ", x, " skipping "))
      return (d)
      
    }
    
    ## Create the searchTaxon column
    message ('Reading ALA data for ', x)
    d[,"searchTaxon"] = x
    d[,"searchTaxon"] = gsub("_ALA_records.RData", "", d[,"searchTaxon"])
    
    ## Choose only the desired columns
    d = d %>%
      select(one_of(ALA.keep))
    
    if(!is.character(d["id"])) {
      
      d["id"] <- as.character(d["id"])
      
    }
    
    ## Then print warnings
    warnings()
    
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


## Total original records for 176 species  1141520 + 1033700 = 2,175,220


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
message('Running TPL taxonomy for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")
ALA.TAXO <- TPL(unique(ALA.TRIM$scientificName), infra = TRUE,
                      corr = TRUE, repeats = 100)  ## to stop it timing out...
sort(names(ALA.TAXO))
saveRDS(ALA.TAXO, paste0(ALA_path, 'ALA_TAXO_', save_run, '.rds'))


## Check the taxonomy by running scientificName through TPL. Then join the GBIF data to the taxonomic check, using 
## "scientificName" as the join field
ALA.TRIM.TAXO <- ALA.TRIM %>%
  left_join(., ALA.TAXO, by = c("scientificName" = "Taxon"))
ALA.TRIM.TAXO$New_binomial = paste(ALA.TRIM.TAXO$New.Genus, ALA.TRIM.TAXO$New.Species, sep = " ") 
names(ALA.TRIM.TAXO)


## Check NAs again
(sum(is.na(ALA.TRIM.TAXO$scientificName)) + dim(subset(ALA.TRIM.TAXO, scientificName == ""))[1])/dim(ALA.TRIM)[1]*100


#########################################################################################################################
## However, the scientificName string and the searchTaxon string are not the same. 
## currently using 'str_detect'
Match.SN = ALA.TRIM.TAXO  %>%
  mutate(Match.SN.ST = 
           str_detect(searchTaxon, New_binomial)) %>%  ## scientificName, New_binomial
  
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
unique(Match.SN$Match.SN.ST)
unique(Match.SN$Taxonomic.status)
unique(Match.SN$New.Taxonomic.status)


#########################################################################################################################
## Incude records where the "scientificName" and the "searchTaxon" match, and where the taxonomic status is 
## accepted, synonym or unresolved


## Also include records where the "scientificName" and the "searchTaxon" don't match, but status is synonym
## Also, the ALA taxonomy is right for some species - catch this with status = "accepted"
## This is the same as the subset of species which are accpeted, but not on our list
match.true  = unique(subset(Match.SN, Match.SN.ST == "TRUE")$scientificName)
match.false = unique(subset(Match.SN, Match.SN.ST == "FALSE" &
                              Taxonomic.status == "Synonym" |
                              Taxonomic.status == "Accepted" )$scientificName)  
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
ALA.CLEAN <- ALA.TRIM.MATCH %>% 
  
  ## Note that these filters are very forgiving...
  ## Unless we include the NAs, very few records are returned!
  filter(!is.na(lon) & !is.na(lat),
         
         ## CULTIVATED == 'CULTIVATED', 
         ## #establishmentMeans!='MANAGED' | is.na(establishmentMeans),
         ## New.Taxonomic.status == 'Accepted'
         
         year >= 1950 & !is.na(year))

## Check
names(ALA.CLEAN)


## How many records were removed by filtering?
message(dim(ALA.TRIM.MATCH)[1] - dim(ALA.CLEAN)[1], " records removed")
message(round((dim(ALA.CLEAN)[1])/dim(ALA.TRIM.MATCH)[1]*100, 2), 
        " % records retained using spatially valid records")






#########################################################################################################################
## 4). REMOVE POINTS OUTSIDE WORLDCLIM LAYERS...
#########################################################################################################################


#########################################################################################################################
## Can use WORLDCIM rasters to get only records where wordlclim data is. 
message('Removing ALA points in the ocean for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")

#str (world.grids.current)


## Now get the XY centroids of the unique 1km * 1km WORLDCLIM blocks where ALA records are found
## Get cell number(s) of WORLDCLIM raster from row and/or column numbers. Cell numbers start at 1 in the upper left corner, 
## and increase from left to right, and then from top to bottom. The last cell number equals the number of raster cells 
xy <- cellFromXY(world.grids.current, ALA.CLEAN[c("lon", "lat")]) %>% 
  
  ## get the unique raster cells
  unique %>% 
  
  ## Get coordinates of the center of raster cells for a row, column, or cell number of WORLDCLIM raster
  xyFromCell(world.grids.current, .)


## For some reason, we need to convert the xy coords to a spatial points data frame, in order to avoid this error:
## 'NAs introduced by coercion to integer range'
xy <- SpatialPointsDataFrame(coords = xy, data = as.data.frame(xy),
                             proj4string = CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))


## Now extract the temperature values for the unique 1km centroids which contain ALA data
class(xy)
z   = raster::extract(world.grids.current[["bio_01"]], xy)
#hist(z, border = NA, col = "orange", breaks = 50, main = "", xlab = "Worldclim Annual temp")


## Then track which values of Z are on land or not
onland = z %>% is.na %>%  `!` # %>% xy[.,]  cells on land or not


## Finally, filter the cleaned ALA data to only those points on land. 
## This is achieved with the final [onland]
ALA.LAND = filter(ALA.CLEAN, cellFromXY(world.grids.current, ALA.CLEAN[c("lon", "lat")]) %in% 
                          unique(cellFromXY(world.grids.current, ALA.CLEAN[c("lon", "lat")]))[onland])


## how many records were on land?
records.ocean = dim(ALA.CLEAN)[1] - dim(ALA.LAND)[1]  ## 91575 records are in the ocean   


## Print the dataframe dimensions to screen
dim(ALA.LAND)
length(unique(ALA.LAND$searchTaxon))


## Add a source column
ALA.LAND$SOURCE = 'ALA'
message(round((dim(ALA.LAND)[1])/dim(ALA.CLEAN)[1]*100, 2), 
        " % records retained using spatially valid records")


## Free some memory
gc()


#########################################################################################################################
## save data
dim(ALA.LAND)
length(unique(ALA.LAND$searchTaxon))


## One of the ALA columns is causing trouble. Reduce the ALA dataset to a minimum set
sort(names(ALA.LAND))
ALA.LAND = ALA.LAND[c("searchTaxon",      "scientificName", "SOURCE", 
                                  "lon", "lat",       "coordinateUncertaintyInMetres", #"geodeticDatum", 
                                  "year", "locality", "country", 
                                  "basisOfRecord",    "institutionCode", "rank",
                                  "Taxonomic.status", "New.Taxonomic.status", "New.Genus", "New.Species")]


#########################################################################################################################
## save data
if(save_data == "TRUE") {
  
  ## save .rds file for the next session
  saveRDS(ALA.LAND, paste0(ALA_path, 'ALA_LAND_', save_run, '.rds'))
  
} else {
  
  message(' skip file saving, not many species analysed')   ##
  
}

## get rid of some memory
gc()






#########################################################################################################################
## OUTSTANDING ALA TASKS:
#########################################################################################################################


## Re-run with updated ALA data


#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################
