#########################################################################################################################
####################################################  FILTER GBIF DATA ################################################## 
#########################################################################################################################


## Aim of this code is to filer the GBIF data to the most reliable recrods. Two sources of uncertainty:

## Spatial 
## Taxonomic


## Which columns overlap between GBIF and ALA/AVH? The problem here is that we don't know what each GBIF field means,

## dataProvider             (ALA), GBIF ?
## basisOfRecord            (ALA), basisOfRecord (GBIF)
## occCultivatedEscapee     (ALA), establishmentMeans (GBIF)
## year
## taxonIdentificationIssue (ALA), GBIF no equivalent
## fatal assertions         (ALA), GBIF no equivalent


## A record with the same lat/long to x decimal places, collected in the same month by the same person is a duplicate


#########################################################################################################################
## LOAD
load("./data/base/HIA_LIST/COMBO/GBIF_TRIM.RData")
dim(GBIF.TRIM)
names(GBIF.TRIM)
length(unique(GBIF.TRIM$searchTaxon))  ## has the list updated with extra species? YES!





#########################################################################################################################
## 1). MARK CULTIVATED RECORDS AND CHECK TAXONOMY
#########################################################################################################################


## Consider whether this should be done before or after merging with ALA? Same applies for the taxonomic check


## Try and make this one command: can these terms be searched for across the whole data frame (i.e. any column)
## Also lot's of Australian species don't have these columns, which might make it tricky to run the clean 
GBIF.TRIM$CULTIVATED <- ifelse(grepl("garden|cultiva",   GBIF.TRIM$locality,           ignore.case = TRUE) | 
                                 grepl("garden|cultiva", GBIF.TRIM$habitat,            ignore.case = TRUE) | 
                                 grepl("garden|cultiva", GBIF.TRIM$eventRemarks,       ignore.case = TRUE) |
                                 grepl("garden|cultiva", GBIF.TRIM$cloc,               ignore.case = TRUE) |
                                 grepl("managed",        GBIF.TRIM$establishmentMeans, ignore.case = TRUE),
                               
                               "CULTIVATED", "UNKNOWN")


## How many records are knocked out by using this definition?
## This is probably a bit strict, in that for some of the fields, garden doesn't = cultivated
GBIF.CULTIVATED = subset(GBIF.TRIM, CULTIVATED == "CULTIVATED")
dim(GBIF.CULTIVATED)[1]/dim(GBIF.TRIM)[1]
View(GBIF.CULTIVATED)


## Also keep the cultivated records:
GBIF.CULTIVATED <- GBIF.TRIM %>% 
  
  ## Note that these filters are very forgiving...
  filter(CULTIVATED == "CULTIVATED")


## Unique(GBIF.UNRESOLVED$New.Taxonomic.status)
dim(GBIF.CULTIVATED)
names(GBIF.CULTIVATED)
unique(GBIF.CULTIVATED$CULTIVATED)
save(GBIF.CULTIVATED, file = paste("./data/base/HIA_LIST/GBIF/GBIF_CULTIVATED.RData"))


#########################################################################################################################
## CHECK TAXONOMY THAT GBIF RETURNED 
## Use "Taxonstand". This also assumes that the ALA data is clean
GBIF.TAXO <- TPL(unique(GBIF.TRIM$scientificName), infra = TRUE,
                 corr = TRUE, repeats = 100)
sort(names(GBIF.TAXO))


# Error in TPLck(sp = d, infra = infra, corr = corr, diffchar = diffchar,  : 
#                  Cannot read TPL website.
# : 'The operation timed out'


# The procedure used for taxonomic standardization is based on function TPLck. A progress bar
# indicates the proportion of taxon names processed so far. In case the TPL website cannot be reached
# temporarily, the function returns an error but repeats trying to match the given taxon multiple times
# (see repeats). If standardization is still not successful, the input taxon is returned in field ’Taxon’
# with NA in all other fields


## Then join the GBIF data to the taxonomic check, using "scientificName" as the join field... 
GBIF.TRIM <- GBIF.TRIM %>%
  left_join(., GBIF.TAXO, by = c("scientificName" = "Taxon"))


## So we can filter by the agreement between "scientificName", and "New.Taxonomic.status"?
## A better way could be to create a new filter for agreement between the search taxon genus and species,
## and the new genus and species. Which of these columns to keep?


#########################################################################################################################
## Now create a column for the agreement between the new genus and the old genus
## First, trim the spaces out
GBIF.TRIM$searchTaxon  = trim.space(GBIF.TRIM$searchTaxon)


## Then combine the genus and species returned by TPL into
GBIF.TRIM$TPL_binomial  = with(GBIF.TRIM, paste(New.Genus, New.Species, sep = " "))


## Now match the searchTaxon with the binomial returned by TPL
GBIF.TRIM$taxo_agree <- ifelse(
  GBIF.TRIM$searchTaxon == GBIF.TRIM$New_binomial, TRUE, FALSE)


## Also keep the managed records:
GBIF.UNRESOLVED <- GBIF.TRIM %>% 
  
  ## Note that these filters are very forgiving...
  ## Unless we include the NAs, very few records are returned!
  filter(New.Taxonomic.status == 'Unresolved')


## Also keep the managed records:
unique(GBIF.UNRESOLVED$New.Taxonomic.status)
dim(GBIF.UNRESOLVED)   ## 1.2 million unresolved records, quite a lot!


## Unique(GBIF.UNRESOLVED$New.Taxonomic.status)
save(GBIF.UNRESOLVED, file = paste("./data/base/HIA_LIST/GBIF/GBIF_UNRESOLVED.RData"))


## Just keep these columns:
## Taxonomic.status, Infraspecific.rank, New.Taxonomic.status, New.ID, New_binomial, taxo_agree
GBIF.TRIM <- GBIF.TIM %>% 
  select(one_of(TPL.keep))
names(GBIF.TRIM)





#########################################################################################################################
## 2). CREATE TABLE OF PRE-CLEAN FLAGS 
#########################################################################################################################


#########################################################################################################################
## Then create a table which counts the number of records meeting each criteria:
## Note that TRUE indicates there is a problem (e.g. if a record has no lat/long, it will = TRUE)
GBIF.PROBLEMS <- with(GBIF.TRIM,
                      
                      table(
                        
                        ## Note this list is incomplete
                        
                        ## No coordinates
                        is.na(lon)|is.na(lat),
                        
                        ## Taxon rank is genus/form?
                        #taxonRank == 'GENUS' & 'FORM',
                        
                        ## Taxonomic status: consider if this is the right filter, it seems to restrictive
                        New.Taxonomic.status == 'Unresolved',
                        
                        ## Cultivated
                        CULTIVATED == 'CULTIVATED',
                        
                        ## Establishment means is "MANAGED" is included in the above
                        #establishmentMeans == 'MANAGED' & !is.na(establishmentMeans),
                        
                        ## Collected before 1950 or na.year
                        year < 1950 & !is.na(year),
                        
                        ## No year
                        is.na(year),
                        
                        ## Coordinate uncertainty is > 100 or is not NA
                        coordinateUncertaintyInMeters > 1000 & 
                          !is.na(coordinateUncertaintyInMeters)
                        
                        ## Add maybe the centre of Australia?
                        ## Lamber centre of Aus: 25.610111, 134.354806
                        
                      )
                      
) %>% 
  
  ## Create a data frame and set the names
  as.data.frame %>%  
  
  setNames(c('NO_COORD',     'TAXON_STATUS', 'CULTIVATED', 
             'PRE_1950',     'NO_YEAR', 
             'COORD_UNCERT', 'COUNT'))


## Print table to screen: in .Rmd file, no need to save
## Note that TRUE indicates there is a problem (e.g. no lat/long = TRUE)
str(GBIF.PROBLEMS)
kable(GBIF.PROBLEMS)


## Quickly check the total record number matches the count of problems
Total.count = sum(GBIF.PROBLEMS$COUNT)
identical(dim(GBIF.TRIM)[1], Total.count)  ## identical matches two objects





#########################################################################################################################
## 3). FILTER RECORDS 
#########################################################################################################################


## Filter the GBIF records using conditions which are not too restrictive
dim(GBIF.TRIM)
GBIF.CLEAN <- GBIF.TRIM %>% 
  
  ## Note that these filters are very forgiving...
  ## Unless we include the NAs, very few records are returned!
  filter(!is.na(lon) & !is.na(lat),
         
         ## CULTIVATED == 'CULTIVATED', 
         ## #establishmentMeans!='MANAGED' | is.na(establishmentMeans),
         
         year >= 1950 & !is.na(year),
         New.Taxonomic.status == 'Accepted')


## The table above gives the details, but worth documenting how many records are knocked out by each filter
## Consider what these filters accomplish. Is it really worth knocking out that many records automatically?
dim(GBIF.CLEAN)

Remaining.percent = dim(GBIF.CLEAN)[1]/Total.count*100
Filters.applied = "NA COORD | < 1950/NA | UNRESOLVED TAXONOMY" ## INCLUDE CULTIVATED RECORDS
Remaining.percent ## 57% of records remain after cleaning 
gc()


## Check
dim(GBIF.CLEAN)
head(GBIF.CLEAN)


## What does this dataframe look like?
names(GBIF.CLEAN)





#########################################################################################################################
## 4). REMOVE POINTS OUTSIDE WORLDCLIM LAYERS...
#########################################################################################################################


## Can use WORLDCIM rasters to get only records where wordlclim data is. 
## At the global scale, there probably is no alterntive to using WORLDCLIM...


## First get one of the BIOCLIM variables
world.temp = raster("//sci-7910/F/data/worldclim/world/0.5/bio/current/bio_01")
plot(world.temp)


## Now get the XY centroids of the unique 1km * 1km WORLDCLIM blocks where GBIF records are found
## Get cell number(s) of WORLDCLIM raster from row and/or column numbers. Cell numbers start at 1 in the upper left corner, 
## and increase from left to right, and then from top to bottom. The last cell number equals the number of raster cells 
xy <- cellFromXY(world.temp, GBIF.CLEAN[c("lon", "lat")]) %>% 
  
  ## get the unique raster cells
  unique %>% 
  
  ## Get coordinates of the center of raster cells for a row, column, or cell number of WORLDCLIM raster
  xyFromCell(world.temp, .)


## Take a look at xy: NA's should be removed...
summary(xy)
str(xy)
points(xy, pch = ".", col = "red")


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
summary(onland)


## Finally, filter the cleaned GBIF data to only those points on land. 
## This is achieved with the final [onland]
GBIF.LAND = filter(GBIF.CLEAN, cellFromXY(world.temp, GBIF.CLEAN[c("lon", "lat")]) %in% 
                     unique(cellFromXY(world.temp, GBIF.CLEAN[c("lon", "lat")]))[onland])


## how many records were on land?
records.ocean = dim(GBIF.CLEAN)[1] - dim(GBIF.LAND)[1]  ## 91575 records are in the ocean   
records.ocean


## Free some memory
gc()


## note that this still leaves lots of points in the Islands. So we need to decide if those are legitimate.
WORLD <- readOGR("./data/base/CONTEXTUAL/TM_WORLD_BORDERS-0.3.shp", layer = "TM_WORLD_BORDERS-0.3")
LAND  <- readOGR("./data/base/CONTEXTUAL/ne_10m_land.shp",          layer = "ne_10m_land")
plot(WORLD)
plot(LAND)



#########################################################################################################################
## save data
save(GBIF.LAND,        file = paste("./data/base/HIA_LIST/GBIF/GBIF_LAND_POINTS.RData"))
gc()


## Now save .RData file for the next session
save.image("STEP_3_GBIF_CLEAN.RData")





#########################################################################################################################
## OUTSTANDING CLEANING TASKS:
#########################################################################################################################


## When should the additional filters be run in? before or after ALA merge?       - Ask John

## GBIF taxonomic errors                                                          - Use taxonstand, check with John, Dave K.

## Keep cultivated records as a separate column/file                              - Have them: check search with Rach, Linda

## Estimate native ranges: as a separate colum too?                               - Ubran polygons for AUS, USA, EU? Ask Dave Kendall

## GBIF duplicates                                                                - Check GBIF issues, but John's SDM code will get rid of more

## GBIF spatial outliers: Ocean, middle of Australia, etc. pp? Duplicated?        - Check each map, email Alex Zika RE geoclean

## Duplicates between GBIF and ALA                                                - See email from CSIRO, check in again with them





#########################################################################################################################
###################################################### TBC ############################################################## 
#########################################################################################################################