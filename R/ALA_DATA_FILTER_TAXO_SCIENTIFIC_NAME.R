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

#  a spot of paranoia
gc()


#########################################################################################################################
## Combine all the taxa into a single dataframe at once
ALA.ALL <- ala.download %>% 
  
  ## Pipe the list into lapply
  lapply(function(x) {
    
    ## Create a character string of each .RData file
    message("Reading ALA data for ", x)
    f <- sprintf(paste0(ALA_path, "%s"), x)
    
    ## Load each file - check if some are already dataframes 
    d <- get(load(f))
    if (length(class(d)) > 1) {
      
      d <- d[["data"]]
      
    } else {
      
      d = d #data.frame("searchTaxon" = c())
      
    }
    
    ## Check if the dataframes have data
    if (nrow(d) <= 2) {
      
      ## If the species has < 2 records, escape the loop
      print (paste ("No ALA records for ", x, " skipping "))
      return (d)
      
    }
    
    #  type standardisation    
    names(d)[names(d) == 'latitude']  <- 'lat'
    names(d)[names(d) == 'longitude'] <- 'lon'
    
    #  standardi[sz]e catnum colname 
    if("catalogueNumber" %in% colnames(d)) {
      message ("Renaming catalogueNumber column to catalogNumber")
      names(d)[names(d) == 'catalogueNumber'] <- 'catalogNumber'
      
    }
    
    if (!is.character(d$catalogNumber)) {
      d$catalogNumber = as.character(d$catalogNumber)
      
    }
    
    #  standardi[sz]e catnum colname 
    if('coordinateUncertaintyinMetres' %in% colnames(d)) {
      message ("Renaming recordID column to id")
      names(d)[names(d) == 'coordinateUncertaintyinMetres'] <- 'coordinateUncertaintyInMetres'
      
    }
    
    #  standardi[sz]e catnum colname 
    if('recordID' %in% colnames(d)) {
      message ("Renaming recordID column to id")
      names(d)[names(d) == 'recordID'] <- 'id'
      
    }
  
    ## Create the searchTaxon column
    message ('Formatting ALA data for ', x)
    d[,"searchTaxon"] = x
    d[,"searchTaxon"] = gsub("_ALA_records.RData", "", d[,"searchTaxon"])
    
    if(!is.character(d["id"])) {
      d["id"] <- as.character(d["id"])
    }
    
    ## Choose only the desired columns
    d = d %>%
      select(one_of(ALA.keep))
    
    ## Then print warnings
    warnings()
    
    ## This is a list of columns in different ALA files which have weird characters
    message ('Formatting numeric ALA data for ', x)
    d[,"coordinateUncertaintyInMetres"] = as.numeric(unlist(d["coordinateUncertaintyInMetres"]))
    d["year"]  = as.numeric(unlist(d["year"]))
    d["month"] = as.numeric(unlist(d["month"]))
    d["id"]    = as.character(unlist(d["id"]))

    return(d)
    
  }) %>%
  
  ## Finally, bind all the rows together
  bind_rows

## Clear the garbage
gc()


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
  dplyr::select(dplyr::one_of(ALA.keep))

dim(ALA.TRIM)
sort(names(ALA.TRIM))


## What are the unique species?
length(unique(ALA.TRIM$searchTaxon))
length(unique(ALA.TRIM$scientificNameOriginal)) 
length(unique(ALA.TRIM$scientificName)) 
length(unique(ALA.TRIM$species))
(sum(is.na(ALA.TRIM$scientificName)) + dim(subset(ALA.TRIM, scientificName == ""))[1])/dim(ALA.TRIM)[1]*100





#########################################################################################################################
## 3). FILTER RECORDS TO THOSE WITH COORDINATES, AND AFTER 1950
#########################################################################################################################


#########################################################################################################################
## Now filter the ALA records using conditions which are not too restrictive
ALA.CLEAN <- ALA.TRIM %>% 
  
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
message(dim(ALA.TRIM)[1] - dim(ALA.CLEAN)[1], " records removed")
message(round((dim(ALA.CLEAN)[1])/dim(ALA.TRIM)[1]*100, 2), 
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
