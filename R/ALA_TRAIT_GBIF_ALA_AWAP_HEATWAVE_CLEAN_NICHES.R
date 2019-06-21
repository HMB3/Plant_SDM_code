#########################################################################################################################
################################################# FLAG SPATIAL OUTLIERS #################################################
#########################################################################################################################


## This code uses the CoordinateCleaner package to Automatcially screen out the dodgy spatial records. See ::
## https://github.com/azizka/CoordinateCleaner


## It then creates a table of the species niche :

## 2). A table with one row for each species, including contextual data and species attributes (niches, traits, etc.)

## These tables are subsequently used to estimate the current global realised niche/climatic tolerance using the best 
## available data, and susequently model the niches using the maxent algorithm.


#########################################################################################################################
## Read in the three data tables
if(read_data == "TRUE") {
  
  ## Load GBIF and ALA data
  COMBO.RASTER.CONVERT = readRDS( paste0(DATA_path, 'COMBO_RASTER_CONVERT_',  save_run, '.rds'))
  
  
} else {
  
  message(' skip file reading, not many species analysed')   ##
  
}


## Check dimensions of the occurrence and inventory data tables.
length(unique(COMBO.RASTER.CONVERT$searchTaxon))
formatC(dim(COMBO.RASTER.CONVERT)[1], format = "e", digits = 2)
names(COMBO.RASTER.CONVERT)


#########################################################################################################################
## Create a unique identifier
COMBO.RASTER.CONVERT$CC.OBS <- 1:nrow(COMBO.RASTER.CONVERT)





#########################################################################################################################
## 1). FLAG GBIF OUTLIERS
#########################################################################################################################


#########################################################################################################################
## Rename the columns to fit the CleanCoordinates format and create a tibble. 
## A Tibble is needed for running the spatial outlier cleaning


## We've already stripped out the records that fall outside
## the worldclim raster boundaries, so the sea test is probably not the most important
## Study area is the globe, but we are only projecting models onto Australia


#########################################################################################################################
## Don't run the outliers test here, it is slower. Also, can't run cleaning on the urban tree inventory data, because this
## removes all the records near capital cities


#########################################################################################################################
## Rename the columns to fit the CleanCoordinates format and create a tibble. 
TIB.GBIF <- COMBO.RASTER.CONVERT %>% dplyr::rename(coord_spp        = searchTaxon,
                                                 decimallongitude = lon, 
                                                 decimallatitude  = lat) %>%
  
  ## The create a tibble for running the spatial outlier cleaning
  timetk::tk_tbl() %>% 
  
  ## Consider the arguments. We've already stripped out the records that fall outside
  ## the worldclim raster boundaries, so the sea test is probably not the most important
  ## Study area is the globe, but we are only projecting models onto Australia
  
  ## The geographic outlier detection is not working here. Try it in setp 6
  ## Also, the duplicates step is not working, flagging too many records
  clean_coordinates(.,
                    verbose         = TRUE,
                    tests = c("capitals",     "centroids", "equal", "gbif", 
                              "institutions", "zeros"), ## duplicates flagged too many
                    
                    capitals_rad    = 10000,  ## remove records within 10km  of capitals
                    centroids_rad   = 5000    ## remove records within 5km of country centroids
                    # outliers_method = "distance", ## Save the other checks for step 6
                    # outliers_td     = 800,
                    # outliers_mtp    = 5
  ) %>%
  
  ## The select the relevant columns and rename
  select(., coord_spp, CC.OBS, .val,  .equ, .zer,  .cap,   
         .cen,    .gbf,   .inst, .summary) 

## Then rename
summary(TIB.GBIF)
names(TIB.GBIF) = c("coord_spp", "CC.OBS",    "coord_val",  "coord_equ",  "coord_zer",  "coord_cap", 
                    "coord_cen", "coord_gbf", "coord_inst", "coord_summary")


## Flagging ~ x%, excluding the spatial outliers. Seems reasonable?
message(round(with(TIB.GBIF, table(coord_summary)/sum(table(coord_summary))*100), 2), " % records removed")


#########################################################################################################################
## Check the order still matches
identical(COMBO.RASTER.CONVERT$CC.OBS, TIB.GBIF$CC.OBS)                                               ## order matches
identical(COMBO.RASTER.CONVERT$searchTaxon, TIB.GBIF$coord_spp)                                       ## order matches


## Is the species column the same as the searchTaxon column?
TEST.GEO   = join(COMBO.RASTER.CONVERT, TIB.GBIF)
identical(COMBO.RASTER.CONVERT$searchTaxon, TEST.GEO$coord_spp)                                       ## order matches 
identical(COMBO.RASTER.CONVERT$CC.OBS,      TEST.GEO$CC.OBS)                                          ## order matches


## Now subset to records that are flagged as outliers
CLEAN.TRUE  = subset(TEST.GEO, coord_summary == "TRUE")
CLEAN.FALSE = subset(TEST.GEO, coord_summary == "FALSE")
table(CLEAN.TRUE$coord_summary) 


## What percentage of records are retained?
length(unique(CLEAN.TRUE$searchTaxon))
message(round(nrow(CLEAN.TRUE)/nrow(TEST.GEO)*100, 2), " % records retained")





#########################################################################################################################
## 2). FLAG SPATIAL OUTLIERS
#########################################################################################################################


#########################################################################################################################
## Create a unique identifier for spatial cleaning. This is used for automated cleaing of the records, and also saving shapefiles
## But this will not be run for all species linearly. So, it probably needs to be a combination of species and number
CLEAN.TRUE$SPOUT.OBS <- 1:nrow(CLEAN.TRUE)
CLEAN.TRUE$SPOUT.OBS <- paste0(CLEAN.TRUE$SPOUT.OBS, "_SPOUT_", CLEAN.TRUE$searchTaxon)
CLEAN.TRUE$SPOUT.OBS <- gsub(" ",     "_",  CLEAN.TRUE$SPOUT.OBS, perl = TRUE)
length(CLEAN.TRUE$SPOUT.OBS);length(unique(CLEAN.TRUE$SPOUT.OBS))


## Check dimensions
dim(CLEAN.TRUE)
length(unique(CLEAN.TRUE$searchTaxon))
length(unique(CLEAN.TRUE$SPOUT.OBS))
unique(CLEAN.TRUE$SOURCE)

                                              



#########################################################################################################################
## 4). CREATE SHAPEFILES TO CHECK OUTLIERS ARCMAP
#########################################################################################################################


## Select columns
GBIF.ALA.CHECK  = select(CLEAN.TRUE, searchTaxon, scientificName, lat, lon, SOURCE, year,
                         country, locality, basisOfRecord, institutionCode)


## Rename the fields so that ArcMap can handle them
GBIF.ALA.CHECK     = dplyr::rename(GBIF.ALA.CHECK,
                                   TAXON     = searchTaxon,
                                   LAT       = lat,
                                   LON       = lon,
                                   SC_NAME   = scientificName,
                                   BASIS     = basisOfRecord,
                                   LOCAL     = locality,
                                   INSTIT    = institutionCode,
                                   COUNTRY   = country,
                                   YEAR      = year)
names(GBIF.ALA.CHECK)


#########################################################################################################################
## Then create SPDF
GBIF.ALA.SPDF    = SpatialPointsDataFrame(coords      = GBIF.ALA.CHECK[c("LON", "LAT")],
                                          data        = GBIF.ALA.CHECK,
                                          proj4string = CRS.WGS.84)


# ## Create list of species
# TAXA <- unique(GBIF.ALA.SPDF$TAXON)
# 
# 
# #########################################################################################################################
# ## Then, loop over the species list and create a shapefile for each 
# # for (i in 1:length(TAXA)) {
# # 
# #   ## Need to check the OBS column matches up - or do we not need this again?
# #   message("Writing shapefile for ", TAXA[i])
# #   tmp <- GBIF.ALA.SPDF[GBIF.ALA.SPDF$TAXON == TAXA[i], ]
# #   writeOGR(tmp, dsn = "./data/ANALYSIS/CLEAN_GBIF", TAXA[i], driver = "ESRI Shapefile", overwrite_layer = TRUE)
# # 
# # }
# 
# 
## Write the whole shapefile out just in case
writeOGR(obj = GBIF.ALA.SPDF, dsn = "./data/ANALYSIS/CLEAN_GBIF",
         layer = paste0("CLEAN_SPDF_", save_run), driver = "ESRI Shapefile")





#########################################################################################################################
## 5). INTERSECT SPECIES RECORDS WITH LOCAL GOV AREAS AND SIGNIFICANT URBAN AREAS
#########################################################################################################################


#########################################################################################################################
## Only run this code if we are saving niches 
if(calc_niche == "TRUE") {
  
  #########################################################################################################################
  ## Create a spatial points object
  ## We want to summarise the niches at 1km, but including all the environmental variables (EG PET, etc,),
  ## Not just those used in the SDM table (i.e. worldclim so far)
  ## CLEAN.INV = CLEAN.INV[CLEAN.INV$searchTaxon %in% GBIF.spp, ]
  NICHE.1KM    <- CLEAN.TRUE
  NICHE.1KM.84 <- SpatialPointsDataFrame(coords      = NICHE.1KM[c("lon", "lat")],
                                         data        = NICHE.1KM,
                                         proj4string = CRS.WGS.84)
  
  
  ## Use a projected, rather than geographic, coordinate system
  ## Not sure why, but this is needed for the spatial overlay step
  AUS.WGS      = spTransform(AUS,       CRS.WGS.84)


  ## Create global niche and Australian niche.
  ## So we need a subset for Australia
  ## The ,] acts like a clip in ArcMap
  NICHE.AUS <-  NICHE.1KM.84[AUS.WGS, ]

  
  
  
  
  #########################################################################################################################
  ## 4). CREATE NICHES FOR PROCESSED TAXA
  #########################################################################################################################
  
  
  #########################################################################################################################
  ## Now summarise the niches. But, figure out a cleaner way of doing this
  message('Estimating global niches for ', length(unique(CLEAN.TRUE$searchTaxon)), ' species across ', 
          length(env.variables), ' climate variables')
  
  
  #########################################################################################################################
  ## Create niche summaries for each environmental condition like this
  ## Create df for worldclim variables
  ENV.DF = NICHE.1KM %>% 
    as.data.frame() %>% 
    dplyr::select(., searchTaxon, one_of(env.variables))
  

  
  #########################################################################################################################
  ## Remove the NA worldclim data :: 99% retained
  ENV.DAT = completeFun(ENV.DF, "PET")
  head(niche_estimate (DF = ENV.DAT, colname = "Annual_mean_temp"))[1:5]  ## Including the q05 and q95
 
   
  message(round(nrow(ENV.DAT)/nrow(ENV.DF)*100, 2), " % records retained")
  message(round(length(unique(ENV.DAT$searchTaxon))/
                  length(unique(ENV.DF$searchTaxon))*100, 2), "% species retained for worldclim")
  
  
  #########################################################################################################################
  message('Estimating environmental niches for ', length(unique(ENV.DAT$searchTaxon)), ' species in the set ', "'", save_run, "'")
  ENV.NICHE <- env.variables %>% 
    
    ## Pipe the list into lapply
    lapply(function(x) {
      
      ## Now, use the niche width function on each colname
      ## Also, need to figure out how to make the aggregating column generic (species, genus, etc.)
      ## currently it only works hard-wired..............................................................................
      niche_estimate(DF = ENV.DAT, colname = x)
      
      ## Would be good to remove the duplicate columns here..............................................................
      
    }) %>% 
    
    ## finally, create one dataframe for all niches
    as.data.frame
  
  ## Remove duplicate Taxon columns and check the output :: would be great to skip these columns when running the function
  ENV.NICHE <- subset(ENV.NICHE, select = -c(searchTaxon.1,  searchTaxon.2,  searchTaxon.3,  searchTaxon.4,
                                             searchTaxon.5,  searchTaxon.6,  searchTaxon.7,  searchTaxon.7,
                                             searchTaxon.8,  searchTaxon.9,  searchTaxon.10, searchTaxon.11,
                                             searchTaxon.12, searchTaxon.13, searchTaxon.14, searchTaxon.15,
                                             searchTaxon.16, searchTaxon.17, searchTaxon.18, searchTaxon.19))
  

  #########################################################################################################################
  ## Add global counts for each species, and record the total number of taxa processed
  ## This is actually redundant - the records will the be the same as the maxent_records
  Global_records = as.data.frame(table(NICHE.1KM$searchTaxon))
  names(Global_records) = c("searchTaxon", "Global_records")
  identical(nrow(Global_records), nrow(ENV.NICHE))

  
  #########################################################################################################################
  ## Join all the tables together - this allows all the data to be used for each vaiable group
  ALL.NICHE = join_all(list(Global_records, ENV.NICHE),  by = 'searchTaxon', type = 'full')
  

  
  
  
  #########################################################################################################################
  ## 5). CALCULATE AREA OF OCCUPANCY
  #########################################################################################################################
  
  
  ## Given a dataframe of georeferenced occurrences of one, or more, taxa, this function provide statistics values 
  ## (Extent of Occurrence, Area of Occupancy, number of locations, number of subpopulations) and provide a preliminary 
  ## conservation status following Criterion B of IUCN.
  
  
  #########################################################################################################################
  ## AREA OF OCCUPANCY (AOO).
  ## For every species in the list: calculate the AOO
  ## x = spp.geo[1]
  spp.geo = as.character(unique(NICHE.1KM$searchTaxon))
  
  GBIF.AOO <- spp.geo %>%
    
    ## Pipe the list into lapply
    lapply(function(x) {
      
      ## Subset the the data frame to calculate area of occupancy according the IUCN.eval
      DF   = subset(NICHE.1KM, searchTaxon == x)[, c("lat", "lon", "searchTaxon")]
      
      message('Calcualting geographic ranges for ', x, ', ', nrow(DF), ' records')
      AOO  = AOO.computing(XY = DF, Cell_size_AOO = 2)  ## Grid size in decimal degrees. Changes the results
      AOO  = as.data.frame(AOO)
      AOO$searchTaxon  <- rownames(AOO)
      rownames(AOO)    <- NULL
      
      ## Remeber to explicitly return the df at the end of loop, so we can bind
      return(AOO)
      
    }) %>%
    
    ## Finally, create one dataframe for all niches
    bind_rows
  
  head(GBIF.AOO)
  
  
  #########################################################################################################################
  ## Now join on the geographic range and glasshouse data
  identical(nrow(GBIF.AOO), length(unique(NICHE.1KM$searchTaxon)))
  ALL.NICHE = join(GBIF.AOO, ALL.NICHE, type = 'full')
  ALL.NICHE = ALL.NICHE[c(2, 1, 3:length(names(ALL.NICHE)))]
  
  
  #########################################################################################################################
  ## save .rds file for the next session
  message('Writing 1km resolution niche data for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")
  saveRDS(ALL.NICHE,              paste0(DATA_path, 'AWAP_NICHE_DATA_',  save_run, '.rds'))
  saveRDS(CLEAN.TRUE,             paste0(DATA_path, 'AWAP_POINT_DATA_',  save_run, '.rds'))
  
} else {
  
  message(' skip file saving, ', length(GBIF.spp), ' species analysed')   ##
  
}




#########################################################################################################################
## OUTSTANDING CLEANING TASKS:
#########################################################################################################################


## Check the automatic cleaning results against manual cleaning results




#########################################################################################################################
#####################################################  TBC ##############################################################
#########################################################################################################################