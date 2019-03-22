#########################################################################################################################
################################################# FLAG SPATIAL OUTLIERS #################################################
#########################################################################################################################


#########################################################################################################################
## This code uses the CoordinateCleaner package to Automatcially screen out the dodgy spatial records. See ::
## https://github.com/azizka/CoordinateCleaner


## It then creates a table of the species niche :

## 1). A table with one row for each species, including contextual data and species attributes (niches, traits, etc.)

## These tables are subsequently used to estimate the current global realised niche/climatic tolerance using the best 
## available data, and susequently model the niches using the maxent algorithm.


## Update line 548 RE niches.............................................................................................


#########################################################################################################################
## Michelle's advice :: use all the global data.
## Does including the global data improve the bad species?
## Re-import Alessandro's inventory data and re-run the models...........................................................


## Print the species run to the screen
message('Cleaning outliers and creating niches for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")


#########################################################################################################################
## Read in the global data
if(read_data == "TRUE") {
  
  ## Read in RDS files from previous step
  COMBO.RASTER.CONVERT = readRDS(paste0(DATA_path, 'COMBO_RASTER_CONVERT_',  save_run, '.rds'))
  
} else {
  
  message(' skip file reading, not many species analysed')   ##
  
}


#########################################################################################################################
## Read in the global data
if(dim(TI.XY.SPP)[1] > 0 & read_data == "TRUE") {
  
  ## Read in RDS files from previous step
  TI.RASTER.CONVERT = readRDS(paste0(DATA_path, 'TI_RASTER_CONVERT_', save_run, '.rds'))
  
} else {
  
  message(' skip file reading, no urban records')   ##
  
}


## Check dimensions of the occurrence and inventory data tables.
length(unique(COMBO.RASTER.CONVERT$searchTaxon))
formatC(dim(COMBO.RASTER.CONVERT)[1], format = "e", digits = 2)
names(COMBO.RASTER.CONVERT)


#########################################################################################################################
## Create a unique identifier. This is used for automated cleaing the records, and also saving shapefiles
COMBO.RASTER.CONVERT$CC.OBS <- 1:nrow(COMBO.RASTER.CONVERT)





#########################################################################################################################
## 1). FLAG GBIF OUTLIERS
#########################################################################################################################


#########################################################################################################################
## Rename the columns to fit the CleanCoordinates format and create a tibble. 
TIB.GBIF <- COMBO.RASTER.CONVERT %>% dplyr::rename(species          = searchTaxon,
                                                   decimallongitude = lon, 
                                                   decimallatitude  = lat) %>%
  
  ## The create a tibble for running the spatial outlier cleaning
  timetk::tk_tbl() %>% 
  
  ## Consider the arguments. We've already stripped out the records that fall outside
  ## the worldclim raster boundaries, so the sea test is probably not the most important
  ## Study area is the globe, but we are only projecting models onto Australia
  
  ## The geographic outlier detection is not working here. Try it in setp 6
  clean_coordinates(.,
                    verbose         = TRUE,
                    tests = c("capitals",     "centroids", "equal", "gbif", 
                              "institutions", "zeros"),
                    
                    capitals_rad    = 10000,  ## remove records within 5km  of capitals
                    centroids_rad   = 5000   ## remove records within 500m of country centroids
                    # outliers_method = "distance", ## remove records > 100km from other records
                    # outliers_td     = 800,
                    # outliers_mtp    = 5
                    ) %>%
  
  ## The select the relevant columns and rename
  select(., species, CC.OBS, .val,  .equ, .zer,  .cap,   
            .cen,    .gbf,   .inst, .summary) 

## Then rename
summary(TIB.GBIF)
names(TIB.GBIF) = c("coord_spp", "CC.OBS",    "coord_val",  "coord_equ",  "coord_zer",  "coord_cap", 
                    "coord_cen", "coord_gbf", "coord_inst", "coord_summary")


## Flagging ~ x%, excluding the spatial outliers. Seems reasonable?
message(round(with(TIB.GBIF, table(coord_summary)/sum(table(coord_summary))*100), 1), " % records removed")


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
unique(CLEAN.TRUE$coord_summary)   





#########################################################################################################################
## 2). PLOT GBIF OUTLIERS TO CHECK
#########################################################################################################################


#########################################################################################################################
## Try plotting the points which are outliers for a subset of species and label them
CLEAN.FALSE.PLOT = SpatialPointsDataFrame(coords      = CLEAN.FALSE[c("lon", "lat")],
                                          data        = CLEAN.FALSE,
                                          proj4string = CRS.WGS.84)


CLEAN.TRUE.PLOT = SpatialPointsDataFrame(coords      = CLEAN.TRUE[c("lon", "lat")],
                                         data        = CLEAN.TRUE,
                                         proj4string = CRS.WGS.84)


## Add land
LAND.84 = LAND %>%
  spTransform(CRS.WGS.84)
AUS.84 = AUS %>%
  spTransform(CRS.WGS.84)


#########################################################################################################################
## Get the first 10 species
## species = plot.taxa[9]
plot.taxa <- as.character(unique(CLEAN.TRUE$searchTaxon))
for (species in plot.taxa) {

  ## Plot a subset of taxa
  CLEAN.TRUE.P.I  = subset(CLEAN.TRUE.PLOT,  searchTaxon == species)
  CLEAN.FALSE.P.I = subset(CLEAN.FALSE.PLOT, searchTaxon == species)
  
  message("plotting occ data for ", species, ", ", 
          nrow(CLEAN.TRUE.P.I), " clean records")
  
  #############################################################
  ## Plot true and false points for the world
  message('Writing map of global coord clean records for ', species)
  png(sprintf("./data/ANALYSIS/CLEAN_GBIF/%s_%s", species, "global_coord_check.png"),
      16, 10, units = 'in', res = 500)
  
  par(mfrow = c(1,2))
  plot(LAND.84, main = paste0("Global points for ", species),
       lwd = 0.01, asp = 1, col = 'grey', bg = 'sky blue')
  
  points(CLEAN.TRUE.P.I,
         pch = ".", cex = 3.3, cex.lab = 3, cex.main = 4, cex.axis = 2, 
         xlab = "", ylab = "", asp = 1,
         col = "blue")
  
  points(CLEAN.FALSE.P.I,
         pch = ".", cex = 3.3, cex.lab = 3, cex.main = 4, cex.axis = 2, 
         xlab = "", ylab = "", asp = 1,
         col = "red")
  
  #############################################################
  ## Plot true and false points for the world
  plot(AUS.84, main = paste0("Australian points for ", species),
       lwd = 0.01, asp = 1, bg = 'sky blue', col = 'grey')
  
  points(CLEAN.TRUE.P.I,
         pch = ".", cex = 3.3, cex.lab = 3, cex.main = 4, cex.axis = 2, 
         xlab = "", ylab = "", asp = 1,
         col = "blue")
  
  points(CLEAN.FALSE.P.I,
         pch = ".", cex = 3.3, cex.lab = 3, cex.main = 4, cex.axis = 2, 
         xlab = "", ylab = "", asp = 1,
         col = "red")
  
  dev.off()
  

}


## What percentage of records are retained?
length(unique(CLEAN.TRUE$searchTaxon))
message(round(nrow(CLEAN.TRUE)/nrow(TEST.GEO)*100, 2), " % records retained")                                               


#########################################################################################################################
## Now bind on the urban tree inventory data. We are assuming this data is clean, after we manually fix the taxonomy
## If you want to clean the inventory records, put a OBS column in here too.
if(dim(TI.XY.SPP)[1] > 0) {
  
  message('Combining Australian inventory data with occurrence data') 
  intersect(names(TI.RASTER.CONVERT), names(CLEAN.TRUE))
  CLEAN.TRUE = bind_rows(CLEAN.TRUE, TI.RASTER.CONVERT)
  
} else {
  
  ## Update with global data
  message('No Australian inventory data for these species')   ##
  CLEAN.TRUE = CLEAN.TRUE
  
}


## How many local municipalities have data for the species analysed
names(CLEAN.TRUE)
unique(CLEAN.TRUE$SOURCE) 
unique(CLEAN.TRUE$INVENTORY)
length(unique(CLEAN.TRUE$INVENTORY))


length(unique(CLEAN.TRUE$searchTaxon))
summary(CLEAN.TRUE$Annual_mean_temp)
summary(CLEAN.TRUE$Annual_mean_temp)


## By how many % does including tree inventories increase the overal number of records?
message("Tree inventory data increases records by ", round(dim(CLEAN.TRUE)[1]/dim(TEST.GEO)[1]*100, 2), " % ")   


#########################################################################################################################
## save data
if(save_data == "TRUE") {
  
  ## save .rds file for the next session
  saveRDS(CLEAN.TRUE, paste0(DATA_path, 'CLEAN_TRUE_', save_run, '.rds'))
  
} else {
  
  message(' skip file saving, not many species analysed')   ##
  
}





#########################################################################################################################
## 2). CREATE SHAPEFILES - IF (!) YOU WANT TO CHECK OUTLIERS ARCMAP
#########################################################################################################################


#########################################################################################################################
## This code could be used to manually clean outlying records, using the OBS column......................................
## Note that for records with catalogue numbers, you can actually search GBIF and ALA for those occurrences. This seems
## like overkill for so many species. But you could make a list of the OBS column, then exclude these from the larger
## dataframe. Seems like you would only do this on the final dataset, filtered to 1km resolution.


## Select columns to send to the shapefile
GBIF.ALA.CHECK  = dplyr::select(CLEAN.TRUE,     CC.OBS, searchTaxon, scientificName, lat, lon, SOURCE, INVENTORY, year, 
                                country, locality, basisOfRecord, institutionCode, 
                                Taxonomic.status, New.Taxonomic.status)


## Rename the fields so that ArcMap can handle them
GBIF.ALA.CHECK  = dplyr::rename(GBIF.ALA.CHECK, 
                                OBS       = CC.OBS,
                                TAXON     = searchTaxon,
                                LAT       = lat,
                                LON       = lon,
                                SC_NAME   = scientificName,
                                BASIS     = basisOfRecord,                
                                LOCAL     = locality,                      
                                INSTIT    = institutionCode,                
                                COUNTRY   = country,                
                                YEAR      = year,
                                TAX_STAT  = Taxonomic.status,
                                NEW_STA   = New.Taxonomic.status)
names(GBIF.ALA.CHECK)


## Then create a SPDF
GBIF.ALA.SPDF    = SpatialPointsDataFrame(coords      = GBIF.ALA.CHECK[c("LON", "LAT")],
                                          data        = GBIF.ALA.CHECK,
                                          proj4string = CRS.WGS.84)


## Create list of species 
TAXA <- unique(GBIF.ALA.SPDF$TAXON)


#########################################################################################################################
## Then, loop over the species list and create a shapefile for each 
## For lots of taxa, this would be needed to 
for (i in 1:length(TAXA)) {

  ## Need to check the OBS column matches up - or do we not need this again?
  message("Writing shapefile for ", TAXA[i])
  tmp <- GBIF.ALA.SPDF[GBIF.ALA.SPDF$TAXON == TAXA[i], ]
  writeOGR(tmp, dsn = "./data/ANALYSIS/CLEAN_GBIF", TAXA[i], driver = "ESRI Shapefile", overwrite_layer = TRUE)

}


## Also write out the whole shapefile, just in case
# writeOGR(obj = GBIF.ALA.SPDF, dsn = "./data/ANALYSIS/CLEAN_GBIF", 
#          layer = paste0("CLEAN_SPDF_", save_run), driver = "ESRI Shapefile")





#########################################################################################################################
## 3). INTERSECT SPECIES RECORDS WITH LOCAL GOV AREAS AND SIGNIFICANT URBAN AREAS
#########################################################################################################################


#########################################################################################################################
## Only run this code if we are saving niches 
if(calc_niche == "TRUE") {
  
  ## We want to know the count of species that occur in 'n' LGAs, across a range of climates. Read in LGA and SUA
  length(unique(CLEAN.TRUE$searchTaxon))
  projection(LGA);projection(AUS);projection(SUA_2016)
  
  
  ## Convert the raster data back into a spdf
  COMBO.RASTER.SP   = SpatialPointsDataFrame(coords      = CLEAN.TRUE[c("lon", "lat")], 
                                             data        = CLEAN.TRUE,
                                             proj4string = CRS.WGS.84)
  
  
  ## Project using a projected, rather than geographic, coordinate system
  ## Not sure why, but this is needed for the spatial overlay step
  LGA.WGS   = spTransform(LGA,       CRS.WGS.84)
  SUA.WGS   = spTransform(SUA_2016 , CRS.WGS.84)
  AUS.WGS   = spTransform(AUS,       CRS.WGS.84)
  LAND.WGS  = spTransform(LAND,      CRS.WGS.84)
  
  
  ## Remove the columns we don't need
  LGA.WGS = LGA.WGS[, c("LGA_CODE16", "LGA_NAME16")] 
  head(SUA.WGS)
  
  
  #########################################################################################################################
  ## Run join between species records and LGAs/SUAs :: Double check they are the same
  ## See the ABS for details :: there are 563 LGAs
  ## http://www.abs.gov.au/ausstats/abs@.nsf/Lookup/by%20Subject/1270.0.55.003~July%202016~Main%20Features~Local%20Government%20Areas%20(LGA)~7
  message('Joining occurence data to SUAs for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")
  
  
  projection(COMBO.RASTER.SP);projection(LGA.WGS);projection(SUA.WGS);projection(AUS.WGS)
  SUA.JOIN      = over(COMBO.RASTER.SP,   SUA.WGS)
  COMBO.SUA.LGA = cbind.data.frame(COMBO.RASTER.SP, SUA.JOIN) 
  
  
  #########################################################################################################################
  ## AGGREGATE THE NUMBER OF SUAs EACH SPECIES IS FOUND IN. NA LGAs ARE OUTSIDE AUS
  SUA.AGG   = tapply(COMBO.SUA.LGA$SUA_NAME16, COMBO.SUA.LGA$searchTaxon, function(x) length(unique(x))) ## group SUA by species name
  SUA.AGG   = as.data.frame(SUA.AGG)
  AUS.AGG   = aggregate(SUA_NAME16 ~ searchTaxon, data = COMBO.SUA.LGA,   function(x) {sum(!is.na(x))}, na.action = NULL)
  SUA.AGG   = cbind.data.frame(AUS.AGG, SUA.AGG)
  names(SUA.AGG) = c("searchTaxon", "AUS_RECORDS", "SUA_COUNT")
  
  
  ## Now create a table of all the SUA's that each species occurrs
  SUA.SPP.COUNT = as.data.frame(table(COMBO.SUA.LGA[["SUA_NAME16"]], COMBO.SUA.LGA[["searchTaxon"]]))
  names(SUA.SPP.COUNT) = c("SUA", "SPECIES", "SUA_RECORDS")
  
  
  ## Save .rds file for the next session
  saveRDS(SUA.SPP.COUNT, paste0(DATA_path, 'SUA_SPP_COUNT_', save_run, '.rds'))
  
  
  ## Check : That's ok, but we want a table of which SUA each species is actually in.
  dim(SUA.AGG)
  head(SUA.AGG)
  
  
  ## Remove duplicate coordinates
  drops <- c("lon.1", "lat.1")
  COMBO.SUA.LGA <- COMBO.SUA.LGA[ , !(names(COMBO.SUA.LGA) %in% drops)]
  names(COMBO.SUA.LGA)
  dim(COMBO.SUA.LGA)
  str(unique(COMBO.SUA.LGA$searchTaxon))
  unique(COMBO.SUA.LGA$SOURCE)
  
  
  
  
  
  #########################################################################################################################
  ## 4). CREATE NICHES FOR SELECTED TAXA
  #########################################################################################################################
  
  
  #########################################################################################################################
  ## Now summarise the niches. But, figure out a cleaner way of doing this
  message('Estimating global niches for ', length(GBIF.spp), ' species across ', length(env.variables), ' climate variables')
  
  
  ## If the occurrence source is ALA, just get that. Only needed for Diana's data
  if(OCC_SOURCE == "ALA") {
    
    ## save .rds file for the next session
    message('Create niches from ', OCC_SOURCE, ' records')
    COMBO.SUA.LGA = COMBO.SUA.LGA[COMBO.SUA.LGA$SOURCE %in% OCC_SOURCE , ]
    message(unique(COMBO.SUA.LGA$SOURCE))
    
  } else {
    
    message('Create niches from all sources')
    message(unique(COMBO.SUA.LGA$SOURCE))
    
  }
  
  
  #########################################################################################################################
  ## Create niche summaries for each environmental condition like this...
  ## Here's what the function will produce :
  NICHE.DF = completeFun(COMBO.SUA.LGA, "PET")
  dim(NICHE.DF)
  head(niche_estimate (DF = NICHE.DF, colname = "Annual_mean_temp"))  ## including the q05 and q95
  
  
  ## So lets use lapply on the "SearchTaxon"
  ## test = run_function_concatenate(list, DF, "DF, colname = x") 
  message('Estimating global niches for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")
  COMBO.NICHE <- env.variables %>% 
    
    ## Pipe the list into lapply
    lapply(function(x) {
      
      ## Now, use the niche width function on each colname
      ## Also, need to figure out how to make the aggregating column generic (species, genus, etc.)
      ## currently it only works hard-wired..............................................................................
      niche_estimate (DF = NICHE.DF, colname = x)
      
      ## would be good to remove the duplicate columns here..............................................................
      
    }) %>% 
    
    ## finally, create one dataframe for all niches
    as.data.frame
  
  
  ## Remove duplicate Taxon columns and check the output :: would be great to skip these columns when running the function
  names(COMBO.NICHE)
  COMBO.NICHE <- subset(COMBO.NICHE, select = -c(searchTaxon.1,  searchTaxon.2,  searchTaxon.3,  searchTaxon.4,
                                                 searchTaxon.5,  searchTaxon.6,  searchTaxon.7,  searchTaxon.7,
                                                 searchTaxon.8,  searchTaxon.9,  searchTaxon.10, searchTaxon.11,
                                                 searchTaxon.12, searchTaxon.13, searchTaxon.14, searchTaxon.15,
                                                 searchTaxon.16, searchTaxon.17, searchTaxon.18, searchTaxon.19))
  
  
  #########################################################################################################################
  ## Add counts for each species, and record the total number of taxa processed
  GLOBAL_RECORDS = as.data.frame(table(COMBO.SUA.LGA$searchTaxon))
  names(GLOBAL_RECORDS) = c("searchTaxon", "GLOBAL_RECORDS")
  identical(nrow(GLOBAL_RECORDS), nrow(COMBO.NICHE))
  
  
  #########################################################################################################################
  ## Add the counts of Australian records for each species to the niche database
  Total.taxa.processed = nrow(COMBO.NICHE)
  COMBO.NICHE  = join(GLOBAL_RECORDS, COMBO.NICHE, type = "right")
  COMBO.NICHE  = join(SUA.AGG, COMBO.NICHE,        type = "right")
  
  
  head(COMBO.NICHE$AUS_RECORDS)
  head(COMBO.NICHE$SUA_COUNT)
  
  
  names(COMBO.NICHE)
  dim(COMBO.NICHE)
  
  

  
  
  #########################################################################################################################
  ## 5). CALCULATE AREA OF OCCUPANCY
  #########################################################################################################################
  
  
  ## Given a dataframe of georeferenced occurrences of one, or more, taxa, this function provide statistics values 
  ## (Extent of Occurrence, Area of Occupancy, number of locations, number of subpopulations) and provide a preliminary 
  ## conservation status following Criterion B of IUCN. A graphical map output is also available.
  
  
  ## Create a species list to estimate the ranges for. Currently we can't estimate ranges for widely distributed species
  ## A workaround for the trait species is to just use the Australian records. But for the full list of species, we will need to 
  ## change this. Otherwise, we could report the global niche, but the Australian range only. 
  ## Change when Gilles Dauby updates......................................................................................
  COMBO.SUA.SP = SpatialPointsDataFrame(coords      = COMBO.SUA.LGA[c("lon", "lat")], 
                                        data        = COMBO.SUA.LGA,
                                        proj4string = CRS.WGS.84)
  
  
  ## Subset the occurrence records to Australia only.
  projection(COMBO.SUA.SP);projection(AUS.WGS)
  AUS.AOO.DAT = COMBO.SUA.SP[AUS.WGS, ]
  
  
  ## This separates the points out 
  plot(LAND)
  points(AUS.AOO.DAT,  col = "red",  pch = ".")
  points(COMBO.SUA.SP, col = "blue", pch = ".")
  
  
  ## The IUCN.eval function needs a data frame, not a spdf
  AUS.AOO.DAT = as.data.frame(AUS.AOO.DAT)
  
  
  #########################################################################################################################
  ## AREA OF OCCUPANCY (AOO)
  ## For every species in the list: calculate the AOO
  ## Change this once Giles Dauby updates..................................................................................
  ## x = spp.geo[1]
  spp.geo = as.character(unique(COMBO.SUA.SP$searchTaxon))
  
  GBIF.AOO <- spp.geo %>%
    
    ## Pipe the list into lapply
    lapply(function(x) {
      
      ## Subset the the data frame to calculate area of occupancy according the IUCN.eval
      ## Is the function very slow for large datasets? YES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DF = subset(AUS.AOO.DAT, searchTaxon == x)[, c("lat", "lon", "searchTaxon")]
      IUCN.eval(DF, Cell_size_AOO = 10, DrawMap = FALSE, country_map = land, SubPop = FALSE) 
      
    }) %>%
    
    ## Finally, create one dataframe for all niches
    bind_rows
  
  names(GBIF.AOO)[names(GBIF.AOO) == 'taxa'] <- 'searchTaxon'
  GBIF.AOO             = dplyr::select(GBIF.AOO, searchTaxon, EOO, AOO, Category_CriteriaB)
  names(GBIF.AOO)      = c("searchTaxon", "AUS_EOO", "AUS_AOO", "ICUN_catb")
  head(GBIF.AOO)
  
  
  #########################################################################################################################
  ## Now join on the geographic range and glasshouse data
  identical(length(GBIF.AOO$AUS_AOO), length(GBIF.spp))
  COMBO.NICHE = join(GBIF.AOO, COMBO.NICHE, type = "right")
  
  
  
  
  
  #########################################################################################################################
  ## 6). JOIN ON CONTEXTUAL DATA
  #########################################################################################################################
  
  
  #########################################################################################################################
  ## Now join the horticultural contextual data onto one or both tables ()
  message('Joining contextual data for raster and niche files', length(GBIF.spp), ' species in the set ', "'", save_run, "'")
  COMBO.RASTER.CONTEXT = CLEAN.TRUE
  names(COMBO.RASTER.CONTEXT)
  
  
  #########################################################################################################################
  ## Now join the hort context data and the tree inventory plantings to the niche
  COMBO.NICHE.CONTEXT = join(CLEAN.GROW, COMBO.NICHE,         type = "right") 
  COMBO.NICHE.CONTEXT = join(TI.LUT,     COMBO.NICHE.CONTEXT, type = "right")
  
  
  ## Check the columns are still there
  head(COMBO.NICHE.CONTEXT$AUS_RECORDS)
  head(COMBO.NICHE.CONTEXT$GLOBAL_RECORDS)
  head(COMBO.NICHE.CONTEXT$Plantings)
  head(COMBO.NICHE.CONTEXT$SUA_COUNT)
  
  
  #########################################################################################################################
  ## Is there a relationship between the number of records in each category
  plot(COMBO.NICHE.CONTEXT$AUS_RECORDS, COMBO.NICHE.CONTEXT$GLOBAL_RECORDS)
  plot(COMBO.NICHE.CONTEXT$AUS_RECORDS, COMBO.NICHE.CONTEXT$Plantings)
  
  
  lm.glob  = lm(COMBO.NICHE.CONTEXT$AUS_RECORDS ~ COMBO.NICHE.CONTEXT$GLOBAL_RECORDS)
  lm.plant = lm(COMBO.NICHE.CONTEXT$AUS_RECORDS ~ COMBO.NICHE.CONTEXT$Plantings)
  
  
  layout(matrix(c(1,1,2,3), 2, 1, byrow = TRUE))
  
  plot(COMBO.NICHE.CONTEXT$GLOBAL_RECORDS, COMBO.NICHE.CONTEXT$AUS_RECORDS, 
       pch = 19, col  = "blue",
       xlab = "Global records", ylab = "Australian records", 
       abline(lm.glob), 
       main = save_run, cex = 2)
  
  legend("topleft", bty = "n", 
         legend = paste("R2 is", format(summary(lm.glob)$adj.r.squared, digits = 4)))
  
  plot(COMBO.NICHE.CONTEXT$AUS_RECORDS, COMBO.NICHE.CONTEXT$Plantings, 
       pch = 19, col  = "blue",
       xlab = "Aus records", ylab = "Plantings", 
       abline(lm.plant), 
       main = save_run, cex = 2)
  
  legend("topleft", bty = "n", 
         legend = paste("R2 is", format(summary(lm.plant)$adj.r.squared, digits = 4)))
  
  #########################################################################################################################
  ## Now combine the SDM output with the niche context data 
  ## CLEAN.GROW needs to be put through GBIF and TPL........................................................................
  COMBO.NICHE.CONTEXT = COMBO.NICHE.CONTEXT [order(COMBO.NICHE.CONTEXT $searchTaxon),]
  
  
  ## Set NA to blank, then sort by no. of growers
  COMBO.NICHE.CONTEXT$Total_growers[is.na(COMBO.NICHE.CONTEXT$Total.growers)] <- 0
  COMBO.NICHE.CONTEXT = COMBO.NICHE.CONTEXT[with(COMBO.NICHE.CONTEXT, rev(order(Total_growers))), ]
  
  
  ## View the data
  dim(COMBO.RASTER.CONTEXT)
  dim(COMBO.NICHE.CONTEXT)

  
  ## Print the dataframe dimensions to screen
  dim(CLEAN.TRUE)
  dim(COMBO.NICHE.CONTEXT)
  dim(COMBO.RASTER.CONTEXT)
  
  length(unique(CLEAN.TRUE$searchTaxon))
  length(COMBO.NICHE.CONTEXT$searchTaxon)
  length(unique(COMBO.RASTER.CONTEXT$searchTaxon))
  unique(COMBO.RASTER.CONTEXT$SOURCE)
  
  
  ## Check the column order for the niche and the record tables
  ## Note that for records with catalogue numbers, you can actually search GBIF and ALA for those occurrences
  head(COMBO.NICHE.CONTEXT)[1:12]
  head(COMBO.RASTER.CONTEXT)[1:16]
  
  
  #########################################################################################################################
  ## save .rds file for the next session
  message('Writing niche and raster data for ', length(GBIF.spp), ' species in the set ', "'", save_run, "'")
  saveRDS(COMBO.NICHE.CONTEXT,    paste0(DATA_path, 'COMBO_NICHE_CONTEXT_',  OCC_SOURCE, '_RECORDS_', save_run, '.rds'))
  saveRDS(COMBO.RASTER.CONTEXT,   paste0(DATA_path, 'COMBO_RASTER_CONTEXT_', OCC_SOURCE, '_RECORDS_', save_run, '.rds'))
  write.csv(COMBO.NICHE.CONTEXT,  paste0(DATA_path, 'COMBO_NICHE_CONTEXT_',  OCC_SOURCE, '_RECORDS_', save_run, '.csv'), row.names = FALSE)
  
  
} else {
  
  message(' skip file saving, ', length(GBIF.spp), ' species analysed')   ##
  
}


#########################################################################################################################
## OUTSTANDING CLEANING TASKS:
#########################################################################################################################


## 1). Check the settings of coordinate cleaner..........................................................................


## 2). Check the automatic cleaning results against manual cleaning results..............................................


## 3). Add a global AOO  ................................................................................................



#########################################################################################################################
#####################################################  TBC ##############################################################
#########################################################################################################################