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
  COMBO.AWAP.CONVERT = readRDS( paste0(DATA_path, 'COMBO_RASTER_CONVERT_',  save_run, '.rds'))
  
  
} else {
  
  message(' skip file reading, not many species analysed')   ##
  
}


## Check dimensions of the occurrence and inventory data tables.
length(unique(COMBO.AWAP.CONVERT$searchTaxon))
formatC(dim(COMBO.AWAP.CONVERT)[1], format = "e", digits = 2)
names(COMBO.AWAP.CONVERT)


#########################################################################################################################
## Create a unique identifier
COMBO.AWAP.CONVERT$CC.OBS <- 1:nrow(COMBO.AWAP.CONVERT)





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
TIB.GBIF <- COMBO.AWAP.CONVERT %>% dplyr::rename(coord_spp        = searchTaxon,
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
identical(COMBO.AWAP.CONVERT$CC.OBS, TIB.GBIF$CC.OBS)                                               ## order matches
identical(COMBO.AWAP.CONVERT$searchTaxon, TIB.GBIF$coord_spp)                                       ## order matches


## Is the species column the same as the searchTaxon column?
TEST.GEO   = join(COMBO.AWAP.CONVERT, TIB.GBIF)
identical(COMBO.AWAP.CONVERT$searchTaxon, TEST.GEO$coord_spp)                                       ## order matches 
identical(COMBO.AWAP.CONVERT$CC.OBS,      TEST.GEO$CC.OBS)                                          ## order matches


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


# #########################################################################################################################
# ## Create a tibble to supply to coordinate cleaner
# SDM.COORDS  <- CLEAN.TRUE %>% 
# 
#   select(searchTaxon, lon, lat, SPOUT.OBS, SOURCE) %>%
#   
#   dplyr::rename(species          = searchTaxon,
#                 decimallongitude = lon, 
#                 decimallatitude  = lat) %>%
#   
#   timetk::tk_tbl()
# 
# 
# #########################################################################################################################
# ## Check how many records each species has
# COMBO.LUT <- SDM.COORDS %>% 
#   as.data.frame() %>%
#   select(species) %>%
#   table() %>%
#   as.data.frame() 
# COMBO.LUT <- setDT(COMBO.LUT, keep.rownames = FALSE)[]
# names(COMBO.LUT) = c("species", "FREQUENCY")
# COMBO.LUT = COMBO.LUT[with(COMBO.LUT, rev(order(FREQUENCY))), ]
# View(COMBO.LUT)
# 
# 
# ## Watch out here - this sorting could cause problems for the order of the data frame once it's stitched back together
# ## If we we use species to join the data back together, will it preserve the order? 
# LUT.100K = as.character(subset(COMBO.LUT, FREQUENCY < 100000)$species)
# LUT.100K = trimws(LUT.100K [order(LUT.100K)])
# length(LUT.100K)
# 
# 
# ## See communications with Alex Zizka
# ## Check the output for patterns - make the settings strict here, as outliers could be bogus after 1km thinning
# ## Generally, species with heaps of records, especially those with clumped/biased records, get more outliers
# 
# 
# ## Create a data frame of species name and spatial outlier
# SPAT.OUT <- LUT.100K  %>%
#   
#   ## pipe the list of species into lapply
#   lapply(function(x) {
#     
#     ## Create the species df by subsetting by species
#     f <- subset(SDM.COORDS, species == x)
#     
#     ## Run the spatial outlier detection
#     message("Running spatial outlier detection for ", nrow(f), " records for ", x)
#     sp.flag <- cc_outl(f,
#                        lon     = "decimallongitude",
#                        lat     = "decimallatitude",
#                        species = "species",
#                        method  = "quantile", #"distance",
#                        mltpl   = 10,
#                        #tdi     = 300,
#                        value   = "flagged",
#                        verbose = "TRUE")
#     
#     ## Now add attache column for species, and the flag for each record
#     d = cbind(searchTaxon = x,
#               SPAT_OUT = sp.flag, f)[c("searchTaxon", "SPAT_OUT", "SPOUT.OBS")]
#     
#     ## Remeber to explicitly return the df at the end of loop, so we can bind
#     return(d)
#     
#   }) %>%
#   
#   ## Finally, bind all the rows together
#   bind_rows
# 
# gc()
# 
# 
# ## These settings produce too few outliers. Try changing the settings..................................................
# print(table(SPAT.OUT$SPAT_OUT, exclude = NULL))
# length(unique(SPAT.OUT$searchTaxon))
# head(SPAT.OUT)
# 
# 
# 
# 
# 
# #########################################################################################################################
# ## 2). FILTER DATA TO THE CLEAN RECORDS
# #########################################################################################################################
# 
# 
# #########################################################################################################################
# ## Join data :: Best to use the 'OBS' column here
# identical(nrow(SDM.COORDS), nrow(SPAT.OUT))
# identical(SDM.DATA.ALL$searchTaxon, SPAT.OUT$searchTaxon)
# length(unique(SPAT.OUT$searchTaxon))
# 
# 
# ## This explicit join is required. Check the species have been analysed in exactly the same order
# SPAT.FLAG = join(as.data.frame(SDM.DATA.ALL), SPAT.OUT, by = c("SPOUT.OBS", "searchTaxon") , type = "left", match = "first")    
# identical(SDM.DATA.ALL$searchTaxon, SPAT.FLAG$searchTaxon)
# 
# 
# ## Check the join is working 
# message('Checking spatial flags for ', length(unique(SPAT.FLAG$searchTaxon)), ' species in the set ', "'", save_run, "'")
# print(table(SPAT.FLAG$SPAT_OUT, exclude = NULL))
# length(unique(SPAT.FLAG$searchTaxon))
# unique(SPAT.FLAG$SOURCE)
# unique(SPAT.FLAG$SPAT_OUT)
# 
# 
# ## Just get the records that were not spatial outliers
# SDM.SPAT.ALL = subset(SPAT.FLAG, SPAT_OUT == "TRUE")
# unique(SDM.SPAT.ALL$SPAT_OUT)   
# unique(SDM.SPAT.ALL$SOURCE) 
# length(unique(SDM.SPAT.ALL$searchTaxon))
# 
# 
# ## What percentage of records are retained? 1-2% seems reasonable
# message(round(nrow(SDM.SPAT.ALL)/nrow(SPAT.FLAG)*100, 2), " % records retained")
# 
# 
# 
# #########################################################################################################################
# ## 3). FILTER DATA TO THE CLEAN RECORDS
# #########################################################################################################################
# 
# 
# #########################################################################################################################
# ## Join data :: exclude the decimal lat/long, check the length 
# identical(dim(TIB.GBIF)[1],dim(FLAGS)[1])#;length(GBIF.SPAT.OUT)
# names(FLAGS)[1]    = c("coord_spp")
# identical(COMBO.AWAP.CONVERT$searchTaxon, FLAGS$coord_spp)                                            ## order matches
# 
# 
# #########################################################################################################################
# ## Is the species column the same as the searchTaxon column?
# TEST.GEO   = cbind(COMBO.AWAP.CONVERT, FLAGS)
# identical(TEST.GEO$searchTaxon, TEST.GEO$coord_spp)                                                     ## order matches
# 
# 
# ## Keep records which passed the GBIF and spatial test
# dim(subset(TEST.GEO, summary == "TRUE")) 
# CLEAN.TRUE = subset(TEST.GEO, summary == "TRUE")
# unique(CLEAN.TRUE$summary)
# 
# 
# ## What percentage of records are retained?
# identical(dim(CLEAN.TRUE)[1], (dim(COMBO.AWAP.CONVERT)[1] - dim(subset(TEST.GEO, summary == "FALSE"))[1]))
# length(unique(CLEAN.TRUE$searchTaxon))
# message(round(dim(CLEAN.TRUE)[1]/dim(TEST.GEO)[1]*100, 2), " % records retained")                                               





#########################################################################################################################
## 4). CREATE SHAPEFILES TO CHECK OUTLIERS ARCMAP
#########################################################################################################################


# ## Select columns
# GBIF.ALA.CHECK  = select(CLEAN.TRUE, searchTaxon, scientificName, lat, lon, SOURCE, year,
#                          country, locality, basisOfRecord, institutionCode)
# 
# 
# ## Rename the fields so that ArcMap can handle them
# GBIF.ALA.CHECK     = dplyr::rename(GBIF.ALA.CHECK,
#                                    TAXON     = searchTaxon,
#                                    LAT       = lat,
#                                    LON       = lon,
#                                    SC_NAME   = scientificName,
#                                    BASIS     = basisOfRecord,
#                                    LOCAL     = locality,
#                                    INSTIT    = institutionCode,
#                                    COUNTRY   = country,
#                                    YEAR      = year)
# names(GBIF.ALA.CHECK)
# 
# 
# #########################################################################################################################
# ## Then create SPDF
# GBIF.ALA.SPDF    = SpatialPointsDataFrame(coords      = GBIF.ALA.CHECK[c("LON", "LAT")],
#                                           data        = GBIF.ALA.CHECK,
#                                           proj4string = CRS.WGS.84)
# 
# 
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
# ## Write the whole shapefile out just in case
# writeOGR(obj = GBIF.ALA.SPDF, dsn = "./data/ANALYSIS/CLEAN_GBIF", 
#          layer = paste0("CLEAN_SPDF_", save_run), driver = "ESRI Shapefile")





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
  ## Here's what the function will produce
  # NICHE.AUS.DF = NICHE.AUS %>% 
  #   as.data.frame() %>% 
  #   dplyr::select(., searchTaxon, one_of(env.variables))
  
  
  ## Create df for worldclim variables...............................................................DIANA CHANGE TO AUS
  ENV.DF = NICHE.1KM %>% 
    as.data.frame() %>% 
    dplyr::select(., searchTaxon, one_of(env.variables))
  
  
  ## Create df for drought variables
  DROUGHT.DF = NICHE.1KM %>% 
    as.data.frame() %>% 
    dplyr::select(., searchTaxon, one_of(drought.variables))
  
  
  ## Create df for heat variables
  HEAT.DF = NICHE.1KM %>% 
    as.data.frame() %>% 
    dplyr::select(., searchTaxon, one_of(heat.variables))
  
  
  ## Create df for radiation variables
  RAD.DF = NICHE.1KM %>% 
    as.data.frame() %>% 
    dplyr::select(., searchTaxon, one_of(rad.variables))
  
  # NICHE.GLO.DF = NICHE.1KM.84 %>% 
  #   as.data.frame() %>% 
  #   dplyr::select(., searchTaxon, one_of(env.variables))
  
  
  #########################################################################################################################
  ## Remove the NA worldclim data :: 99% retained
  ENV.DAT = completeFun(ENV.DF, "PET")
  head(niche_estimate (DF = ENV.DAT, colname = "Annual_mean_temp"))[1:5]  ## Including the q05 and q95
  
  message(round(nrow(ENV.DAT)/nrow(ENV.DF)*100, 2), " % records retained")
  message(round(length(unique(ENV.DAT$searchTaxon))/
                  length(unique(ENV.DF$searchTaxon))*100, 2), "% species retained for worldclim")
  
  ## Remove the NA drought data :: 60% retained
  DROUGHT.DAT = completeFun(DROUGHT.DF, "Drought_freq_extr")
  head(niche_estimate (DF = DROUGHT.DAT, colname = "Drought_freq_extr"))[1:5]  ## Including the q05 and q95
  
  message(round(nrow(DROUGHT.DAT)/nrow(DROUGHT.DF)*100, 2), " % records retained")
  message(round(length(unique(DROUGHT.DAT$searchTaxon))/
                  length(unique(ENV.DF$searchTaxon))*100, 2), "% species retained for drought")
  
  
  ## Remove the NA heat data :: 35% retained
  HEAT.DAT = completeFun(HEAT.DF, "HW_CUM_ALL")
  head(niche_estimate (DF = HEAT.DAT, colname = "HW_CUM_ALL"))[1:5]  ## Including the q05 and q95
  
  message(round(nrow(HEAT.DAT)/nrow(HEAT.DF)*100, 2), " % records retained")
  message(round(length(unique(HEAT.DAT$searchTaxon))/
                  length(unique(ENV.DF$searchTaxon))*100, 2), "% species retained for heat")
  
  ## Remove the NA radiation data :: 62.98 % retained
  RAD.DAT = completeFun(RAD.DF, "mean_monthly_par")
  head(niche_estimate (DF = RAD.DAT, colname = "mean_monthly_par"))[1:5]  ## Including the q05 and q95
  
  message(round(nrow(RAD.DAT)/nrow(RAD.DF)*100, 2), " % records retained for radiation")
  message(round(length(unique(RAD.DAT$searchTaxon))/
                  length(unique(ENV.DF$searchTaxon))*100, 2), "% species retained for radiation")
  
  
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
  message('Estimating drought niches for ', length(unique(DROUGHT.DAT$searchTaxon)), ' species in the set ', "'", save_run, "'")
  DROUGHT.NICHE <- drought.variables %>% 
    
    ## Pipe the list into lapply
    lapply(function(x) {
      
      ## Now, use the niche width function on each colname
      ## Also, need to figure out how to make the aggregating column generic (species, genus, etc.)
      ## currently it only works hard-wired..............................................................................
      niche_estimate(DF = DROUGHT.DAT, colname = x)
      
      ## Would be good to remove the duplicate columns here..............................................................
      
    }) %>% 
    
    ## finally, create one dataframe for all niches
    as.data.frame
  
  ## Remove duplicate Taxon columns and check the output :: would be great to skip these columns when running the function
  DROUGHT.NICHE <- subset(DROUGHT.NICHE, select = -c(searchTaxon.1,  searchTaxon.2,  searchTaxon.3,  searchTaxon.4,
                                                     searchTaxon.5,  searchTaxon.6))
  
  #########################################################################################################################
  message('Estimating drought niches for ', length(unique(HEAT.DAT$searchTaxon)), ' species in the set ', "'", save_run, "'")
  HEAT.NICHE <- heat.variables %>% 
    
    ## Pipe the list into lapply
    lapply(function(x) {
      
      ## Now, use the niche width function on each colname
      ## Also, need to figure out how to make the aggregating column generic (species, genus, etc.)
      ## currently it only works hard-wired..............................................................................
      niche_estimate(DF = HEAT.DAT, colname = x)
      
      ## Would be good to remove the duplicate columns here..............................................................
      
    }) %>% 
    
    ## finally, create one dataframe for all niches
    as.data.frame
  
  ## Remove duplicate Taxon columns and check the output :: would be great to skip these columns when running the function
  HEAT.NICHE <- subset(HEAT.NICHE, select = -c(searchTaxon.1,  searchTaxon.2,  searchTaxon.3,  searchTaxon.4,
                                               searchTaxon.5,  searchTaxon.6, searchTaxon.7))
  
  
  #########################################################################################################################
  message('Estimating drought niches for ', length(unique(HEAT.DAT$searchTaxon)), ' species in the set ', "'", save_run, "'")
  RAD.NICHE <- rad.variables %>% 
    
    ## Pipe the list into lapply
    lapply(function(x) {
      
      ## Now, use the niche width function on each colname
      ## Also, need to figure out how to make the aggregating column generic (species, genus, etc.)
      ## currently it only works hard-wired..............................................................................
      niche_estimate(DF = RAD.DAT, colname = x)
      
      ## Would be good to remove the duplicate columns here..............................................................
      
    }) %>% 
    
    ## finally, create one dataframe for all niches
    as.data.frame
  
  ## Remove duplicate Taxon columns and check the output :: would be great to skip these columns when running the function
  RAD.NICHE <- subset(RAD.NICHE, select = -c(searchTaxon.1))
  


  #########################################################################################################################
  ## Add global counts for each species, and record the total number of taxa processed
  ## This is actually redundant - the records will the be the same as the maxent_records
  Global_records = as.data.frame(table(NICHE.1KM$searchTaxon))
  names(Global_records) = c("searchTaxon", "Global_records")
  identical(nrow(Global_records), nrow(ENV.NICHE))
  
  
  ## Add the count of Australian records - this is not necessarily the same as maxent_records 
  # Aus_records = as.data.frame(table(NICHE.AUS.DF$searchTaxon))
  # names(Aus_records) = c("searchTaxon", "Aus_records")
  # identical(nrow(Aus_records), nrow(AUS.NICHE))
  
  
  #########################################################################################################################
  ## Join all the tables together - this allows all the data to be used for each vaiable group
  ALL.NICHE = join_all(list(Global_records, ENV.NICHE, DROUGHT.NICHE, HEAT.NICHE, RAD.NICHE),  by = 'searchTaxon', type = 'full')
  

  
  
  
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