#########################################################################################################################
## GET SPATIAL OUTLIERS AND DUPLICATES?
#########################################################################################################################


## Load raw GBIF data
load("./data/base/HIA_LIST/GBIF/GBIF_TRIM.RData")
dim(GBIF.TRIM)
length(unique(GBIF.TRIM$searchTaxon))  ## has the list update with extra species? YES!


#########################################################################################################################
## Try 'Geoclean' : this provides several different tests to clean datasets with geographic coordinates...
# GBIF.TRIM.GEO = dplyr::rename(GBIF.TRIM, identifier = searchTaxon, 
#                               XCOOR = lat, YCOOR = lon)
# 
# ## 
# GBIF.GeoClean.table = GeoClean(GBIF.TRIM.GEO, #countrycentroid = TRUE, 
#                                outp = 'detailed') ## one column for each check
# GBIF.GeoClean       = GeoClean(GBIF.TRIM.GEO, outp = 'summary')
# save(GBIF.GeoClean.table, file = paste("./data/base/HIA_LIST/GBIF/GBIF_GEOCLEAN_TABLE.RData"))
# head(GBIF.GeoClean.table)                                                                ##  Too restrictive!
# unique(GBIF.GeoClean)                                                                    ##  FALSE = suspicious coordinates
# 
# 
# ## Check the output
# length(GBIF.GeoClean)
# length(GBIF.GeoClean[GBIF.GeoClean == TRUE])
# length(GBIF.GeoClean[GBIF.GeoClean == FALSE])                                            ##  FALSE = suspicious coordinates


#########################################################################################################################
## And 'duplicated'...returns a logical vector of which rows of a table are duplicates of a row with smaller subscripts.
# GBIF.dups <- duplicated(GBIF.TRIM)
# unique(GBIF.dups)                                                                ##  TRUE  = Duplicated
# 
# 
# ## Check the output
# length(GBIF.dups)
# length(GBIF.dups[GBIF.dups == TRUE])
# length(GBIF.dups[GBIF.dups == FALSE])


#############################################################################################################################
## Try PPP? This is too restrictive. Needs to be run within each species, not across whole dataset
# sp.n      = "Syzygium floribundum"
# test      = subset(GBIF.RASTER.CONTEXT, searchTaxon == sp.n)[, c("lon", "lat")]

# x    <- test$lon ; y<-test$lat
# w    <- ripras(x, y)
# wp   <- ppp(x,y, window = w)
# dupv <- duplicated.ppp(wp)
# 
# 
# ## Check
# length(dupv[dupv == TRUE])          ## Run on each species, not the whole dataset
# length(dupv[dupv == FALSE])         ## Check the maps: is it getting rid of good data?
# 
# 
# ## Assing to
# x2   <- x[which(dupv == FALSE)] 
# y2   <- y[which(dupv == FALSE)]
# 
# 
# ## coordinates of points with no duplicates
# x2<-x[which(dupv==FALSE)] ; y2<-y[which(dupv == FALSE)]


#############################################################################################################################
## Now add these two columns to the intial GBIF file
# GBIF.TRIM$GEOCLEAN    = GBIF.GeoClean
# GBIF.TRIM$DUPLICATED  = GBIF.dups
# names(GBIF.TRIM)