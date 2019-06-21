#########################################################################################################################
## Background data proportional sampling
## Currently, the background data doesn't have exactly the same bias as the analysis data. But if we re-run all the data 
## creation with the new urban data, we could just include this in the background data? Then include the source column


## Spatial pattern of BG data needs to match 
## Binary raster, rasterize ALA and GBIF :: grid cell. Far more presences, due to spatial scale 
## Which is less intensive computational?
## Tree inventory data rarefy this too. But could -
## Maybe thinning just the inventory data, not ALA.
## Biases vs. 

## Just how much difference does autocorrelation make?
## IDW to occurrences and background

## Add inventory data to GB points - rarefy to 5km, 10km..................................................................
## Run 5 species with BG using all data, vs rarefy inventory data.........................................................
## Binary file for ALA/GBIF


## Different data sources :: targetted background combines all the data together
## Retain the true bias of occurrences.
## Sensitivity analysis for this :: raregy 5km vs. 10km..................................................................
## estimate density surface across ALA and GBIF, use this to rarefy?
## matching occ to inventory data loses info...

#########################################################################################################################
## MESS map sampling

## Try running a few species manually to test. Did the code changes affect the MESS map creation? EG the names of rasters,
## etc. Could be a minor change that affects the processing within the loop

## Also,  try running 5 species on separate cores - this mimics what the cluster analysis will do. Are the errors simply
## Environment issues?


# Converting raster MESS maps to polygons under ac85bi30 scenario for Eucalyptus_camaldulensis
# H:\green_cities_sdm
# Eucalyptus_camaldulensis failed
# [[1]]
# [1] TRUE

#Either there was an error with the request, or the servers may be down (try again later). If this problem persists please #notify the ALA4R maintainers by lodging an issue at https://github.com/AtlasOfLivingAustralia/ALA4R/issues/ or emailing #support@ala.org.au

#Error in ras_dist(x = k, lat = lat, lon = lon, ras = ras, weights = TRUE) : 
# object 'ras' not found


# Warning messages:
# 1: In ogrFIDs(dsn = dsn, layer = layer) : no features found
# 2: In value[[3L]](cond) : Eucalyptus_camaldulensis: no features found


# >             novel_current_poly <- polygonizer(sprintf('%s/%s%s.tif',   MESS_dir, species, "_current_novel"))
# H:\green_cities_sdm
# Error in readOGR(dirname(outshape), layer = basename(outshape), verbose = !quietish) : 
  # no features found
# In addition: Warning message:
# In ogrFIDs(dsn = dsn, layer = layer) : no features found


## Why don't the files list?
#Running summary of SDM predictions within SUAs for Acer_palmatum using SUA_CODE16 shapefile
#Calcualting mean of 2030 GCMs for Acer_palmatum
#doing Acer_palmatum | Logistic > 0.1924 for 2030
#Error in .local(.Object, ...) :


#Creating polygon list under mc85bi30 scenario for Pyrus_calleryana
#Error in RGEOSBinTopoFunc(spgeom1, spgeom2, byid, id, drop_lower_td, unaryUnion_if_byid_false,  : 
#  TopologyException: Input geom 1 is invalid: Ring Self-intersection at or near point 316091.36989912001 -1277934.1729915701 at 316091.36989912001 -1277934.1729915701
#In addition: There were 28 warnings (use warnings() to see them)


#Error in `names<-`(`*tmp*`, value = grid_names) : 
#  incorrect number of layer names

#Creating polygon list under ac85bi30 scenario for Callistemon_viminalis
#Create MESS panel maps for Callistemon_viminalis under ac85bi30 scenario
#Error in UseMethod("depth") : 
#  no applicable method for 'depth' applied to an object of class "NULL"


#Creating MESS directory for Prunus_cerasifera
#Creating mess maps of each current environmental predictor for Prunus_cerasifera
#Error in UseMethod("depth") : 
#  no applicable method for 'depth' applied to an object of class "NULL"
#In addition: There were 50 or more warnings (use warnings() to see the first 50)


#########################################################################################################################
## Test exotics BS, INV
GBIF.spp      = popular.test ## your list of species
GBIF.spp      = as_utf8(GBIF.spp, normalize = TRUE)  ## Check the species names have the right characters
save_run      = "HIA_TEST_POP"                    ## a variable to append the run name to the output files
read_data     = "FALSE"                              ## only needed to load the SDM table in
save_data     = "TRUE"                              ## only needed to load the SDM table in
check_maps    = "TRUE"                               ## Create maps, shapefiles and histograms of each speices?

GBIF_path     = "./data/base/HIA_LIST/GBIF/"  ## The path where GBIF data is stored
ALA_path      = "./data/base/HIA_LIST/ALA/"   ## The path where ALA data is stored place
DATA_path     = "./data/ANALYSIS/"            ## The path where the final data from/for analyses are stored 
OCC_SOURCE    = c("ALA", "GBIF", "INVENTORY") ## Data source :: for trees, ALL. For non trees, ALA/GBIF, for occ and bg
calc_niche    = "TRUE"

maxent_path   = './output/maxent/HIA_TEST_INV/'    ## The directory where maxent files are saved               
maxent_dir    = 'output/maxent/HIA_TEST_INV'       ## Another version of the path needed to run maxent loop

bs_path       = './output/maxent/HIA_TEST_INV_BS/' ## Backwards selection directory
bs_dir        = 'output/maxent/HIA_TEST_INV_BS'    ## the directory to harvest results : BS dir?
results_dir   = bs_dir   


#########################################################################################################################
## Test exotics BS, INV
GBIF.spp      = popular.test ## your list of species
GBIF.spp      = as_utf8(GBIF.spp, normalize = TRUE)  ## Check the species names have the right characters
save_run      = "HIA_TEST_ALA_BS"                    ## a variable to append the run name to the output files
read_data     = "FALSE"                              ## only needed to load the SDM table in
save_data     = "TRUE"                              ## only needed to load the SDM table in
check_maps    = "TRUE"                               ## Create maps, shapefiles and histograms of each speices?

GBIF_path     = "./data/base/HIA_LIST/GBIF/"         ## The path where GBIF data is stored
ALA_path      = "./data/base/HIA_LIST/ALA/"          ## The path where ALA data is stored place
DATA_path     = "./data/ANALYSIS/"                   ## The path where the final data from/for analyses are stored 
OCC_SOURCE    = c("ALA", "GBIF")                     ## Data source :: for trees, ALL. For non trees, ALA/GBIF, for occ and bg
calc_niche    = "TRUE"

maxent_path   = './output/maxent/HIA_TEST_ALA/'      ## The directory where maxent files are saved               
maxent_dir    = 'output/maxent/HIA_TEST_ALA'         ## Another version of the path needed to run maxent loop

bs_path       = './output/maxent/HIA_TEST_ALA_BS/'   ## Backwards selection directory
bs_dir        = 'output/maxent/HIA_TEST_ALA_BS'      ## the directory to harvest results : BS dir?
results_dir   = bs_dir   



#########################################################################################################################
## Test exotics BS, ALA
GBIF.spp      = inv.test ## your list of species
GBIF.spp      = as_utf8(GBIF.spp, normalize = TRUE)  ## Check the species names have the right characters
save_run      = "HIA_TEST_ALA_BS"                    ## a variable to append the run name to the output files
read_data     = "FALSE"                              ## only needed to load the SDM table in
save_data     = "TRUE"                              ## only needed to load the SDM table in
check_maps    = "TRUE"                               ## Create maps, shapefiles and histograms of each speices?

GBIF_path     = "./data/base/HIA_LIST/GBIF/"         ## The path where GBIF data is stored
ALA_path      = "./data/base/HIA_LIST/ALA/"          ## The path where ALA data is stored place
DATA_path     = "./data/ANALYSIS/"                   ## The path where the final data from/for analyses are stored 
OCC_SOURCE    = c("ALA", "GBIF")                     ## Data source :: for trees, ALL. For non trees, ALA/GBIF, for occ and bg
calc_niche    = "TRUE"

maxent_path   = './output/maxent/HIA_TEST_ALA/'      ## The directory where maxent files are saved               
maxent_dir    = 'output/maxent/HIA_TEST_ALA'         ## Another version of the path needed to run maxent loop

bs_path       = './output/maxent/HIA_TEST_ALA_BS/'   ## Backwards selection directory
bs_dir        = 'output/maxent/HIA_TEST_ALA_BS'      ## the directory to harvest results : BS dir?
results_dir   = bs_dir   


#########################################################################################################################
## Test exotics BS, INV
GBIF.spp      = inv.test ## your list of species
GBIF.spp      = as_utf8(GBIF.spp, normalize = TRUE)  ## Check the species names have the right characters
save_run      = "HIA_TEST_INV_BS"                    ## a variable to append the run name to the output files
read_data     = "FALSE"                              ## only needed to load the SDM table in
save_data     = "FALSE"                              ## only needed to load the SDM table in
check_maps    = "TRUE"                               ## Create maps, shapefiles and histograms of each speices?

GBIF_path     = "./data/base/HIA_LIST/GBIF/"  ## The path where GBIF data is stored
ALA_path      = "./data/base/HIA_LIST/ALA/"   ## The path where ALA data is stored place
DATA_path     = "./data/ANALYSIS/"            ## The path where the final data from/for analyses are stored 
OCC_SOURCE    = c("ALA", "GBIF", "INVENTORY") ## Data source :: for trees, ALL. For non trees, ALA/GBIF, for occ and bg
calc_niche    = "TRUE"

maxent_path   = './output/maxent/HIA_TEST_INV/'    ## The directory where maxent files are saved               
maxent_dir    = 'output/maxent/HIA_TEST_INV'       ## Another version of the path needed to run maxent loop

bs_path       = './output/maxent/HIA_TEST_INV_BS/' ## Backwards selection directory
bs_dir        = 'output/maxent/HIA_TEST_INV_BS'    ## the directory to harvest results : BS dir?
results_dir   = bs_dir   







#########################################################################################################################
## Test exotics BS, ALA
GBIF.spp      = hollow.test.spp ## your list of species
GBIF.spp      = as_utf8(GBIF.spp, normalize = TRUE)  ## Check the species names have the right characters
save_run      = "HOLLOW_10_ALA_BS"                    ## a variable to append the run name to the output files
read_data     = "FALSE"                              ## only needed to load the SDM table in
save_data     = "FALSE"                              ## only needed to load the SDM table in
check_maps    = "TRUE"                               ## Create maps, shapefiles and histograms of each speices?

GBIF_path     = "./data/base/HIA_LIST/GBIF/"  ## The path where GBIF data is stored
ALA_path      = "./data/base/HIA_LIST/ALA/"   ## The path where ALA data is stored place
DATA_path     = "./data/ANALYSIS/"            ## The path where the final data from/for analyses are stored 
OCC_SOURCE    = c("ALA", "GBIF")              ## Data source :: for trees, ALL. For non trees, ALA/GBIF, for occ and bg
calc_niche    = "TRUE"

maxent_path   = './output/maxent/HOLLOW_10_ALA/'    ## The directory where maxent files are saved               
maxent_dir    = 'output/maxent/HOLLOW_10_ALA'       ## Another version of the path needed to run maxent loop

bs_path       = './output/maxent/HOLLOW_10_ALA_BS/' ## Backwards selection directory
bs_dir        = 'output/maxent/HOLLOW_10_ALA_BS'    ## the directory to harvest results : BS dir?
results_dir   = bs_dir   





















