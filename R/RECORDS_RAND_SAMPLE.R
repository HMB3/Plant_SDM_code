#########################################################################################################################
################################################# THIN RECORDS ##########################################################
#########################################################################################################################


## This code spatially thins species occurrence records to help address problems associated with spatial sampling biases. 
## Ideally, thinning removes the fewest records necessary to substantially reduce the effects of sampling bias, while 
## simultaneously retaining the greatest amount of useful information. 

## Use the SDM dataset - one record per 1km cell - rather than the original dataset......................................

## Use the stratified function: https://www.rdocumentation.org/packages/fifer/versions/1.0/topics/stratified
## spThin didn't really work https://cran.r-project.org/web/packages/spThin/vignettes/spThin_vignette.html


## Create lists
source('./R/HIA_LIST_MATCHING.R')


#########################################################################################################################
## 1). ADD COLUMNS FOR AUS RAINFALL AND STATE CATEGORIES
#########################################################################################################################


#########################################################################################################################
## Load GBIF data and rain shapefile
BIAS.DATA.ALL        = readRDS("./data/base/HIA_LIST/COMBO/SDM_DATA_CLEAN_052018.rds")        
SPP.BIAS             = intersect(SPP.BIAS, SUA.spp)    ## just re-run the models for species on the list
SPP_BIAS             = gsub(" ", "_", SPP.BIAS)


## Project the SDM data into WGS
BIAS.DATA.ALL <- spTransform(BIAS.DATA.ALL, CRS("+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
projection(BIAS.DATA.ALL)


## Get the coordinates
BIAS.COORDS = coordinates(BIAS.DATA.ALL)
class(BIAS.COORDS)
tail(BIAS.COORDS)
tail(BIAS.DATA.ALL)[, c(1:3)]   ## indices match


## Bind the coordinates to the SDM table
BIAS.COORDS   = as.data.frame(BIAS.COORDS)
BIAS.DATA.DF = cbind(as.data.frame(BIAS.DATA.ALL), BIAS.COORDS)
BIAS.DATA.DF = BIAS.DATA.DF[, c(1:22)]
names(BIAS.DATA.DF)
dim(BIAS.DATA.DF)


#########################################################################################################################
## Intersect BOM with rainfall data
BIAS.DATA.SP   = SpatialPointsDataFrame(coords      = BIAS.DATA.DF[c("lon", "lat")], 
                                        data        = BIAS.DATA.DF,
                                        proj4string = CRS.WGS.84)


## Project data :: this takes a long time for rain, lots of features
AUS.WGS    = spTransform(aus, CRS.WGS.84)
RAIN.WGS   = spTransform(AUS_RAIN, CRS.WGS.84)
AUS.STATE  = AUS.WGS[, c("name")]
names(AUS.STATE)[names(AUS.STATE) == 'name'] <- 'AUS_STATE'
AUS.RAIN   = RAIN.WGS[, c("AUS_RN_ZN")] 


## Check projections and plot
projection(BIAS.DATA.SP);projection(AUS.STATE);projection(AUS.RAIN)
BIAS.DATA.CHECK  <- BIAS.DATA.SP[BIAS.DATA.SP$searchTaxon == "Acacia floribunda",]
plot(BIAS.DATA.CHECK)
plot(RAIN.WGS, main = "Average Australian rainfall, 1961-1990")
plot(AUS.STATE, main = "Australian states")


#########################################################################################################################
## Run join between species records and Australian states...
STATE.JOIN        = over(BIAS.DATA.SP, AUS.STATE)                   
BIAS.STATE        = cbind.data.frame(BIAS.DATA.SP, STATE.JOIN)
names(BIAS.STATE)
unique(BIAS.STATE$AUS_STATE)
saveRDS(BIAS.STATE, file = paste("./data/base/HIA_LIST/COMBO/BIAS_STATE.rds"))


## Recreate sp dataframe
BIAS.STATE.SP   = SpatialPointsDataFrame(coords       = BIAS.STATE[c("lon", "lat")], 
                                         data        = BIAS.STATE,
                                         proj4string = CRS.WGS.84)


#########################################################################################################################
## Now join RAIN categories
RAIN.JOIN          = over(BIAS.STATE.SP, AUS.RAIN)
BIAS.STATE.RAIN    = cbind.data.frame(BIAS.STATE.SP, RAIN.JOIN)
names(BIAS.STATE.RAIN)
unique(BIAS.STATE$AUS_STATE)


## Rename and remove
BIAS.STATE.RAIN = subset(BIAS.STATE.RAIN, select = -c(lon.1,  lat.1, optional, optional.1, lon.2, lat.2))
names(BIAS.STATE.RAIN)
unique(BIAS.STATE.RAIN$AUS_STATE)
unique(BIAS.STATE.RAIN$AUS_RN_ZN)


BIAS.STATE.SP   = SpatialPointsDataFrame(coords       = BIAS.STATE[c("lon", "lat")], 
                                         data        = BIAS.STATE,
                                         proj4string = CRS.WGS.84)


## Save data
saveRDS(BIAS.STATE.RAIN, file = paste("./data/base/HIA_LIST/COMBO/BIAS_STATE_RAIN.rds"))
#BIAS.STATE.RAIN = readRDS("./data/base/HIA_LIST/COMBO/BIAS_STATE_RAIN.rds")





#########################################################################################################################
## 2). CREATE A RANDOM SUB-SAMPLE OF X% OF SPP RECORDS, STRATIFIED BY RAINFALL, THE RUN SDMs AGAIN FOR SELECTED TAXA
#########################################################################################################################


# Shawn :: The Biodiverse approach would snap the coordinates to the centre of each cell.  It would be better to select 
# one point at random within each cell, which should be pretty simple in R.  

# The random sample for a large starting number would probably be reasonable. With either of the above approaches, you also 
# have the problem that you might lose some of your environmental variability, in which case it might be better to thin 
# in environmental space as well as geospatial.  Maybe a stratified thinning would be sensible?


## Use the stratified function: https://www.rdocumentation.org/packages/fifer/versions/1.0/topics/stratified

## For each species in the list :: 

## subset the data to just those species records

## subset the data to just those species records in NSW, and outside NSW

## Subsample within all rainfall bands - take 50% of total records

## Then combine the stratified sample with those records outside NSW


#########################################################################################################################
## Now create the variables needed to access current environmental conditions + their names in the functions
sdm.predictors <- c("Annual_mean_temp",    "Mean_diurnal_range",  "Isothermality",      "Temp_seasonality",  
                    "Max_temp_warm_month", "Min_temp_cold_month", "Temp_annual_range",  
                    "Mean_temp_wet_qu",    "Mean_temp_dry_qu",    
                    "Mean_temp_warm_qu",   "Mean_temp_cold_qu",   
                    
                    "Annual_precip",       "Precip_wet_month",   "Precip_dry_month",    "Precip_seasonality",  
                    "Precip_wet_qu",       "Precip_dry_qu",      
                    "Precip_warm_qu",      "Precip_col_qu")


## A-priori worldclim predictors
sdm.select     <- c("Annual_mean_temp", "Temp_seasonality",    "Max_temp_warm_month", "Min_temp_cold_month",
                    "Annual_precip",    "Precip_seasonality",  
                    "Precip_wet_month", "Precip_dry_month")      


## Create a raster stack of current environmental conditions if needed
## This code is probably not needed... 
i  <- match(sdm.predictors, sdm.predictors)
ff <- file.path('./data/base/worldclim/world/0.5/bio/current',
                sprintf('bio_%02d.tif', i))


## Name the grids :: these should be indentical
env.grids.current = stack(sub('0.5', '1km', ff))
names(env.grids.current) <- sdm.predictors[i]
identical(names(env.grids.current),sdm.predictors)


#########################################################################################################################
## Check data subsets
unique(BIAS.STATE.RAIN$AUS_RN_ZN)
unique(BIAS.STATE.RAIN$AUS_STATE)


#########################################################################################################################
## Loop over all the species with bias: spp = SPP.BIAS[1]
lapply(SPP.BIAS, function(spp){ 
  
  ## Skip the species if the directory already exists, before the loop
  ## Create a separate directory for these subsampled species. Then we can compare the folder output with the originals
  outdir <- 'output/maxent/SPP_BIAS'
  if(dir.exists(file.path(outdir, gsub(' ', '_', spp)))) {
    message('Skipping ', spp, ' - already run.')
    invisible(return(NULL))
    
  }
  
  ## Print the taxa being processed to screen
  if(spp %in% BIAS.STATE.RAIN$searchTaxon) {
    message('Doing ', spp) 
    
    ## Do the stratification by rainfall, etc, in here...................................................................
    
    ## First, subset the records to only the taxa being processed
    occurrence <- subset(BIAS.STATE.RAIN, searchTaxon == spp)
    
    ## Now subset to the species data to those records inside and outside NSW
    ## Could change this to a generic "state" argument
    occurrence.nsw <- subset(occurrence, AUS_STATE == "New South Wales")
    occurrence.out <- subset(occurrence, AUS_STATE != "New South Wales")
    
    unique(occurrence.nsw$AUS_STATE)
    unique(occurrence.out$AUS_STATE)
    
    ## Now take a 50% random sample from all rainfall groups, just for the points within NSW
    ## Set a large starting number R? How large is large?
    set.seed(123458)
    strat.nsw = stratified(occurrence.nsw, "AUS_RN_ZN", .5)
    dim(strat.nsw)[1]/dim(occurrence.nsw)[1]*100
    
    ## Now combine the two data frames :: about 50% of the original data remains
    occ.strat = as.data.frame(bind_rows(strat.nsw, occurrence.out))
    names(occ.strat)
    dim(occ.strat)[1]/dim(occurrence)[1]*100
    
    ## Now get the background points. These can come from any spp, other than the modelled species.
    ## Also should the background points come from NSW too?
    background <- subset(BIAS.STATE.RAIN, searchTaxon != spp)
    
    ## Re-convert to a spatial points df
    background    = SpatialPointsDataFrame(coords      = background[c("lon", "lat")], 
                                           data        = background,
                                           proj4string = CRS.WGS.84)
    
    occ.strat    = SpatialPointsDataFrame(coords      = occ.strat[c("lon", "lat")], 
                                          data        = occ.strat,
                                          proj4string = CRS.WGS.84)
    
    ## Re-project
    background <- spTransform(background,  CRS("+init=ESRI:54009 +proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs84=0,0,0"))
    occ.strat  <- spTransform(occ.strat,   CRS("+init=ESRI:54009 +proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs84=0,0,0"))
    #writeOGR(obj = SPP.SAMPLE, dsn = "./data/base/HIA_LIST/BIAS", layer = "SPP_RAND_SAMPLE", driver = "ESRI Shapefile")
    
    ## Finally fit the models using FIT_MAXENT_TARG_BG. Also use tryCatch to skip any exceptions
    tryCatch(
      FIT_MAXENT_TARG_BG(occ                     = occ.strat, 
                         bg                      = background, 
                         sdm.predictors          = sdm.select, 
                         name                    = spp, 
                         outdir                  = outdir, 
                         template.raster,
                         min_n                   = 20,     ## This should be higher...
                         max_bg_size             = 100000, 
                         Koppen                  = Koppen_1975,
                         background_buffer_width = 200000,
                         shapefiles              = TRUE,
                         features                = 'lpq',
                         replicates              = 5,
                         responsecurves          = TRUE),
      
      ## https://stackoverflow.com/questions/19394886/trycatch-in-r-not-working-properly
      #function(e) message('Species skipped ', spp)) ## skip any species for which the function fails
      error = function(cond) {
        
        message(paste('Species skipped ', spp))
        
      })
    
  } else {
    
    message(spp, ' skipped - no data.')         ## This condition ignores species which have no data...
    
  }  
  
})





#########################################################################################################################
## Simple check of the stratification -


## Compare the occurrence points that are output by the SDM code from the BIA folder with the previous folder





#########################################################################################################################
## 3). CREATE MAPS FOR SIX GCMs
#########################################################################################################################


#########################################################################################################################
## Create a list of GCM scenarios, which are used to create maps of habitat suitability 


## Eight of the 40 CMIP5 models assessed in this project have been selected for use in provision of application-ready data. 
## This facilitates efficient exploration of climate projections for Australia.
## https://www.climatechangeinaustralia.gov.au/en/support-and-guidance/faqs/eight-climate-models-data/

## The full list is:
# scen = c("ip85bi50", "mc85bi50", "mg85bi50", "mi85bi50", "mp85bi50", 
#          "mr85bi50", "no85bi50", "ac85bi50", "bc85bi50", "cc85bi50", 
#          "cn85bi50", "gf85bi50", "gs85bi50", "hd85bi50", "he85bi50",
#          "hg85bi50", "in85bi50")


## We can only get the bioclim variables for 6 of these projections......................................................
## Create a lookup table of GCMs using the WORLDCLIM website
h <- read_html('http://www.worldclim.org/cmip5_30s') 
gcms <- h %>% 
  html_node('table') %>% 
  html_table(header = TRUE) %>% 
  filter(rcp85 != '')

id.50 <- h %>% 
  html_nodes(xpath = "//a[text()='bi']/@href") %>% 
  as.character %>% 
  grep('85bi50', ., value = TRUE) %>% 
  basename %>% 
  sub('\\.zip.', '', .)

id.70 <- h %>% 
  html_nodes(xpath = "//a[text()='bi']/@href") %>% 
  as.character %>% 
  grep('85bi70', ., value = TRUE) %>% 
  basename %>% 
  sub('\\.zip.', '', .)


## Work around because the 2030 data comes from CC in Aus, not worldclim
id.30 = gsub("50", "30", id.50)


## Create the IDs
gcms.30 <- cbind(gcms, id.30)
gcms.30$GCM = sub(" \\(#\\)", "", gcms$GCM)

gcms.50 <- cbind(gcms, id.50)
gcms.50$GCM = sub(" \\(#\\)", "", gcms$GCM)  ## sub replaces first instance in a string, gsub = global

gcms.70 <- cbind(gcms, id.70)
gcms.70$GCM = sub(" \\(#\\)", "", gcms$GCM) 


## Now create the scenario lists across which to loop
gcms.50 ; gcms.70 ; gcms.30


## Just get the 6 models picked by CSIRO for Australia, for 2030, 2050 and 2070
scen_2030 = c("mc85bi30", "no85bi30", "ac85bi30", "cc85bi30", "gf85bi30", "hg85bi30")
scen_2050 = c("mc85bi50", "no85bi50", "ac85bi50", "cc85bi50", "gf85bi50", "hg85bi50")
scen_2070 = c("mc85bi70", "no85bi70", "ac85bi70", "cc85bi70", "gf85bi70", "hg85bi70")


## Now divide the current environmental grids by 10
env.grids.current <- stack(
  file.path('./data/base/worldclim/aus/1km/bio/current',   ## ./data/base/worldclim/aus/1km/bio
            sprintf('bio_%02d.tif', 1:19)))

for(i in 1:11) {
  
  ## simple loop
  message(i)
  env.grids.current[[i]] <- env.grids.current[[i]]/10
  
}


#########################################################################################################################
## Name the environmental grids to be used in the mapping code :: this is all 19, because we are using the directory structure
grid.names = c('Annual_mean_temp',    'Mean_diurnal_range',  'Isothermality',      'Temp_seasonality', 
               'Max_temp_warm_month', 'Min_temp_cold_month', 'Temp_annual_range',  'Mean_temp_wet_qu',
               'Mean_temp_dry_qu',    'Mean_temp_warm_qu',   'Mean_temp_cold_qu',  'Annual_precip',
               'Precip_wet_month',    'Precip_dry_month',    'Precip_seasonality', 'Precip_wet_qu',
               'Precip_dry_qu',       'Precip_warm_qu',      'Precip_col_qu')


#########################################################################################################################
## Create 2030 maps BIAS_REV = sort(SPP_BIAS, decreasing = TRUE)
env.grids.2030 = tryCatch(project_maxent_grids(scen_list     = scen_2030,
                                               species_list  = SPP_BIAS,
                                               maxent_path   = "./output/maxent/SPP_BIAS/",
                                               climate_path  = "./data/base/worldclim/aus/1km/bio",
                                               grid_names    = grid.names,
                                               time_slice    = 30,
                                               current_grids = env.grids.current),
                          
                          ## Will this work outside a loop?
                          error = function(cond) {
                            
                            message(paste('Species skipped - check', spp))
                            
                          })





#########################################################################################################################
## 3). COMPARE MAXENT RESULTS WITH AND WITHOUT SAMPLING
#########################################################################################################################


#########################################################################################################################
## First, read in the list of files for the current models, and specify the file path
BIAS.set.var             = "./output/maxent/SPP_BIAS/"


## Create an object for the maxent settings :: using the same variable for every model
BIAS.selection.settings = "Set_variables"  
records_setting          = "COORD_CLEAN"


## Create a file list for each model run
BIAS.tables = list.files(BIAS.set.var)                 ## Chagne this for each variable selection strategy
BIAS_path   = BIAS.set.var                             ## Chagne this for each variable selection strategy
length(BIAS.tables)                                    ## Should match the number of taxa tested


#########################################################################################################################
## Could turn this into a function, and loop over a list of subfolders...
## Then pipe the table list into lapply
## BIAS.RESULTS = bind_BIAS_tables(table_list)
BIAS.RESULTS <- BIAS.tables[c(1:length(BIAS.tables))] %>%         ## currently x species with enough records
  
  ## pipe the list into lapply
  lapply(function(x) {
    
    ## Create the character string
    f <- paste0(BIAS_path, x, "/full/maxentResults.csv")
    
    ## Load each .RData file
    d <- read.csv(f)
    
    ##
    #m <- readRDS(sprintf('%s/%s/full/model.rds', BIAS_path, x))
    m <- readRDS(sprintf('%s/%s/full/maxent_fitted.rds', BIAS_path, x))
    m <- m$me_full
    number_var = length(m@lambdas) - 4                                   ## (the last 4 slots of lambdas are not variables) 
    
    ## Now add a model column
    d = cbind(searchTaxon = x, 
              Settings    = model.selection.settings, 
              Number_var  = number_var, 
              Records     = records_setting, d)  ## see step 7, make a variable for multiple runs
    dim(d)
    
    ## Remove path gunk, and species
    d$Species    = NULL
    d
    
  }) %>%
  
  ## Finally, bind all the rows together
  bind_rows


#########################################################################################################################
## Create a TSS file list for each model run
TSS.tables = list.files(TSS.set.var, pattern = 'species_omission\\.csv$', full.names = TRUE, recursive = TRUE)
TSS_path   = TSS.set.var                            
length(TSS.tables)                                   


max_tss <- sapply(TSS.tables, function(f) {
  
  d <- read.csv(f)
  i <- which.min(d$Test.omission + d$Fractional.area)
  
  c(max_tss = 1 - min(d$Test.omission + d$Fractional.area),
    thr     = d$Corresponding.logistic.value[i])

  
})


## How can we process max TSS?
max_tss = as.data.frame(max_tss)
class(max_tss)
head(max_tss)
BIAS.RESULTS = c(BIAS.RESULTS, max_tss)
BIAS.RESULTS = BIAS.RESULTS[, c(1:9, 65, 10:64)]
names(BIAS.RESULTS)


## This is a summary of BIAS output for current conditions
## Also which species have AUC < 0.7?
dim(BIAS.RESULTS)
head(BIAS.RESULTS, 20)[1:9]
dim(subset(BIAS.RESULTS, Training.AUC < 0.7))


#########################################################################################################################
## Calculate TSS






#########################################################################################################################
#####################################################  TBC ##############################################################
#########################################################################################################################