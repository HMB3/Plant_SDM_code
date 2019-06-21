Update CALCULATE_1KM_NICHES.R

- Calculate environmental ranges using the 1km SDM data, but also the extra variables like PET, etc.

fatal: Unable to create 'H:/green_cities_sdm/.git/index.lock': File exists.

Another git process seems to be running in this repository, e.g.
an editor opened by 'git commit'. Please make sure all processes
are terminated then try again. If it still fails, a git process
may have crashed in this repository earlier:
remove the file manually to continue.



## Test the code on one species at a time
GBIF.spp = HOLLOW.SPP[1]  ## 20, 40, 60, 80, 100


## Step 7 :: Run maxent on a table of all species, using targetted background selection, then backwards selection
source('./R/7_RUN_MAXENT.R', echo = TRUE)


## Step 8 :: Create habitat suitability maps for each species using six GCMs and three time slices (2030/50/70). 
## Then summarise maxent results and estimate species presences in significant urban areas under climate change
## Then run MESS maps. Or, the MESS maps could be run before the combination step
source('./R/8_MAP_SDM_COMBINE.R', echo = TRUE)



# ImportError: No module named site
# Error in ogrInfo(dsn = dsn, layer = layer, encoding = encoding, use_iconv = use_iconv,  : 
  # Cannot open data source
  
#H:\green_cities_sdm\data\Koppen_10000m_WGS84.tif