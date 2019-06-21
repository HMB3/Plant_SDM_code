#########################################################################################################################
######################################### DOWNLOAD CLIMATE DATA ######################################################### 
#########################################################################################################################


## This code can be use to download climate data, eg from the Chelsa website:

## https://www.wsl.ch/lud/chelsa/data/


## It'll be easy to loop through CHELSA urls to download. You can identify the target urls by looking at the page source, 
## right-clicking one of the Download buttons and clicking Inspect. Look for the associated onClick action.


## E.g. for the baseline precip climatologies, urls are like this: 
## https://www.wsl.ch/lud/chelsa/data/climatologies/prec/CHELSA_prec_01_land.7z
## Can help with loops if you need. Just make sure you use mode='wb' as an argument to download.file. 
 
## You can also download in batches by using method='libcurl', in which case you need to provide a vector of source urls 
## and corresponding vector of destinati## on filepaths. The following should work (for the precipitation month-by-year time series):


## Crete a list of year/month names to loop over when downloading
ym <- do.call(paste, c(expand.grid(1979:2013, sprintf('%02d', 1:12)), sep='_'))
chelsa.baseline.precip <- sprintf('https://www.wsl.ch/lud/chelsa/data/timeseries/prec/CHELSA_prec_%s_V1.2.1.tif', ym)
chelsa.baseline.temp   <- sprintf('https://www.wsl.ch/lud/chelsa/data/timeseries/temp/CHELSA_temp_%s_V1.2.1.tif', ym)



## Check the format?
chelsa.baseline.temp[1:5]
chelsa.baseline.precip[1:5]


## Create the directories needed to store the data
dir.create('data/base/CHELSA/precip', recursive = TRUE)
dir.create('data/base/CHELSA/temp',   recursive = TRUE)


## Loop over the list of year/months and download the CHELSA data using the 
download.file(chelsa.baseline.precip, file.path('data/base/CHELSA/precip', 
                                                basename(chelsa.baseline.precip)), mode = 'wb', method = 'libcurl')





