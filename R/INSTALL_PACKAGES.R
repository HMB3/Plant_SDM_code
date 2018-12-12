#########################################################################################################################
###################################### INSTALL MULTIPLE PACKAGES AT ONCE  ############################################### 
#########################################################################################################################


## This code installs a list of packages. On this computer, they currently sit here ::
## .libPaths() 
## C:\Users\user\Documents\R\win-library


#########################################################################################################################
## RE-INSTALL R
#########################################################################################################################


#########################################################################################################################
## 1).  Before you upgrade, build a temp file with all of your old packages.
tmp <- installed.packages()
installedpkgs <- as.vector(tmp[is.na(tmp[,"Priority"]), 1])
save(installedpkgs, file="installed_old.rda")

## 2). Update R here too: https://cran.r-project.org/web/packages/installr/index.html
# installing/loading the package:
# if(!require(installr)) {
#   install.packages("installr"); require(installr)} ## load / install+load installr

## 3. Once you’ve got the new version up and running, reload the saved packages and re-install them from CRAN.
tmp <- installed.packages()
installedpkgs.new <- as.vector(tmp[is.na(tmp[,"Priority"]), 1])
missing <- setdiff(installedpkgs, installedpkgs.new)
install.packages(missing)
update.packages()





#########################################################################################################################
## INSTALL R PACKAGES ON CRAN
#########################################################################################################################


#########################################################################################################################
## create a big list of all the packages needed for a project
## can't do the special github ones like this though...
packages <- c('ggplot2',    'plyr',      'reshape2',  'RColorBrewer',  'scales',      'grid',
              'raster',     'spdep',     'rgdal',     'GISTools',      'data.table',  'dtplyr',
              'Hmisc',      'Cairo',     'lattice',   'gtools',        'ggplot2',     'sjPlot',
              'gridExtra',  'grid',      'sqldf',     'knitr',         'mgcv',        'MuMIn',
              'mvabund',    'asbio',     'mvtnorm',   'testthat',      'rgl',         'R.matlab',
              'ENMeval',    'lavaan',    'semPlot',   'rgdal',         'sp',          'formula.tools',
              'dismo',      'mctest',    'rJava',     'ENMeval',       'SDMTools',    'ALA4R',
              'maxent',     'devtools',  'knitr',     'yaml',          'htmltools',   'microbenchmark',
              'knitr',      'yaml',      'caTools',   'bitops',        'rmarkdown',   'speciesgeocodeR',                  
              'bitops',     'rmarkdown', 'cluster',   'gsubfn',        'functional',  'splitstackshape',
              'EML',        'taxize',    'geonames',  'rWBclimate',    'rfigshare',   'tidyverse',
              'jsonlite',   'zoom',      'bigmemory', 'installr',      'sfsmisc',     'Rtools',
              'RCurl',      'httr',      'scrubr',    'red',           'spatialEco',  'rnaturalearth',
              'gdalUtils',  'rvest',     'ConR',      'rvest',         'rasterVis',   'latticeExtra',
              'ff',         'things',    'raster',    'rgdal',         'data.table',  'RColorBrewer',
              'sp',         'rgeos',     'gdalUtils', 'rmaxent',       'dplyr',       'rasterVis',
              'readr',      'readr',     'parallel',  'dismo',         'tidyr',       'envirem',
              'Taxonstand', 'biomod2',   'vegan',     'session',       'rapportools', 'BiodiversityR',
              'spThin',     'digest',    'CoordinateCleaner',          'PerformanceAnalytics', 'timetk',
              'dat', 'rcreds', 'pkgmaker', 'utf8') # class(packages)

p <- c('ff',    'things',    'raster',        'dismo',        'sp',           'latticeExtra', 'data.table',
       'rgdal', 'rgeos',     'gdalUtils',     'rmaxent',      'readr',        'plyr',         'dplyr',
       'tidyr', 'readr',     'rnaturalearth', 'rasterVis',    'RColorBrewer', 'latticeExtra', 'parallel',
       'ALA4R', 'stringr',   'Taxonstand',    'CoordinateCleaner', 'gsubfn',  'PerformanceAnalytics',
       'rvest', 'magrittr',  'devtools',      'ggplot2',      'reshape2',     'rmarkdown', 'flexdashboard', 'shiny', 'rgbif',
       'ENMeval', 'tibble',  'ncdf4',         'Cairo', 'taxonlookup', 'kgc', 'maptools', 'DataCombine', 'mgcv', 'rsq')


#########################################################################################################################
## This function from the internet that takes a list of packages, installs and loads them:
## See https://gist.github.com/stevenworthington/3178163 
ipak <- function(pkg){
  
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  
  sapply(pkg, require, character.only = TRUE)
  
}


## Run the function on the list.....
ipak(p)





#########################################################################################################################
## INSTALL PACKAGES NOT ON CRAN
#########################################################################################################################


## Also to create pdf/html documents, you need to install a latex program. EG MiKtex for windows.
## easiest to use the installr function, see:
## https://stackoverflow.com/questions/24239420/tex-package-not-installing-in-r-version-3-1-0
## installr::installr() and pick MikTeX (at least).
## install.packages("climates",,"http://rforge.net/",type="source")


## To rsemove all the packages
## rownames(installed.packages())
## remove.packages()


#########################################################################################################################
## Install github packages that can't be downloaded from CRAN
devtools::install_github("ropensci/rgbif")
library(devtools)
install_github('johnbaums/rmaxent')
install_github('johnbaums/jagstools')
install_github('johnbaums/hues')
install_github('johnbaums/trees')
install_github('johnbaums/things')
install_github("ropenscilabs/datastorr")
install_github("wcornwell/taxonlookup")
install_github('KarelMokany/AdaptR')
install_github("danlwarren/ENMTools")
devtools::install_github("richfitz/datastorr")
devtools::install_github("traitecoevo/baad.data")


#########################################################################################################################
## To run maxent in R, you need Java installed on your computer, and the main version of Java needs to match the R Version
## See these instructions:
## https://www.rdocumentation.org/packages/dismo/versions/1.1-4/topics/maxent
## check that R maxent is installed
# list.files(system.file("java", package = "dismo"))
# ENMevaluate



## Check the settings again
## https://www.r-statistics.com/2013/03/updating-r-from-r-on-windows-using-the-installr-package/





#########################################################################################################################
###################################### INSTALL MULTIPLE PACKAGES AT ONCE  ############################################### 
#########################################################################################################################