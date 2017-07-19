#########################################################################################################################
###################################### INSTALL MULTIPLE PACKAGES AT ONCE  ############################################### 
#########################################################################################################################


#########################################################################################################################
## create a big list of all the packages needed for a project
## can't do the special github ones like this though...
packages <- c("ggplot2",   "plyr",   "reshape2", "RColorBrewer",  "scales",     "grid",
              "raster",    "spdep",  "rgdal",    "GISTools",      "data.table", "dtplyr",
              "Hmisc",     "Cairo",  "lattice",  "gtools",        "ggplot2",    "sjPlot",
              "gridExtra", "grid",   "sqldf",    "knitr",         "mgcv",       "MuMIn",
              "mvabund",   "asbio",  "mvtnorm",  "testthat",      "rgl",        "R.matlab",
              "ENMeval",   "lavaan", "semPlot",  "rgdal",         "sp",         "formula.tools",
              "dismo",     "mctest", "rJava",    "ENMeval",       "SDMTools",   "ALA4R",
              "statisticalModeling", 
              "knitr", "yaml", "htmltools", "caTools", "bitops", "rmarkdown") # class(packages)


## also to create pdf/html documents, you need to install a latex program. EG MiKtex for windows.
## easiest to use the installr function, see:
## https://stackoverflow.com/questions/24239420/tex-package-not-installing-in-r-version-3-1-0
## installr::installr() and pick MikTeX (at least).
## install.packages("climates",,"http://rforge.net/",type="source")


#########################################################################################################################
## now use a function from the internet that takes a list of packages, installs and loads them
ipak <- function(pkg){
  
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  
  sapply(pkg, require, character.only = TRUE)
  
}


## run the function on the list.....
ipak(packages)


## check that R maxent is installed
list.files(system.file("java", package = "dismo"))
ENMevaluate





#########################################################################################################################
###################################### INSTALL MULTIPLE PACKAGES AT ONCE  ############################################### 
#########################################################################################################################