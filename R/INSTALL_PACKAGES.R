#########################################################################################################################
###################################### INSTALL MULTIPLE PACKAGES AT ONCE  ############################################### 
#########################################################################################################################


## create a mega list of all the packages you could ever need?
## can't do the special github ones like this though...
packages <- c("ggplot2",   "plyr",   "reshape2", "RColorBrewer",  "scales",     "grid",
              "raster",    "spdep",  "rgdal",    "GISTools",      "data.table", "dtplyr",
              "Hmisc",     "Cairo",  "lattice",  "gtools",        "ggplot2",    "sjPlot",
              "gridExtra", "grid",   "sqldf",    "knitr",         "mgcv",       "MuMIn",
              "mvabund",   "asbio",  "mvtnorm",  "testthat",      "rgl",        "R.matlab",
              "ENMeval",   "lavaan", "semPlot",  "rgdal",         "sp",         "formula.tools",
              "dismo",     "mctest", "rJava",    "ENMeval",       "SDMTools",   "ALA4R",
              "statisticalModeling") # class(packages)


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


## now commit these changes...why can't I commit them?



