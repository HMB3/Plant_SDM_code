#########################################################################################################################
##################################################  ALA DATA ############################################################ 
#########################################################################################################################


# John: In case my comment doesn't come through: Occurrence records for all vascular plants (native and non-native) in Australia,
# to the extent that they exist in OEH (Bionet) and the AVH hub of ALA.

 
# OEH data not otherwise cleaned other than removing records with coordinate uncertainty > 1000 m. ALA data have had records
# removed if they: were duplicates (based on ALA duplicate flags); had invalid geodetic datums; were flagged as cultivated/
# escapees (occCultivatedEscapee); were recorded earlier than 1950; had taxonomic identification issues (taxonIdentificationIssue)
# or were otherwise not "taxonomically kosher" (taxonomicKosher); had basis of obs other than PreservedSpecimen (basisOfRecord);
# had coordinate uncertainty > 1000 m (coordinateUncertaintyInMetres); were outliers for one or more environmental layers
# (outlierForLayer); were "detected outliers" (detectedOutlier); had unrecognised (nameNotRecognised) or invalid (invalidScientificName)
# names; or had geospatial_kosher=false. The source code for the construction of this background file is in
# c:/projects/maxent_in_R/download_all_avh.R (but I can't find the original bionet extract provided by OEH, nor the code used to
# read it in.. I suspect it might be the series of txt files in \\sci-6944.mqauth.uni.mq.edu.au\C\data\OEH\OEH_Atlas_Licensed\
# native_flora.7z).

 
# I believe the OEH data are licensed and use would not be permitted beyond the refugia project. I suggest we recreate these data
# using a new bionet extract (will need to contact OEH to ask for it) and redownload the ALA data as well. Alternatively, just use
# ALA alone, but for our NSW work we complemented it with OEH data because it is more complete for NSW, and this better matched
# our focal species' data.


# Also, please don't pass those data on any further. OEH recently contacted us about unlicensed use of their data. They're happy to 
# provide licenses, but don't want data floating around being used without their permission and attribution.
# If you're working with any threatened species (probably not?) then you may need to request as-held records for those species, since 
# they degrade the coordinates otherwise (as in the dataset I shared).



#########################################################################################################################
## 1). READ IN AVH AND MATCH TO HIA SPECIES LIST
#########################################################################################################################


## Read in ALA point data and shapefiles: need readRDS
AVH.OEH.VASC = readRDS("./data/base/HIA_LIST/ALA/SPECIES/background_all_plants_oeh_avh_with_ibra_dodgy_removed.rds", refhook = NULL)
SUA          = readOGR("./data/base/CONTEXTUAL/SUA_2011_AUST.shp", layer = "SUA_2011_AUST")
LGA          = readOGR("./data/base/CONTEXTUAL/LGA_2016_AUST.shp", layer = "LGA_2016_AUST")


## plot to check
plot(AVH.OEH.VASC)
plot(SUA)
plot(LGA)


## What is this file?
str(AVH.OEH.VASC)
names(AVH.OEH.VASC)
head(AVH.OEH.VASC)


## Project the ALA data into WGS84
projection(AVH.OEH.VASC)
projection(LGA)
projection(SUA)

CRS.new <- CRS("+init=epsg:4326") # EPSG:3577
AVH.WGS = spTransform(AVH.OEH.VASC, CRS.new)
LGA.WGS = spTransform(LGA, CRS.new)
SUA.WGS = spTransform(SUA, CRS.new)

projection(AVH.WGS)
projection(LGA.WGS)
projection(LGA.WGS)


## write AVH shapefile
writeOGR(AVH.OEH.VASC, ".", "AVH_OEH_VASC", driver = "ESRI Shapefile")


## Get the coordinates
AVH.XY          = as.data.frame(coordinates(AVH.WGS))
head(AVH.XY)
AVH.OEH.VASC.DF = as.data.frame(AVH.OEH.VASC)
AVH.OEH.VASC.XY = cbind(AVH.XY, AVH.OEH.VASC.DF)
AVH.OEH.VASC.XY[12] <- NULL
AVH.OEH.VASC.XY[13] <- NULL
AVH.OEH.VASC.XY[12] <- NULL


## check the output
names(AVH.OEH.VASC.XY)
str(AVH.OEH.VASC.XY)


#########################################################################################################################
## Now match the species list
HIA.FIN = unique(HIA.SPP$Binomial)
AVH.SPP = unique(AVH.OEH.VASC.XY$scientificname)


## Intersection of the HIA list and the AVH list
HIA.AVH.OVERLAP = intersect(HIA.FIN, AVH.SPP)
HIA.AVH.DIFF    = setdiff(HIA.FIN, AVH.SPP)


## Read in Kate's data
## set wd to where kates files are, eg LGA folder
file_names <- list.files(pattern = ".csv") ## "./data/base/HIA_LIST/LGA/") # where you have your files
your_data_frame <- do.call(bind_rows, lapply(file_names, read.csv))
## change WD back to the top...





#########################################################################################################################
## 2). INTERSECT ALA WITH LGA DATA
#########################################################################################################################


## point.in.poly Intersects point and polygon feature classes and adds polygon attributes to points
## Check columns
names(LGA.WGS)
names(SUA.WGS)
AVH.LGA <- point.in.poly(AVH.WGS, LGA.WGS)
AVH.SUA <- point.in.poly(AVH.WGS, SUA.WGS)


## save
save(AVH.LGA , file = paste("./data/base/HIA_LIST/ALA/SPECIES/AVH.LGA .RData", sep = ""))
save(AVH.SUA , file = paste("./data/base/HIA_LIST/ALA/SPECIES/AVH.SUA .RData", sep = ""))


#########################################################################################################################
## OUTSTANDING ALA TASKS:
#########################################################################################################################


## Get permission to use the ALA data.

## Get the file which Kate from WSU made: spatial join of ALA and LGAs?

## Filter the records that Paul needs...

## When do we join the ALA data to others? 


#########################################################################################################################
#################################################  END OF ALA CODE ###################################################### 
#########################################################################################################################