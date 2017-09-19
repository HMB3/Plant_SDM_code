#########################################################################################################################
#######################################  COMBINE ALL SPECIES DATA FRAMES INTO ONE ####################################### 
#########################################################################################################################


# John: In case my comment doesn't come through: Occurrence records for all vascular plants (native and non-native) in Australia,
# to the extent that they exist in OEH (Bionet) and the AVH hub of ALA.
# 
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
# 
# 
# I believe the OEH data are licensed and use would not be permitted beyond the refugia project. I suggest we recreate these data
# using a new bionet extract (will need to contact OEH to ask for it) and redownload the ALA data as well. Alternatively, just use
# ALA alone, but for our NSW work we complemented it with OEH data because it is more complete for NSW, and this better matched
# our focal species' data.



#########################################################################################################################
## 1). READ IN AVH
#########################################################################################################################


## Read in ALA data: need readRDS
AVH.OEH.VASC = readRDS("./data/base/HIA_LIST/ALA/SPECIES/background_all_plants_oeh_avh_with_ibra_dodgy_removed.rds", refhook = NULL)
write.csv(AVH.OEH.VASC, "./data/base/HIA_LIST/ALA/SPECIES/AVH_OEH_VASC.csv", row.names = FALSE)


## What is this file?
str(AVH.OEH.VASC)
names(AVH.OEH.VASC)
head(AVH.OEH.VASC)


## Now when do we join on the others? And do we blend the AVH and GBIF? Is there a simple way of combining the two sources? 
## Consider that GBIF has data for both sources. So are we topping up the native ranges with the AVH. So it will be important 
## to get rid of the duplicates. 
names(GBIF.RASTER.CONTEXT)








