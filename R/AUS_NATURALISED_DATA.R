#########################################################################################################################
################################################ APC NATIVE/NATURALISED STATUS  #########################################
#########################################################################################################################


# The purpose of this code is to use the Australian Plant Census taxon distribution data to compute for each record the 
# native and/or naturalised status, based on the Australian state & territory data.

# The data source is the Australian Plant Census, which can be freely downloaded (currently ~100,000 records). Relevant 
# fields are canonicalName and taxonDistribution. The latter contains distribution data indicating native range on 
# a state-by-state basis 

  
# NSW, ACT, Vic, Tas, SA, Vic (naturalised) WA, NT, Qld, NSW, Vic (extinct, see https://github.com/snubian/native_range)

# and sometimes with additional qualifying text in parentheses to indicate naturalised or extinct status (among others).
# These are occasionally ambiguous or messy, using words such as doubtfully, sparingly and possibly, or question marks:
# Note that numerous offshore territories are included in distribution data, referred to here generally as regions:


# Concept

# The idea is to split each taxon's taxonDistribution text into separate regions, and then set TRUE/FALSE flags to indicate 
# the distribution status based on the qualifying text (in parentheses):

# native
# naturalised
# native and naturalised (i.e. both native and naturalised populations exist)
# extinct
# unresolved

# A taxon is assumed to be native unless there is some indication that it is naturalised. Because of uncertainty/ambiguity 
# there are also flags set for potentially naturalised, potentially extinct and potentially native and naturalised. The 
# unresolved flag is set to TRUE in cases where status could not be determined.

# Each record in the output represents a combination of taxon and region. If a taxon originally had distribution data for 
# three regions (e.g. "qld, nsw, vic"), then there will be three records for that taxon in the output, one for each region. 
# The distribution flags apply to the given region.





#########################################################################################################################
## 1). Read in naturalised data
#########################################################################################################################


## Read naturalised list in
load("./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT.RData")
load("./data/base/HIA_LIST/COMBO/COMBO_RASTER_CONTEXT.RData")
source('./R/HIA_LIST_MATCHING.R')
source('./R/HIA_CLEAN_MATCHING.R')


##
APC.NAT   = read.csv("./data/base/TRAITS/apc_taxon_distribution.csv", stringsAsFactors = FALSE)
View(APC.NAT)
unique(APC.NAT$taxonDistribution)


## how many species?
NAT.SPP = unique(APC.NAT$canonicalName)
str(NAT.SPP)


## How many species are on the target list?
## The trial species
test.spp = sort(unique(c(renee.full$Species, 
                         "Betula pendula", "Fraxinus excelsior", "Quercus robur", "Fagus sylvatica",
                         Manuel.test)))


## Combine with 45 from the main list
HIA.SAMPLE = head(COMBO.NICHE.CONTEXT, 53)[, c("searchTaxon")]
test.spp   = sort(unique(c(test.spp, HIA.SAMPLE)))


#########################################################################################################################
## How many species are already on our list?
intersect(HIA.SPP$Binomial, NAT.SPP)  ## 380 native species
setdiff(HIA.SPP$Binomial, NAT.SPP)

intersect(test.spp, NAT.SPP)          ## 84
setdiff(test.spp, NAT.SPP)            ## 21


#########################################################################################################################
## Try filtering the records: is there a way of using these filters to generate the native lists, rather than per row?
# Don't use the `taxonDistribution` column - this is the raw and often messy data that was 
# originally in the APC dataset. Most of what I tried to do was to parse that field and put the information contained 
# in there into some sort of standard structure. I only put that column in the final output for reference actually, 
# so one could see what the raw data was.
APC.NAT.DIST = APC.NAT[!duplicated(APC.NAT$canonicalName), ][, c("canonicalName", "regionName", "native", "naturalised")]
names(APC.NAT.DIST)





#########################################################################################################################
## 2). JOIN STATE DATA TO RECORDS
#########################################################################################################################


#########################################################################################################################
## First, restrict the dataset to just the target species
APC.RECORDS = COMBO.RASTER.CONTEXT[COMBO.RASTER.CONTEXT$searchTaxon %in% intersect(test.spp, NAT.SPP), ]
APC.RECORDS = APC.RECORDS[, c("searchTaxon", "lon", "lat")]
unique(APC.RECORDS$searchTaxon)


## Convert the records to a spatial points data frame
APC.RECORDS = SpatialPointsDataFrame(coords = APC.RECORDS[c("lon", "lat")], 
                                     data   = APC.RECORDS,
                                     proj4string = CRS("+init=epsg:4326"))


## Read in the state boundaries and the world
LAND      = readOGR("./data/base/CONTEXTUAL/ne_10m_land.shp",            layer = "ne_10m_land")
AUS.STATE = readOGR("./data/base/CONTEXTUAL/STE11aAust_BUFFER_1KM.shp",  layer = "STE11aAust_BUFFER_1KM")
CRS.new   = CRS("+init=epsg:4326") # EPSG:3577
AUS.WGS   = spTransform(AUS.STATE, CRS.new)
identical(projection(AUS.WGS), projection(APC.RECORDS))
AUS.WGS$STATE_NAME


#########################################################################################################################
## Run the intersect between the APC records and the AUS states. The states and records line up
APC.JOIN = over(APC.RECORDS, AUS.WGS)[c("STATE_NAME")]
dim(APC.JOIN)
head(APC.JOIN, 30)


## Now join this data back on to the original
APC.GEO  = cbind.data.frame(APC.RECORDS, APC.JOIN) 
drops  <- c("lon.1",  "lat.1",  "optional")
APC.GEO = APC.GEO[ , !(names(APC.GEO) %in% drops)]
head(APC.GEO, 30)


## Save the intersection which took ages to create...
save(APC.GEO, file = paste("./data/base/TRAITS/APC_GEO.RData", sep = ""))
load("./data/base/TRAITS/APC_GEO.RData")





#########################################################################################################################
## 3). MERGE WITH STU / RACH'S NATURALISED DATA
#########################################################################################################################


## Just using the taxondistribution column
APC.NAT.DIST = dplyr::rename(APC.NAT.DIST, searchTaxon = canonicalName)
names(APC.NAT.DIST);names(APC.RECORDS)
dim(APC.NAT.DIST)[1];dim(APC.RECORDS)[1]


#########################################################################################################################
## Join the data
## So the join code seems to be working. However, 87% of the records are NA. Not sure why this is occurring...


## How could these be matched to the taxondistribution column?
class(APC.GEO$STATE_NAME);class(APC.NAT.DIST$regionName)
APC.NAT.DIST$regionName = as.factor(APC.NAT.DIST$regionName)

str(APC.GEO$STATE_NAME);str(APC.NAT.DIST$regionName)
unique(APC.GEO$STATE_NAME);unique(APC.NAT.DIST$regionName)


## Rename state factors to match Stu's names
APC.GEO$STATE_NAME = revalue(APC.GEO$STATE_NAME, c("Australian Capital Territory" = "act", 
                                                   "New South Wales"              = "nsw",
                                                   "Northern Territory"           = "nt",
                                                   "Queensland"                   = "qld",
                                                   "South Australia"              = "sa",
                                                   "Tasmania"                     = "tas",
                                                   "Victoria"                     = "vic",
                                                   "Western Australia"            = "wa"))

## Don't I need to make sure that all the state categories are the same in both data sets?
APC.NAT.DIST$regionName = revalue(APC.NAT.DIST$regionName, c("?nsw" = "nsw",
                                                             "?nt"  = "nt",
                                                             "?qld" = "qld",
                                                             "?sa"  = "sa",
                                                             "?tas" = "tas",
                                                             "?vic" = "vic",
                                                             "?wa"  = "wa",
                                                             "chi"  = "Other Territories",
                                                             "ni"   = "Other Territories",
                                                             "csi"  = "Other Territories",
                                                             "hi"   = "Other Territories",
                                                             "lhi"  = "Other Territories",
                                                             "coi"  = "Other Territories",
                                                             "mi"   = "Other Territories",
                                                             "ar"   = "Other Territories",
                                                             "?lhi" = "Other Territories",
                                                             "?chi" = "Other Territories"))


## Check the species
str(unique(APC.RECORDS$searchTaxon));str(unique(APC.NAT.DIST$searchTaxon)) ## How many species in each?
APC.GEO$STATE_NAME       = as.character(APC.GEO$STATE_NAME)
APC.NAT.DIST$regionName  = as.character(APC.NAT.DIST$regionName)


## What is the difference between these columns...
setdiff(unique(APC.NAT.DIST$regionName), unique(APC.GEO$STATE_NAME))
intersect(unique(APC.NAT.DIST$regionName), unique(APC.GEO$STATE_NAME))
setdiff(unique(APC.GEO$STATE_NAME), unique(APC.NAT.DIST$regionName))


#########################################################################################################################
## If you have dataset A with `searchTaxon`, `lon`, `lat` and `STATE_NAME`, and dataset B being `canonicalName`, 
## `regionName` and `native` etc, then merge A and B by a combination of taxon name and state. 
## This will in effect append the native/naturalised/etc flags to your lonlat data, so you should can then determine 
## whether a particular occurrence is in the native range or not. I hope this makes sense!
names(APC.GEO);names(APC.NAT.DIST)
str(APC.GEO);str(APC.NAT.DIST)


## Why does this code create so many NA records? Is it to do with the column names?
## When common ids have different names, use by.x and by.y to match them. R will keep the name of the first dataset (by.x) 
GBIF.APC <- merge(APC.GEO, APC.NAT.DIST, by.x = c("searchTaxon", "STATE_NAME"), by.y = c("searchTaxon", "regionName"), all.x = TRUE)
dim(GBIF.APC);dim(APC.GEO);dim(APC.NAT.DIST)  ## dimensions are ok?


## Check the join fields: not sure why "regionName" disappears?
head(GBIF.APC)
setdiff(unique(GBIF.APC$STATE_NAME), unique(GBIF.APC$STATE_NAME))


#########################################################################################################################
## Consider the exceptions: 
## Records outside AUS will have "NA" for naturalised, but which problems are to do with the merge, column naming?

## Other Territories	150.6586	-35.1311	NA	NA
## qld	142.3	-10.6	NA	NA
## nsw	147.927347	-30.14320403	TRUE	FALSE
## nt	132.816667	-23.71666702	NA	NA
## sa	132.221771	-30.912699	NA	NA
GBIF.APC.NA = GBIF.APC[rowSums(is.na(GBIF.APC)) > 0,]
dim(GBIF.APC[rowSums(is.na(GBIF.APC)) > 0,])[1] /dim(GBIF.APC)[1] ## 87% of records are NA naturalised?


## Get the non-NA rows
GBIF.APC.VALID       = subset(GBIF.APC, (!is.na(GBIF.APC$native))) 
GBIF.APC.STATE.NA    = subset(GBIF.APC, (is.na(GBIF.APC$STATE_NAME))) 
GBIF.APC.STATE.VALID = subset(GBIF.APC, (!is.na(GBIF.APC$STATE_NAME))) 


## What are the percentages?
dim(GBIF.APC.NA)[1]/dim(GBIF.APC)[1]
dim(GBIF.APC.VALID)[1]/dim(GBIF.APC)[1]
dim(GBIF.APC.STATE.NA)[1]/dim(GBIF.APC)[1]
dim(GBIF.APC.STATE.VALID)[1]/dim(GBIF.APC)[1]


## Plot them
plot(AUS.STATE)
points(GBIF.APC.NA[c("lon", "lat")],  pch = ".", col = "red", cex = 0.5)
points(GBIF.APC.VALID [c("lon", "lat")],  pch = ".", col = "blue")


## Could there be a mis-match between the state where the record was located, and where it is supposed to be naturalised?
unique(GBIF.APC.NA$STATE_NAME)


## Why are there so many NA records, and why are they concentrated in NSW?
head(GBIF.APC.NA)
str(unique(GBIF.APC.NA$searchTaxon));str(unique(GBIF.APC$searchTaxon))





#########################################################################################################################
## 4). CREATE A FIELD FOR OUTSIDE AUS
#########################################################################################################################


## Not sure if this will do the job...
GBIF.APC$OUT_SIDE_AUS = ifelse(is.na(GBIF.APC$STATE_NAME), "TRUE", "FALSE")
unique(GBIF.APC$OUT_SIDE_AUS)


#########################################################################################################################
## Write to file
#save(GBIF.APC, file = paste("./data/base/TRAITS/GBIF_APC.RData", sep = ""))
save.image("GBIF_APC_NATIVE_RANGE.RData")
#load("./data/base/TRAITS/GBIF_APC_NATIVE_RANGE.RData")
write.csv(GBIF.APC, "./data/base/TRAITS/GBIF_APC_NATIVE_RANGE.csv", row.names = FALSE)






#########################################################################################################################
## OUTSTANDING NATURALISED TASKS:
#########################################################################################################################


## Consider the exceptions...







#########################################################################################################################
#########################################################  TBC ########################################################## 
#########################################################################################################################