#########################################################################################################################
################################################ APC NATIVE/NATURALISED STATUS  #########################################
#########################################################################################################################


# The purpose of this code is to use the Australian Plant Census taxon distribution data to compute for each taxon the 
# native and/or naturalised status in each Australian state & territory.The data source is the Australian Plant Census, 
# which can be freely downloaded (currently ~100,000 records). Relevant fields are canonicalName and taxonDistribution. 
# The latter contains distribution data indicating native range on a state-by-state basis 

  
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

# In order to retrieve, for example, all taxa native in NSW, you could do something like:
   
# x <- data[data$regionName == "nsw" & data$native, ]
# x <- data %>% filter(regionName == "nsw", native)





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
intersect(HIA.SPP$Binomial, NAT.SPP)  ## 380
setdiff(HIA.SPP$Binomial, NAT.SPP)

intersect(test.spp, NAT.SPP)          ## 84
setdiff(test.spp, NAT.SPP)            ## 21


#########################################################################################################################
## Try filtering the records: is there a way of using these filters to generate the native lists, rather than per row?
APC.NAT.DIST = APC.NAT[!duplicated(APC.NAT$canonicalName), ][, c("canonicalName", "taxonDistribution")]





#########################################################################################################################
## 2). JOIN STATE DATA TO RECORDS
#########################################################################################################################


#########################################################################################################################
## First, restrict the dataset to just the target species
APC.RECORDS = COMBO.RASTER.CONTEXT[COMBO.RASTER.CONTEXT$searchTaxon %in% intersect(test.spp, NAT.SPP), ]
APC.RECORDS = APC.RECORDS[, c("searchTaxon", "lon", "lat")]
unique(APC.RECORDS$searchTaxon)


## We want to know the count of species that occur in n LGAs, across a range of climates. Read in LGA and SUA
APC.RECORDS = SpatialPointsDataFrame(coords = APC.RECORDS[c("lon", "lat")], 
                                     data   = APC.RECORDS,
                                     proj4string = CRS("+init=epsg:4326"))


## Read in the state boundaries
AUS.STATE = readOGR("./data/base/CONTEXTUAL/STE11aAust.shp", layer = "STE11aAust")
CRS.new   = CRS("+init=epsg:4326") # EPSG:3577
AUS.WGS   = spTransform(AUS.STATE, CRS.new)
identical(projection(AUS.WGS), projection(APC.RECORDS))
AUS.WGS$STATE_NAME


#########################################################################################################################
## This should be a simple intersect
APC.JOIN = over(APC.RECORDS, AUS.WGS)[c("STATE_NAME")]
head(APC.JOIN, 30)


## Now join this data back on to the original
APC.GEO  = cbind.data.frame(APC.RECORDS, APC.JOIN) 
drops  <- c("lon.1",  "lat.1",  "optional")
APC.GEO = APC.GEO[ , !(names(APC.GEO) %in% drops)]
head(APC.GEO, 30)





#########################################################################################################################
## 3). CHANGE STU's DATA to SUIT
#########################################################################################################################


## Just using the taxondistribution column
APC.NAT.DIST = dplyr::rename(APC.NAT.DIST, searchTaxon = canonicalName)
names(APC.NAT.DIST)
names(APC.RECORDS)
dim(APC.NAT.DIST)
dim(APC.RECORDS)



#########################################################################################################################
## Join the data
unique(APC.RECORDS$searchTaxon)
unique(APC.NAT.DIST$searchTaxon)

GBIF.APC = merge(APC.GEO, APC.NAT.DIST, by = "searchTaxon", all = FALSE)
dim(GBIF.APC)
head(GBIF.APC)

unique(GBIF.APC$taxonDistribution)
unique(GBIF.APC$STATE_NAME)


## Rename state factors to match Stu's names
GBIF.APC$STATE_NAME = revalue(TEST$STATE_NAME, c("Australian Capital Territory" = "act", 
                                                 "New South Wales"              = "nsw",
                                                 "Northern Territory"           = "nt",
                                                 "Queensland"                   = "qld",
                                                 "South Australia"              = "sa",
                                                 "Tasmania"                     = "tas",
                                                 "Victoria"                     = "vic",
                                                 "Western Australia"            = "wa"))


## How could these be matched to the taxondistribution column?
unique(GBIF.APC$STATE_NAME)
head(GBIF.APC$STATE_NAME, 20)


## We just want one row for each record - taxonDistribution - to link to the state occurrence.
## EG If STATE_NAME = NSW, and the plant is naturalised in NSW, another column would say "naturalised".
## So we would need to scrape the taxonDistribution column for which state's in is naturlaised in.

## This might need another column with "naturalised" or the like
View(GBIF.APC)
unique(GBIF.APC$taxonDistribution)

## Some kind of if/else statement? 
## frame$twohouses <- ifelse(frame$data>=2, 2, 1)
#GBIF.APC$native_region      = 
#GBIF.APC$naturalised_region =



#########################################################################################################################
## Write to file
write.csv(GBIF.APC, "./data/base/TRAITS/GBIF_APC_NATIVE_RANGE.csv", row.names = FALSE)


#########################################################################################################################
#########################################################  TBC ########################################################## 
#########################################################################################################################