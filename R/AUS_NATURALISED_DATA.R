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

intersect(test.spp, NAT.SPP)  ## 84
setdiff(test.spp, NAT.SPP)    ## 21





#########################################################################################################################
## 2). JOIN STATE DATA TO RECORDS
#########################################################################################################################


## First, restrict the dataset to just the target species
APC.RECORDS = COMBO.RASTER.CONTEXT[COMBO.RASTER.CONTEXT$searchTaxon %in% intersect(test.spp, NAT.SPP), ]
APC.RECORDS = head(RECORDS.APC, 500)[, c("searchTaxon", "lon", "lat")]
unique(APC.TEST$searchTaxon)


## We want to know the count of species that occur in n LGAs, across a range of climates. Read in LGA and SUA
APC.TEST.SP   = SpatialPointsDataFrame(coords = APC.RECORDS[c("lon", "lat")], 
                                       data   = APC.RECORDS,
                                       proj4string = CRS("+init=epsg:4326"))


## Read in the state boundaries
AUS.STATE = readOGR("./data/base/CONTEXTUAL/STE11aAust.shp", layer = "STE11aAust")


## Project
CRS.new  <- CRS("+init=epsg:4326") # EPSG:3577
AUS.WGS  = spTransform(AUS.STATE, CRS.new)


projection(AUS.WGS)
projection(APC.TEST.SP)


## Now try intersecting the recrods with the states
names(AUS.WGS)
AUS.WGS$STATE_NAME


## This should be a simple intersect
APC.JOIN   = over(APC.TEST.SP, AUS.WGS)   
head(APC.JOIN, 30)


## Now join this data back on to the original
APC.RECORDS  = cbind.data.frame(APC.RECORDS, APC.JOIN) 
unique(APC.RECORDS$STATE_NAME)





#########################################################################################################################
## 3). CHANGE STU's DATA to SUIT
#########################################################################################################################


## Which columns to use?
## We just want one row for each record - taxonDistribution - to say link to the state occurrence.
## If location state = NSW and the plant is naturalised in NSW, then this record would have "naturalised" associated with the
## lat/long. So we would need to scrape the taxonDistribution column for which state it's naturlaised in.

## This might need another column with "naturalised" or the like

names(APC.NAT)
unique(APC.NAT$regionName)


## Subset and rename
APC.NAT = APC.NAT[, c("canonicalName", "taxonDistribution", "regionName",
                       "native", "naturalised")]

APC.NAT = dplyr::rename(APC.NAT, 
                        native_region      = native,
                        naturalised_region = naturalised,
                        searchTaxon        = canonicalName)

## Check the names..
names(APC.NAT)
names(APC.RECORDS)
dim(APC.NAT)
dim(APC.RECORDS)
unique(APC.NAT$taxonDistribution)


#########################################################################################################################
## Join the data
unique(APC.RECORDS$searchTaxon)
unique(APC.NAT$searchTaxon)



TEST = merge(APC.RECORDS, APC.NAT, by = "searchTaxon", all = FALSE)
dim(TEST)
head(TEST)


unique(TEST$region)
unique(TEST$STATE_NAME)

## rename factor?
TEST$STATE_NAME = revalue(TEST$STATE_NAME, c("Australian Capital Territory" = "act", 
                                             "New South Wales"              = "nsw",
                                             "Northern Territory"           = "nt",
                                             "Queensland"                   = "qld",
                                             "South Australia"              = "sa",
                                             "Tasmania"                     = "tas",
                                             "Victoria"                     = "vic",
                                             "Western Australia"            = "wa"))
unique(TEST$STATE_NAME)
head(TEST$STATE_NAME, 20)






