#########################################################################################################################
########################################### SUMMARISE GBIF DATARECORDS ################################################## 
#########################################################################################################################


#########################################################################################################################
## 1). SUMMARISE SPECIES NUMERICALLY
#########################################################################################################################


#########################################################################################################################
## Now create a master summary of the recrods so far. What do we need to know?
## These numbers could be variables that update each time code is re-run with changes to the variables
## Also a table could be made for each source (COMBO, ALA, Council, etc.), But ideally it is just one table
load("./data/base/HIA_LIST/COMBO/COMBO_RASTER_CONTEXT.RData")   ## All the environmental data, one row for each record
load("./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT.RData")    ## The niches for each variable, one row for each species
load("./data/base/HIA_LIST/GBIF/GBIF_TRIM.RData")


## The table above gives the details, but worth documenting how many records are knocked out by each
Total.records        = dim(GBIF.TRIM)[1]
Total.taxa.processed = dim(COMBO.NICHE.CONTEXT)[1]
Remaining.records    = dim(COMBO.RASTER.CONTEXT)[1] 
Remaining.percent    = Remaining.records/Total.records * 100
Filters.applied      = "NA COORD | MANAGED/NA | < 1950/NA"

## Total number of number of taxa (ie. all variables are dynamic)
## Total number of records (uncleaned, and broken down by source)
## % records retained after filtering 
## 6-figure summary across all taxa (min, max, median)
## number of all taxa with > n records (e.g. +50, +100, +1000, +10,000) - assuming it doesn't matter above certain level
## number of top 200 taxa with > n records
summary(COMBO.NICHE.CONTEXT$COMBO.count)


## How many species where knocked out by using filters?
load("./data/base/HIA_LIST/GBIF/GBIF_PROBLEMS.RData")
kable(GBIF.PROBLEMS)   ## doesn't exist


#########################################################################################################################
## Which species have more than n records, are on the top 200, etc?
## first create column headings for final results table  
## these match data frame created by bivariate lm function
n.samples = 1
table.columns = c(  
  "Dataset",
  "Taxa processed",
  "Rank",
  "Total records",
  "Filters applied",
  "Cleaned records",
  "% Records retained",
  "Min",  
  "Max",
  "Median",
  "All taxa count 50+",
  "All taxa count 100+",
  "All taxa count 1000+",
  "All taxa count 10,000+",
  "Top taxa count 50+",
  "Top taxa count 100+",
  "Top taxa count 1000+",
  "Top taxa count 10,000+"
)


## Create big matrices with n rows = n samples
table.length                    = length(table.columns)
COMBO.RECORD.SUMMARY            = matrix(0, n.samples, table.length)
colnames(COMBO.RECORD.SUMMARY)  = table.columns 
COMBO.RECORD.SUMMARY            = as.data.frame(COMBO.RECORD.SUMMARY)


#########################################################################################################################
## Now put the data in. Could use "melt", etc?
COMBO.RECORD.SUMMARY[1, "Dataset"]                 = "COMBO_occurrence"                          ## H1.ms[1, "Intercept"]
COMBO.RECORD.SUMMARY[1, "Taxa processed"]          = Total.taxa.processed
COMBO.RECORD.SUMMARY[1, "Rank"]                    = "Species"
COMBO.RECORD.SUMMARY[1, "Total records"]           = Total.records
COMBO.RECORD.SUMMARY[1, "Filters applied"]         = Filters.applied
COMBO.RECORD.SUMMARY[1, "Cleaned records"]         = Remaining.records
COMBO.RECORD.SUMMARY[1, "% Records retained"]      = Remaining.percent

COMBO.RECORD.SUMMARY[1, "Min"]                     = summary(COMBO.NICHE.CONTEXT$COMBO.count)[1]
COMBO.RECORD.SUMMARY[1, "Max"]                     = summary(COMBO.NICHE.CONTEXT$COMBO.count)[6]
COMBO.RECORD.SUMMARY[1, "Median"]                  = summary(COMBO.NICHE.CONTEXT$COMBO.count)[4]

COMBO.RECORD.SUMMARY[1, "All taxa count 50+"]      = sum(COMBO.NICHE.CONTEXT$COMBO.count >= 50,   na.rm = TRUE)
COMBO.RECORD.SUMMARY[1, "All taxa count 100+"]     = sum(COMBO.NICHE.CONTEXT$COMBO.count >= 100,  na.rm = TRUE)
COMBO.RECORD.SUMMARY[1, "All taxa count 1000+"]    = sum(COMBO.NICHE.CONTEXT$COMBO.count >= 1000, na.rm = TRUE)
COMBO.RECORD.SUMMARY[1, "All taxa count 10,000+"]  = sum(COMBO.NICHE.CONTEXT$COMBO.count >= 10000, na.rm = TRUE)

COMBO.RECORD.SUMMARY[1, "Top taxa count 50+"]      = sum(COMBO.NICHE.CONTEXT$COMBO.count >= 50    & COMBO.NICHE.CONTEXT$Top_200 == "TRUE", na.rm = TRUE)
COMBO.RECORD.SUMMARY[1, "Top taxa count 100+"]     = sum(COMBO.NICHE.CONTEXT$COMBO.count >= 100   & COMBO.NICHE.CONTEXT$Top_200 == "TRUE", na.rm = TRUE)
COMBO.RECORD.SUMMARY[1, "Top taxa count 1000+"]    = sum(COMBO.NICHE.CONTEXT$COMBO.count >= 1000  & COMBO.NICHE.CONTEXT$Top_200 == "TRUE", na.rm = TRUE)
COMBO.RECORD.SUMMARY[1, "Top taxa count 10,000+"]  = sum(COMBO.NICHE.CONTEXT$COMBO.count >= 10000 & COMBO.NICHE.CONTEXT$Top_200 == "TRUE", na.rm = TRUE)


#########################################################################################################################
## Visualise this
str(COMBO.RECORD.SUMMARY)
colnames(COMBO.RECORD.SUMMARY)
rownames(COMBO.RECORD.SUMMARY) = "VALUE"

COMBO.RECORD.SUMMARY = as.data.frame(t(COMBO.RECORD.SUMMARY)) ## do melting, etc later
kable(COMBO.RECORD.SUMMARY)

## Save data:
save(COMBO.RECORD.SUMMARY,  file = paste("./data/base/HIA_LIST/COMBO/COMBO_RECORD_SUMMARY.RData",  sep = ""))


#########################################################################################################################
## Now visualise the distribution of records
histogram(COMBO.NICHE.CONTEXT$COMBO.count,
          breaks = 50, border = NA, col = "grey",
          xlab = "Number of COMBO records per species", 
          main = "Distribution of COMBO records for HIA species")





#########################################################################################################################
## 2). QUERY SPECIES TO FIND NEW ONES: : WHICH COULD BE MODELLED?
#########################################################################################################################


## Find infrequently sold spp, big environmental & geographic range, but could have similar traits to popular species
summary(COMBO.NICHE.CONTEXT$AREA_OCCUPANCY)
summary(COMBO.NICHE.CONTEXT$Annual_mean_temp_range)
summary(COMBO.NICHE.CONTEXT$Number.of.growers)
summary(COMBO.NICHE.CONTEXT$COMBO.count)


## Rare species we can't model?
rare.spp = subset(COMBO.NICHE.CONTEXT, COMBO.count < 100)[ c("searchTaxon", "Top_200", "Origin")]


## Potential new species
new.spp = subset(COMBO.NICHE.CONTEXT, AREA_OCCUPANCY > 7000 & 
                   Annual_mean_temp_range > 18 & 
                   Number.of.growers < 25)[ c("searchTaxon", "Top_200", "Origin")]


## Other queries?





#########################################################################################################################
## OUTSTANDING NICHE TASKS:
#########################################################################################################################


## Check on species which seem to have been knocked out: E.G Fagus sylvatica. Individual filter doesn't knock them all out  

## Convert WORLDCLIM values back into decimals            - multiply by 10

## Check geographic range: doesn't look right for some species. Calc extent of occurrnece as well 

## Return species EG:                                     -

## Find infrequently sold spp., big environmental & geographic range, but could have similar traits to popular species

## Find rarest species (are there popular species with not many records?)

## Ask Linda and others about the kind of queries they want to run...


#########################################################################################################################
#####################################################  END ############################################################## 
#########################################################################################################################