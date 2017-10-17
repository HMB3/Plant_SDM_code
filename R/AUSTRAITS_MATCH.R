#########################################################################################################################
############################################# READ IN AUSTRAITS DATA #################################################### 
#########################################################################################################################


#########################################################################################################################
## In this file, we are mathing the AUS TRAITS data to the HIA species niches from GBIF  


## Consider the HIA brief again:

# The first module will focus on fifty plant species identified in the project’s Target Species List, and will develop maps 
# that demonstrate each species’ suitability to both current and future climates across Australia.
# 
# These maps will be used to demonstrate how well or poorly a particular species will be able to tolerate future conditions 
# in urban centres across Australia as the climate changes, based on our current understanding of species’ climatic 
# requirements.

# Our research might demonstrate, for example, that a particular species of tree is already at the very limit of its 
# ability to cope with heat, and that the only suitable place to plant this species in the future will be in cool-climate 
# or more temperate locations. This kind of information would be very useful to a council seeking to avoid investing in 
# tree species for street planting that are unlikely to cope with higher temperatures.

# We will also use information from national herbaria and other sources to quantify each species’ climatic limits - the 
# warmest, coldest, driest or wettest conditions they can cope with. This information will then be tested through the 
# Planting Successes and Failures module of the research programme to ensure that the Interactive Plant Features Tool 
# matches the right plant in the right region with an eye on the future.

# We will also be working with growers, nurseries, landscape architects and many others to capture their recordings of 
# major plant traits including:


# Growth rate and form
# height
# canopy density
# Ground cover
# Longevity
# Seasonality
# Water quality
# Allergenicity
# Air and water quality influences and urban temperatures
# Insect resistance
# Ornamental and amenity features
# and biodiversity impacts.


## NEED TO DECIDE WHICH TRAITS TO USE...................................................................................



#########################################################################################################################
## 1). READ IN AUSTRAITS DATA
#########################################################################################################################


#########################################################################################################################
# Austraits has at least one trait measured for 321 of these 'species' (90 different traits across all species (same 
# definitions as in TRY database) across roughly 17,000 observations (i.e. some species have many traits measured, sometimes 
# multiple measures for the one trait). To do the merge I used the genus and species name only - leaving out all cultivar, 
# var or subsp info. Please don't pass this data on to anyone else, including those in the project without checking with me 
# (co-authors need to give express permission to use it at this point; it will be open access next year)
HIA.AUST = read.csv("./data/base/TRAITS/HIA_austraits.csv", stringsAsFactors = FALSE)


## So, start by deciding which trait source is the most Authoritative (E.G. RSBG), then get unique trait values for that
# AUST.UNIQUE = HIA.AUST[!duplicated(HIA.AUST$trait_name), ]
# str(unique(AUST.UNIQUE$taxon))


## The problem here is that we can't get unique values for either taxon or trait, because this knocks out data. So, can
## we get unique traits within each species


## Just get the needed columns
AUST.LOOKUP <- 
  HIA.AUST %>% 
  select(taxon, study, trait_name, value, unit, value_type) %>%
  as.data.frame()

AUST.JOIN <- 
  HIA.AUST %>% 
  select(taxon, trait_name, value) %>%
  as.data.frame()

View(AUST.LOOKUP)


## 
test.spread = spread(AUST.JOIN, key = trait_name, value = value) ## Key = the column we want to spread (e.g. trait name), value = variable value   

##
library(reshape2)
dcast(melt(AUST.JOIN), taxon ~ Time + variable)


##
library(dplyr)
library(tidyr)
df %>%
  group_by(start_end, id) %>%
  mutate(ind = row_number()) %>%
  spread(start_end, date) %>% 
  select(start, end)


##
View(test.spread)
