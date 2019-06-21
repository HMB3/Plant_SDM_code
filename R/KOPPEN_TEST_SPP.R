#########################################################################################################################
############################################# KOPPEN TEST SPP ########################################################### 
#########################################################################################################################


#########################################################################################################################
## This code combines the table of Koppen * SPP records, and uses Renee's species list to choose 10 species to trial
## The Koppen approach, outlined below :: 


## We could then group species according to their distribution in Koppen zones. For each of these major groupings, we 
## would then come up with a list of variables to use in the models, based on our understanding of the climate of those 
## zones. So, for example, we would be able to identify a group of species that are primarily Montane. We would then decide 
## on which variables to use for these. A second group of species might be primarily ari##semi-arid, and so on.
## So the point of overlaying the occurrences with the Koppen is simply to identify what zones the species falls in. We 
## would do that by calculating the % of a species' records in each zone. I have done similar things before, and I guarantee 
## it still wont be straightforward with grouping species, but it will be a start.


#########################################################################################################################
## To save time, load in previous data
source('./R/HIA_LIST_MATCHING.R')
load("./data/base/HIA_LIST/COMBO/KOPPEN_BY_SPP.RData")


## Intersect this with the test spp
dim(KOPPEN.CAST)
names(KOPPEN.CAST)
View(KOPPEN.CAST)


## Just look at the experimental species
mod.compare.spp = trimws(sort(unique(c(renee.full$Species, 
                                       MQ.glasshouse$Species))))


## Check the contextual data for just the test species
intersect(mod.compare.spp, COMBO.ALL$searchTaxon)
MOD.COMP.NICHE = COMBO.ALL[COMBO.ALL$searchTaxon %in% intersect(mod.compare.spp, COMBO.ALL$searchTaxon), ]
View(MOD.COMP.NICHE[, c(2:18)][with(MOD.COMP.NICHE[, c(2:18)], rev(order(COMBO.count))), ])
View(MOD.COMP.NICHE[, c(2:18)][with(MOD.COMP.NICHE[, c(2:18)], rev(order(Number.of.growers))), ])


## So lets choose :
## 5 natives: 2 with many records, 3 with few
## 5 exotics: 2 with many records, 3 with few


## Create a list of just the species to compare model outputs
mod.comp.spp = c('Lomandra longifolia', 'Dianella caerulea',  'Backhousia citriodora',   'Eucalyptus erythrocorys',  'Ficus brachypoda',
                 'Cordyline australis', 'Murraya paniculata', 'Calodendrum capense',     'Liriope muscari',          'Ulmus parvifolia')

##
MOD.COMP = COMBO.ALL[COMBO.ALL$searchTaxon %in% mod.comp.spp, ]
View(MOD.COMP[, c(2:18)][with(MOD.COMP[, c(2:18)], rev(order(COMBO.count))), ])
MOD.SPP = MOD.COMP[, c("searchTaxon", "Origin", "Plant.type", "Number.of.growers", "Top_200")]
  

## Note that the number of records is a bit different (by one or two records). This is a versioning problem and not important, but should be the same
KOPPEN.MOD.SPP = join(MOD.SPP, KOPPEN.CAST)
View(KOPPEN.MOD.SPP)


## save
write.csv(KOPPEN.MOD.SPP, "./data/base/HIA_LIST/COMBO/KOPPEN_TEST_SPP.csv", row.names = FALSE)


