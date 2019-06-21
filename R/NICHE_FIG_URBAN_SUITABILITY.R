#########################################################################################################################
############################################# ESTIMATE SPECIES NICHES ################################################### 
#########################################################################################################################



#########################################################################################################################
## This code creates a figure of one narrowly distributed plant that is widely cultivated, and one widely distributed plant
## that is not rarely cultivated


## Read in data
COMBO.NICHE.CONTEXT  = readRDS('data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT_APRIL_2018_COORD_CLEAN.rds')
COMBO.RASTER.CONTEXT = readRDS('data/base/HIA_LIST/COMBO/CLEAN_ONLY_HIA_SPP.rds')
dim(COMBO.NICHE.CONTEXT)
names(COMBO.NICHE.CONTEXT)


## So we are looking for two species to compare :
## One species which is narrowly distributed in climates space, but widely sold
## Ones species which is widely distributed in climates space, but rarely sold


## What is the distribution of raingfall and temperature niches
summary(COMBO.NICHE.CONTEXT$Annual_mean_temp_q95_q05)
summary(COMBO.NICHE.CONTEXT$Annual_precip_q95_q05)


## Wide unsold?
wide_unsold     = subset(COMBO.NICHE.CONTEXT, Annual_mean_temp_q95_q05 >= 15 & Annual_precip_q95_q05 >= 1500 & Number.of.growers <= 5)
wide.unsold.spp = wide_unsold$searchTaxon
View(wide_unsold)


## Now try subsetting the records to just the wide_unsold
REC.WIDE.UNSOLD = COMBO.RASTER.CONTEXT[COMBO.RASTER.CONTEXT$searchTaxon %in% wide.unsold.spp, ]


## Create polygon of land surface
LAND       = readRDS("F:/green_cities_sdm/data/base/CONTEXTUAL/LAND_world.rds")
aus = readRDS("F:/green_cities_sdm/data/base/CONTEXTUAL/aus_states.rds") %>%
  spTransform(ALB.CONICAL)

## Sambucus canadensis



#########################################################################################################################
## Loop over the species list and plot the occurrence data for each to check the data bias
# TAXA = as.list(sort(unique(intersect(SDM.DATA.ALL$searchTaxon, spp.mile))))
# 
# for (i in 1:length(TAXA)) {
# 
#   ## Create points for each species
#   spp.points <- SDM.DATA.ALL[SDM.DATA.ALL$searchTaxon == TAXA[i], ] %>%
#     spTransform(CRS('+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'))
# 
#   ## Print to file
#   save_name = gsub(' ', '_', TAXA[i])
#   save_dir  = "output/maxent/summary"
#   png(sprintf('%s/%s_%s.png', save_dir,
#               save_name, "Australian_points"),
#       3236, 2000, units = 'px', res = 300)
# 
#   ## set margins
#   par(mar   = c(3, 3, 5, 3),  ## b, l, t, r
#       #mgp   = c(9.8, 2.5, 0),
#       oma   = c(1.5, 1.5, 1.5, 1.5))
# 
#   ## Plot just the Australian points
#   plot(aus, main = TAXA[i])
#   points(spp.points, col = "red", cex = .3, pch = 19)
# 
#   ## Finish the device
#   dev.off()
# 
# }


#########################################################################################################################
## Update
write.csv(wide_unsold, "./data/base/HIA_LIST/COMBO/WIDE_DIST_UNSOLD_SPP.csv", row.names = FALSE)






#########################################################################################################################
#####################################################  TBC ##############################################################
#########################################################################################################################
