#########################################################################################################################
########################################  GLASSHOUSE VS GLOBAL NICHES ################################################### 
#########################################################################################################################


#########################################################################################################################
## This code quantifies the correlations between the global environmental niches and their corresponding temperature
## values in the glasshouse. Extremes are thought to be more closelyt related than average values...


#########################################################################################################################
## 1). COMBINE NICHE DATA WITH GLASSHOUSE DATA
#########################################################################################################################


#########################################################################################################################
## Read in niche data
TRAIT.NICHE.CONTEXT  = readRDS("./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT_APRIL_2018_COORD_CLEAN_NEW_SPP.rds")
TRAIT.RASTER.CONTEXT = readRDS("./data/base/HIA_LIST/COMBO/CLEAN_ONLY_HIA_SPP.rds")


##
dim(TRAIT.NICHE.CONTEXT)
dim(TRAIT.RASTER.CONTEXT)


## Read in heat risk data
HEAT.RISK  = read.csv("./data/base/HIA_LIST/RENEE/MOD3_HEAT_RANKS_072018.csv", stringsAsFactors = FALSE)
TRAIT.SPP  = read.csv("./data/base/HIA_LIST/RENEE/MOD3_TRAIT_SPP.csv", stringsAsFactors = FALSE)
RANK.TRAIT = read.csv("./data/base/HIA_LIST/RENEE/RankingTraits_Control_latest.csv", stringsAsFactors = FALSE)
HEAT.DAYS = raster("./data/base/AWAP/mean_days_35_degrees.tif")
plot(HEAT.DAYS, main = "DAYS > 35 deg")


#########################################################################################################################
## Read in temperature-max_monthly-mean
temp_max_monthly <- raster("./data/base/ANUCLIM/ANUClimate_v1-0_temperature-max_monthly-mean_0-01deg_1976-2005.nc")
class(temp_max_monthly);projection(temp_max_monthly)


##
PET = raster("./data/base/worldclim/world/1km/pet_he_yr1.tif")
AI  = raster("./data/base/worldclim/world/1km/AI_annual/ai_yr")
AI  = AI/10000


## Now find the match between the trait species and the trait species... 
colnames(HEAT.RISK)[colnames(HEAT.RISK)=="Species"]   <- "searchTaxon"
colnames(TRAIT.SPP)[colnames(TRAIT.SPP)=="Species"]   <- "searchTaxon"
colnames(RANK.TRAIT)[colnames(RANK.TRAIT)=="Species"] <- "searchTaxon"
TRAIT.SPP = merge(TRAIT.SPP, HEAT.RISK, by = "searchTaxon")


## How many experimental species do we have climate data for?
setdiff(TRAIT.SPP$searchTaxon,    TRAIT.NICHE.CONTEXT$searchTaxon)
intersect(TRAIT.SPP$searchTaxon,  TRAIT.NICHE.CONTEXT$searchTaxon)
intersect(RANK.TRAIT$searchTaxon, TRAIT.NICHE.CONTEXT$searchTaxon)


## Merge traits and species
HEAT.NICHE.AUS = merge(TRAIT.NICHE.CONTEXT,  RANK.TRAIT, by = "searchTaxon")
HEAT.NICHE.AUS = completeFun(HEAT.NICHE.AUS, "Tcrit_C")
HEAT.NICHE.AUS = subset(HEAT.NICHE.AUS, Native_Woody == "yes")
unique(HEAT.NICHE.AUS$Native_Woody)
dim(HEAT.NICHE.AUS)
NICHE.SPP = HEAT.NICHE.AUS[, c("searchTaxon", "OsmPot_MPa", "TLP_MPa", "Tcrit_C", "Desicc_Max", "FvFm_HW", "TLP_Tcrit")]



dim(HEAT.NICHE.AUS)
names(HEAT.NICHE.AUS)


## But how many species have good records? This will change a bit with Alessandro's data
COMBO.NICHE.GLASSHOUSE  = TRAIT.NICHE.CONTEXT[TRAIT.NICHE.CONTEXT$searchTaxon %in% unique(RANK.TRAIT$searchTaxon), ] 
TRAIT.RASTER.CONTEXT    = TRAIT.RASTER.CONTEXT[TRAIT.RASTER.CONTEXT$searchTaxon %in% unique(RANK.TRAIT$searchTaxon), ] 

GLASSHOUSE.COUNT         = head(COMBO.NICHE.GLASSHOUSE, 36)[, c("searchTaxon", "COMBO.count", "AUS_RECORDS", "Origin", "Plant.type")]
GLASSHOUSE.COUNT         = GLASSHOUSE.COUNT[with(GLASSHOUSE.COUNT, order(COMBO.count)), ]
View(GLASSHOUSE.COUNT)
summary(GLASSHOUSE.COUNT$searchTaxon)

#saveRDS(TRAIT.RASTER.CONTEXT, file = paste("./data/base/HIA_LIST/COMBO/TRAIT_WORLDCLIM_DATA.rds"))



#########################################################################################################################
## Create a spatial points df of the worldclim data, for just the experimental species
TRAIT.POINTS   = SpatialPointsDataFrame(coords      = TRAIT.RASTER.CONTEXT [c("lon", "lat")], 
                                        data        = TRAIT.RASTER.CONTEXT ,
                                        proj4string = CRS.MOL)
dim(TRAIT.POINTS)


## Now extract just the heat days data
HEAT.RASTER <- extract(HEAT.DAYS, TRAIT.POINTS) %>% 
  cbind(TRAIT.RASTER.CONTEXT , .)
colnames(HEAT.RASTER )[90] <- "HEAT_DAYS"
colnames(HEAT.RASTER )[90]


## Now extract just the mean temperature
TEMP.RASTER <- extract(temp_max_monthly, TRAIT.POINTS) %>% 
  cbind(TRAIT.RASTER.CONTEXT , .)
colnames(TEMP.RASTER )[90] <- "Mnt_mn_daily_max_tmp"
colnames(TEMP.RASTER )[90]


## Now extract just the mean temperature
AI.RASTER <- extract(AI, TRAIT.POINTS) %>% 
  cbind(TRAIT.RASTER.CONTEXT , .)
colnames(AI.RASTER )[90] <- "Aridity_index"
colnames(AI.RASTER )[90]



#########################################################################################################################
## Loop over the species list and plot the occurrence data for each to check the data bias
TAXA = as.list(unique(RANK.TRAIT$searchTaxon))

for (i in 1:length(TAXA)) {

  ## Create points for each species
  spp.points <- TRAIT.POINTS[TRAIT.POINTS$searchTaxon == TAXA[i], ] %>%
    spTransform(CRS('+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'))

  ## Print to file
  save_name = gsub(' ', '_', TAXA[i])
  save_dir  = "data/base/HIA_LIST/RENEE"
  png(sprintf('%s/%s_%s.png', save_dir,
              save_name, "Australian_points"),
      3236, 2000, units = 'px', res = 300)

  ## set margins
  par(mar   = c(3, 3, 5, 3),  ## b, l, t, r
      #mgp   = c(9.8, 2.5, 0),
      oma   = c(1.5, 1.5, 1.5, 1.5))

  ## Plot just the Australian points
  plot(aus, main = TAXA[i])
  points(spp.points, col = "red", cex = .3, pch = 19)

  ## Finish the device
  dev.off()

}



#########################################################################################################################
## And estimate the heatwave niche
HEAT.DF = completeFun(HEAT.RASTER, "HEAT_DAYS")
dim(HEAT.DF)
HEAT.NICHE = niche_estimate(DF = HEAT.DF, colname = "HEAT_DAYS")
dim(HEAT.NICHE) ## Just the Australian species, because we extracted for only an AUS raster


## Estimate the Monthly 1976-2005 mean daily maximum temperature niche
TEMP.DF = completeFun(TEMP.RASTER, "Mnt_mn_daily_max_tmp")
dim(TEMP.DF)
TEMP.NICHE = niche_estimate(DF = TEMP.DF, colname = "Mnt_mn_daily_max_tmp")
dim(TEMP.NICHE) ## Just the Australian species, because we extracted for only an AUS raster


## Estimate the Monthly 1976-2005 mean daily maximum temperature niche
AI.DF = completeFun(AI.RASTER, "Aridity_index")
dim(AI.DF)
AI.NICHE = niche_estimate(DF = AI.DF, colname = "Aridity_index")
dim(AI.NICHE) ## Just the Australian species, because we extracted for only an AUS raster


## Now merge the full niche df with the heatwave niche df
HEAT.NICHE.RISK = merge(HEAT.NICHE.AUS,  HEAT.NICHE, by = "searchTaxon")
HEAT.NICHE.RISK = merge(HEAT.NICHE.RISK, TEMP.NICHE, by = "searchTaxon")
HEAT.NICHE.RISK = merge(HEAT.NICHE.RISK, AI.NICHE,   by = "searchTaxon")

dim(HEAT.NICHE.RISK);dim(HEAT.NICHE.AUS)
names(HEAT.NICHE.RISK)


## Renee :: test how the relationship changes if we use the max_q95 MAT throughout the species range, or the northern-edge MAT


#########################################################################################################################
## Get target cols and rename
## Worldclim variables
OSMP.WORLDCL.PLOT = HEAT.NICHE.RISK[, c("OsmPot_MPa", "Annual_mean_temp_max",  "Max_temp_warm_month_max", 
                                       "Max_temp_warm_month_q95_q05", "Temp_annual_range_mean")]

TLP.WORLDCL.PLOT  = HEAT.NICHE.RISK[, c("TLP_MPa", "Annual_mean_temp_max",  "Annual_mean_temp_min", 
                                       "Max_temp_warm_month_max", "Max_temp_warm_month_min")]

TCR.WORLDCL.PLOT  = HEAT.NICHE.RISK[, c("Tcrit_C", 
                                        "Annual_mean_temp_mean",    "Annual_mean_temp_max",    "Annual_mean_temp_min",
                                        "Max_temp_warm_month_mean", "Max_temp_warm_month_max", "Max_temp_warm_month_min")]

DSM.WORLDCL.PLOT  = HEAT.NICHE.RISK[, c("Desicc_Max", "Annual_mean_temp_max",  "Max_temp_warm_month_max", 
                                       "Max_temp_warm_month_q95_q05", "Temp_annual_range_mean")]


## Rename
names(OSMP.WORLDCL.PLOT) = c("OSMPO", "TEMP_MX","WARM_MX","TEMP_95_5", "TEMP_RG");
names(TLP.WORLDCL.PLOT)  = c("TLP", "MAT_MEAN",  "MAT_MX","MAT_MN", "WARM_MX", "WARM_MX");
names(TCR.WORLDCL.PLOT)  = c("TCR",  "MAT_MEAN", "MAT_MX","MAT_MN", "WARM_MN", "WARM_MX", "WARM_MX");
names(DSM.WORLDCL.PLOT)  = c("DSM",   "TEMP_MX","WARM_MX","TEMP_95_5", "TEMP_RG")
summary(TCR.WORLDCL.PLOT)


# #########################################################################################################################
# ## Heat risk variables
# OSMP.HEAT.PLOT = HEAT.NICHE.RISK[, c("OsmPot_MPa", "HEAT_DAYS_max",  "HEAT_DAYS_q95", 
#                                     "HEAT_DAYS_q95_q05", "HEAT_DAYS_range")]
# 
# TLP.HEAT.PLOT  = HEAT.NICHE.RISK[, c("TLP_MPa", "HEAT_DAYS_max",  "HEAT_DAYS_q95", 
#                                     "HEAT_DAYS_q95_q05", "HEAT_DAYS_range")]
# 
# TCR.HEAT.PLOT  = HEAT.NICHE.RISK[, c("Tcrit_C", "HEAT_DAYS_max",  "HEAT_DAYS_q95", 
#                                     "HEAT_DAYS_q95_q05", "HEAT_DAYS_range")]
# 
# DSM.HEAT.PLOT  = HEAT.NICHE.RISK[, c("Desicc_Max", "HEAT_DAYS_max",  "HEAT_DAYS_q95", 
#                                     "HEAT_DAYS_q95_q05", "HEAT_DAYS_range")]
# 
# 
# ## Rename
# names(OSMP.HEAT.PLOT) = c("OSMPO", "HEAT_MX", "HEAT_95", "HEAT_95_5", "HEAT_RG");
# names(TLP.HEAT.PLOT)  = c("TLP",   "HEAT_MX", "HEAT_95", "HEAT_95_5", "HEAT_RG");
# names(TCR.HEAT.PLOT)  = c("Tcrit", "HEAT_MX", "HEAT_95", "HEAT_95_5", "HEAT_RG");
# names(DSM.HEAT.PLOT)  = c("DSM",   "HEAT_MX", "HEAT_95", "HEAT_95_5", "HEAT_RG")
# summary(OSMP.HEAT.PLOT)





#########################################################################################################################
## 2). PLOT GLASSHOUSE TRAITS VS WORLDCLIM RASTER DATA
#########################################################################################################################


#########################################################################################################################
## Plot leaf turgor loss point vs the worldclim niches
CairoPNG(width = 8090, height = 8090, 
         file = "./output/figures/Glasshouse_niches/TLP_v_WORLDCLIM.png", 
         canvas = "white", bg = "white", units = "px", dpi = 600)

par(mar = c(8, 8, 6, 4), 
    mgp = c(6, 2, 0))

par(lwd = 2)

pairs(TLP.WORLDCL.PLOT,  
      lower.panel = panel.cor,
      #diag.panel  = panel.hist,
      upper.panel = panel.smooth, #function(...) smoothScatter(..., nrpoints = 0, add = TRUE))
      cex.labels  = 1.5, cex.axis = 2, font.labels = 2,
      main = "Leaf turgor loss point (TLP) vs worldclim")

## finish the device
dev.off()


#########################################################################################################################
## Plot leaf critical temperature (Tcrit) vs the worldclim niches
CairoPNG(width = 8090, height = 8090, 
         file = "./output/figures/Glasshouse_niches/OSMOPOT_v_WORLDCLIM.png", 
         canvas = "white", bg = "white", units = "px", dpi = 600)

par(mar = c(8, 8, 6, 4), 
    mgp = c(6, 2, 0))

par(lwd = 2)

pairs(TCR.WORLDCL.PLOT,  
      lower.panel = panel.cor,
      #diag.panel  = panel.hist,
      upper.panel = panel.smooth, #function(...) smoothScatter(..., nrpoints = 0, add = TRUE))
      cex.labels  = 1.5, cex.axis = 2, font.labels = 2,
      main = "Leaf critical temperature (Tcrit) vs worldclim")

## finish the device
dev.off()






#########################################################################################################################
## 3). PLOT GLASSHOUSE TRAITS VS HEAT DAYS RASTER DATA
#########################################################################################################################


#########################################################################################################################
## Plot leaf turgor loss point vs the worldclim niches
CairoPNG(width = 8090, height = 8090, 
         file = "./output/figures/Glasshouse_niches/TLP_v_WORLDCLIM.png", 
         canvas = "white", bg = "white", units = "px", dpi = 600)

par(mar = c(8, 8, 6, 4), 
    mgp = c(6, 2, 0))

par(lwd = 2)

pairs(TLP.HEAT.PLOT,  
      lower.panel = panel.cor,
      #diag.panel  = panel.hist,
      upper.panel = panel.smooth, #function(...) smoothScatter(..., nrpoints = 0, add = TRUE))
      cex.labels  = 1.5, cex.axis = 2, font.labels = 2,
      main = "Leaf turgor loss point (TLP) vs Heat days")

## finish the device
dev.off()


#########################################################################################################################
## Plot leaf critical temperature (Tcrit) vs the worldclim niches
CairoPNG(width = 8090, height = 8090, 
         file = "./output/figures/Glasshouse_niches/OSMOPOT_v_WORLDCLIM.png", 
         canvas = "white", bg = "white", units = "px", dpi = 600)

par(mar = c(8, 8, 6, 4), 
    mgp = c(6, 2, 0))

par(lwd = 2)

pairs(TCR.HEAT.PLOT,  
      lower.panel = panel.cor,
      #diag.panel  = panel.hist,
      upper.panel = panel.smooth, #function(...) smoothScatter(..., nrpoints = 0, add = TRUE))
      cex.labels  = 1.5, cex.axis = 2, font.labels = 2,
      main = "Leaf critical temperature (Tcrit) vs heat days")

## finish the device
dev.off()



#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################
