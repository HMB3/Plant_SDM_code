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
CLEAN.NICHE.CONTEXT  = readRDS("./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT_APRIL_2018_COORD_CLEAN.rds")
COMBO.RASTER.CONTEXT = readRDS("./data/base/HIA_LIST/COMBO/CLEAN_ONLY_HIA_SPP.rds")


## Read in heat risk data
HEAT.RISK = read.csv("./data/base/HIA_LIST/RENEE/MOD3_HEAT_RANKS_072018.csv", stringsAsFactors = FALSE)
HEAT.DAYS = raster("./data/base/AWAP/mean_days_35_degrees.tif")
plot(HEAT.DAYS, main = "DAYS > 35 deg")


## Why don't they match up exactly?
colnames(HEAT.RISK)[colnames(HEAT.RISK)=="Species"] <- "searchTaxon"
HEAT.NICHE.AUS = merge(COMBO.NICHE.CONTEXT, HEAT.RISK, by = "searchTaxon")
HEAT.NICHE.AUS = completeFun(HEAT.NICHE.AUS , "Tcrit_C")

dim(HEAT.NICHE.AUS)
names(HEAT.NICHE.AUS)


## But how many species have good records?
COMBO.RASTER.GLASSHOUSE  = COMBO.RASTER.CONTEXT[COMBO.RASTER.CONTEXT$searchTaxon %in% unique(HEAT.RISK$searchTaxon), ] 
COMBO.NICHE.GLASSHOUSE   = COMBO.NICHE.CONTEXT[COMBO.NICHE.CONTEXT$searchTaxon %in% unique(HEAT.RISK$searchTaxon), ] 
GLASSHOUSE.COUNT         = head(COMBO.NICHE.GLASSHOUSE, 36)[, c("searchTaxon", "COMBO.count", "AUS_RECORDS" )]
GLASSHOUSE.COUNT         = GLASSHOUSE.COUNT[with(GLASSHOUSE.COUNT, order(COMBO.count)), ]
View(GLASSHOUSE.COUNT)


#########################################################################################################################
## Create a spatial points df of the worldclim data, for just the experimental species
COMBO.POINTS   = SpatialPointsDataFrame(coords      = COMBO.RASTER.CONTEXT[c("lon", "lat")], 
                                        data        = COMBO.RASTER.CONTEXT,
                                        proj4string = CRS.MOL)
dim(COMBO.POINTS)


## Now extract just the heat days data
HEAT.RASTER <- extract(HEAT.DAYS, COMBO.POINTS) %>% 
  cbind(COMBO.RASTER.CONTEXT, .)
colnames(HEAT.RASTER )[90] <- "HEAT_DAYS"
colnames(HEAT.RASTER )[90]


## And estimate the heatwave niche
HEAT.DF = completeFun(HEAT.RASTER, "HEAT_DAYS")
dim(HEAT.DF)
HEAT.NICHE = niche_estimate(DF = HEAT.DF, colname = "HEAT_DAYS")
dim(HEAT.NICHE) ## Just the Australian species, because we extracted for only an AUS raster


## Now merge the full niche df with the heatwave niche df
HEAT.NICHE.RISK = merge(HEAT.NICHE.AUS, HEAT.NICHE, by = "searchTaxon")
dim(HEAT.NICHE.RISK);dim(HEAT.NICHE.AUS)


#########################################################################################################################
## Get target cols and rename
## Worldclim variables
OSMP.WORLDCL.PLOT = HEAT.NICHE.AUS[, c("OsmPot_MPa", "Annual_mean_temp_max",  "Max_temp_warm_month_max", 
                                        "Max_temp_warm_month_q95_q05", "Temp_annual_range_mean")]

TLP.WORLDCL.PLOT  = HEAT.NICHE.AUS[, c("TLP_MPa", "Annual_mean_temp_max",  "Max_temp_warm_month_max", 
                                       "Max_temp_warm_month_q95_q05", "Temp_annual_range_mean")]

TCR.WORLDCL.PLOT  = HEAT.NICHE.AUS[, c("Tcrit_C", "Annual_mean_temp_max",  "Max_temp_warm_month_max", 
                                       "Max_temp_warm_month_q95_q05", "Temp_annual_range_mean")]

DSM.WORLDCL.PLOT  = HEAT.NICHE.AUS[, c("Desicc_Max", "Annual_mean_temp_max",  "Max_temp_warm_month_max", 
                                       "Max_temp_warm_month_q95_q05", "Temp_annual_range_mean")]


## rename
names(OSMP.WORLDCL.PLOT) = c("OSMPO", "TEMP_MX","WARM_MX","TEMP_95_5", "TEMP_RG");
names(TLP.WORLDCL.PLOT)  = c("TLP",   "TEMP_MX","WARM_MX","TEMP_95_5", "TEMP_RG");
names(TCR.WORLDCL.PLOT)  = c("Tcrit", "TEMP_MX","WARM_MX","TEMP_95_5", "TEMP_RG");
names(DSM.WORLDCL.PLOT)  = c("DSM",   "TEMP_MX","WARM_MX","TEMP_95_5", "TEMP_RG")
summary(TCR.WORLDCL.PLOT)


#########################################################################################################################
## Heat risk variables
OSMP.HEAT.PLOT = HEAT.NICHE.RISK[, c("OsmPot_MPa", "HEAT_DAYS_max",  "HEAT_DAYS_q95", 
                                    "HEAT_DAYS_q95_q05", "HEAT_DAYS_range")]

TLP.HEAT.PLOT  = HEAT.NICHE.RISK[, c("TLP_MPa", "HEAT_DAYS_max",  "HEAT_DAYS_q95", 
                                    "HEAT_DAYS_q95_q05", "HEAT_DAYS_range")]

TCR.HEAT.PLOT  = HEAT.NICHE.RISK[, c("Tcrit_C", "HEAT_DAYS_max",  "HEAT_DAYS_q95", 
                                    "HEAT_DAYS_q95_q05", "HEAT_DAYS_range")]

DSM.HEAT.PLOT  = HEAT.NICHE.RISK[, c("Desicc_Max", "HEAT_DAYS_max",  "HEAT_DAYS_q95", 
                                    "HEAT_DAYS_q95_q05", "HEAT_DAYS_range")]


## rename
names(OSMP.HEAT.PLOT) = c("OSMPO", "HEAT_MX","HEAT_95","HEAT_95_5", "HEAT_RG");
names(TLP.HEAT.PLOT)  = c("TLP",   "HEAT_MX","HEAT_95","HEAT_95_5", "HEAT_RG");
names(TCR.HEAT.PLOT)  = c("Tcrit",   "HEAT_MX","HEAT_95","HEAT_95_5", "HEAT_RG");
names(DSM.HEAT.PLOT)  = c("DSM",   "HEAT_MX","HEAT_95","HEAT_95_5", "HEAT_RG")
summary(OSMP.HEAT.PLOT)





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
