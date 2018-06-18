#########################################################################################################################
########################################  GLASSHOUSE VS GLOBAL NICHES ################################################### 
#########################################################################################################################


#########################################################################################################################
## This code quantifies the correlations between the global environmental niches and their corresponding temperature
## values in the glasshouse. Extremes are thought to be more closelyt related than average values...


#########################################################################################################################
## Read in niche data
CLEAN.NICHE.CONTEXT  = readRDS("./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT_APRIL_2018_COORD_CLEAN.rds")
COMBO.RASTER.CONTEXT = readRDS("./data/base/HIA_LIST/COMBO/CLEAN_ONLY_HIA_SPP.rds")


## Read in heat risk data
HEAT.RISK = read.csv("./data/base/HIA_LIST/RENEE/MOD3_HEAT_RANKS.csv", stringsAsFactors = FALSE)
HEAT.DAYS = raster("./data/base/AWAP/mean_days_35_degrees.tif")
plot(HEAT.DAYS, main = "DAYS > 35 deg")


## Why don't they match up exactly?
HEAT.RISK.NICHE = merge(COMBO.NICHE.CONTEXT, HEAT.RISK, by = "searchTaxon")
dim(HEAT.RISK.NICHE)
names(HEAT.RISK.NICHE)


#########################################################################################################################
## Create a spatial points df
COMBO.POINTS   = SpatialPointsDataFrame(coords      = COMBO.RASTER.CONTEXT[c("lon", "lat")], 
                                        data        = COMBO.RASTER.CONTEXT,
                                        proj4string = CRS.MOL)


## Now extract just the heat days data
HEAT.RASTER <- extract(HEAT.DAYS, COMBO.POINTS) %>% 
  cbind(COMBO.RASTER.CONTEXT, .)
colnames(HEAT.RASTER )[90] <- "HEAT_DAYS"
colnames(HEAT.RASTER )[90]


## Estimate the heatwave niche
HEAT.DF = completeFun(HEAT.RASTER, "HEAT_DAYS")
dim(HEAT.DF)
HEAT.NICHE = niche_estimate(DF = HEAT.DF, colname = "HEAT_DAYS")
dim(HEAT.NICHE) ## Just the Australian species...


## Now merge the full niche df with the heatwave niche df
dim(HEAT.NICHE);dim(HEAT.RISK.NICHE)
HEAT.NICHE.RISK = merge(HEAT.RISK.NICHE, HEAT.NICHE, by = "searchTaxon")


## Get target cols and rename
## Worldclim variables
RISK.WORLDCL.PLOT = HEAT.NICHE.RISK[, c("heat_diff.mean", "Annual_mean_temp_max",  "Max_temp_warm_month_max", 
                                        "Max_temp_warm_month_q95_q05", "Temp_annual_range_mean")] 
names(RISK.WORLDCL.PLOT) = c("HEAT_DIF","TEMP_MX","WARM_MX","TEMP_95_5", "TEMP_RG")
summary(RISK.WORLDCL.PLOT)


## Update this
RISK.DAYS.PLOT = HEAT.NICHE.RISK[, c("heat_diff.mean", "HEAT_DAYS_max",  "HEAT_DAYS_95", 
                                     "HEAT_DAYS_q95_q05", "HEAT_DAYS_range")] 
names(RISK.DAYS.PLOT) = c("HEAT_DIF","HEAT_MX","HEAT_95","HEAT_95_5", "HEAT_RG")
summary(RISK.DAYS.PLOT)


#########################################################################################################################
## Plot the glasshouse data vs the worldclim niches
CairoPNG(width = 8090, height = 8090, 
         file = "./output/figures/Glasshouse_niches/HEAT_DIFF_v_WORLDCLIM.png", 
         canvas = "white", bg = "white", units = "px", dpi = 600)

par(mar = c(8, 8, 6, 4), 
    mgp = c(6, 2, 0))

par(lwd = 2)

pairs(HEAT.WORLDCL.PLOT,  
      lower.panel = panel.cor,
      #diag.panel  = panel.hist,
      upper.panel = panel.smooth, #function(...) smoothScatter(..., nrpoints = 0, add = TRUE))
      cex.labels  = 1.5, cex.axis = 2, font.labels = 2)

## finish the device
dev.off()


#########################################################################################################################
## Plot the glasshouse data vs the worldclim niches
CairoPNG(width = 8090, height = 8090, 
         file = "./output/figures/Glasshouse_niches/HEAT_DIFF_v_HEAT_DAYS.png", 
         canvas = "white", bg = "white", units = "px", dpi = 600)

par(mar = c(8, 8, 6, 4), 
    mgp = c(6, 2, 0))

par(lwd = 2)

pairs(RISK.DAYS.PLOT,  
      lower.panel = panel.cor,
      #diag.panel  = panel.hist,
      upper.panel = panel.smooth, #function(...) smoothScatter(..., nrpoints = 0, add = TRUE))
      cex.labels  = 1.5, cex.axis = 2, font.labels = 2)

## finish the device
dev.off()



#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################
