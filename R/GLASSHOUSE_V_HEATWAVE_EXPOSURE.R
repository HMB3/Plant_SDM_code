#########################################################################################################################
########################################  GLASSHOUSE VS GLOBAL NICHES ################################################### 
#########################################################################################################################


#########################################################################################################################
## This code quantifies the correlations between the global environmental niches and their corresponding 
## values in the glasshouse. Extremes are thought to be more closely related than average values...


#########################################################################################################################
## 1). COMBINE NICHE DATA WITH GLASSHOUSE DATA
#########################################################################################################################


## Pick up from here.....................................................................................................


#########################################################################################################################
## Read in niche data
TRAIT.NICHE.CONTEXT  = readRDS(paste0('data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT_HEATWAVE',  save_run, '.rds'))
TRAIT.RASTER.CONTEXT = readRDS(paste0('data/base/HIA_LIST/COMBO/COMBO_AWAP_CONVERT_',  save_run, '.rds'))


##
dim(TRAIT.NICHE.CONTEXT)
dim(TRAIT.RASTER.CONTEXT)


## Now find the match between the trait species and the trait species... 
length(intersect(TRAIT.SPP$searchTaxon,    TRAIT.NICHE.CONTEXT$searchTaxon))


#########################################################################################################################
## Just get the species with data for Tcrit?
TRAIT.NICHE.AUS = completeFun(TRAIT.SPP, "Tcrit_C")
TRAIT.NICHE.WOODY = subset(TRAIT.NICHE.AUS, Native_Woody == "yes")
unique(TRAIT.NICHE.WOODY$Native_Woody)
dim(TRAIT.NICHE.WOODY)
NICHE.SPP = TRAIT.NICHE.WOODY 


## But how many species have good records? This will change a bit with Alessandro's data
COMBO.NICHE.GLASSHOUSE  = TRAIT.NICHE.CONTEXT[TRAIT.NICHE.CONTEXT$searchTaxon %in% unique(TRAIT.SPP$searchTaxon), ] 
TRAIT.RASTER.CONTEXT    = TRAIT.RASTER.CONTEXT[TRAIT.RASTER.CONTEXT$searchTaxon %in% unique(TRAIT.SPP$searchTaxon), ] 
length(unique(TRAIT.RASTER.CONTEXT$searchTaxon))


#########################################################################################################################
## Create a spatial points df of the worldclim data, for just the experimental species
TRAIT.POINTS   = SpatialPointsDataFrame(coords      = TRAIT.RASTER.CONTEXT [c("lon", "lat")], 
                                        data        = TRAIT.RASTER.CONTEXT ,
                                        proj4string = CRS.WGS.84)

TRAIT.POINTS  <- spTransform(TRAIT.POINTS, CRS("+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))




#########################################################################################################################
## 4). ESTIMATE HEAT EXPOSURE NICHES
#########################################################################################################################


#########################################################################################################################
## merge the full niche df with the heatwave niche df
TRAIT.NICHE.RISK = merge(TRAIT.NICHE.WOODY ,   TRAIT.NICHE.CONTEXT,  by = "searchTaxon")
TRAIT.NICHE.AUS = merge(TRAIT.NICHE.AUS,       TRAIT.NICHE.CONTEXT,  by = "searchTaxon")
dim(TRAIT.NICHE.RISK)[1]
names(TRAIT.NICHE.RISK)


## Check Diana's calculation
plot(TRAIT.NICHE.RISK$Annual_mean_temp_median, TRAIT.NICHE.RISK$Tcrit_C)

#########################################################################################################################
## TLP vs drought
TRAIT.NICHE.RISK$searchTaxon
length(TRAIT.NICHE.RISK$searchTaxon)
TCRIT.HWN    = TRAIT.NICHE.RISK[, c("Tcrit_C", 
                                    "HWN_q95",
                                    "HWN_max", 
                                    "HWN_median",
                                    "HWN_mode")] 

TCRIT.HWA    = TRAIT.NICHE.RISK[, c("Tcrit_C", 
                                    "HWA_q95",
                                    "HWA_max", 
                                    "HWA_median",
                                    "HWA_mode")] 

TCRIT.HWF    = TRAIT.NICHE.RISK[, c("Tcrit_C", 
                                    "HWF_q95",
                                    "HWF_max", 
                                    "HWF_median",
                                    "HWF_mode")]

TCRIT.HCH    = TRAIT.NICHE.RISK[, c("Tcrit_C", 
                                    "HW_CUM_HOT_q95",
                                    "HW_CUM_HOT_max", 
                                    "HW_CUM_HOT_median",
                                    "HW_CUM_HOT_mode")] 

## Rename
names(TCRIT.HWN)     = c("TLP",   "95%",  "MAX",  "MEDIAN", "MODE")
names(TCRIT.HWA)     = c("TLP",   "95%",  "MAX",  "MEDIAN", "MODE")
names(TCRIT.HWF)     = c("Tcrit", "95%",  "MAX",  "MEDIAN", "MODE")
names(TCRIT.HCH)     = c("Tcrit", "95%",  "MAX",  "MEDIAN", "MODE")
summary(TCRIT.HWN)





#########################################################################################################################
## 3). PLOT GLASSHOUSE TRAITS VS HEAT EXPOSURE NICHES
#########################################################################################################################


#########################################################################################################################
## Plot Turgor loss point vs. relative drought intensity
##  ............................................................................................


#########################################################################################################################
## Plot TCRIT vs the average heatwave intensity
plot(awap.heatwave[["HWN"]],        main = "No. indiv heatwaves in a season (AWAP, 1961-19901)")

# CairoPNG(width = 8090, height = 8090, 
#          file = "./output/figures/Glasshouse_niches/TLP_v_AWAP_DROUGHT_INTENSITY.png", 
#          canvas = "white", bg = "white", units = "px", dpi = 600)
# 
par(mar = c(8, 8, 6, 4),
    mgp = c(6, 2, 0))

par(lwd = 2)

pairs(TCRIT.HWN,  
      lower.panel = panel.cor,
      #diag.panel  = panel.hist,
      upper.panel = panel.smooth, #function(...) smoothScatter(..., nrpoints = 0, add = TRUE))
      cex.labels  = 1.5, cex.axis = 2, font.labels = 2,
      main = "Tcrit vs No. indiv heatwaves (AWAP 1961-1990)")


#########################################################################################################################
## Plot Tcrit vs hottest day (1961-1990)
plot(awap.heatwave[["HWA"]], main = "Hottest day of hottest heatwave (AWAP 1961-1990)")

pairs(TCRIT.HWA,  
      lower.panel = panel.cor,
      #diag.panel  = panel.hist,
      upper.panel = panel.smooth, #function(...) smoothScatter(..., nrpoints = 0, add = TRUE))
      cex.labels  = 1.5, cex.axis = 2, font.labels = 2,
      main = "Tcrit vs Hottest day of hottest heatwave (AWAP 1961-1990)")


#########################################################################################################################
## Plot Tcrit vs number of heatwave days in a season (1961-1990)
plot(awap.heatwave[["HWF"]], main = "Number of heatwave days in a season (AWAP 1961-1990)")

pairs(TCRIT.HWF,  
      lower.panel = panel.cor,
      #diag.panel  = panel.hist,
      upper.panel = panel.smooth, #function(...) smoothScatter(..., nrpoints = 0, add = TRUE))
      cex.labels  = 1.5, cex.axis = 2, font.labels = 2,
      main = "Tcrit vs No. heatwave days in season (AWAP 1961-1990)")


#########################################################################################################################
## Plot Tcrit vs number of heatwave days in a season (1961-1990)
plot(awap.heatwave[["HW_CUM_HOT"]], main = "Cumulative heat during the hottest heatwave (AWAP 1961-1990)")

pairs(TCRIT.HCH,  
      lower.panel = panel.cor,
      #diag.panel  = panel.hist,
      upper.panel = panel.smooth, #function(...) smoothScatter(..., nrpoints = 0, add = TRUE))
      cex.labels  = 1.5, cex.axis = 2, font.labels = 2,
      main = "Tcrit vs Cumulative heat during the hottest heatwave (AWAP 1961-1990)")





#########################################################################################################################
## 4). RUN GAMS
#########################################################################################################################


#########################################################################################################################
## Run GAMs
TLP.MAX.GAM = gam(TLP ~ s (MAX, k = 5), 
                  data = TLP.DR.REL.INT, 
                  #family = tw(theta = NULL, link = "log", a = 1.01, b = 1.99), 
                  method = "REML")
summary(TLP.MAX.GAM)[["dev.expl"]]    ## Not much going on there....

TLP.MAX.TEST  = data.frame(MAX = seq(min(TLP.DR.REL.INT[,"MAX"]),
                                     max(TLP.DR.REL.INT[,"MAX"]), 
                                     length = length(TLP.DR.REL.INT[["TLP"]])))

PRED.TLP.MAX = predict(TLP.MAX.GAM, newdata = TLP.MAX.TEST, type ='response')


########################################################
## RAIN
CairoPNG(width = 10000, height = 16180, 
         file = "./output/figures/figure_3/CHAP2_FIG4_UPDATE_TRANSPARENT.png", 
         canvas = "white", bg = "white", units = "px", dpi = 600)

par(mfrow = c(4, 2),
    mar   = c(13.5, 16, 4, 4.8), 
    mgp   = c(11.8, 3, 0),
    oma   = c(1, 1, 1, 1))


## PANEL 1
par(font.axis = 1)

plot(TLP.DR.REL.INT[,"MAX"], TLP.DR.REL.INT[,"TLP"], 
     col = alpha("blue", 0.3), pch = 19, cex = 2, 
     #cex.axis = 5, cex.lab = 6,
     las = 1, xlab = "Max reltive drought (% 1900-2011)", ylab = "TLP")
box(lwd=3)

lines(TLP.MAX.TEST$MAX, PRED.TLP.MAX, col = "orange",  lwd = 8)
















# OCC.TAXA = as.list(unique(TRAIT.RASTER.CONTEXT$searchTaxon))
# 
# for (i in 1:length(PLOT.DROUGHT)) {
#   
#   
#   # CairoPNG(width = 8090, height = 8090, 
#   #          file = "./output/figures/Glasshouse_niches/TLP_v_AWAP_DROUGHT_INTENSITY.png", 
#   #          canvas = "white", bg = "white", units = "px", dpi = 600)
#   # 
#   # par(mar = c(8, 8, 6, 4), 
#   #     mgp = c(6, 2, 0))
#   
#   par(lwd = 2)
#   
#   pairs(TLP.DR.INT,  
#         lower.panel = panel.cor,
#         #diag.panel  = panel.hist,
#         upper.panel = panel.smooth, #function(...) smoothScatter(..., nrpoints = 0, add = TRUE))
#         cex.labels  = 1.5, cex.axis = 2, font.labels = 2,
#         main = "Leaf turgor loss point (TLP) vs Max drought intensity (AWAP 1900-2011)")
#   
#   ## finish the device
#   #dev.off()
#   
#   ## Then plot the GBIF histogram?
#   ## Plot just the Australian points
#   plot(aus, main = OCC.TAXA[i])
#   points(spp.points, col = "red", cex = .5, pch = 19)
#   
#   
#   ## Then plot the AWAP histogram
#   
#   
#   ## Finish the device
#   # dev.off()
#   
# }



#########################################################################################################################
## Save the R object
save.image("GLASS_HOUSE_VS_NICHE.RData")


#########################################################################################################################
## Focus on TLP vs "Drought_max_rel_int_extr"


#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################