#########################################################################################################################
########################################  GLASSHOUSE VS GLOBAL NICHES ################################################### 
#########################################################################################################################


#########################################################################################################################
## This code quantifies the correlations between the global environmental niches and their corresponding 
## values in the glasshouse. Extremes are thought to be more closely related than average values...


## 1). Renee's work - physiology data

## 2). Hybrid manuscript: Tcrit = f range size + MAT + extreme climate + GLS (HIE + MQ) + ecology traits

## 3). Renee's multivariate analysis? 

## 4). Diana's tcrit vs. lat, temp, radiation etc.


## Calculate AOO for all species - check what's going wrong
## Send Alessandro a shapefile of Diana's species

## Build the table ::
## Email Anna and Sarah
## Get the list of species with provenances and calculate the niches separately for this  





#########################################################################################################################
## 1). READ NICHE AND POINT DATA FOR GLASSHOUSE SPECIES
#########################################################################################################################


## Try reading in diana's species too
DIANA.NICHE.CONTEXT  = readRDS(paste0(DATA_path, 'AWAP_NICHE_DATA_TCRIT_SPP.rds'))
DIANA.RASTER.CONTEXT  = readRDS(paste0(DATA_path, 'AWAP_POINT_DATA_TCRIT_SPP.rds'))


#########################################################################################################################
## Read in niche data
TRAIT.NICHE.CONTEXT  = readRDS(paste0(DATA_path, 'AWAP_NICHE_DATA_ALL_TRAIT_SPP.rds'))
TRAIT.RASTER.CONTEXT = readRDS(paste0(DATA_path, 'AWAP_POINT_DATA_ALL_TRAIT_SPP.rds'))


##  Check
dim(DIANA.NICHE.CONTEXT)
dim(DIANA.RASTER.CONTEXT)

dim(TRAIT.NICHE.CONTEXT)
dim(TRAIT.RASTER.CONTEXT)


#########################################################################################################################
## Just get the species with data for Tcrit?
#TRAIT.NICHE.WOODY = subset(TRAIT.NICHE.CONTEXT, Native_Woody == "yes")
# unique(TRAIT.NICHE.WOODY$Native_Woody)
# dim(TRAIT.NICHE.WOODY)
# NICHE.SPP = TRAIT.NICHE.WOODY 


## But how many species have good records? This will change a bit with Alessandro's data
# COMBO.NICHE.GLASSHOUSE = TRAIT.NICHE.CONTEXT[TRAIT.NICHE.CONTEXT$searchTaxon %in% unique(TRAIT.SPP$searchTaxon), ] 
# TRAIT.RASTER.CONTEXT   = TRAIT.RASTER.CONTEXT[TRAIT.RASTER.CONTEXT$searchTaxon %in% unique(TRAIT.SPP$searchTaxon), ] 
# length(unique(TRAIT.RASTER.CONTEXT$searchTaxon))


#########################################################################################################################
## Create a spatial points df of the worldclim data, for just the experimental species
TRAIT.POINTS   = SpatialPointsDataFrame(coords      = TRAIT.RASTER.CONTEXT[c("lon", "lat")], 
                                        data        = TRAIT.RASTER.CONTEXT,
                                        proj4string = CRS.WGS.84)

DIANA.POINTS   = SpatialPointsDataFrame(coords      = DIANA.RASTER.CONTEXT[c("lon", "lat")], 
                                        data        = DIANA.RASTER.CONTEXT,
                                        proj4string = CRS.WGS.84)


#########################################################################################################################
## Now save the points as shapefile to map in Arc
TRAIT.SPDF   = TRAIT.POINTS[c("searchTaxon", "lon", "lat")]
projection(TRAIT.SPDF );names(TRAIT.SPDF)


## Save the shapefile, to be subsampled in ArcMap
writeOGR(obj = TRAIT.POINTS, dsn = "./data/base/HIA_LIST/COMBO", layer = "GLASSHOUSE_TRAIT_SPDF ", driver = "ESRI Shapefile")
writeOGR(obj = DIANA.POINTS, dsn = "./data/base/HIA_LIST/COMBO", layer = "TCRIT_TRAIT_SPDF ", driver = "ESRI Shapefile")





#########################################################################################################################
## 3). CREATE HISTOGRAMS
#########################################################################################################################


##############################################################################################
## histograms of temperature and rainfall
# HIST.TAXA = (as.list(unique(TRAIT.RASTER.CONTEXT$searchTaxon)))
# HIST.TAXA = losers
# names(TRAIT.RASTER.CONTEXT)
# 
# 
# ## Print the histograms to screen
# Print_global_histogram(taxa.list    = HIST.TAXA, 
#                        DF           = TRAIT.RASTER.CONTEXT,  ## 33 is a problem: Cupianopsis anacardiodes
#                        env.var.1    = "Max_temp_warm_month",   
#                        env.col.1    = "orange",  
#                        env.units.1  = "°C",
#                        env.var.2    = "Annual_precip",   
#                        env.col.2    = "blue",     
#                        env.units.2  = "°C")
# 
# 
# Boxplot_GBIF_records(taxa.list = HIST.TAXA,       DF = TRAIT.RASTER.CONTEXT,
#                      env.1 = "Max_temp_warm_month",   env.col.1 = "orange",     env.units.1 = "°C",
#                      env.2 = "Annual_precip",         env.col.2 = "blue",       env.units.2 = "°mm")





#########################################################################################################################
## 4). ESTIMATE HEAT EXPOSURE NICHES
#########################################################################################################################


#########################################################################################################################
## Merge the full niche df with the heatwave niche df
TRAIT.NICHE.DROUGHT = join(DROUGHT.TRAIT, TRAIT.NICHE.CONTEXT)
TRAIT.NICHE.TCRIT   = join(TCRIT,         DIANA.NICHE.CONTEXT)


## Write out 
write.csv(TRAIT.NICHE.DROUGHT, paste0(DATA_path, 'TRAIT_NICHE_DROUGHT.csv'), row.names = FALSE)


TRAIT.NICHE.DROUGHT = na.omit(TRAIT.NICHE.DROUGHT)
TRAIT.NICHE.TCRIT   = na.omit(TRAIT.NICHE.TCRIT)


#########################################################################################################################
## TLP vs drought
TLP.DR.INT        = TRAIT.NICHE.DROUGHT[, c("TLP_MPa", 
                                            "Drought_max_int_extr_q95",
                                            "Drought_max_int_extr_max", 
                                            "Drought_max_int_extr_median",
                                            "Drought_max_int_extr_mode")]
TLP.DR.INT = na.omit()


TLP.DR.REL.INT    = TRAIT.NICHE.DROUGHT[, c("TLP_MPa",                         ## This is the one 
                                            "Drought_max_rel_int_extr_q95",
                                            "Drought_max_rel_int_extr_max", 
                                            "Drought_max_rel_int_extr_median",
                                            "Drought_max_rel_int_extr_mode")] 

TCRIT.DR.REL.INT  = TRAIT.NICHE.TCRIT[, c("Tcrit", 
                                          "Drought_max_rel_int_extr_q95",
                                          "Drought_max_rel_int_extr_max", 
                                          "Drought_max_rel_int_extr_median",
                                          "Drought_max_rel_int_extr_mode")] 

## Rename
names(TLP.DR.INT)         = c("TLP",     "95%",  "MAX",  "MEDIAN", "MODE")
names(TLP.DR.REL.INT)     = c("TLP",     "95%",  "MAX",  "MEDIAN", "MODE")
names(TCRIT.DR.REL.INT )  = c("Tcrit",   "95%",  "MAX",  "MEDIAN", "MODE")
summary(TLP.DR.REL.INT)
dim(TLP.DR.REL.INT)




#########################################################################################################################
## 3). PLOT GLASSHOUSE TRAITS VS HEAT EXPOSURE NICHES
#########################################################################################################################


#########################################################################################################################
## Plot Turgor loss point vs. relative drought intensity
##  ............................................................................................


#########################################################################################################################
## Plot leaf turgor loss point vs the AUS Max drought intesntiy - MAT (1950-2000)
plot(awap.extreme[["Drought_max_int_extr"]], main = "Maximum drought intensity (AWAP, mm 1900-2011)")
points(TRAIT.POINTS, add = "TRUE")


# CairoPNG(width = 8090, height = 8090, 
#          file = "./output/figures/Glasshouse_niches/TLP_v_AWAP_DROUGHT_INTENSITY.png", 
#          canvas = "white", bg = "white", units = "px", dpi = 600)
# 
par(mar = c(8, 8, 6, 4),
    mgp = c(6, 2, 0))

par(lwd = 2)

pairs(TLP.DR.INT,  
      lower.panel = panel.cor,
      #diag.panel  = panel.hist,
      upper.panel = panel.smooth, #function(...) smoothScatter(..., nrpoints = 0, add = TRUE))
      cex.labels  = 1.5, cex.axis = 2, font.labels = 2,
      main = "Leaf turgor loss point (TLP) vs Max drought intensity (AWAP 1900-2011)")


#########################################################################################################################
## Plot leaf turgor loss point vs the AUS relative drought intesntiy - MAT (1950-2000)
plot(awap.extreme[["Drought_max_rel_int_extr"]], main = "Relative drought intensity (AWAP, % 1900-2011)")

pairs(TLP.DR.REL.INT,  
      lower.panel = panel.cor,
      #diag.panel  = panel.hist,
      upper.panel = panel.smooth, #function(...) smoothScatter(..., nrpoints = 0, add = TRUE))
      cex.labels  = 1.5, cex.axis = 2, font.labels = 2,
      main = "TLP vs Max relative drought intensity (AWAP 1900-2011)")


pairs(TCRIT.DR.REL.INT,  
      lower.panel = panel.cor,
      #diag.panel  = panel.hist,
      upper.panel = panel.smooth, #function(...) smoothScatter(..., nrpoints = 0, add = TRUE))
      cex.labels  = 1.5, cex.axis = 2, font.labels = 2,
      main = "Tcrit vs Max relative drought intensity (AWAP 1900-2011)")

## finish the device
#dev.off()





#########################################################################################################################
## 4). RUN GAMS AND PLOT FOR INDIVIDUAL RELATIONSHIPS?
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
# CairoPNG(width = 10000, height = 16180, 
#          file = "./output/figures/figure_3/CHAP2_FIG4_UPDATE_TRANSPARENT.png", 
#          canvas = "white", bg = "white", units = "px", dpi = 600)
# 
# par(mfrow = c(4, 2),
#     mar   = c(13.5, 16, 4, 4.8), 
#     mgp   = c(11.8, 3, 0),
#     oma   = c(1, 1, 1, 1))


## PANEL 1
par(font.axis = 1)

plot(TLP.DR.REL.INT[,"MAX"], TLP.DR.REL.INT[,"TLP"], 
     col = alpha("blue", 0.3), pch = 19, cex = 2, 
     #cex.axis = 5, cex.lab = 6,
     las = 1, xlab = "Max reltive drought (% 1900-2011)", ylab = "TLP")

box(lwd = 3)

lines(TLP.MAX.TEST$MAX, PRED.TLP.MAX, col = "orange",  lwd = 8)

legend("bottomright", bty = "n", cex = 2, pt.cex = 2, 
       text.col = "orange", legend = paste("DE =", format(summary(TLP.MAX.GAM)$dev.expl, digits = 3)))





#########################################################################################################################
## Focus on TLP vs "Drought_max_rel_int_extr"


#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################