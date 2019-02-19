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


## Calculate AOO for all species
## Send Alessandro a shapefile of Diana's species

## Build the table ::
## Email Anna and Sarah
## Get the list of species with provenances and calculate the niches separately for this  





#########################################################################################################################
## 1). COMBINE NICHE DATA WITH GLASSHOUSE DATA
#########################################################################################################################

## Try reading in diana's species too
DIANA.NICHE.CONTEXT  = readRDS(paste0(DATA_path, 'COMBO_NICHE_CONTEXT_',  OCC_SOURCE, '_RECORDS_', 'DIANA_SPP', '.rds'))
DIAN.RASTER.CONTEXT  = readRDS(paste0(DATA_path, 'COMBO_RASTER_CONTEXT_', OCC_SOURCE, '_RECORDS_', 'DIANA_SPP', '.rds'))

#########################################################################################################################
## Read in niche data
## TRAIT.NICHE.CONTEXT  = COMBO.NICHE.CONTEXT
## TRAIT.RASTER.CONTEXT = COMBO.RASTER.CONTEXT
TRAIT.NICHE.CONTEXT  = readRDS(paste0(DATA_path, 'COMBO_NICHE_CONTEXT_',  OCC_SOURCE, '_RECORDS_', save_run, '.rds'))
TRAIT.RASTER.CONTEXT = readRDS(paste0(DATA_path, 'COMBO_RASTER_CONTEXT_', OCC_SOURCE, '_RECORDS_', save_run, '.rds'))


##
dim(TRAIT.NICHE.CONTEXT)
dim(TRAIT.RASTER.CONTEXT)


## Now find the match between the trait species and the trait species... 
length(intersect(TRAIT.SPP$searchTaxon, TRAIT.NICHE.CONTEXT$searchTaxon))


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
## 2). PLOT OCCURRENCE DATA
#########################################################################################################################


#########################################################################################################################
## Combine Australian and global maps with histograms of GBIF and AWAP data
## Loop over the species list and plot the occurrence data for each to check the data bias
OCC.TAXA = as.list(unique(TRAIT.RASTER.CONTEXT$searchTaxon))

for (i in 1:length(OCC.TAXA)) {

  ## Create points for each species
  spp.points <- TRAIT.POINTS[TRAIT.POINTS$searchTaxon == OCC.TAXA[i], ] #%>%
    #spTransform(CRS('+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'))

  ## Print to file
  # save_name = gsub(' ', '_', TAXA[i])
  # save_dir  = "data/base/HIA_LIST/RENEE"
  # png(sprintf('%s/%s_%s.png', save_dir,
  #             save_name, "Australian_points"),
  #     3236, 2000, units = 'px', res = 300)

  ## set margins
  # par(mar   = c(3, 3, 5, 3),  ## b, l, t, r
  #     #mgp   = c(9.8, 2.5, 0),
  #     oma   = c(1.5, 1.5, 1.5, 1.5))

  ## Plot just the Australian points
  # plot(LAND, main = OCC.TAXA[i])
  # points(spp.points, col = "red", cex = .5, pch = 19)
  # 
  ## Then plot the GBIF histogram?
  
  ## Plot just the Australian points
  plot(aus, main = OCC.TAXA[i])
  points(spp.points, col = "red", cex = .5, pch = 19)
  
  
  ## Then plot the AWAP histogram

  
  ## Finish the device
  # dev.off()

}


#########################################################################################################################
## Now save the points as shapefile to map in Arc
TRAIT.SPDF   = TRAIT.POINTS[c("searchTaxon", "lon", "lat")]
projection(TRAIT.SPDF );names(TRAIT.SPDF)


## Save the shapefile, to be subsampled in ArcMap
writeOGR(obj = TRAIT.SPDF , dsn = "./data/base/HIA_LIST/COMBO", layer = "GLASSHOUSE_TRAIT_SPDF ", driver = "ESRI Shapefile")





#########################################################################################################################
## 3). CREATE HISTOGRAMS
#########################################################################################################################


##############################################################################################
## histograms of temperature and rainfall
HIST.TAXA = (as.list(unique(TRAIT.RASTER.CONTEXT$searchTaxon)))
HIST.TAXA = losers
names(TRAIT.RASTER.CONTEXT)


## Print the histograms to screen
Print_global_histogram(taxa.list    = HIST.TAXA, 
                       DF           = TRAIT.RASTER.CONTEXT,  ## 33 is a problem: Cupianopsis anacardiodes
                       env.var.1    = "Max_temp_warm_month",   
                       env.col.1    = "orange",  
                       env.units.1  = "°C",
                       env.var.2    = "Annual_precip",   
                       env.col.2    = "blue",     
                       env.units.2  = "°C")


Boxplot_GBIF_records(taxa.list = HIST.TAXA,       DF = TRAIT.RASTER.CONTEXT,
                     env.1 = "Max_temp_warm_month",   env.col.1 = "orange",     env.units.1 = "°C",
                     env.2 = "Annual_precip",         env.col.2 = "blue",       env.units.2 = "°mm")





#########################################################################################################################
## 4). ESTIMATE HEAT EXPOSURE NICHES
#########################################################################################################################


#########################################################################################################################
## merge the full niche df with the heatwave niche df
TRAIT.NICHE.RISK = merge(TRAIT.NICHE.WOODY ,   TRAIT.NICHE.CONTEXT,  by = "searchTaxon")
TRAIT.NICHE.AUS  = merge(TRAIT.NICHE.AUS,      TRAIT.NICHE.CONTEXT,  by = "searchTaxon")
dim(TRAIT.NICHE.RISK)[1]
names(TRAIT.NICHE.RISK)


## Check Diana's calculation
plot(TRAIT.NICHE.RISK$Annual_mean_temp_median, TRAIT.NICHE.RISK$Tcrit_C)


#########################################################################################################################
## TLP vs drought
TLP.DR.INT    = TRAIT.NICHE.RISK[, c("TLP_MPa", 
                                     "Drought_max_int_extr_q95",
                                     "Drought_max_int_extr_max", 
                                     "Drought_max_int_extr_median",
                                     "Drought_max_int_extr_mode")] 

TLP.DR.REL.INT    = TRAIT.NICHE.RISK[, c("TLP_MPa", 
                                         "Drought_max_rel_int_extr_q95",
                                         "Drought_max_rel_int_extr_max", 
                                         "Drought_max_rel_int_extr_median",
                                         "Drought_max_rel_int_extr_mode")] 

TCRIT.DR.REL.INT    = TRAIT.NICHE.RISK[, c("Tcrit_C", 
                                           "Drought_max_rel_int_extr_q95",
                                           "Drought_max_rel_int_extr_max", 
                                           "Drought_max_rel_int_extr_median",
                                           "Drought_max_rel_int_extr_mode")] 

## Rename
names(TLP.DR.INT)         = c("TLP",     "95%",  "MAX",  "MEDIAN", "MODE")
names(TLP.DR.REL.INT)     = c("TLP",     "95%",  "MAX",  "MEDIAN", "MODE")
names(TCRIT.DR.REL.INT )  = c("Tcrit",   "95%",  "MAX",  "MEDIAN", "MODE")
summary(TLP.DR.REL.INT)





#########################################################################################################################
## 3). PLOT GLASSHOUSE TRAITS VS HEAT EXPOSURE NICHES
#########################################################################################################################


#########################################################################################################################
## Plot Turgor loss point vs. relative drought intensity
##  ............................................................................................


#########################################################################################################################
## Plot leaf turgor loss point vs the AUS Max drought intesntiy - MAT (1950-2000)
plot(awap.extreme[["Drought_max_int_extr"]], main = "Maximum drought intensity (AWAP, mm 1900-2011)")


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
## 4). RUN GAMS AND PLOT FOR INDIVIDUAL RELATIONSHIPS
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
## Focus on TLP vs "Drought_max_rel_int_extr"


#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################