#########################################################################################################################
########################################  GLASSHOUSE VS GLOBAL NICHES ################################################### 
#########################################################################################################################


#########################################################################################################################
## This code quantifies the correlations between the global environmental niches and their corresponding 
## values in the glasshouse. Extremes are thought to be more closely related than average values...


#########################################################################################################################
## 1). COMBINE NICHE DATA WITH GLASSHOUSE DATA
#########################################################################################################################


#########################################################################################################################
## Read in niche data
TRAIT.NICHE.CONTEXT  = readRDS("./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT_TRAIT_SPP.rds")
TRAIT.RASTER.CONTEXT = readRDS("./data/base/HIA_LIST/COMBO/COMBO_AWAP_CONVERT_TRAIT_SPP.rds")


##
View(TRAIT.NICHE.CONTEXT)
dim(TRAIT.RASTER.CONTEXT)


#########################################################################################################################
## Read in the trait and raster data
awap.extreme = stack(list.files(as.character('./data/base/AWAP'), pattern = '.nc$', full.names = TRUE))
names(awap.extreme)
names(awap.extreme) = c("Drought_freq_extr", "Drought_max_dur_extr", "Drought_max_int_extr", "Drought_max_rel_int_extr",
                        "Drought_mean_dur_extr", "Drought_mean_int_extr", "Drought_mean_rel_int_extr")
names(awap.extreme)
projection(awap.extreme)


# PET        = raster("./data/base/worldclim/world/1km/pet_he_yr1.tif")
# AI         = raster("./data/base/worldclim/world/1km/AI_annual/ai_yr")
# AI         = AI/10000
# temp_max_monthly <- raster("./data/base/ANUCLIM/ANUClimate_v1-0_temperature-max_monthly-mean_0-01deg_1976-2005.nc")
# class(temp_max_monthly);projection(temp_max_monthly)


## Now find the match between the trait species and the trait species... 
intersect(TRAIT.SPP$searchTaxon,    TRAIT.NICHE.CONTEXT$searchTaxon)


#########################################################################################################################
## Just get the species with data for Tcrit?
TRAIT.NICHE.AUS = completeFun(TRAIT.SPP, "Tcrit_C")
TRAIT.NICHE.AUS = subset(TRAIT.NICHE.AUS, Native_Woody == "yes")
unique(TRAIT.NICHE.AUS$Native_Woody)
dim(TRAIT.NICHE.AUS)
NICHE.SPP = TRAIT.NICHE.AUS[, c("searchTaxon", "OsmPot_MPa", "TLP_MPa", "Tcrit_C", "Desicc_Max", "FvFm_HW", "TLP_Tcrit")]


## But how many species have good records? This will change a bit with Alessandro's data
COMBO.NICHE.GLASSHOUSE  = TRAIT.NICHE.CONTEXT[TRAIT.NICHE.CONTEXT$searchTaxon %in% unique(TRAIT.SPP$searchTaxon), ] 
TRAIT.RASTER.CONTEXT    = TRAIT.RASTER.CONTEXT[TRAIT.RASTER.CONTEXT$searchTaxon %in% unique(TRAIT.SPP$searchTaxon), ] 
length(unique(TRAIT.RASTER.CONTEXT$searchTaxon))


# GLASSHOUSE.COUNT        = head(COMBO.NICHE.GLASSHOUSE, 36)[, c("searchTaxon", "COMBO.count", "AUS_RECORDS", "Origin", "Plant.type")]
# GLASSHOUSE.COUNT        = GLASSHOUSE.COUNT[with(GLASSHOUSE.COUNT, order(COMBO.count)), ]
# View(GLASSHOUSE.COUNT)
# summary(GLASSHOUSE.COUNT$searchTaxon)

#saveRDS(TRAIT.RASTER.CONTEXT, file = paste("./data/base/HIA_LIST/COMBO/TRAIT_WORLDCLIM_DATA.rds"))
#TRAIT.RASTER.CONTEXT = readRDS("./data/base/HIA_LIST/COMBO/TRAIT_WORLDCLIM_DATA.rds")



#########################################################################################################################
## Create a spatial points df of the worldclim data, for just the experimental species
TRAIT.POINTS   = SpatialPointsDataFrame(coords      = TRAIT.RASTER.CONTEXT [c("lon", "lat")], 
                                        data        = TRAIT.RASTER.CONTEXT ,
                                        proj4string = CRS.WGS.84)

TRAIT.POINTS  <- spTransform(TRAIT.POINTS, CRS("+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# crs(DEG.MAT)  <- "+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
# crs(MAX.TMAX) <- "+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
# DEG.MAT      <- projectRaster(DEG.MAT,    CRS("+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# MAX.TMAX     <- projectRaster(MAX.TMAX,   CRS("+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# projection(TRAIT.POINTS);projection(DEG.MAT);projection(MAX.TMAX)


## Now raster::extract just the max temp - MAT
DEG.MAT.RASTER <- raster::extract(DEG.MAT, TRAIT.POINTS) %>% 
  cbind(TRAIT.RASTER.CONTEXT , .)
names(DEG.MAT.RASTER)[names(DEG.MAT.RASTER) == "."] <- "Diff_max_MAT"
colnames(DEG.MAT.RASTER)[90]


## Now raster::extract just the heat Max tmax
MAX.TMAX.RASTER <- raster::extract(MAX.TMAX, TRAIT.POINTS) %>% 
  cbind(DEG.MAT.RASTER , .)
names(MAX.TMAX.RASTER)[names(MAX.TMAX.RASTER) == "."] <- "Max_tmax"
colnames(MAX.TMAX.RASTER )[91]





#########################################################################################################################
## 2). PLOT OCCURRENCE DATA
#########################################################################################################################


## Combine Australian and global maps with histograms of GBIF and AWAP data


#########################################################################################################################
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
  plot(LAND, main = OCC.TAXA[i])
  points(spp.points, col = "red", cex = .5, pch = 19)
  
  
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
HIST.TAXA = as.list(unique(MAX.TMAX.RASTER$searchTaxon))
names(MAX.TMAX.RASTER)


## Print the histograms to screen
Print_global_histogram(taxa.list    = HIST.TAXA, 
                       DF           = MAX.TMAX.RASTER,  ## 33 is a problem: Cupianopsis anacardiodes
                       env.var.1    = "Max_temp_warm_month",   
                       env.col.1    = "orange",  
                       env.units.1  = "°C",
                       env.var.2    = "Max_tmax",   
                       env.col.2    = "red",     
                       env.units.2  = "°C")


## Save the histograms to file?
histogram_GBIF_records(taxa.list = HIST.TAXA[1:2], DF = MAX.TMAX.RASTER,
                       env.var.1 = "Max_temp_warm_month",   env.col.1 = "orange",     env.units.1 = "°C",
                       env.var.2 = "Max_tmax",              env.col.2 = "firebrick1", env.units.2 = "°C")



#########################################################################################################################
## 4). ESTIMATE HEAT EXPOSURE NICHES
#########################################################################################################################


#########################################################################################################################
## And estimate the heatwave niche
DEG.MAT.DF = completeFun(DEG.MAT.RASTER, "Diff_max_MAT")
dim(DEG.MAT.DF)
DEG.MAT.NICHE = niche_estimate(DF = DEG.MAT.DF, colname = "Diff_max_MAT")
dim(DEG.MAT.NICHE) ## Just the Australian species, because we raster::extracted for only an AUS raster


## Estimate the Monthly 1976-2005 mean daily maximum temperature niche
MAX.TMAX.DF = completeFun(MAX.TMAX.RASTER, "Max_tmax")
dim(MAX.TMAX.DF)
MAX.TMAX.NICHE = niche_estimate(DF = MAX.TMAX.DF, colname = "Max_tmax")
dim(MAX.TMAX.NICHE) ## Just the Australian species, because we raster::extracted for only an AUS raster


## Estimate the Monthly 1976-2005 mean daily maximum temperature niche
# AI.DF = completeFun(AI.RASTER, "Aridity_index")
# dim(AI.DF)
# AI.NICHE = niche_estimate(DF = AI.DF, colname = "Aridity_index")
# dim(AI.NICHE) ## Just the Australian species, because we raster::extracted for only an AUS raster


## Now merge the full niche df with the heatwave niche df
TRAIT.NICHE.RISK = merge(TRAIT.NICHE.AUS,   TRAIT.NICHE.CONTEXT,  by = "searchTaxon")
#TRAIT.NICHE.RISK = merge(TRAIT.NICHE.RISK,  MAX.TMAX.NICHE, by = "searchTaxon")
dim(TRAIT.NICHE.RISK)
names(TRAIT.NICHE.RISK)


#########################################################################################################################
## Heat risk variables
TLP.DR.FREQ    = TRAIT.NICHE.RISK[, c("TLP_MPa", "Drought_freq_extr_max",  "Drought_freq_extr_q95",
                                      "Drought_freq_extr_q95_q05", "Drought_freq_extr_range")]

TLP.DR.DUR     = TRAIT.NICHE.RISK[, c("TLP_MPa", "Drought_max_dur_extr_max",  "Drought_max_dur_extr_q95",
                                      "Drought_max_dur_extr_q95_q05", "Drought_max_dur_extr_range")]

TLP.DR.INT    = TRAIT.NICHE.RISK[, c("TLP_MPa", "Drought_max_int_extr_max",  "Drought_max_int_extr_q95",
                                     "Drought_max_int_extr_q95_q05", "Drought_max_int_extr_range")]

TLP.DR.REL.INT    = TRAIT.NICHE.RISK[, c("TLP_MPa", "Drought_max_rel_int_extr_max",  "Drought_max_rel_int_extr_q95",
                                         "Drought_max_rel_int_extr_q95_q05", "Drought_max_rel_int_extr_range")]

TLP.DR.REL.INT    = TRAIT.NICHE.RISK[, c("TLP_MPa", "Drought_max_rel_int_extr_max",  "Drought_max_rel_int_extr_q95",
                                         "Drought_max_rel_int_extr_q95_q05", "Drought_max_rel_int_extr_range")]




PLOT.DROUGHT    = c("TLP.DR.FREQ", "TLP.DR.DUR", "TLP.DR.INT")


## Rename
names(TLP.DR.FREQ) = c("TLP",   "DRF_MX",  "DRF_95",   "DRF_95_5",   "DRF_RG");
names(TLP.DR.DUR)  = c("TLP",   "DRDU_MX", "DRDU_95",  "DRDU_95_5",  "DRDU_RG");
names(TLP.DR.INT)  = c("TLP",   "DRIN_MX", "DRIN_95",  "DRIN_95_5",  "DRIN_RG")


## Do the values make sense?
summary(TLP.DR.FREQ);
summary(TLP.DR.DUR);
summary(TLP.DR.INT)





#########################################################################################################################
## 3). PLOT GLASSHOUSE TRAITS VS HEAT EXPOSURE NICHES
#########################################################################################################################


#########################################################################################################################
## Plot Turgor loss point vs. 
plot(awap.extreme[["Drought_freq_extr"]], main = "Maximum drought frequency (AWAP, 1900-2011)")


CairoPNG(width = 8090, height = 8090, 
         file = "./output/figures/Glasshouse_niches/TLP_v_AWAP_DROUGHT_FREQ.png", 
         canvas = "white", bg = "white", units = "px", dpi = 600)

par(mar = c(8, 8, 6, 4), 
    mgp = c(6, 2, 0))

par(lwd = 2)

pairs(TLP.DR.FREQ,  
      lower.panel = panel.cor,
      #diag.panel  = panel.hist,
      upper.panel = panel.smooth, #function(...) smoothScatter(..., nrpoints = 0, add = TRUE))
      cex.labels  = 1.5, cex.axis = 2, font.labels = 2,
      main = "Turgor Loss Point (TLP) vs Drought frequency (1900-2011)")

## finish the device
dev.off()


#########################################################################################################################
## Plot leaf turgor loss point vs the AUS max temp - MAT (1950-2000)
plot(awap.extreme[["Drought_max_dur_extr"]], main = "Maximum drought duration (AWAP, 1900-2011)")

CairoPNG(width = 8090, height = 8090, 
         file = "./output/figures/Glasshouse_niches/TLP_v_AWAP_DROUGHT_DURATION.png", 
         canvas = "white", bg = "white", units = "px", dpi = 600)

par(mar = c(8, 8, 6, 4), 
    mgp = c(6, 2, 0))

par(lwd = 2)

pairs(TLP.DR.DUR,  
      lower.panel = panel.cor,
      #diag.panel  = panel.hist,
      upper.panel = panel.smooth, #function(...) smoothScatter(..., nrpoints = 0, add = TRUE))
      cex.labels  = 1.5, cex.axis = 2, font.labels = 2,
      main = "Leaf turgor loss point (TLP) vs Max drought duration (AWAP 1900-2011)")

## finish the device
dev.off()




#########################################################################################################################
## Plot leaf turgor loss point vs the AUS max temp - MAT (1950-2000)
plot(awap.extreme[["Drought_max_int_extr"]], main = "Maximum drought intensity (AWAP, 1900-2011)")
plot(awap.extreme[["Drought_max_rel_int_extr"]], main = "Maximum drought intensity (AWAP, 1900-2011)")

CairoPNG(width = 8090, height = 8090, 
         file = "./output/figures/Glasshouse_niches/TLP_v_AWAP_DROUGHT_INTENSITY.png", 
         canvas = "white", bg = "white", units = "px", dpi = 600)

par(mar = c(8, 8, 6, 4), 
    mgp = c(6, 2, 0))

par(lwd = 2)

pairs(TLP.DR.INT,  
      lower.panel = panel.cor,
      #diag.panel  = panel.hist,
      upper.panel = panel.smooth, #function(...) smoothScatter(..., nrpoints = 0, add = TRUE))
      cex.labels  = 1.5, cex.axis = 2, font.labels = 2,
      main = "Leaf turgor loss point (TLP) vs Max drought intensity (AWAP 1900-2011)")

pairs(TLP.DR.REL.INT ,  
      lower.panel = panel.cor,
      #diag.panel  = panel.hist,
      upper.panel = panel.smooth, #function(...) smoothScatter(..., nrpoints = 0, add = TRUE))
      cex.labels  = 1.5, cex.axis = 2, font.labels = 2,
      main = "Leaf turgor loss point (TLP) vs Max drought relative intensity (AWAP 1900-2011)")


## finish the device
dev.off()





#########################################################################################################################
## 4). LOOP TO PLOT VARIABLES
#########################################################################################################################


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