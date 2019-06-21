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
DAYS.35    = raster("./data/base/AWAP/mean_days_35_degrees.tif")
DEG.MAT    = raster("./data/base/AWAP/degrees_above_MAT.tif")
MAX.TMAX   = raster("./data/base/AWAP/tmaxmax_19502000.tif")
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
crs(DEG.MAT)  <- "+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
crs(MAX.TMAX) <- "+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
# DEG.MAT      <- projectRaster(DEG.MAT,    CRS("+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# MAX.TMAX     <- projectRaster(MAX.TMAX,   CRS("+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
projection(TRAIT.POINTS);projection(DEG.MAT);projection(MAX.TMAX)


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


## Now raster::extract just the mean temperature
# ANUCLIM.TEMP.RASTER <- raster::extract(temp_max_monthly, TRAIT.POINTS) %>% 
#   cbind(TRAIT.RASTER.CONTEXT , .)
# colnames(TEMP.RASTER )[90] <- "Mnt_mn_daily_max_tmp"
# colnames(TEMP.RASTER )[90]


## Now raster::extract just the mean temperature
# AI.RASTER <- raster::extract(AI, TRAIT.POINTS) %>% 
#   cbind(TRAIT.RASTER.CONTEXT , .)
# colnames(AI.RASTER )[90] <- "Aridity_index"
# colnames(AI.RASTER )[90]



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
TRAIT.NICHE.RISK = merge(TRAIT.NICHE.AUS,   DEG.MAT.NICHE,  by = "searchTaxon")
TRAIT.NICHE.RISK = merge(TRAIT.NICHE.RISK,  MAX.TMAX.NICHE, by = "searchTaxon")
dim(TRAIT.NICHE.RISK)
names(TRAIT.NICHE.RISK)


#########################################################################################################################
## Heat risk variables
OSMP.HEAT.PLOT = TRAIT.NICHE.RISK[, c("OsmPot_MPa", "Diff_max_MAT_max",  "Diff_max_MAT_q95",
                                      "Diff_max_MAT_q95_q05", "Diff_max_MAT_range")]

OSMP.TMAX.PLOT = TRAIT.NICHE.RISK[, c("OsmPot_MPa", "Max_tmax_max",  "Max_tmax_q95",
                                      "Max_tmax_q95_q05", "Max_tmax_range")]

TLP.HEAT.PLOT  = TRAIT.NICHE.RISK[, c("TLP_MPa", "Diff_max_MAT_max",  "Diff_max_MAT_q95",
                                    "Diff_max_MAT_q95_q05", "Diff_max_MAT_range")]

TLP.TMAX.PLOT  = TRAIT.NICHE.RISK[, c("TLP_MPa", "Max_tmax_max",  "Max_tmax_q95",
                                      "Max_tmax_q95_q05", "Max_tmax_range")]

TCR.HEAT.PLOT  = TRAIT.NICHE.RISK[, c("Tcrit_C", "Diff_max_MAT_max",  "Diff_max_MAT_q95",
                                    "Diff_max_MAT_q95_q05", "Diff_max_MAT_range")]

TCR.TMAX.PLOT  = TRAIT.NICHE.RISK[, c("Tcrit_C", "Max_tmax_max",  "Max_tmax_q95",
                                      "Max_tmax_q95_q05", "Max_tmax_range")]

DSM.HEAT.PLOT  = TRAIT.NICHE.RISK[, c("Desicc_Max", "Diff_max_MAT_max",  "Diff_max_MAT_q95",
                                    "Diff_max_MAT_q95_q05", "Diff_max_MAT_range")]

DSM.TMAX.PLOT  = TRAIT.NICHE.RISK[, c("Desicc_Max", "Max_tmax_max",  "Max_tmax_q95",
                                      "Max_tmax_q95_q05", "Max_tmax_range")]


## Rename
names(OSMP.HEAT.PLOT) = c("OSMPO", "HEAT_MX", "HEAT_95", "HEAT_95_5", "HEAT_RG");
names(OSMP.TMAX.PLOT) = c("OSMPO", "TMAX_MX", "TMAX_95", "TMAX_95_5", "TMAX_RG");
names(TLP.HEAT.PLOT)  = c("TLP",   "HEAT_MX", "HEAT_95", "HEAT_95_5", "HEAT_RG");
names(TLP.TMAX.PLOT)  = c("TLP",   "TMAX_MX", "TMAX_95", "TMAX_95_5", "TMAX_RG");
names(TCR.HEAT.PLOT)  = c("Tcrit", "HEAT_MX", "HEAT_95", "HEAT_95_5", "HEAT_RG");
names(TCR.TMAX.PLOT)  = c("Tcrit", "TMAX_MX", "TMAX_95", "TMAX_95_5", "TMAX_RG");
names(DSM.HEAT.PLOT)  = c("DSM",   "HEAT_MX", "HEAT_95", "HEAT_95_5", "HEAT_RG")
names(DSM.TMAX.PLOT)  = c("DSM",   "TMAX_MX", "TMAX_95", "TMAX_95_5", "TMAX_RG")


## Do the values make sense?
summary(TCR.HEAT.PLOT)
summary(TCR.TMAX.PLOT)





#########################################################################################################################
## 3). PLOT GLASSHOUSE TRAITS VS HEAT EXPOSURE NICHES
#########################################################################################################################


## Plot rasters to check : looks ok?
plot(DEG.MAT,  main = "AUS max temp - MAT (1950-2000)")


#########################################################################################################################
## Plot leaf critical temperature (Tcrit) vs the AUS max temp - MAT (1950-2000)
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
      main = "Leaf critical temperature (Tcrit) vs AUS max temp - MAT (1950-2000)")

## finish the device
dev.off()


#########################################################################################################################
## Plot leaf turgor loss point vs the AUS max temp - MAT (1950-2000)
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
      main = "Leaf turgor loss point (TLP) vs AUS max temp - MAT (1950-2000)")

## finish the device
dev.off()


## Plot max of max
plot(MAX.TMAX, main = "AUS Max tmax (1950-2000)")


#########################################################################################################################
## Plot leaf critical temperature (Tcrit) vs the AUS Max tmax (1950-2000)
CairoPNG(width = 8090, height = 8090, 
         file = "./output/figures/Glasshouse_niches/OSMOPOT_v_WORLDCLIM.png", 
         canvas = "white", bg = "white", units = "px", dpi = 600)

par(mar = c(8, 8, 6, 4), 
    mgp = c(6, 2, 0))

par(lwd = 2)

pairs(TCR.TMAX.PLOT,  
      lower.panel = panel.cor,
      #diag.panel  = panel.hist,
      upper.panel = panel.smooth, #function(...) smoothScatter(..., nrpoints = 0, add = TRUE))
      cex.labels  = 1.5, cex.axis = 2, font.labels = 2,
      main = "Leaf critical temperature (Tcrit) vs AUS Max tmax (1950-2000)")

## finish the device
dev.off()


#########################################################################################################################
## Plot leaf turgor loss point vs the AUS Max tmax (1950-2000)
CairoPNG(width = 8090, height = 8090, 
         file = "./output/figures/Glasshouse_niches/TLP_v_WORLDCLIM.png", 
         canvas = "white", bg = "white", units = "px", dpi = 600)

par(mar = c(8, 8, 6, 4), 
    mgp = c(6, 2, 0))

par(lwd = 2)

pairs(TLP.TMAX.PLOT,  
      lower.panel = panel.cor,
      #diag.panel  = panel.hist,
      upper.panel = panel.smooth, #function(...) smoothScatter(..., nrpoints = 0, add = TRUE))
      cex.labels  = 1.5, cex.axis = 2, font.labels = 2,
      main = "Leaf turgor loss point (TLP) vs AUS Max tmax (1950-2000)")

## finish the device
dev.off()



#########################################################################################################################
## Save the R object
save.image("GLASS_HOUSE_VS_NICHE.RData")





#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################