#########################################################################################################################
################################# CREATE GAIN/LOSS BARPLOTS FOR SUAS ####################################################
#########################################################################################################################


## This code creates a barplot of the % of cells in each SUA that are gained/lost/stable.
## So the Y axis is the % change in suitable habitat. Error bars would be for the species.


#########################################################################################################################
## 1). AGGREGATE SPECIES CELLS LOST/GAINED BY SPATIAL UNIT (e.g. SUA)
#########################################################################################################################


#########################################################################################################################
## Read in data and aggregate by SUA.
SUA.LGS           = SUA.PLOT.GOOD
PERIOD            = 70
SUA.LGS.PERIOD    = SUA.LGS[SUA.LGS$PERIOD == PERIOD, ]
unique(SUA.LGS.PERIOD$PERIOD)

SUA.LGS.PERIOD    = SUA.LGS.PERIOD[c("SPECIES", "SUA", "LOST", "GAINED", "STABLE", "NEVER", "CELL_COUNT", "PERIOD")]
SPP.LGS.PERIOD    = SUA.LGS.PERIOD


#########################################################################################################################
## First Calcualte the L/G/S % using dis-aggregated data
SPP.LGS.PERIOD$LOSS_CHANGE   = SPP.LGS.PERIOD$LOST/(SPP.LGS.PERIOD$CELL_COUNT)   * 100
SPP.LGS.PERIOD$GAIN_CHANGE   = SPP.LGS.PERIOD$GAINED/(SPP.LGS.PERIOD$CELL_COUNT) * 100
SPP.LGS.PERIOD$STABLE_CHANGE = SPP.LGS.PERIOD$STABLE/(SPP.LGS.PERIOD$CELL_COUNT) * 100
SPP.LGS.PERIOD$NEVER_CHANGE  = SPP.LGS.PERIOD$NEVER/(SPP.LGS.PERIOD$CELL_COUNT)  * 100

SPP.LGS.PERIOD$SUA           = as.factor(SPP.LGS.PERIOD$SUA)
SPP.LGS.PERIOD$SPECIES       = as.factor(SPP.LGS.PERIOD$SPECIES)
SPP.LGS.PERIOD.M             = melt(SPP.LGS.PERIOD)
SPP.LGS.PERIOD.M             = subset(SPP.LGS.PERIOD.M, variable %in% c("LOSS_CHANGE", "GAIN_CHANGE", "STABLE_CHANGE"))


#########################################################################################################################
## Then aggregate the data by SUA - alculating the variance might not make sense
SUA.LGS.CHANGE.PERIOD = aggregate(. ~ SUA, data = SPP.LGS.PERIOD, sum, na.rm = TRUE)
SUA.LGS.CHANGE.PERIOD = SUA.LGS.CHANGE.PERIOD[c("SUA", "LOST", "GAINED", "STABLE", "NEVER", "CELL_COUNT", "PERIOD")]


## Calculate species change
SUA.LGS.CHANGE.PERIOD$LOSS_CHANGE   = SUA.LGS.CHANGE.PERIOD$LOST/(SUA.LGS.CHANGE.PERIOD$CELL_COUNT)   * 100
SUA.LGS.CHANGE.PERIOD$GAIN_CHANGE   = SUA.LGS.CHANGE.PERIOD$GAINED/(SUA.LGS.CHANGE.PERIOD$CELL_COUNT) * 100
SUA.LGS.CHANGE.PERIOD$STABLE_CHANGE = SUA.LGS.CHANGE.PERIOD$STABLE/(SUA.LGS.CHANGE.PERIOD$CELL_COUNT) * 100
SUA.LGS.CHANGE.PERIOD$NEVER_CHANGE  = SUA.LGS.CHANGE.PERIOD$NEVER/(SUA.LGS.CHANGE.PERIOD$CELL_COUNT)  * 100


## Just check the aggregate works
## We are effectively taking the average % of cells lost/gained/stable in each SUA across all species
View(SPP.LGS.PERIOD)
View(SUA.LGS.CHANGE.PERIOD)
mean(subset(SPP.LGS.PERIOD, SUA == "Adelaide")$LOSS_CHANGE)
SUA.LGS.CHANGE.PERIOD.M             = melt(SUA.LGS.CHANGE.PERIOD)
SUA.VAR                             = subset(SUA.LGS.CHANGE.PERIOD.M, variable %in% c("LOSS_CHANGE", "GAIN_CHANGE", "STABLE_CHANGE"))





#########################################################################################################################
## 2). PLOT LOSS/GAIN PATTERNS FOR EACH CAPITAL CITY
#########################################################################################################################


#########################################################################################################################
## Now create a stacked barplot for each capital city
capital.cities = c("Adelaide", "Perth", "Brisbane", "Sydney", "Canberra - Queanbeyan", "Melbourne", "Hobart")
city           = capital.cities[1]


## Make all the figures in the MS colorblind safe. The "Paried" color scheme is aparently safe.
## No. 1 = light blue, 4 = green, and 8 = dark orange 
SUA.plot.cols = brewer.pal(12, "Paired")
city = capital.cities[1]


for (city in capital.cities) { 
  
  ## Subset bigger plot to create error bars.
  ## Unfortunately the error bars are too big
  plot.SUA = subset(SUA.VAR, SUA == city)
  plot.L   = subset(SUA.LGS.CHANGE.PERIOD.M, SUA.LGS.CHANGE.PERIOD.M$SUA == city & variable == "LOSS_CHANGE")
  plot.G   = subset(SUA.LGS.CHANGE.PERIOD.M, SUA.LGS.CHANGE.PERIOD.M$SUA == city & variable == "GAIN_CHANGE")
  plot.S   = subset(SUA.LGS.CHANGE.PERIOD.M, SUA.LGS.CHANGE.PERIOD.M$SUA == city & variable == "STABLE_CHANGE")
  
  L.SD     = aggregate(. ~ SUA, data = plot.L, sd, na.rm = TRUE)[c("value")]
  G.SD     = aggregate(. ~ SUA, data = plot.G, sd, na.rm = TRUE)[c("value")]
  S.SD     = aggregate(. ~ SUA, data = plot.S, sd, na.rm = TRUE)[c("value")]
  
  ## Now create the barplot - 
  plot.SUA$SD = as.numeric(c(L.SD, G.SD, S.SD))
  png(sprintf("./output/figures/FIG_1/%s_LGS_BAR_PLOT_20%s.png", city, PERIOD),
      6, 6, units = 'in', res = 500)
  
  barplot = ggplot(plot.SUA, aes(x = SUA, y = value)) +
    geom_bar(aes(fill = as.factor(variable)), stat = "identity", position = "dodge") +
    xlab("") + 
    ylab("% Change") +
    
    scale_fill_manual(values = rev(colorRampPalette(c(SUA.plot.cols[1], SUA.plot.cols[4], SUA.plot.cols[8]))(3))) +
    
    ## Error bars are huge, either miscalcualte or misleading
    # geom_text(aes(label = variable), vjust = 0.5) +
    # geom_errorbar(aes(x = SUA, ymin = value - SD, ymax = value + SD), 
    #               width = 0.4, colour = "orange", alpha = 0.9, size = 1.3) +
    
    ## Add themes
    theme(axis.title.x     = element_blank(),
          axis.text.x      = element_blank(),
          axis.ticks.x     = element_blank(),
          axis.text.y      = element_text(vjust = 0.5, size = 35),
          axis.title.y     = element_text(face = "bold", colour = "black", size = 40),
          panel.background = element_blank(),
          panel.border     = element_rect(colour = "black", fill = NA, size = 3),
          legend.position  = "none")
  
  ## Print the plot and close the device 
  print(barplot)
  dev.off()
  
}





#########################################################################################################################
## 3). PLOT CLIMOGRAPH FOR EACH CAPITAL CITY
#########################################################################################################################


## Try using climogragh
# names(MaquehueTemuco)
# 
# 
# ## Get climograph data
# Syd.r  = get_historical(stationid = "066062", type = "rain") 
# Syd.tm = get_historical(stationid = "066062", type = "min")
# Syd.tx = get_historical(stationid = "066062", type = "max")
# 
# 
# Syd.r  = transform(Syd.r,  Index = as.numeric(paste0(Year, Month, Day)))#  sep = "-"))
# Syd.tm = transform(Syd.tm, Index = paste0(Year, Month, Day))#  sep = "-"))
# Syd.tx = transform(Syd.tx, Index = paste0(Year, Month, Day))#  sep = "-"))
# 
# Syd.r  = Syd.r[c("Index",  "Rainfall")]
# Syd.tm = Syd.tm[c("Index", "Min_temperature")]
# Syd.tx = Syd.tx[c("Index", "Max_temperature")]
# 
# 
# Syd.clim = Reduce(function(x,y) merge(x = x, y = y, by = "Index", all = FALSE), 
#                   list(Syd.r, Syd.tm, Syd.tx))
# Syd.clim$Index <- as.Date(as.character(Syd.clim$Index), format="%m/%d/%Y") 
# Syd.clim = as.zoo(Syd.clim)
# 
# #Syd.clim = Syd.clim[complete.cases(Syd.clim), ]
# 
# 
# t = data(MaquehueTemuco)
# pcp <- MaquehueTemuco[, 1]
# tmx <- MaquehueTemuco[, 2]
# tmn <- MaquehueTemuco[, 3]
# 
# pcp <- Syd.clim[, 2]
# tmx <- Syd.clim[, 4]
# tmn <- Syd.clim[, 3]
# 
# 
# ## Plotting the climograph
# m <- climograph(pcp=pcp, tmx=tmx, tmn=tmn, na.rm=TRUE)
# m <- climograph(pcp = Syd.clim[, 2], tmx = Syd.clim[, 4], tmn = Syd.clim[, 3], na.rm = TRUE)
# 
# 
# ## 
# for (city in capital.cities) { 
#   
#   ## Subset bigger plot to create error bars.
#   ## Unfortunately the error bars are too big
#   plot.SUA = subset(SUA.VAR, SUA == city)
#   plot.L   = subset(SUA.LGS.CHANGE.PERIOD.M, SUA.LGS.CHANGE.PERIOD.M$SUA == city & variable == "LOSS_CHANGE")
#   plot.G   = subset(SUA.LGS.CHANGE.PERIOD.M, SUA.LGS.CHANGE.PERIOD.M$SUA == city & variable == "GAIN_CHANGE")
#   plot.S   = subset(SUA.LGS.CHANGE.PERIOD.M, SUA.LGS.CHANGE.PERIOD.M$SUA == city & variable == "STABLE_CHANGE")
#   
#   L.SD     = aggregate(. ~ SUA, data = plot.L, sd, na.rm = TRUE)[c("value")]
#   G.SD     = aggregate(. ~ SUA, data = plot.G, sd, na.rm = TRUE)[c("value")]
#   S.SD     = aggregate(. ~ SUA, data = plot.S, sd, na.rm = TRUE)[c("value")]
#   
#   ## Now create the barplot - 
#   plot.SUA$SD = as.numeric(c(L.SD, G.SD, S.SD))
#   png(sprintf("./output/figures/FIG_1/%s_LGS_BAR_PLOT_20%s.png", city, PERIOD),
#       6, 6, units = 'in', res = 500)
#   
#   barplot = ggplot(plot.SUA, aes(x = SUA, y = value)) +
#     geom_bar(aes(fill = as.factor(variable)), stat = "identity", position = "dodge") +
#     xlab("") + 
#     ylab("% Change") +
#     
#     scale_fill_manual(values = rev(colorRampPalette(c(SUA.plot.cols[1], SUA.plot.cols[4], SUA.plot.cols[8]))(3))) +
#     
#     ## Error bars are huge, either miscalcualte or misleading
#     # geom_text(aes(label = variable), vjust = 0.5) +
#     # geom_errorbar(aes(x = SUA, ymin = value - SD, ymax = value + SD), 
#     #               width = 0.4, colour = "orange", alpha = 0.9, size = 1.3) +
#     
#     ## Add themes
#     theme(axis.title.x     = element_blank(),
#           axis.text.x      = element_blank(),
#           axis.ticks.x     = element_blank(),
#           axis.text.y      = element_text(vjust = 0.5, size = 35),
#           axis.title.y     = element_text(face = "bold", colour = "black", size = 40),
#           panel.background = element_blank(),
#           panel.border     = element_rect(colour = "black", fill = NA, size = 3),
#           legend.position  = "none")
#   
#   ## Print the plot and close the device 
#   print(barplot)
#   dev.off()
#   
# }

