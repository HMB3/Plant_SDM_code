#########################################################################################################################
## SUBSET DATA TO MAKE PLOTTING EASIER
#########################################################################################################################


#########################################################################################################################
## Create plots in the format needed for ggplot bar  
SUA.PLOT.GOOD.30 = subset(SUA.PLOT.GOOD, PERIOD == 30)
SUA.PLOT.GOOD.50 = subset(SUA.PLOT.GOOD, PERIOD == 50)
SUA.PLOT.GOOD.70 = subset(SUA.PLOT.GOOD, PERIOD == 70)
unique(SUA.PLOT.GOOD.30$PERIOD);unique(SUA.PLOT.GOOD.50$PERIOD);unique(SUA.PLOT.GOOD.70$PERIOD)
dim(SUA.PLOT.GOOD.30);dim(SUA.PLOT.GOOD.50);dim(SUA.PLOT.GOOD.70)


## Melt the table into the right format : but does this mean the different categories are mutually exclusive?
## Species can only fall in the categories gain, loss, stable, never, in each SUA. So yes, they are exclusive.
SUA.PLOT.30          = table(SUA.PLOT.GOOD.30$SUA, SUA.PLOT.GOOD.30$GAIN_LOSS)
SUA.PLOT.50          = table(SUA.PLOT.GOOD.50$SUA, SUA.PLOT.GOOD.50$GAIN_LOSS)
SUA.PLOT.70          = table(SUA.PLOT.GOOD.70$SUA, SUA.PLOT.GOOD.70$GAIN_LOSS)
SUA.PLOT.30.M        = melt(SUA.PLOT.30)
SUA.PLOT.50.M        = melt(SUA.PLOT.50)
SUA.PLOT.70.M        = melt(SUA.PLOT.70)

names(SUA.PLOT.30.M) = c("SUA", "AREA_CHANGE", "SPECIES_COUNT")
names(SUA.PLOT.50.M) = c("SUA", "AREA_CHANGE", "SPECIES_COUNT")
names(SUA.PLOT.70.M) = c("SUA", "AREA_CHANGE", "SPECIES_COUNT")
SUA.PLOT.30.M$PERIOD = 30 
SUA.PLOT.50.M$PERIOD = 50 
SUA.PLOT.70.M$PERIOD = 70 


SUA.PLOT.30.M        = subset(SUA.PLOT.30.M, AREA_CHANGE != "NEVER")# & AREA_CHANGE != "STABLE")
SUA.PLOT.50.M        = subset(SUA.PLOT.50.M, AREA_CHANGE != "NEVER")# & AREA_CHANGE != "STABLE")
SUA.PLOT.70.M        = subset(SUA.PLOT.70.M, AREA_CHANGE != "NEVER")# & AREA_CHANGE != "STABLE") 

SUA.30.M.LOSS        = subset(SUA.PLOT.30.M, AREA_CHANGE %in% c("LOSS"))
SUA.30.M.GAIN        = subset(SUA.PLOT.30.M, AREA_CHANGE %in% c("GAIN"))
SUA.30.M.STABLE      = subset(SUA.PLOT.30.M, AREA_CHANGE %in% c("STABLE"))

SUA.50.M.LOSS        = subset(SUA.PLOT.50.M, AREA_CHANGE %in% c("LOSS"))
SUA.50.M.GAIN        = subset(SUA.PLOT.50.M, AREA_CHANGE %in% c("GAIN"))
SUA.50.M.STABLE      = subset(SUA.PLOT.50.M, AREA_CHANGE %in% c("STABLE"))

SUA.70.M.LOSS        = subset(SUA.PLOT.70.M, AREA_CHANGE %in% c("LOSS"))
SUA.70.M.GAIN        = subset(SUA.PLOT.70.M, AREA_CHANGE %in% c("GAIN"))
SUA.70.M.STABLE      = subset(SUA.PLOT.70.M, AREA_CHANGE %in% c("STABLE"))

head(SUA.30.M.LOSS)
head(SUA.30.M.GAIN)
head(SUA.30.M.STABLE)


#########################################################################################################################
## Attach the climate
SUA.CLIM      = SUA.PREDICT[!duplicated(SUA.PLOT.GOOD[,c('SUA')]),][c("SUA", "CURRENT_MAT", "CURRENT_MAP", 
                                                                      "CURRENT_PET", "CURRENT_AI", "CURRENT_MAXT", 
                                                                      "AREASQKM16", "ClimateZ", "POP_2017")]


## Find a more efficient way to join everything on to the subsets
SUA.PLOT.30.M = join(SUA.PLOT.30.M, SUA.CLIM)
SUA.PLOT.50.M = join(SUA.PLOT.50.M, SUA.CLIM)
SUA.PLOT.70.M = join(SUA.PLOT.70.M, SUA.CLIM)

SUA.30.M.LOSS   = join(SUA.30.M.LOSS,   SUA.CLIM)
SUA.30.M.GAIN   = join(SUA.30.M.GAIN,   SUA.CLIM)
SUA.30.M.STABLE = join(SUA.30.M.STABLE, SUA.CLIM)

SUA.70.M.LOSS   = join(SUA.70.M.LOSS,   SUA.CLIM)
SUA.70.M.GAIN   = join(SUA.70.M.GAIN,   SUA.CLIM)
SUA.70.M.STABLE = join(SUA.70.M.STABLE, SUA.CLIM)



#########################################################################################################################
## Use ggplot to plot the percentage of species inside an LGA, which is being lost or gained
## Calculate the gain/loss for counts of all species?

## Option 1).
## Gain % = species gained / (stable + lost)   * 100
## Lost % = species lost   / (stable + lost)   * 100

## Calculate the gain/loss for for only species that are recorded within each SUA?

## Option 2).
## Gain % = species gained / (stable + lost + gained)   * 100
## Lost % = species lost   / (stable + lost + gained)   * 100



#########################################################################################################################
## Run different variations
#"./R/SUA_GAIN_LOSS_ALL_SUAs.R"       ## 1).
#"./R/SUA_%100_PLOT_LARGER_SUAs.R"    ## 2).
#"./R/SUA_RAIN_PLOT.R"                ## 2), using rain not temp





#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################