#########################################################################################################################
## CREATE A PANEL OF SCATTERPLOTS FOR PREDICTED SPECIES GAINS AND LOSSES
#########################################################################################################################


#########################################################################################################################
## Create the gain and loss variables for plotting
SUA.GAIN.2070                = subset(SUA.PLOT.70.M, AREA_CHANGE %in% c("GAIN"))
SUA.LOSS.2070                = subset(SUA.PLOT.70.M, AREA_CHANGE %in% c("LOSS"))

SUA.GAIN.2070$SPECIES_GAIN   = SUA.GAIN.2070$SPECIES_COUNT/(SUA.70.M.LOSS$SPECIES_COUNT + 
                                                              SUA.70.M.STABLE$SPECIES_COUNT + 
                                                              SUA.70.M.GAIN$SPECIES_COUNT)

SUA.LOSS.2070$SPECIES_LOSS   = -SUA.LOSS.2070$SPECIES_COUNT/(SUA.70.M.LOSS$SPECIES_COUNT + 
                                                               SUA.70.M.STABLE$SPECIES_COUNT + 
                                                               SUA.70.M.GAIN$SPECIES_COUNT)

SUA.GAIN.2070                = SUA.GAIN.2070[c("SPECIES_GAIN", SUA_ORDER)]
SUA.LOSS.2070                = SUA.LOSS.2070[c("SPECIES_LOSS", SUA_ORDER)]


#########################################################################################################################
## Run GAMs of species gains vs MAT for 2070 : eval(parse(text = SUA_ORDER))
GAIN.2070.TEST  = data.frame(x = seq(min(SUA.GAIN.2070[, SUA_ORDER]),
                                     max(SUA.GAIN.2070[, SUA_ORDER]), 
                                     length = length(SUA.GAIN.2070[["SPECIES_GAIN"]])))
colnames(GAIN.2070.TEST) = SUA_ORDER

PRED.GAIN.2070 = predict(GAIN.2070.GAM, newdata = GAIN.2070.TEST, type ='response')


#########################################################################################################################
## Run GAMs of species losses vs MAT for 2070
LOSS.2070.GAM = gam(SPECIES_LOSS ~ s (CURRENT_MAT, k = 5), 
                    data = SUA.LOSS.2070, 
                    method = "REML")
summary(LOSS.2070.GAM)[["dev.expl"]]               

LOSS.2070.TEST  = data.frame(CURRENT_MAT = seq(min(SUA.LOSS.2070[,"CURRENT_MAT"]),
                                               max(SUA.LOSS.2070[,"CURRENT_MAT"]), 
                                               length = length(SUA.LOSS.2070[["SPECIES_LOSS"]])))

PRED.LOSS.2070 = predict(LOSS.2070.GAM, newdata = LOSS.2070.TEST, type ='response')


#########################################################################################################################
## Create PNG
CairoPNG(width = 13000, height = 10000, 
         file = sprintf('output/figures/SUA_percent/SUA_2070_SPP_GAIN_vs_MAT_GAM_%s_%s.png', SUAs, SUA_SPP),
         canvas = "white", bg = "white", units = "px", dpi = 600)

## Add mfrow
par(mar   = c(13.5, 16, 4, 4.8), 
    mgp   = c(11.8, 3, 0),
    oma   = c(1, 1, 1, 1))


#################################################################
## PANEL 1 :: GAMs of species gains vs MAT for 2070
par(font.axis = 1)

plot(SUA.GAIN.2070[,"CURRENT_MAT"], SUA.GAIN.2070[,"SPECIES_GAIN"], 
     col = alpha("blue", 0.3), pch = 19, cex = 4, 
     cex.axis = 3, cex.lab = 5,
     las = 1, ylab = "Species gained (% current - 2070)", xlab = "Current MAT of SUA (1960-1990)")

box(lwd = 3)

lines(GAIN.2070.TEST$CURRENT_MAT, PRED.GAIN.2070, col = "orange",  lwd = 8)

legend("topright", bty = "n", cex = 5, pt.cex = 5, 
       text.col = "orange", legend = paste0("DE = ", format(summary(GAIN.2070.GAM)$dev.expl *100, digits = 3), "%"))

dev.off()	   


#################################################################
## PANEL 2 :: GAMs of species losses vs MAT for 2070
CairoPNG(width = 13000, height = 10000, 
         file = sprintf('output/figures/SUA_percent/SUA_2070_SPP_LOSS_vs_MAT_GAM_%s_%s.png', SUAs, SUA_SPP),
         canvas = "white", bg = "white", units = "px", dpi = 600)

## Add mfrow
par(mar   = c(13.5, 16, 4, 4.8), 
    mgp   = c(11.8, 3, 0),
    oma   = c(1, 1, 1, 1))

par(font.axis = 1)

plot(SUA.LOSS.2070[,"CURRENT_MAT"], SUA.LOSS.2070[,"SPECIES_LOSS"], 
     col = alpha("blue", 0.3), pch = 19, cex = 4, 
     cex.axis = 3, cex.lab = 5,
     las = 1, ylab = "Species lost (% current - 2070)", xlab = "Current MAT of SUA (1960-1990)")

box(lwd = 3)

lines(LOSS.2070.TEST$CURRENT_MAT, PRED.LOSS.2070, col = "orange",  lwd = 8)

legend("topleft", bty = "n", cex = 5, pt.cex = 5, 
       text.col = "orange", legend = paste0("DE = ", format(summary(LOSS.2070.GAM)$dev.expl *100, digits = 3), "%"))

dev.off()	