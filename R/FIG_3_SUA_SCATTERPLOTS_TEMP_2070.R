#########################################################################################################################
## CREATE A PANEL OF SCATTERPLOTS FOR PREDICTED SPECIES GAINS AND LOSSES
#########################################################################################################################


## This code creates scatterplots of species gained/lost in each of the 101 SUAs
## Check that barplots and scatterplots are working with the same number of SUAs?


#########################################################################################################################
## Create the gain and loss variables for plotting.
SUA.GAIN.2070                   = subset(SUA.PLOT.70.M, AREA_CHANGE %in% c("GAIN"))
SUA.LOSS.2070                   = subset(SUA.PLOT.70.M, AREA_CHANGE %in% c("LOSS"))
SUA.STABLE.2070                 = subset(SUA.PLOT.70.M, AREA_CHANGE %in% c("STABLE"))

SUA.GAIN.2070.CAP               = subset(SUA.GAIN.2070,   AREASQKM16 > 200 & POP_2017 > 80000)
SUA.LOSS.2070.CAP               = subset(SUA.LOSS.2070,   AREASQKM16 > 200 & POP_2017 > 80000)
SUA.STABLE.2070.CAP             = subset(SUA.STABLE.2070, AREASQKM16 > 200 & POP_2017 > 80000)


## How many SUA's are temperate
length(unique(SUA.GAIN.2070$SUA))


#########################################################################################################################
## All SUAs
SUA.GAIN.2070$SPECIES_GAIN      = SUA.GAIN.2070$SPECIES_COUNT/(SUA.LOSS.2070$SPECIES_COUNT + 
                                                                 SUA.STABLE.2070$SPECIES_COUNT + 
                                                                 SUA.GAIN.2070$SPECIES_COUNT)
SUA.GAIN.2070$SPECIES_GAIN      = SUA.GAIN.2070$SPECIES_GAIN*100

SUA.LOSS.2070$SPECIES_LOSS      = SUA.LOSS.2070$SPECIES_COUNT/(SUA.LOSS.2070$SPECIES_COUNT + 
                                                                 SUA.STABLE.2070$SPECIES_COUNT + 
                                                                 SUA.GAIN.2070$SPECIES_COUNT)
SUA.LOSS.2070$SPECIES_LOSS      = SUA.LOSS.2070$SPECIES_LOSS*100


#########################################################################################################################
## Big SUAs
SUA.GAIN.2070.CAP$SPECIES_GAIN  = SUA.GAIN.2070.CAP$SPECIES_COUNT/(SUA.LOSS.2070.CAP$SPECIES_COUNT + 
                                                                     SUA.STABLE.2070.CAP$SPECIES_COUNT + 
                                                                     SUA.GAIN.2070.CAP$SPECIES_COUNT)
SUA.GAIN.2070.CAP$SPECIES_GAIN  = SUA.GAIN.2070.CAP$SPECIES_GAIN*100

SUA.LOSS.2070.CAP$SPECIES_LOSS  = SUA.LOSS.2070.CAP$SPECIES_COUNT/(SUA.LOSS.2070.CAP$SPECIES_COUNT + 
                                                                     SUA.STABLE.2070.CAP$SPECIES_COUNT + 
                                                                     SUA.GAIN.2070.CAP$SPECIES_COUNT)
SUA.LOSS.2070.CAP$SPECIES_LOSS   = SUA.LOSS.2070.CAP$SPECIES_LOSS*100





#########################################################################################################################
## 1). GAMS of GAIN/LOSS
#########################################################################################################################


#########################################################################################################################
## Repeat from here
SUA.GAIN.2070.GAM                = SUA.GAIN.2070[c("SPECIES_GAIN", "CURRENT_MAT")]
SUA.LOSS.2070.GAM                = SUA.LOSS.2070[c("SPECIES_LOSS", "CURRENT_MAT")]

SUA.GAIN.2070.GAM.CAP            = SUA.GAIN.2070.CAP[c("SPECIES_GAIN", "CURRENT_MAT")]
SUA.LOSS.2070.GAM.CAP            = SUA.LOSS.2070.CAP[c("SPECIES_LOSS", "CURRENT_MAT")]


#########################################################################################################################
## Remove SUA's with < 2 species as their count
message('Remove SUAs with a species count < 2')
SUA.GAIN.2070     = subset(SUA.GAIN.2070, SPECIES_GAIN >= 2)
SUA.LOSS.2070     = subset(SUA.LOSS.2070, SPECIES_LOSS >= 2)

SUA.GAIN.2070.GAM = subset(SUA.GAIN.2070.GAM, SPECIES_GAIN >= 2)
SUA.LOSS.2070.GAM = subset(SUA.LOSS.2070.GAM, SPECIES_LOSS >= 2)

SUA.GAIN.2070.CAP = subset(SUA.GAIN.2070.CAP, SPECIES_GAIN >= 2)
SUA.LOSS.2070.CAP = subset(SUA.LOSS.2070.CAP, SPECIES_LOSS >= 2)

SUA.GAIN.2070.GAM.CAP = subset(SUA.GAIN.2070.GAM.CAP, SPECIES_GAIN >= 2)
SUA.LOSS.2070.GAM.CAP = subset(SUA.LOSS.2070.GAM.CAP, SPECIES_LOSS >= 2)


  
#########################################################################################################################
## Run GAMs of species gains vs MAT for 2070, ALL SUA
GAIN.2070.GAM = gam(SPECIES_GAIN ~ s (CURRENT_MAT, k = 5), 
                    data = SUA.GAIN.2070.GAM, 
                    method = "REML")
summary(GAIN.2070.GAM)[["dev.expl"]]               

GAIN.2070.TEST  = data.frame(CURRENT_MAT = seq(min(SUA.GAIN.2070[,"CURRENT_MAT"]),
                                               max(SUA.GAIN.2070[,"CURRENT_MAT"]), 
                                               length = length(SUA.GAIN.2070[["SPECIES_GAIN"]])))
                             
                            # ClimateZ = SUA.GAIN.2070[,"ClimateZ"])

PRED.GAIN.2070 = predict(GAIN.2070.GAM, newdata = GAIN.2070.TEST, type ='response')


#########################################################################################################################
## Run GAMs of species gains vs MAT for 2070, BIG SUA
GAIN.2070.GAM.CAP = gam(SPECIES_GAIN ~ s (CURRENT_MAT, k = 5), 
                    data = SUA.GAIN.2070.GAM.CAP, 
                    method = "REML")
summary(GAIN.2070.GAM.CAP)[["dev.expl"]]               

GAIN.2070.TEST.CAP  = data.frame(CURRENT_MAT = seq(min(SUA.GAIN.2070.CAP[,"CURRENT_MAT"]),
                                               max(SUA.GAIN.2070.CAP[,"CURRENT_MAT"]), 
                                               length = length(SUA.GAIN.2070.CAP[["SPECIES_GAIN"]])))

PRED.GAIN.2070.CAP = predict(GAIN.2070.GAM.CAP, newdata = GAIN.2070.TEST.CAP, type ='response')


#########################################################################################################################
## Run GAMs of species losses vs MAT for 2070
LOSS.2070.GAM = gam(SPECIES_LOSS ~ s (CURRENT_MAT, k = 5), 
                    data = SUA.LOSS.2070.GAM, 
                    method = "REML")
summary(LOSS.2070.GAM)[["dev.expl"]]               

LOSS.2070.TEST  = data.frame(CURRENT_MAT = seq(min(SUA.LOSS.2070[,"CURRENT_MAT"]),
                                               max(SUA.LOSS.2070[,"CURRENT_MAT"]), 
                                               length = length(SUA.LOSS.2070[["SPECIES_LOSS"]])))

PRED.LOSS.2070 = predict(LOSS.2070.GAM, newdata = LOSS.2070.TEST, type ='response')


#########################################################################################################################
## Run GAMs of species losses vs MAT for 2070, BIG SUAs
LOSS.2070.GAM.CAP = gam(SPECIES_LOSS ~ s (CURRENT_MAT, k = 5), 
                    data = SUA.LOSS.2070.GAM.CAP, 
                    method = "REML")
summary(LOSS.2070.GAM.CAP)[["dev.expl"]]               

LOSS.2070.TEST.CAP  = data.frame(CURRENT_MAT = seq(min(SUA.LOSS.2070.CAP[,"CURRENT_MAT"]),
                                               max(SUA.LOSS.2070.CAP[,"CURRENT_MAT"]), 
                                               length = length(SUA.LOSS.2070.CAP[["SPECIES_LOSS"]])))

PRED.LOSS.2070.CAP = predict(LOSS.2070.GAM.CAP, newdata = LOSS.2070.TEST.CAP, type ='response')





#########################################################################################################################
## 2). SCATTERPLOT FOR GAMS of GAIN/LOSS, WITH COLORS FOR KOPPEN ZONE
#########################################################################################################################


#########################################################################################################################
## First, create a data frame of colbrewer colors :: display.brewer.all()
sua.col <-
   with(SUA.GAIN.2070,
        data.frame(ClimateZ = unique(SUA.GAIN.2070$ClimateZ),
                   color = I(brewer.pal(length(unique(SUA.GAIN.2070$ClimateZ)), 
                                        name = 'Set2'))))

SUA.GAIN.2070     <- merge(SUA.GAIN.2070,     sua.col)
SUA.GAIN.2070.CAP <- merge(SUA.GAIN.2070.CAP, sua.col)

SUA.LOSS.2070     <- merge(SUA.LOSS.2070,     sua.col)
SUA.LOSS.2070.CAP <- merge(SUA.LOSS.2070.CAP, sua.col)


#########################################################################################################################
## Now try saving the datframe so Alessandro can re-produce
# SUA.GL.PLOT = bind_rows(SUA.GAIN.2030, SUA.LOSS.2030, SUA.GAIN.2070, SUA.LOSS.2070)
# saveRDS(SUA.GL.PLOT, file =  paste0(DATA_path, 'SUA_GAIN_LOSS_PLOT_',  save_run, '.rds'))
# write.csv(SUA.GAIN.2030[c("SUA", "ClimateZ")], 
#           file =  paste0(DATA_path, 'SUA_TEMPERATE.csv'), row.names = FALSE)


#########################################################################################################################
## Create PNG
CairoPNG(width = 18000, height = 16000, 
         file = sprintf('output/figures/FIG_2/SUA_2070_SPP_GL_PANEL_COLS_vs_%s_GAM_%s_%s_%s.png',  
                        'temp', SUAs, SUA_SPP, KOP_ZONE),
         canvas = "white", bg = "white", units = "px", dpi = 600)


## Add mfrow
par(mfrow = c(2, 2),
    mar   = c(13.5, 16, 4, 4.8), 
    mgp   = c(11.8, 3, 0),
    oma   = c(1, 1, 1, 1),
    xpd   = TRUE)


#################################################################
## PANEL 1 :: GAMs of species gains vs MAT for 2070, ALL SUAs
par(font.axis = 1)

plot(SUA.GAIN.2070[,"CURRENT_MAT"], SUA.GAIN.2070[,"SPECIES_GAIN"], 
     col = SUA.GAIN.2070$color, pch = 19, cex = 6,
     cex.axis = 5, cex.lab = 5,
     las = 1, ylab = "Species gained (%)", # (% current - 2070)",
     xlab = "")

box(lwd = 3)

lines(GAIN.2070.TEST$CURRENT_MAT, PRED.GAIN.2070, col = "blue",  lwd = 10)

legend("topright", bty = "n", cex = 5, pt.cex = 5, 
       text.col = "black", 
       legend = paste0("DE = ", format(summary(GAIN.2070.GAM)$dev.expl *100, digits = 3), "%"))


#################################################################
## PANEL 2 :: GAMs of species gains vs MAT for 2070, BIG SUAs
par(font.axis = 1)

plot(SUA.GAIN.2070.CAP[,"CURRENT_MAT"], SUA.GAIN.2070.CAP[,"SPECIES_GAIN"], 
     col = SUA.GAIN.2070.CAP$color, pch = 19, cex = 6,
     cex.axis = 5, cex.lab = 5,
     las = 1,
     ylab = "", 
     xlab = "")

box(lwd = 3)

lines(GAIN.2070.TEST.CAP$CURRENT_MAT, PRED.GAIN.2070.CAP, col = "blue",  lwd = 10)

legend("topright", bty = "n", cex = 5, pt.cex = 5, 
       text.col = "black", 
       legend = paste0("DE = ", format(summary(GAIN.2070.GAM.CAP)$dev.expl *100, digits = 3), "%"))


#################################################################
## PANEL 3 :: GAMs of species gains vs MAT for 2070 ALL SUAs
par(font.axis = 1)

plot(SUA.LOSS.2070[,"CURRENT_MAT"], SUA.LOSS.2070[,"SPECIES_LOSS"], 
     col = SUA.LOSS.2070$color, pch = 19, cex = 6,
     cex.axis = 5, cex.lab = 5,
     las = 1, 
     ylab = "Species lost (%)", # (% current - 2070)", 
     xlab = "Current MAT of SUA (1960-1990)")

box(lwd = 3)

lines(LOSS.2070.TEST$CURRENT_MAT, PRED.LOSS.2070, col = "blue",  lwd = 10)

legend("bottomright", bty = "n", cex = 5, pt.cex = 5, 
       text.col = "black", 
       legend = paste0("DE = ", format(summary(LOSS.2070.GAM)$dev.expl *100, digits = 3), "%"))


#################################################################
## PANEL 4:: GAMs of species gains vs. MAT for 2070 BIG SUAs
par(font.axis = 1)

plot(SUA.LOSS.2070.CAP[,"CURRENT_MAT"], SUA.LOSS.2070.CAP[,"SPECIES_LOSS"], 
     col = SUA.LOSS.2070.CAP$color, pch = 19, cex = 6, 
     cex.axis = 4, cex.lab = 5,
     las = 1, 
     ylab = "", 
     xlab = "Current MAT of SUA (1960-1990)")

box(lwd = 3)

lines(LOSS.2070.TEST.CAP$CURRENT_MAT, PRED.LOSS.2070.CAP, col = "blue",  lwd = 10)

legend("bottomright", bty = "n", cex = 5, pt.cex = 5, 
       text.col = "black", 
       legend = paste0("DE = ", format(summary(LOSS.2070.GAM.CAP)$dev.expl *100, digits = 3), "%"))

## Finish the device
dev.off()





#########################################################################################################################
## 4). CREATE LEGEND
#########################################################################################################################


#########################################################################################################################
## Create PNG
CairoPNG(width = 10000, height = 16000,
         file = sprintf('output/figures/FIG_2/SUA_KOPPEN_LEGEND.png'),
         canvas = "white", bg = "white", units = "px", dpi = 600)


## Add mfrow
par(mar   = c(13.5, 16, 4, 4.8), 
    mgp   = c(11.8, 3, 0),
    oma   = c(1, 1, 1, 1))


#################################################################
## PANEL 1 :: GAMs of species gains vs MAT for 2070, ALL SUAs
par(font.axis = 1, xpd   = TRUE)

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)


##
legend(x = "center", 
       legend = as.character(sua.col$ClimateZ),
       col = sua.col$color, 
       bty = 'n', xjust = 1,
       pch = 16, pt.cex = 8, cex = 8,
       title = "Koppen Zone")


## finish the device
dev.off()





#########################################################################################################################
## 5). PLOT RESIDUALS FROM GAMs
#########################################################################################################################


#########################################################################################################################
## GAIN.2070.GAM, plot residuals: heteroscedasticity present
plot(fitted(GAIN.2070.GAM), resid(GAIN.2070.GAM), 
     pch = 19, col = "blue", cex = 1.2)

abline(h = 0, col = "orange", lwd = 4)


#########################################################################################################################
## GAIN.2070.GAM, plot residuals: heteroscedasticity present
plot(fitted(GAIN.2070.GAM.CAP), resid(GAIN.2070.GAM.CAP), 
     pch = 19, col = "blue", cex = 1.2)

abline(h = 0, col = "orange", lwd = 4)


#########################################################################################################################
## GAIN.2070.GAM, plot residuals: heteroscedasticity present
plot(fitted(LOSS.2070.GAM), resid(LOSS.2070.GAM), 
     pch = 19, col = "blue", cex = 1.2)

abline(h = 0, col = "orange", lwd = 4)


#########################################################################################################################
## GAIN.2070.GAM, plot residuals: heteroscedasticity present
plot(fitted(LOSS.2070.GAM.CAP), resid(LOSS.2070.GAM.CAP), 
     pch = 19, col = "blue", cex = 1.2)

abline(h = 0, col = "orange", lwd = 4)





#########################################################################################################################
## 5). NLS of GAIN/LOSS VS. TEMP OF SUA
#########################################################################################################################


#########################################################################################################################
## fit non-linear least squares regression (nls) model for Species gains
## you gotta play a bit to find the right k and z start values
GAIN.2070.NLS = nls(SPECIES_GAIN ~ z * exp(-k * CURRENT_MAT), 
                    data = SUA.GAIN.2070, trace = TRUE, start = list(k = 0.01, z = 50), model = TRUE)
plot(GAIN.2070.NLS, pch = 19)


## Try exponential and quadratic models
GAIN.2070.EXP   <- lm(log(SPECIES_GAIN) ~ CURRENT_MAT, data = SUA.GAIN.2070)
MAT2            <- SUA.GAIN.2070$CURRENT_MAT^2
quadratic.model <- lm(SPECIES_GAIN ~ CURRENT_MAT + MAT2, data = SUA.GAIN.2070)


## Check
plot(GAIN.2070.EXP) 
plot(quadratic.model) 

summary(GAIN.2070.NLS)
summary(GAIN.2070.EXP)


#calculate pseudo-R2
timevalues <- seq(0, 30, 1)
Counts.exponential2 <- exp(predict(GAIN.2070.EXP, list(CURRENT_MAT = timevalues)))
plot(SUA.GAIN.2070$CURRENT_MAT, SUA.GAIN.2070$SPECIES_GAIN, pch = 16)
lines(timevalues, Counts.exponential2, lwd = 2, col = "blue", xlab = "Time (s)", ylab = "Counts")





## plot points based on their “type”, we can change point pch/color based on SUAs Koppen if 
## desidered, otherwise remove “, type="n“ from previous command
plot(SUA.GAIN.2070$CURRENT_MAT, SUA.GAIN.2070$SPECIES_GAIN, 
     ylab = "Species gained (% current - 2070)", 
     xlab = "Current MAT of SUA (1960-1990)", type = "n")


points(SUA.GAIN.2070$CURRENT_MAT[type == "ClimateZ"], SUA.GAIN.2070$SPECIES_GAIN[type == "ClimateZ"], pch = 15)
points(MmSpecies[type=="HCR"], Mass6L365[type=="HCR"], pch=17)
points(MmSpecies[type=="LCP"], Mass6L365[type=="LCP"], pch=16)



#plot the NLS model
curve(coef(GAIN.2070.NLS)[2]*exp(-coef(GAIN.2070.NLS)[1]*x),0,6,add=T)



#add model as text to scatterplot. You need to change parameters from summary(nls)
# text(4,98,"y=77.11*exp(-0.29*x)", cex=0.7)
# text(3.8,90,"Pseudo R^2 = 0.54",  cex=0.7)



#########################################################################################################################
## 2). SCATTERPLOT FOR GAMS of GAIN/LOSS,
#########################################################################################################################


# #########################################################################################################################
# ## Create PNG
# CairoPNG(width = 18000, height = 16000, 
#          file = sprintf('output/figures/SUA_percent/SUA_2070_SPP_GL_PANEL_vs_%s_GAM_%s_%s.png',  
#                         'temp', SUAs, SUA_SPP),
#          canvas = "white", bg = "white", units = "px", dpi = 600)
# 
# ## Add mfrow
# par(mfrow = c(2, 2),
#     mar   = c(13.5, 16, 4, 4.8), 
#     mgp   = c(11.8, 3, 0),
#     oma   = c(1, 1, 1, 1))
# 
# 
# #################################################################
# ## PANEL 1 :: GAMs of species gains vs MAT for 2070, ALL SUAs
# par(font.axis = 1)
# 
# plot(SUA.GAIN.2070[,"CURRENT_MAT"], SUA.GAIN.2070[,"SPECIES_GAIN"], 
#      col = alpha("blue", 0.3), pch = 19, cex = 6, 
#      cex.axis = 5, cex.lab = 5,
#      las = 1, ylab = "Species gained (%)", # (% current - 2070)",
#      xlab = "")
# 
# box(lwd = 3)
# 
# lines(GAIN.2070.TEST$CURRENT_MAT, PRED.GAIN.2070, col = "orange",  lwd = 10)
# 
# legend("topright", bty = "n", cex = 5, pt.cex = 5, 
#        text.col = "black", 
#        legend = paste0("DE = ", 
#                        format(summary(GAIN.2070.GAM)$dev.expl *100, digits = 3), "%"))
# 
# 
# #################################################################
# ## PANEL 2 :: GAMs of species gains vs MAT for 2070, BIG SUAs
# par(font.axis = 1)
# 
# plot(SUA.GAIN.2070.CAP[,"CURRENT_MAT"], SUA.GAIN.2070.CAP[,"SPECIES_GAIN"], 
#      col = alpha("blue", 0.3), pch = 19, cex = 6, 
#      cex.axis = 5, cex.lab = 5,
#      las = 1,
#      ylab = "", 
#      xlab = "")
# 
# box(lwd = 3)
# 
# lines(GAIN.2070.TEST.CAP$CURRENT_MAT, PRED.GAIN.2070.CAP, col = "orange",  lwd = 10)
# 
# legend("topright", bty = "n", cex = 5, pt.cex = 5, 
#        text.col = "black", 
#        legend = paste0("DE = ", 
#                        format(summary(GAIN.2070.GAM.CAP)$dev.expl *100, digits = 3), "%"))
# 
# 
# #################################################################
# ## PANEL 3 :: GAMs of species gains vs MAT for 2070 ALL SUAs
# par(font.axis = 1)
# 
# plot(SUA.LOSS.2070[,"CURRENT_MAT"], SUA.LOSS.2070[,"SPECIES_LOSS"], 
#      col = alpha("blue", 0.3), pch = 19, cex = 6, 
#      cex.axis = 5, cex.lab = 5,
#      las = 1, 
#      ylab = "Species lost (%)", # (% current - 2070)", 
#      xlab = "Current MAT of SUA (1960-1990)")
# 
# box(lwd = 3)
# 
# lines(LOSS.2070.TEST$CURRENT_MAT, PRED.LOSS.2070, col = "orange",  lwd = 10)
# 
# legend("bottomright", bty = "n", cex = 5, pt.cex = 5, 
#        text.col = "black", 
#        legend = paste0("DE = ", format(summary(LOSS.2070.GAM)$dev.expl *100, digits = 3), "%"))
# 
# 
# #################################################################
# ## PANEL 4:: GAMs of species gains vs MAT for 2070 BIG SUAs
# par(font.axis = 1)
# 
# plot(SUA.LOSS.2070.CAP[,"CURRENT_MAT"], SUA.LOSS.2070.CAP[,"SPECIES_LOSS"], 
#      col = alpha("blue", 0.3), pch = 19, cex = 6, 
#      cex.axis = 5, cex.lab = 5,
#      las = 1, 
#      ylab = "", 
#      xlab = "Current MAT of SUA (1960-1990)")
# 
# box(lwd = 3)
# 
# lines(LOSS.2070.TEST.CAP$CURRENT_MAT, PRED.LOSS.2070.CAP, col = "orange",  lwd = 10)
# 
# legend("bottomright", bty = "n", cex = 5, pt.cex = 5, 
#        text.col = "black", 
#        legend = paste0("DE = ", format(summary(LOSS.2070.GAM.CAP)$dev.expl *100, digits = 3), "%"))
# 
# ## Finish the device
# dev.off()




#########################################################################################################################
## OUTSTANDING PLOT TASKS:
#########################################################################################################################