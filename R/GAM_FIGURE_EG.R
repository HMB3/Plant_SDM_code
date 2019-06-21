#########################################################################################################################
#############################################  EXAMPLE OF GAM FIGURE #################################################### 
#########################################################################################################################


## This code can be used to run and plot GAMs for one explanatory variable, while holding all other explanatory variables 
## constant. EG Tcrit vs geographic range, but controlling for soil and FPAR


## Partial regression plots would also be handy ::
## http://zevross.com/blog/2014/09/15/recreate-the-gam-partial-regression-smooth-plots-from-r-package-mgcv-with-a-little-style/

  
## Read data
library(mgcv)
library(plyr)
library(Cairo)
library(cairoDevice)
library(RColorBrewer) 
library(scales)





####################################################################################################################
## 1). CREATE SCATTERPLOTS
####################################################################################################################


## Read data in a table of all the variables needed to run the analysis
COMBO.NICHE = read.csv("./data/ANALYSIS/COMBO_NICHE_CONTEXT_ALA_RECORDS_DIANA_SPP.csv", stringsAsFactors = FALSE)
str(COMBO.NICHE)





####################################################################################################################
## 2). TCRIT GAMs
####################################################################################################################



## Run GAMS of Tcrit vs. different environmental variables..........................................................


#########################################################################################################
## FIT GAM OF TRCIT VS GEOGRAPHIC RANGE
fit1av = gam(Tcrit ~ s (rain_av, k = 5) + s(temp_av, k = 5) + s(species_count, k = 5)  + s(mean_simpson, k = 5), 
             data = COMBO.NICHE, family = tw(theta = NULL, link = "log", a = 1.01, b = 1.99), method = "REML")

testdata1av = data.frame(mean_simpson = mean(fit1av$model$mean_simpson),
						  species_count = mean(fit1av$model$species_count),
                         temp_av = mean(fit1av$model$temp_av),
                         rain_av = seq(min(COMBO.NICHE[,"rain_av"]),
                                       max(COMBO.NICHE[,"rain_av"]), length = length(COMBO.NICHE$gpp_av)))

pred_GPP1av = predict(fit1av, newdata = testdata1av, type ='response')

fitRav = gam(gpp_av ~ s (rain_av, k = 5), data = COMBO.NICHE, 
             family = tw(theta = NULL, link = "log", a = 1.01, b = 1.99), method = "REML")

testdataRav = data.frame(rain_av = seq(min(COMBO.NICHE[,"rain_av"]),
                                       max(COMBO.NICHE[,"rain_av"]), length = length(COMBO.NICHE$gpp_av)))

pred_GPPRav = predict(fitRav, newdata = testdataRav, type ='response')


########################################################
## TEMP

fit2av = gam(gpp_av ~ s (rain_av, k = 5) + s(temp_av, k = 5) + s(species_count, k = 5)  + s(mean_simpson, k = 5), 
             data = COMBO.NICHE, family = tw(theta = NULL, link = "log", a = 1.01, b = 1.99), method = "REML")

testdata2av = data.frame(mean_simpson = mean(fit2av$model$mean_simpson),
                         species_count = mean(fit2av$model$species_count),
						 rain_av = mean(fit2av$model$rain_av),
                         temp_av = seq(min(COMBO.NICHE[,"temp_av"]),
                                       max(COMBO.NICHE[,"temp_av"]), length = length(COMBO.NICHE$gpp_av)))

pred_GPP2av = predict(fit2av, newdata = testdata2av, type ='response')

fitTav = gam(gpp_av ~ s (temp_av, k = 5), data = COMBO.NICHE, 
             family = tw(theta = NULL, link = "log", a = 1.01, b = 1.99), method = "REML")

testdataTav = data.frame(temp_av = seq(min(COMBO.NICHE[,"temp_av"]),
                                       max(COMBO.NICHE[,"temp_av"]), length = length(COMBO.NICHE$gpp_av)))

pred_GPPTav = predict(fitTav, newdata = testdataTav, type ='response')


########################################################
## species_count

fit3av = gam(gpp_av ~ s (rain_av, k = 5) + s(temp_av, k = 5) + s(species_count, k = 5)  + s(mean_simpson, k = 5), 
             data = COMBO.NICHE, family = tw(theta = NULL, link = "log", a = 1.01, b = 1.99), method = "REML")

testdata3av = data.frame(mean_simpson = mean(fit3av$model$mean_simpson),
						 rain_av = mean(fit3av$model$rain_av),
                         temp_av = mean(fit3av$model$temp_av),
                         species_count = seq(min(COMBO.NICHE[,"species_count"]),
                                             max(COMBO.NICHE[,"species_count"]), length = length(COMBO.NICHE$gpp_av)))

pred_GPP3av = predict(fit3av, newdata = testdata3av, type ='response')

fitAav = gam(gpp_av ~ s (species_count, k = 5), data = COMBO.NICHE, 
             family = tw(theta = NULL, link = "log", a = 1.01, b = 1.99), method = "REML")

testdataAav = data.frame(species_count = seq(min(COMBO.NICHE[,"species_count"]),
                                             max(COMBO.NICHE[,"species_count"]), length = length(COMBO.NICHE$gpp_av)))

pred_GPPAav = predict(fitAav, newdata = testdataAav, type ='response')


########################################################
## BETA

fit4av = gam(gpp_av ~ s (rain_av, k = 5) + s(temp_av, k = 5) + s(species_count, k = 5)  + s(mean_simpson, k = 5), 
             data = COMBO.NICHE, family = tw(theta = NULL, link = "log", a = 1.01, b = 1.99), method = "REML")

testdata4av = data.frame(rain_av = mean(fit4av$model$rain_av),
                         temp_av = mean(fit4av$model$temp_av),
                         species_count = mean(fit4av$model$species_count),
						 mean_simpson = seq(min(COMBO.NICHE[,"mean_simpson"]),
                                            max(COMBO.NICHE[,"mean_simpson"]), length = length(COMBO.NICHE$gpp_av)))

pred_GPP4av = predict(fit4av, newdata = testdata4av, type ='response')


fitBav = gam(gpp_av ~ s (mean_simpson, k = 5), data = COMBO.NICHE, 
             family = tw(theta = NULL, link = "log", a = 1.01, b = 1.99), method = "REML")

testdataBav = data.frame(mean_simpson = seq(min(COMBO.NICHE[,"mean_simpson"]),
                                            max(COMBO.NICHE[,"mean_simpson"]), length = length(COMBO.NICHE$gpp_av)))

pred_GPPBav = predict(fitBav, newdata = testdataBav, type ='response')




####################################################################################################################
## CV GAMs
####################################################################################################################

########################################################
## RAIN
fit1cv = gam(gpp_cv ~ s (rain_cv, k = 5) + s(temp_cv, k = 5) + s(species_count, k = 5)  + s(mean_simpson, k = 5), 
             data = COMBO.NICHE, family = tw(theta = NULL, link = "log", a = 1.01, b = 1.99), method = "REML")

testdata1cv = data.frame(mean_simpson = mean(fit1cv$model$mean_simpson),
						 species_count = mean(fit1cv$model$species_count),
                         temp_cv = mean(fit1cv$model$temp_cv),
                         rain_cv = seq(min(COMBO.NICHE[,"rain_cv"]),
                                       max(COMBO.NICHE[,"rain_cv"]), length = length(COMBO.NICHE$gpp_cv)))

pred_GPP1cv = predict(fit1cv, newdata = testdata1cv, type ='response')

fitRcv = gam(gpp_cv ~ s (rain_cv, k = 5), data = COMBO.NICHE, 
             family = tw(theta = NULL, link = "log", a = 1.01, b = 1.99), method = "REML")

testdataRcv = data.frame(rain_cv = seq(min(COMBO.NICHE[,"rain_cv"]),
                                       max(COMBO.NICHE[,"rain_cv"]), length = length(COMBO.NICHE$gpp_cv)))

pred_GPPRcv = predict(fitRcv, newdata = testdataRcv, type ='response')


########################################################
## TEMP

fit2cv = gam(gpp_cv ~ s (rain_cv, k = 5) + s(temp_cv, k = 5) + s(species_count, k = 5)  + s(mean_simpson, k = 5), 
             data = COMBO.NICHE, family = tw(theta = NULL, link = "log", a = 1.01, b = 1.99), method = "REML")

testdata2cv = data.frame(mean_simpson = mean(fit2cv$model$mean_simpson),
                         species_count = mean(fit2cv$model$species_count),
						 rain_cv = mean(fit2cv$model$rain_cv),
                         temp_cv = seq(min(COMBO.NICHE[,"temp_cv"]),
                                       max(COMBO.NICHE[,"temp_cv"]), length = length(COMBO.NICHE$gpp_cv)))

pred_GPP2cv = predict(fit2cv, newdata = testdata2cv, type ='response')

fitTcv = gam(gpp_cv ~ s (temp_cv, k = 5), data = COMBO.NICHE, 
             family = tw(theta = NULL, link = "log", a = 1.01, b = 1.99), method = "REML")

testdataTcv = data.frame(temp_cv = seq(min(COMBO.NICHE[,"temp_cv"]),
                                       max(COMBO.NICHE[,"temp_cv"]), length = length(COMBO.NICHE$gpp_cv)))

pred_GPPTcv = predict(fitTcv, newdata = testdataTcv, type ='response')


########################################################
## species_count

fit3cv = gam(gpp_cv ~ s (rain_cv, k = 5) + s(temp_cv, k = 5) + s(species_count, k = 5)  + s(mean_simpson, k = 5), 
             data = COMBO.NICHE, family = tw(theta = NULL, link = "log", a = 1.01, b = 1.99), method = "REML")

testdata3cv = data.frame(mean_simpson = mean(fit2cv$model$mean_simpson),
						 rain_cv = mean(fit3cv$model$rain_cv),
                         temp_cv = mean(fit3cv$model$temp_cv),
                         species_count = seq(min(COMBO.NICHE[,"species_count"]),
                                             max(COMBO.NICHE[,"species_count"]), length = length(COMBO.NICHE$gpp_cv)))

pred_GPP3cv = predict(fit3cv, newdata = testdata3cv, type ='response')

fitAcv = gam(gpp_cv ~ s (species_count, k = 5), data = COMBO.NICHE, 
             family = tw(theta = NULL, link = "log", a = 1.01, b = 1.99), method = "REML")

testdataAcv = data.frame(species_count = seq(min(COMBO.NICHE[,"species_count"]),
                                             max(COMBO.NICHE[,"species_count"]), length = length(COMBO.NICHE$gpp_cv)))

pred_GPPAcv = predict(fitAcv, newdata = testdataAcv, type ='response')


########################################################
## BETA
## WHY IS THIS SHORTER THAN OTHER VARIABLES????

fit4cv = gam(gpp_cv ~ s (rain_cv, k = 5) + s(temp_cv, k = 5) + s(species_count, k = 5)  + s(mean_simpson, k = 5), 
             data = COMBO.NICHE, family = tw(theta = NULL, link = "log", a = 1.01, b = 1.99), method = "REML")

testdata4cv = data.frame(rain_cv = mean(fit4cv$model$rain_cv),
                         temp_cv = mean(fit4cv$model$temp_cv),
                         species_count = mean(fit4cv$model$species_count),
						 mean_simpson = seq(min(COMBO.NICHE[,"mean_simpson"]),
                                            max(COMBO.NICHE[,"mean_simpson"]), length = length(COMBO.NICHE$gpp_cv)))

pred_GPP4cv = predict(fit4cv, newdata = testdata4cv, type ='response')


fitBcv = gam(gpp_cv ~ s (mean_simpson, k = 5), data = COMBO.NICHE, 
             family = tw(theta = NULL, link = "log", a = 1.01, b = 1.99), method = "REML")

testdataBcv = data.frame(mean_simpson = seq(min(COMBO.NICHE[,"mean_simpson"]),
                                            max(COMBO.NICHE[,"mean_simpson"]), length = length(COMBO.NICHE$gpp_cv)))

pred_GPPBcv = predict(fitBcv, newdata = testdataBcv, type ='response')


## END MODEL CODE


#####################################################################################################################
## create one panel for AV and CV plots


## vary colour of plot using RColorBrewer, Default is probably PuBu
## ArcMap beta diversity colour scheme is PuBuGn

# colramp=colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))

## GEB needs to be smaller than 4000*4000 pixels.
## height = 4000, width = (10000/16180)*4000

CairoPNG(width = 10000, height = 16180, 
         file = "./output/figures/figure_3/CHAP2_FIG4_UPDATE_TRANSPARENT.png", 
         canvas = "white", bg = "white", units = "px", dpi = 600)


par(mfrow = c(4, 2),
    mar   = c(13.5, 16, 4, 4.8), 
    mgp   = c(11.8, 3, 0),
    oma   = c(1, 1, 1, 1))


## PANEL 1
par(font.axis = 1)

plot(COMBO.NICHE[,"rain_av"], COMBO.NICHE[,"gpp_av"], 
              col = alpha("blue", 0.3), pch = 19, cex = 2, 
              cex.axis = 5, cex.lab = 6,
              las = 1, xlab = "Monthly rainfall (mm)", ylab = "GPP")
box(lwd=3)

lines(testdataRav$rain_av, pred_GPPRav, col = "green",  lwd = 8)
lines(testdata1av$rain_av, pred_GPP1av, col = "orange", lwd = 8)

## PANEL 2
par(font.axis = 1) #par(font.lab = 2)

plot(COMBO.NICHE[,"rain_cv"], COMBO.NICHE[,"gpp_cv"], 
              col = alpha("blue", 0.3), pch = 19, cex = 2, 
              cex.axis = 5, cex.lab = 6,
              las = 1, xlab = "CV rain", ylab = "CV GPP")
box(lwd=3)

lines(testdataRcv$rain_cv, pred_GPPRcv, col = "green", lwd = 8)
lines(testdata1cv$rain_cv, pred_GPP1cv, col = "orange", lwd = 8)

## PANEL 3
par(font.axis = 1) #par(font.lab = 2)
plot(COMBO.NICHE[,"temp_av"], COMBO.NICHE[,"gpp_av"], 
              col = alpha("blue", 0.3), pch = 19, cex = 2, 
              cex.axis = 5, cex.lab = 6, 
              las = 1, xlab = "Monthly max temp (°C)", ylab = "GPP")
box(lwd=3)

lines(testdataTav$temp_av,pred_GPPTav,col = "green", lwd = 8)
lines(testdata2av$temp_av, pred_GPP2av, col = "orange", lwd = 8)

## PANEL 4
par(font.axis = 1) #par(font.lab = 2)

plot(COMBO.NICHE[,"temp_cv"], COMBO.NICHE[,"gpp_cv"], 
              col = alpha("blue", 0.3), pch = 19, cex = 2, 
              cex.axis = 5, cex.lab = 6, 
              las = 1, xlab = "CV temp", ylab = "CV GPP")

box(lwd=3)

lines(testdataTcv$temp_cv,pred_GPPTcv,col = "green", lwd = 8)
lines(testdata2cv$temp_cv, pred_GPP2cv, col = "orange", lwd = 8)

## PANEL 5
par(font.axis = 1) #par(font.lab = 2)

plot(COMBO.NICHE[,"species_count"], COMBO.NICHE[,"gpp_av"], 
              col = alpha("blue", 0.3), pch = 19, cex = 2, 
              cex.axis = 5, cex.lab = 6, 
              las = 1, xlab = "Alpha diversity", ylab = "GPP")

box(lwd=3)

lines(testdataAav$species_count, pred_GPPAav,col = "green", lwd = 8)
lines(testdata3av$species_count, pred_GPP3av, col = "orange", lwd = 8)


## PANEL 6
par(font.axis = 1) #par(font.lab = 2)

plot(COMBO.NICHE[,"species_count"], COMBO.NICHE[,"gpp_cv"], 
              col = alpha("blue", 0.3), pch = 19, cex = 2, 
              cex.axis = 5, cex.lab = 6, 
              las = 1, xlab = "Alpha diversity", ylab = "CV GPP")

box(lwd=3)

lines(testdataAcv$species_count, pred_GPPAcv,col = "green", lwd = 8)
lines(testdata3cv$species_count, pred_GPP3cv, col = "orange", lwd = 8)

## PANEL 7
par(font.axis = 1) #par(font.lab = 2)

plot(COMBO.NICHE[,"mean_simpson"], COMBO.NICHE[,"gpp_av"], 
              col = alpha("blue", 0.3), pch = 19, cex = 2, 
              cex.axis = 5, cex.lab = 6, 
              las = 1, xlab = "Beta diversity", ylab = "GPP")

box(lwd=3)

lines(testdataBav$mean_simpson, pred_GPPBav,col = "green", lwd = 8)
lines(testdata4av$mean_simpson, pred_GPP4av, col = "orange", lwd = 8)


## PANEL 8
par(font.axis = 1) #par(font.lab = 2)

plot(COMBO.NICHE[,"mean_simpson"], COMBO.NICHE[,"gpp_cv"], 
              col = alpha("blue", 0.3), pch = 19, cex = 2, 
              cex.axis = 5, cex.lab = 6, 
              las = 1, xlab = "Beta diversity", ylab = "CV GPP")

box(lwd=3)

lines(testdataBcv$mean_simpson, pred_GPPBcv,col = "green", lwd = 8)
lines(testdata4cv$mean_simpson, pred_GPP4cv, col = "orange", lwd = 8)


## END OF PANEL CODE
dev.off() ## CREATE LOOP for this






#####################################################################################################################
## REST OF FIGURE IN GIS, or COMBINE R SPATIAL with PLOTTING
#####################################################################################################################