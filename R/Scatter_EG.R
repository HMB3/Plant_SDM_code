library(mgcv)
library(plyr)
library(Cairo)
library(cairoDevice)
library(RColorBrewer) 
library(scales)



########################################################
## RAIN
fit1av = gam(gpp_av ~ s (rain_av, k = 5), 
             data = R1, family = tw(theta = NULL, link = "log", a = 1.01, b = 1.99), method = "REML")

testdata1av = data.frame(
  # mean_simpson = mean(fit1av$model$mean_simpson),
  #                        species_count = mean(fit1av$model$species_count),
  #                        temp_av = mean(fit1av$model$temp_av),
                         rain_av = seq(min(R1[,"rain_av"]),
                                       max(R1[,"rain_av"]), length = length(R1$gpp_av)))

pred_GPP1av = predict(fit1av, newdata = testdata1av, type ='response')



#####################################################################################################################
## now plot the GPP ~ rainfall by itself for the manuscript
## plot negative exponential nls
CairoPNG(width  = 10000, height = 13000, ## 13000, 10000
         file   = "./output/figs/fig_3/FIG3_A_GPP_rain_NLS_B_GPP_delta_new_rain.png", 
         canvas = "white", bg = "white", units = "px", dpi = 600)

par(mfrow = c(2, 1),            ## c(nrows, ncols)
    mar   = c(11.8, 15, 2.5, 4),  ## b, l, t, r
    mgp   = c(9.8, 2.5, 0),
    oma   = c(1.5, 1.5, 1.5, 1.5))

par(font.axis = 1)#, font.lab = 2, lwd = 2)
par(bg = NA)

plot(SITES.REMOVE$RAIN_0112, SITES.REMOVE$TOTAL_AV_REPLACE, 
     pch  = 19, col = "blue", #cex = 1.2,
     xlab = "",
     ylab = expression(paste("GPP (gC m"^"-2", "month"^"-1", ")", sep = " ")),
     cex = 2.5, cex.axis = 4, cex.lab = 5, las = 1)

lines(sort(SITES.REMOVE$RAIN_0112), fitted(nls.exp.raw)
      [order(SITES.REMOVE$RAIN_0112)], col= 'orange', type = 'l', lwd = 12)

abline(v = 139, col = "grey", lwd = 8, lty = 99)

legend("bottomright", bty = "n", cex = 5, pt.cex = 5,
       text.col = "orange", 
       legend = paste("DE =", format(EP.exp.raw, digits = 3)))

box(lwd = 3)











#####################################################################################################################
## Now do a pairs plot of all the analysis variables for figure 3
#detach("package:dplyr")
#source("./R/SITE_ENW_REGRESSION_FUNCTIONS.R")  ## might need to load the latest version of the pairs function
#library(plyr)
SITES.PLOTS = SITES.TREES
SITES.PLOTS = rename(SITES.PLOTS,   c("TOTAL_AV_REPLACE"  = "GPP",
                                      "RAIN_0112"         = "RAIN",
                                      "RAIN_NW_MN"        = "RNWMN",
                                      "TEMP_0112"         = "TEMP",
                                      "TEMP_NW_MN"        = "TNWMN",
                                      "NICHE_AREA_MN"     = "NAMN",
                                      
                                      "RAIN_NW_MD"        = "RNW",
                                      "TEMP_NW_MD"        = "TNW"))


####################################################
## GPP magnitude vs. median variables, median niches
CairoPNG(width = 8090, height = 8090, 
         file = "./output/figs/fig_3/FIG3_C_median_pairs_plot_NW_first.png", 
         canvas = "white", bg = "white", units = "px", dpi = 600)

par(mar = c(8, 8, 6, 4), 
    mgp = c(6, 2, 0))

par(lwd = 2)

pairs(SITES.PLOTS[,c("GPP", "RNW", "RAIN",
                     "TNW", "TEMP")],  
      lower.panel = panel.cor,
      #diag.panel  = panel.hist,
      upper.panel = panel.smooth, #function(...) smoothScatter(..., nrpoints = 0, add = TRUE))
      cex.labels  = 5.5, cex.axis = 3, font.labels = 2)

## finish the device
dev.off()


####################################################
## GPP magnitude vs. RNW and TNW
CairoPNG(width = 8090, height = 8090, 
         file = "./output/figs/fig_3/FIG3_RNW_vs_TNW.png", 
         canvas = "white", bg = "white", units = "px", dpi = 600)

par(mar = c(8, 8, 6, 4), 
    mgp = c(6, 2, 0))

par(lwd = 2)

pairs(SITES.PLOTS[,c("GPP", "RNW", "TNW")],  
      lower.panel = panel.cor,
      #diag.panel  = panel.hist,
      upper.panel = panel.smooth, #function(...) smoothScatter(..., nrpoints = 0, add = TRUE))
      cex.labels  = 5.5, cex.axis = 3, font.labels = 2)

## finish the device
dev.off()



