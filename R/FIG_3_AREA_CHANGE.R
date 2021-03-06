#########################################################################################################################
############################################### FIG 3 STOTEN ############################################################ 
#########################################################################################################################


#########################################################################################################################
## Code for Producing Fig.3 STOTEN Manuscript - Burley et al., 2019
## Which Plant Where Project - Macquarie University Sydney, Australia
## Ossola Alessandro, 29.01.2019


#########################################################################################################################
## 1). SUBSET
#########################################################################################################################


## Read data
unique(SUA.PLOT$AREA_CHANGE)
table(SUA.PLOT$PERIOD, SUA.PLOT$AREA_CHANGE)
unique(SUA.PLOT$ClimateZ)


## Subset data
allSUAsGAIN <- subset(SUA.PLOT, AREA_CHANGE=="GAIN") 
allSUAsLOSS <- subset(SUA.PLOT, AREA_CHANGE=="LOSS")

bigSUAsGAIN <- subset(SUA.PLOT, AREASQKM16>200 & POP_2017>80000 & AREA_CHANGE=="GAIN")
bigSUAsLOSS <- subset(SUA.PLOT, AREASQKM16>200 & POP_2017>80000 & AREA_CHANGE=="LOSS")


## To calculate GAM deviance
allSUAsGAIN2070 <- subset(SUA.PLOT, PERIOD==2070 & AREA_CHANGE=="GAIN")
allSUAsLOSS2070 <- subset(SUA.PLOT, PERIOD==2070 & AREA_CHANGE=="LOSS")
allSUAsGAIN2030 <- subset(SUA.PLOT, PERIOD==2030 & AREA_CHANGE=="GAIN")
allSUAsLOSS2030 <- subset(SUA.PLOT, PERIOD==2030 & AREA_CHANGE=="LOSS")

bigSUAsGAIN2070 <- subset(SUA.PLOT, PERIOD==2070 & AREA_CHANGE=="GAIN" & AREASQKM16>200 & POP_2017>80000)
bigSUAsLOSS2070 <- subset(SUA.PLOT, PERIOD==2070 & AREA_CHANGE=="LOSS" & AREASQKM16>200 & POP_2017>80000)
bigSUAsGAIN2030 <- subset(SUA.PLOT, PERIOD==2030 & AREA_CHANGE=="GAIN" & AREASQKM16>200 & POP_2017>80000)
bigSUAsLOSS2030 <- subset(SUA.PLOT, PERIOD==2030 & AREA_CHANGE=="LOSS" & AREASQKM16>200 & POP_2017>80000)





#########################################################################################################################
## 2). RUN GAMs FOR EACH PANEL
#########################################################################################################################


#########################################################################################################################
## Fig 3A 
## Run GAMs of species gains vs MAT for 2070, ALL SUA
ALL.GAIN.2070.GAM = gam(SPECIES_GAIN ~ s (CURRENT_MAT, k = 3), 
                    data = allSUAsGAIN2070, 
                    method = "REML")
summary(ALL.GAIN.2070.GAM)[["dev.expl"]][1] 


######################################################
## Run GAMs of species gains vs MAT for 2030, ALL SUA
ALL.GAIN.2030.GAM = gam(SPECIES_GAIN ~ s (CURRENT_MAT, k = 3), 
                    data = allSUAsGAIN2030, 
                    method = "REML")
summary(ALL.GAIN.2030.GAM)[["dev.expl"]][1]  


#########################################################################################################################
## 
dev.new(width=17, height=13)

fig3A<-ggplot(allSUAsGAIN, aes(x=CURRENT_MAT, y=SPECIES_GAIN, color=factor(PERIOD))) +
  
       geom_point(size = 6) + theme_bw() +  theme(axis.title.x=element_text(margin = margin(t = 20)), axis.title.y=element_text(margin = margin(r = 20)),
       axis.text.x=element_text(margin = margin(t = 10)), axis.text.y=element_text(margin = margin(r = 10)), legend.box.margin=margin(10,10,10,10)) +
  
       geom_smooth(method="gam", formula = y ~ s(x, k = 3), se=FALSE, fullrange=TRUE, size = 2) +
       xlab(bquote('')) + ylim(0,70) + ylab(bquote('Range gain (%)')) + ggtitle("All SUAs") + theme(text = element_text(size=30)) + 
       scale_color_manual(values=c("darkgrey", "black"), name = "Period") + 
  
       annotate(geom="text", size=8, x=20, y=65, label=paste("Deviance (2030) =", signif(summary(ALL.GAIN.2030.GAM)[["dev.expl"]][1]*100, digits = 3),"%")) +
       annotate(geom="text", size=8, x=20, y=62, label=paste("Deviance (2070) =", signif(summary(ALL.GAIN.2070.GAM)[["dev.expl"]][1]*100, digits = 3),"%")) +
       theme(legend.position="none")  #this remove legend on the right

fig3A





#########################################################################################################################
### Fig 3B 
## Run GAMs of species gains vs MAT for 2070, BIG SUA
BIG.GAIN.2070.GAM = gam(SPECIES_GAIN ~ s (CURRENT_MAT, k = 3), 
                    data = bigSUAsGAIN2070, 
                    method = "REML")
summary(BIG.GAIN.2070.GAM)[["dev.expl"]][1]


######################################################
## Run GAMs of species gains vs MAT for 2030, BIG SUA
BIG.GAIN.2030.GAM = gam(SPECIES_GAIN ~ s (CURRENT_MAT, k = 3), 
                    data = bigSUAsGAIN2030, 
                    method = "REML")
summary(BIG.GAIN.2030.GAM)[["dev.expl"]][1]


#########################################################################################################################
dev.new(width=17, height=13)

fig3B<-ggplot(bigSUAsGAIN, aes(x=CURRENT_MAT, y=SPECIES_GAIN, color=factor(PERIOD))) +
  
       geom_point(size = 6) + theme_bw() +  theme(axis.title.x=element_text(margin = margin(t = 20)), axis.title.y=element_text(margin = margin(r = 20)),
       axis.text.x=element_text(margin = margin(t = 10)), axis.text.y=element_text(margin = margin(r = 10)), legend.box.margin=margin(10,10,10,10)) +
       geom_smooth(method="gam", formula = y ~ s(x, k = 3), se=FALSE, fullrange=TRUE, size = 2) +
  
       xlab(bquote('')) + ylim(0,70) + 
       ylab(bquote('')) + ggtitle("Largest SUAs") + theme(text = element_text(size=30)) + 
       scale_color_manual(values=c("darkgrey", "black"), name = "Period") +
  
       annotate(geom="text", size=8, x=20, y=65, label=paste("Deviance (2030) =", signif(summary(BIG.GAIN.2030.GAM)[["dev.expl"]][1]*100, digits = 3),"%")) +
       annotate(geom="text", size=8, x=20, y=62, label=paste("Deviance (2070) =", signif(summary(BIG.GAIN.2070.GAM)[["dev.expl"]][1]*100, digits = 3),"%")) 
       #theme(legend.position="none")  #this remove legend on the right
       #
fig3B





#########################################################################################################################
## Fig 3C 
## Run GAMs of species losses vs MAT for 2070, ALL SUA
ALL.LOSS.2070.GAM = gam(SPECIES_LOSS ~ s (CURRENT_MAT, k = 3), 
                    data = allSUAsLOSS2070, 
                    method = "REML")
summary(ALL.LOSS.2070.GAM)[["dev.expl"]][1] 


######################################################
## Run GAMs of species losses vs MAT for 2030, ALL SUA
ALL.LOSS.2030.GAM = gam(SPECIES_LOSS ~ s (CURRENT_MAT, k = 3), 
                    data = allSUAsLOSS2030, 
                    method = "REML")
summary(ALL.LOSS.2030.GAM)[["dev.expl"]][1]


#########################################################################################################################
dev.new(width=17, height=13)

fig3C<-ggplot(allSUAsLOSS, aes(x=CURRENT_MAT, y=SPECIES_LOSS, color=factor(PERIOD))) +
  
       geom_point(size = 6) + theme_bw() +  theme(axis.title.x=element_text(margin = margin(t = 20)), axis.title.y=element_text(margin = margin(r = 20)),
       axis.text.x=element_text(margin = margin(t = 10)), axis.text.y=element_text(margin = margin(r = 10)), legend.box.margin=margin(10,10,10,10)) +
  
       geom_smooth(method="gam", formula = y ~ s(x, k = 3), se=FALSE, fullrange=TRUE, size = 2) +
       xlab(bquote('Current MAT of SUA (1960-1990)')) + ylim(0,80) + ylab(bquote('Range loss (%)')) + theme(text = element_text(size=30)) + 
       scale_color_manual(values=c("darkgrey", "black"), name = "Period") + 
  
       annotate(geom="text", size=8, x=20, y=7, label=paste("Deviance (2030) =", signif(summary(ALL.LOSS.2030.GAM)[["dev.expl"]][1]*100, digits = 3),"%")) +
       annotate(geom="text", size=8, x=20, y=4, label=paste("Deviance (2070) =", signif(summary(ALL.LOSS.2070.GAM)[["dev.expl"]][1]*100, digits = 3),"%")) +
       theme(legend.position="none")  #this remove legend on the right

fig3C





#########################################################################################################################
## Fig 3D 
## Run GAMs of species losses vs MAT for 2070, BIG SUA
BIG.LOSS.2070.GAM = gam(SPECIES_LOSS ~ s (CURRENT_MAT, k = 3), 
                    data = bigSUAsLOSS2070, 
                    method = "REML")
summary(BIG.LOSS.2070.GAM)[["dev.expl"]][1] 


######################################################
## Run GAMs of species losses vs MAT for 2030, BIG SUA
BIG.LOSS.2030.GAM = gam(SPECIES_LOSS ~ s (CURRENT_MAT, k = 3), 
                    data = bigSUAsLOSS2030, 
                    method = "REML")
summary(BIG.LOSS.2030.GAM)[["dev.expl"]][1]  


#########################################################################################################################
dev.new(width=17, height=13)

fig3D<-ggplot(bigSUAsLOSS, aes(x=CURRENT_MAT, y=SPECIES_LOSS, color=factor(PERIOD))) +
  
       geom_point(size = 6) + theme_bw() +  theme(axis.title.x=element_text(margin = margin(t = 20)), axis.title.y=element_text(margin = margin(r = 20)),
       axis.text.x=element_text(margin = margin(t = 10)), axis.text.y=element_text(margin = margin(r = 10)), legend.box.margin=margin(10,10,10,10)) +
  
       geom_smooth(method="gam", formula = y ~ s(x, k = 3), se=FALSE, fullrange=TRUE, size = 2) +
       xlab(bquote('Current MAT of SUA (1960-1990)')) + ylim(0,80) + 
       ylab(bquote('')) + theme(text = element_text(size=30)) + 
  
       scale_color_manual(values=c("darkgrey", "black"), name = "Period") + 
       annotate(geom="text", size=8, x=20, y=7, label=paste("Deviance (2030) =", signif(summary(BIG.LOSS.2030.GAM)[["dev.expl"]][1]*100, digits = 3),"%")) +
       annotate(geom="text", size=8, x=20, y=4, label=paste("Deviance (2070) =", signif(summary(BIG.LOSS.2070.GAM)[["dev.expl"]][1]*100, digits = 3),"%")) + 
       theme(legend.position="none")  #this remove legend on the right

fig3D





#########################################################################################################################
## 3). CREATE A PANEL FOR ALL FIGURES
#########################################################################################################################


#########################################################################################################################
dev.new(width = 40, height = 30)
fig3<- grid_arrange_shared_legend(fig3A, fig3B, fig3C, fig3D, ncol = 2, nrow = 2, position = "right")
fig3


## 
ggsave(filename = "./output/figures/FIG_3/FIG_3_PERCENT_SUA_SPECIES.png", width = 20, height = 18, plot = fig2, dpi = 300)





#########################################################################################################################
## 4). CHECK WHETHER  POPULATION oR AREA MAKES A DIFFERENCE
#########################################################################################################################


## These patterns of change were consistent regardless of the area and population of the SUAs (Fig. 3)?
## How do we test that? 

#########################################################################################################################
## First get rid of the NA pop rows, then log population and area
allSUAsLOSS.AOV          = completeFun(allSUAsLOSS, "POP_2017")
allSUAsLOSS.AOV$LOG_POP  = log(allSUAsLOSS.AOV$POP_2017)
allSUAsLOSS.AOV$LOG_AREA = log(allSUAsLOSS.AOV$AREASQKM16)


## For the gains
allSUAsGAIN.AOV          = completeFun(allSUAsGAIN, "POP_2017")
allSUAsGAIN.AOV$LOG_POP  = log(allSUAsGAIN.AOV$POP_2017)
allSUAsGAIN.AOV$LOG_AREA = log(allSUAsGAIN.AOV$AREASQKM16)

bigSUAsGAIN.AOV <- subset(allSUAsGAIN.AOV, AREASQKM16>200 & POP_2017>80000)
bigSUAsLOSS.AOV <- subset(allSUAsLOSS.AOV, AREASQKM16>200 & POP_2017>80000)


#########################################################################################################################
## To calculate GAM deviance
allSUAsGAIN2070.AOV <- subset(allSUAsGAIN.AOV, PERIOD==2070 & AREA_CHANGE=="GAIN")
allSUAsLOSS2070.AOV <- subset(allSUAsLOSS.AOV, PERIOD==2070 & AREA_CHANGE=="LOSS")
allSUAsGAIN2030.AOV <- subset(allSUAsGAIN.AOV, PERIOD==2030 & AREA_CHANGE=="GAIN")
allSUAsLOSS2030.AOV <- subset(allSUAsLOSS.AOV, PERIOD==2030 & AREA_CHANGE=="LOSS")

bigSUAsGAIN2070.AOV <- subset(bigSUAsGAIN.AOV, PERIOD==2070 & AREA_CHANGE=="GAIN" & AREASQKM16>200 & POP_2017>80000)
bigSUAsLOSS2070.AOV <- subset(bigSUAsLOSS.AOV, PERIOD==2070 & AREA_CHANGE=="LOSS" & AREASQKM16>200 & POP_2017>80000)
bigSUAsGAIN2030.AOV <- subset(bigSUAsGAIN.AOV, PERIOD==2030 & AREA_CHANGE=="GAIN" & AREASQKM16>200 & POP_2017>80000)
bigSUAsLOSS2030.AOV <- subset(bigSUAsLOSS.AOV, PERIOD==2030 & AREA_CHANGE=="LOSS" & AREASQKM16>200 & POP_2017>80000)


#########################################################################################################################
## Now plot the difference
allSUAsLOSS.PLOT =  allSUAsLOSS.AOV[c("SPECIES_LOSS", "CURRENT_MAT", 
                                      "CURRENT_MAP",  "CURRENT_MAXT", 
                                      "LOG_AREA",     "LOG_POP")]
names(allSUAsLOSS.PLOT) = c("SUA_LOSS", "MAT","MAP",  "MAXT", "LOG_AREA", "LOG_POP")

allSUAsGAIN.PLOT =  allSUAsGAIN.AOV[c("SPECIES_GAIN", "CURRENT_MAT", 
                                      "CURRENT_MAP",  "CURRENT_MAXT", 
                                      "LOG_AREA",     "LOG_POP")]
names(allSUAsGAIN.PLOT) = c("SUA_GAIN", "MAT","MAP",  "MAXT", "LOG_AREA", "LOG_POP")


#########################################################################################################################
## Some strong relationships with MAXT
par(mar = c(8, 8, 6, 4),
    mgp = c(6, 2, 0))

par(lwd = 2)

pairs(allSUAsLOSS.PLOT,  
      lower.panel = panel.cor,
      upper.panel = panel.smooth,
      cex.labels  = 1.5, cex.axis = 2, font.labels = 2,
      main = "Range lost (%) vs SUA climate (1960-1990)")

#########################################################################################################################
## Some strong relationships with MAXT
par(mar = c(8, 8, 6, 4),
    mgp = c(6, 2, 0))

par(lwd = 2)

pairs(allSUAsGAIN.PLOT,  
      lower.panel = panel.cor,
      upper.panel = panel.smooth,
      cex.labels  = 1.5, cex.axis = 2, font.labels = 2,
      main = "Range gained (%) vs SUA climate (1960-1990)")



#########################################################################################################################
## Run GAMs of species gains vs MAT for 2070, ALL SUA
ALL.GAIN.2070.POP = gam(SPECIES_GAIN ~ s(CURRENT_MAT, k = 3) + s(LOG_AREA, k = 3) + s(LOG_POP, k = 3),
                        data = allSUAsGAIN2070.AOV, 
                        method = "REML")
print(summary(ALL.GAIN.2070.POP))


######################################################
## Run GAMs of species gains vs MAT for 2030, ALL SUA
ALL.GAIN.2030.POP = gam(SPECIES_GAIN ~ s (CURRENT_MAT, k = 3)+ s(AREASQKM16, k = 3) + s(POP_2017, k = 3), 
                        data = allSUAsGAIN2030, 
                        method = "REML")
print(summary(ALL.GAIN.2030.POP))


#########################################################################################################################
### Fig 3B 
## Run GAMs of species gains vs MAT for 2070, BIG SUA
BIG.GAIN.2070.POP = gam(SPECIES_GAIN ~ s (CURRENT_MAT, k = 3)+ s(AREASQKM16, k = 3) + s(POP_2017, k = 3), 
                        data = bigSUAsGAIN2070, 
                        method = "REML")
print(summary(BIG.GAIN.2070.POP))


######################################################
## Run GAMs of species gains vs MAT for 2030, BIG SUA
BIG.GAIN.2030.POP = gam(SPECIES_GAIN ~ s (CURRENT_MAT, k = 3)+ s(AREASQKM16, k = 3) + s(POP_2017, k = 3),
                        data = bigSUAsGAIN2030, 
                        method = "REML")
print(summary(BIG.GAIN.2070.POP))


#########################################################################################################################
## Fig 3C 
## Run GAMs of species losses vs MAT for 2070, ALL SUA
ALL.LOSS.2070.POP = gam(SPECIES_LOSS ~ s (CURRENT_MAT, k = 3)+ s(AREASQKM16, k = 3) + s(POP_2017, k = 3),
                        data = allSUAsLOSS2070, 
                        method = "REML")
print(summary(ALL.LOSS.2070.POP))


######################################################
## Run GAMs of species losses vs MAT for 2030, ALL SUA
ALL.LOSS.2030.POP = gam(SPECIES_LOSS ~ s (CURRENT_MAT, k = 3)+ s(AREASQKM16, k = 3) + s(POP_2017, k = 3),
                        data = allSUAsLOSS2030, 
                        method = "REML")
print(summary(ALL.LOSS.2030.POP))


#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################