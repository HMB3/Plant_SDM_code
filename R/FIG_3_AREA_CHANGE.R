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
SUA.PLOT  = readRDS(paste0(DATA_path, 'SUA_GAIN_LOSS_', save_run, '.rds'))
unique(SUA.PLOT$AREA_CHANGE)
unique(SUA.PLOT$PERIOD)
unique(SUA.PLOT$ClimateZ)


## Subset data
allSUAsGAIN <- subset(SUA.PLOT, AREA_CHANGE=="GAIN") 
allSUAsLOSS <- subset(SUA.PLOT, AREA_CHANGE=="LOSS")

bigSUAsGAIN <- subset(SUA.PLOT, AREASQKM16>200 & POP_2017>80000 & AREA_CHANGE=="GAIN")
bigSUAsLOSS <- subset(SUA.PLOT, AREASQKM16>200 & POP_2017>80000 & AREA_CHANGE=="LOSS")


## To calculate GAM deviance
allSUAsGAIN2070 <- subset(SUA.PLOT, PERIOD=="2070" & AREA_CHANGE=="GAIN")
allSUAsLOSS2070 <- subset(SUA.PLOT, PERIOD=="2070" & AREA_CHANGE=="LOSS")
allSUAsGAIN2030 <- subset(SUA.PLOT, PERIOD=="2030" & AREA_CHANGE=="GAIN")
allSUAsLOSS2030 <- subset(SUA.PLOT, PERIOD=="2030" & AREA_CHANGE=="LOSS")

bigSUAsGAIN2070 <- subset(SUA.PLOT, PERIOD=="2070" & AREA_CHANGE=="GAIN" & AREASQKM16>200 & POP_2017>80000)
bigSUAsLOSS2070 <- subset(SUA.PLOT, PERIOD=="2070" & AREA_CHANGE=="LOSS" & AREASQKM16>200 & POP_2017>80000)
bigSUAsGAIN2030 <- subset(SUA.PLOT, PERIOD=="2030" & AREA_CHANGE=="GAIN" & AREASQKM16>200 & POP_2017>80000)
bigSUAsLOSS2030 <- subset(SUA.PLOT, PERIOD=="2030" & AREA_CHANGE=="LOSS" & AREASQKM16>200 & POP_2017>80000)





#########################################################################################################################
## 2). RUN GAMs
#########################################################################################################################


####### Figure 3 #####
### make subpolots ###
SUA.plot.cols = brewer.pal(12, "Paired")



### Fig 3A ####################################################################################
#####################################################
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


######################################################
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
ggsave(filename="C:/Users/MQ20174608/Desktop/fig3A.pdf", width=10, height=10, plot=fig3A, dpi=300)





### Fig 3B ####################################################################################
#####################################################
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


######################################################

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
ggsave(filename="C:/Users/MQ20174608/Desktop/fig3B.pdf", width=10, height=10, plot=fig3B, dpi=300)





### Fig 3C ####################################################################################
#####################################################
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


######################################################
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
ggsave(filename="C:/Users/MQ20174608/Desktop/fig3C.pdf", width=10, height=10, plot=fig3C, dpi=300)





### Fig 3D ####################################################################################
#####################################################
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


######################################################
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
ggsave(filename="C:/Users/MQ20174608/Desktop/fig3D.pdf", width=10, height=10, plot=fig3D, dpi=300)





#########################################################################################################################
## 3). CREATE A PANEL FOR ALL FIGURES
#########################################################################################################################


#########################################################################################################################
dev.new(width = 40, height = 30)
fig3<- grid_arrange_shared_legend(fig3A, fig3B, fig3C, fig3D, ncol = 2, nrow = 2, position="right")
fig3

ggsave(filename="C:/Users/MQ20174608/Desktop/fig3.pdf", width=20, height=18, plot=fig3, dpi=300)





#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################