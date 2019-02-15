#########################################################################################################################
############################################### FIG 2 STOTEN ############################################################ 
#########################################################################################################################


# Code for Producing Fig.2 STOTEN Manuscript - Burley et al., 2019
# Which Plant Where Project - Macquarie University Sydney, Australia
# Ossola Alessandro, 29.01.2019


## The default link function for mgcv::gam is Gaussian, and that can be problematic when the response is a percentage/proportion, due to the [0, 1] bounds
## try using the 'beta' distribution in ggplot 
## https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/Beta.html
## https://ggplot2.tidyverse.org/reference/geom_smooth.html


## EG
# ALL.GAIN.2070.GAM = gam(Species.gain/100 ~ s (CURRENT_MAT, k = 3), 
#                         data = allSUAs2070, 
#                         method = "REML", family = betar(link = "logit"))
# summary(ALL.GAIN.2070.GAM)[["dev.expl"]][1]




#########################################################################################################################
## 1). READ DATA
#########################################################################################################################


## prepare data
turnover =  read.csv("./data/ANALYSIS/Fig2.csv", sep = ",", header=T, na.string = "NA")
allSUAs <- subset(turnover, LARGE_SUA=="n") 
bigSUAs <- subset(turnover, LARGE_SUA=="y") 


## to calculate GAM deviance
allSUAs2070 <- subset(allSUAs, PERIOD=="2070")
allSUAs2030 <- subset(allSUAs, PERIOD=="2030")
bigSUAs2070 <- subset(bigSUAs, PERIOD=="2070")
bigSUAs2030 <- subset(bigSUAs, PERIOD=="2030")





#########################################################################################################################
## 2). RUN GAMS FOR EACH PANEL
#########################################################################################################################


#########################################################################################################################
## Fig 2A 
## Run GAMs of species gains vs MAT for 2070, ALL SUA
ALL.GAIN.2070.GAM = gam(Species.gain ~ s (CURRENT_MAT, k = 3), 
                        data = allSUAs2070, 
                        method = "REML")
summary(ALL.GAIN.2070.GAM)[["dev.expl"]][1] 


######################################################
## Run GAMs of species gains vs MAT for 2030, ALL SUA
ALL.GAIN.2030.GAM = gam(Species.gain ~ s (CURRENT_MAT, k = 3), 
                        data = allSUAs2030, 
                        method = "REML")
summary(ALL.GAIN.2030.GAM)[["dev.expl"]][1]  


#########################################################################################################################
## 
dev.new(width = 17, height = 13)


## Plot the gams 
fig2A<-ggplot(allSUAs, aes(x=CURRENT_MAT, y=Species.gain, color=factor(PERIOD))) +
  
  geom_point(size = 6) + theme_bw() +  theme(axis.title.x=element_text(margin = margin(t = 20)), axis.title.y=element_text(margin = margin(r = 20)),
                                             axis.text.x=element_text(margin = margin(t = 10)), axis.text.y=element_text(margin = margin(r = 10)), 
                                             legend.box.margin=margin(10,10,10,10)) +
  
  geom_smooth(method="gam", formula = y ~ s(x, k = 3), se=FALSE, fullrange=TRUE, size = 2) +
  xlab(bquote('')) + ylim(0,70) + ylab(bquote('New species gained (%)')) + ggtitle("All SUAs") + theme(text = element_text(size=30)) + 
  
  scale_color_manual(values=c("darkgrey", "black"), name = "Period") + 
  annotate(geom="text", size=8, x=20, y=65, label=paste("Deviance (2030) =", signif(summary(ALL.GAIN.2030.GAM)[["dev.expl"]][1]*100, digits = 3),"%")) +
  annotate(geom="text", size=8, x=20, y=62, label=paste("Deviance (2070) =", signif(summary(ALL.GAIN.2070.GAM)[["dev.expl"]][1]*100, digits = 3),"%")) +
  
  theme(legend.position="none")  #this remove legend on the right

fig2A






#########################################################################################################################
## Fig 2B 
## Run GAMs of species gains vs MAT for 2070, BIG SUA
BIG.GAIN.2070.GAM = gam(Species.gain ~ s (CURRENT_MAT, k = 3), 
                        data = bigSUAs2070, 
                        method = "REML")
summary(BIG.GAIN.2070.GAM)[["dev.expl"]][1] 


######################################################
## Run GAMs of species gains vs MAT for 2030, BIG SUA
BIG.GAIN.2030.GAM = gam(Species.gain ~ s (CURRENT_MAT, k = 3), 
                        data = bigSUAs2030, 
                        method = "REML")
summary(BIG.GAIN.2030.GAM)[["dev.expl"]][1]  


#########################################################################################################################
dev.new(width=17, height=13)

##
fig2B<-ggplot(bigSUAs, aes(x=CURRENT_MAT, y=Species.gain, color=factor(PERIOD))) +
  
  geom_point(size = 6) + theme_bw() +  theme(axis.title.x=element_text(margin = margin(t = 20)), axis.title.y=element_text(margin = margin(r = 20)),
                                             axis.text.x=element_text(margin = margin(t = 10)), axis.text.y=element_text(margin = margin(r = 10)), 
                                             legend.box.margin=margin(10,10,10,10)) +
  
  geom_smooth(method="gam", formula = y ~ s(x, k = 3), se=FALSE, fullrange=TRUE, size = 2) +
  xlab(bquote('')) + ylim(0,70) + 
  ylab(bquote('')) + ggtitle("Largest SUAs") + theme(text = element_text(size=30)) + 
  
  scale_color_manual(values=c("darkgrey", "black"), name = "Period") + 
  annotate(geom="text", size=8, x=20, y=65, label=paste("Deviance (2030) =", signif(summary(BIG.GAIN.2030.GAM)[["dev.expl"]][1]*100, digits = 3),"%")) +
  annotate(geom="text", size=8, x=20, y=62, label=paste("Deviance (2070) =", signif(summary(BIG.GAIN.2070.GAM)[["dev.expl"]][1]*100, digits = 3),"%")) 
#theme(legend.position="none")  #this remove legend on the right

fig2B






#########################################################################################################################
## Fig 2C 
## Run GAMs of species losses vs MAT for 2070, ALL SUA
ALL.LOSS.2070.GAM = gam(Species.loss ~ s (CURRENT_MAT, k = 3), 
                        data = allSUAs2070, 
                        method = "REML")
summary(ALL.LOSS.2070.GAM)[["dev.expl"]][1]


######################################################
## Run GAMs of species losses vs MAT for 2030, ALL SUA
ALL.LOSS.2030.GAM = gam(Species.loss ~ s (CURRENT_MAT, k = 3), 
                        data = allSUAs2030, 
                        method = "REML")
summary(ALL.LOSS.2030.GAM)[["dev.expl"]][1] 


#########################################################################################################################
## 
dev.new(width=17, height=13)
fig2C<-ggplot(allSUAs, aes(x=CURRENT_MAT, y=Species.loss, color=factor(PERIOD))) +
  
  geom_point(size = 6) + theme_bw() +  theme(axis.title.x=element_text(margin = margin(t = 20)), axis.title.y=element_text(margin = margin(r = 20)),
                                             axis.text.x=element_text(margin = margin(t = 10)), axis.text.y=element_text(margin = margin(r = 10)), 
                                             legend.box.margin=margin(10,10,10,10)) +
  
  geom_smooth(method="gam", formula = y ~ s(x, k = 3), se=FALSE, fullrange=TRUE, size = 2) +
  xlab(bquote('Current MAT of SUA (1960-1990)')) + ylim(0,80) + ylab(bquote('Species lost (%)')) + theme(text = element_text(size=30)) + 
  scale_color_manual(values=c("darkgrey", "black"), name = "Period") + 
  
  annotate(geom="text", size=8, x=20, y=7, label=paste("Deviance (2030) =", signif(summary(ALL.LOSS.2030.GAM)[["dev.expl"]][1]*100, digits = 3),"%")) +
  annotate(geom="text", size=8, x=20, y=4, label=paste("Deviance (2070) =", signif(summary(ALL.LOSS.2070.GAM)[["dev.expl"]][1]*100, digits = 3),"%")) +
  theme(legend.position="none")  #this remove legend on the right

fig2C





#########################################################################################################################
### Fig 2D 
## Run GAMs of species losses vs MAT for 2070, BIG SUA
BIG.LOSS.2070.GAM = gam(Species.loss ~ s (CURRENT_MAT, k = 3), 
                        data = bigSUAs2070, 
                        method = "REML")
summary(BIG.LOSS.2070.GAM)[["dev.expl"]][1]    


######################################################
## Run GAMs of species losses vs MAT for 2030, BIG SUA
BIG.LOSS.2030.GAM = gam(Species.loss ~ s (CURRENT_MAT, k = 3), 
                        data = bigSUAs2030, 
                        method = "REML")
summary(BIG.LOSS.2030.GAM)[["dev.expl"]][1]  


#########################################################################################################################
##
dev.new(width=17, height=13)

fig2D<-ggplot(bigSUAs, aes(x=CURRENT_MAT, y=Species.loss, color=factor(PERIOD))) +
  
  geom_point(size = 6) + theme_bw() +  theme(axis.title.x=element_text(margin = margin(t = 20)), axis.title.y=element_text(margin = margin(r = 20)),
                                             axis.text.x=element_text(margin = margin(t = 10)), axis.text.y=element_text(margin = margin(r = 10)), 
                                             legend.box.margin=margin(10,10,10,10)) +
  
  geom_smooth(method="gam", formula = y ~ s(x, k = 3), se=FALSE, fullrange=TRUE, size = 2) +
  xlab(bquote('Current MAT of SUA (1960-1990)')) + ylim(0,80) +
  ylab(bquote('')) + theme(text = element_text(size=30)) + 
  
  scale_color_manual(values=c("darkgrey", "black"), name = "Period") + 
  annotate(geom="text", size=8, x=20, y=7, label=paste("Deviance (2030) =", signif(summary(BIG.LOSS.2030.GAM)[["dev.expl"]][1]*100, digits = 3),"%")) +
  annotate(geom="text", size=8, x=20, y=4, label=paste("Deviance (2070) =", signif(summary(BIG.LOSS.2070.GAM)[["dev.expl"]][1]*100, digits = 3),"%")) + 
  theme(legend.position="none")  #this remove legend on the right

fig2D






#########################################################################################################################
## 3). CREATE A PANEL FOR ALL FIGURES
#########################################################################################################################


#########################################################################################################################
## Use a function Alessandro created to stack the GPPLOTs
dev.new(width = 40, height = 30)
fig2 <- grid_arrange_shared_legend(fig2A, fig2B, fig2C, fig2D, ncol = 2, nrow = 2, position = "right")
fig2

ggsave(filename = "./output/figures/FIG_2/FIG_2_NEW_SUA_SPECIES.png", width = 20, height = 18, plot = fig2, dpi = 300)




#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################
