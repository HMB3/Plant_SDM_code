#########################################################################################################################
## If we are just analysing the biggest SUAs, restrict them. If not, 
if(SUAs    == "LARGE_SUAs") {
  
  ## Save basic results and SUA results to file
  message('Create barplots for the largest SUAs')
  SUA.PLOT.30.M   = subset(SUA.PLOT.30.M,   AREASQKM16 > 200 & POP_2017 > 80000)
  SUA.30.M.LOSS   = subset(SUA.30.M.LOSS,   AREASQKM16 > 200 & POP_2017 > 80000)
  SUA.30.M.GAIN   = subset(SUA.30.M.GAIN,   AREASQKM16 > 200 & POP_2017 > 80000)
  SUA.30.M.STABLE = subset(SUA.30.M.STABLE, AREASQKM16 > 200 & POP_2017 > 80000)
  
  SUA.PLOT.70.M   = subset(SUA.PLOT.70.M,   AREASQKM16 > 200 & POP_2017 > 80000)
  SUA.70.M.LOSS   = subset(SUA.70.M.LOSS,   AREASQKM16 > 200 & POP_2017 > 80000)
  SUA.70.M.GAIN   = subset(SUA.70.M.GAIN,   AREASQKM16 > 200 & POP_2017 > 80000)
  SUA.70.M.STABLE = subset(SUA.70.M.STABLE, AREASQKM16 > 200 & POP_2017 > 80000)
  
} else {
  
  message('Create barplots for the smallest SUAs')   ##
  
}
  




#########################################################################################################################
## Use ggplot to plot the percentage of species inside an LGA, which is being lost or gained
## How to calculate the gain/loss? Could do :
## Final - Now / Now *100. OR
## Gain % = gain/(stable + lost)   * 100
## Lost % = lost/(stable + lost)   * 100


#########################################################################################################################
## Create PNG output for all SUAs for 2030, ordered by mean annual temperature
png(sprintf('output/figures/SUA_percent/SUA_BAR_PLOT_SUM_100_%s_%s_%s_%s.png', 2030, SUA_ORDER_1, SUAs, SUA_SPP),      
    10, 8, units = 'in', res = 500)

## 2030
ggplot(SUA.PLOT.30.M,  aes(x = reorder(SUA, eval(parse(text = SUA_ORDER_1))), fill = AREA_CHANGE)) + 
  
  ## Plot species being lost from each SUA
  geom_bar(data        = subset(SUA.PLOT.30.M, AREA_CHANGE %in% c("LOSS")),
           
           ## Create % lost: the species lost, divided by those that were both lost + those always there
           aes(y = -(SPECIES_COUNT/(SUA.30.M.LOSS$SPECIES_COUNT + SUA.30.M.STABLE$SPECIES_COUNT + SUA.30.M.GAIN$SPECIES_COUNT))),
           position = "stack", stat = "identity") +
  
  ## Plot species being gained inside each SUA
  geom_bar(data  = subset(SUA.PLOT.30.M, AREA_CHANGE %in% c("GAIN")), 
           
           ## Create % gained: the species gained, divided by those that were both lost + those always there
           aes(y = (SPECIES_COUNT/(SUA.30.M.LOSS$SPECIES_COUNT + SUA.30.M.STABLE$SPECIES_COUNT + SUA.30.M.GAIN$SPECIES_COUNT))),
           position = "stack", stat = "identity") +
  
  scale_y_continuous(labels = scales::percent) +
  
  scale_fill_manual(values = rev(colorRampPalette(c('brown1', 'seagreen3'))(2))) +
  
  labs(#title = "Predicted native species gain/loss (number) within SUAs to 2030", 
    x = paste0("SUA by increasing ", SUA_ORDER_1), y = "Species change")  +
  
  ## Format axes
  theme(axis.title.x     = element_text(face = "bold", colour = "black", size = 15),
        axis.text.x      = element_text(angle = 90, vjust = 0.5, size = 8),
        axis.title.y     = element_text(face = "bold", colour = "black", size = 15),
        axis.text.y      = element_text(vjust = 0.5, size = 12),
        title            = element_text(face = "bold", colour = "black", size = 15),
        legend.title     = element_text(face = "bold", colour = "black", size = 12),
        legend.text      = element_text(face = "bold", size = 12),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border     = element_rect(colour = "black", fill = NA, size = 2)) +
  
  guides(fill = guide_legend(title = "Current vs. 2030"))


dev.off()


#########################################################################################################################
## Create PNG output for all SUAs for 2070, ordered by mean annual temperature
png(sprintf('output/figures/SUA_percent/SUA_BAR_PLOT_SUM_100_%s_%s_%s_%s.png', 2070, SUA_ORDER_1, SUAs, SUA_SPP),      
    10, 8, units = 'in', res = 500)

## 2030
ggplot(SUA.PLOT.70.M,  aes(x = reorder(SUA, eval(parse(text = SUA_ORDER_1))), fill = AREA_CHANGE)) + 
  
  ## Plot species being lost from each SUA
  geom_bar(data        = subset(SUA.PLOT.70.M, AREA_CHANGE %in% c("LOSS")),
           
           ## Create % lost: the species lost, divided by those that were both lost + those always there
           aes(y = -(SPECIES_COUNT/(SUA.70.M.LOSS$SPECIES_COUNT + SUA.70.M.STABLE$SPECIES_COUNT + SUA.70.M.GAIN$SPECIES_COUNT))),
           position = "stack", stat = "identity") +
  
  ## Plot species being gained inside each SUA
  geom_bar(data  = subset(SUA.PLOT.70.M, AREA_CHANGE %in% c("GAIN")), 
           
           ## Create % gained: the species gained, divided by those that were both lost + those always there
           aes(y = (SPECIES_COUNT/(SUA.70.M.LOSS$SPECIES_COUNT + SUA.70.M.STABLE$SPECIES_COUNT + SUA.70.M.GAIN$SPECIES_COUNT))),
           position = "stack", stat = "identity") +
  
  scale_y_continuous(labels = scales::percent) +
  
  scale_fill_manual(values = rev(colorRampPalette(c('brown1', 'seagreen3'))(2))) +
  
  labs(#title = "Predicted native species gain/loss (number) within SUAs to 2070", 
    x = paste0("SUA by increasing ", SUA_ORDER_1), y = "Species change")  +
  
  ## Format axes
  theme(axis.title.x     = element_text(face = "bold", colour = "black", size = 15),
        axis.text.x      = element_text(angle = 90, vjust = 0.5, size = 8),
        axis.title.y     = element_text(face = "bold", colour = "black", size = 15),
        axis.text.y      = element_text(vjust = 0.5, size = 12),
        title            = element_text(face = "bold", colour = "black", size = 15),
        legend.title     = element_text(face = "bold", colour = "black", size = 12),
        legend.text      = element_text(face = "bold", size = 12),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border     = element_rect(colour = "black", fill = NA, size = 2)) +
  
  guides(fill = guide_legend(title = "Current vs. 2070"))

dev.off()





#########################################################################################################################
## CREATE A PANEL OF SCATTERPLOTS FOR PREDICTED SPECIES GAINS AND LOSSES
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
GAIN.2070.GAM = gam(SPECIES_GAIN ~ s (CURRENT_MAT, k = 5), 
                    data = SUA.GAIN.2070, 
                    method = "REML")
summary(GAIN.2070.GAM)[["dev.expl"]]               


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





#########################################################################################################################
## OUTSTANDING PLOT TASKS:
#########################################################################################################################