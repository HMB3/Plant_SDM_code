#########################################################################################################################
################################# SUMMARISE SPEICES GAINS BY SUA ######################################################## 
#########################################################################################################################


## This code creates barplots of species gained/lost in each of the 101 SUAs.


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
  
  message('Create barplots for all the SUAs')   ##
  
}
  

#########################################################################################################################
## If only analysing temperate SUAs, remove the others
if(KOP_ZONE == "TEMPERATE") {
  
  ## Save basic results and SUA results to file
  message('Analyse only temperate SUAs')
  non_temperate = c("Am", "Aw", "BSh", "BSk", "BWh")
  temperate     = c("Cfa", "Cfb", "Csa", "Csb", "Cwa")
  
  SUA.PLOT.70.M   = SUA.PLOT.70.M[SUA.PLOT.70.M$ClimateZ %in% temperate , ]
  
  SUA.PLOT.30.M   = SUA.PLOT.30.M[SUA.PLOT.30.M$ClimateZ %in% temperate , ]
  SUA.30.M.LOSS   = SUA.30.M.LOSS[SUA.30.M.LOSS$ClimateZ %in% temperate , ]
  SUA.30.M.GAIN   = SUA.30.M.GAIN[SUA.30.M.GAIN$ClimateZ %in% temperate , ]
  SUA.30.M.STABLE = SUA.30.M.STABLE[SUA.30.M.STABLE$ClimateZ %in% temperate , ]
  
  SUA.PLOT.70.M   = SUA.PLOT.70.M[SUA.PLOT.70.M$ClimateZ %in% temperate , ]
  SUA.70.M.LOSS   = SUA.70.M.LOSS[SUA.70.M.LOSS$ClimateZ %in% temperate , ]
  SUA.70.M.GAIN   = SUA.70.M.GAIN[SUA.70.M.GAIN$ClimateZ %in% temperate , ]
  SUA.70.M.STABLE = SUA.70.M.STABLE[SUA.70.M.STABLE$ClimateZ %in% temperate , ]
  
} else {
  
  message('Analyse all SUAs') 
  
}



#########################################################################################################################
## Use ggplot to plot the percentage of species inside an LGA, which is being lost or gained
## How to calculate the gain/loss? Could do :
## Final - Now / Now *100. OR
## Gain % = gain/(stable + lost)   * 100
## Lost % = lost/(stable + lost)   * 100


## Also, consider if we should exclude SUAs with < 2 species. This uses a somewhat different column from the data used
## to make the scatterplots, but check if it's analagous.


#########################################################################################################################
## Create PNG output for all SUAs for 2030, ordered by mean annual temperature
png(sprintf('output/figures/SUA_percent/SUA_BAR_PLOT_SUM_100_%s_%s_%s_%s.png', 2030, SUA_ORDER, SUAs, SUA_SPP),      
    10, 8, units = 'in', res = 500)

## 2030
ggplot(SUA.PLOT.30.M,  aes(x = reorder(SUA, eval(parse(text = SUA_ORDER))), fill = AREA_CHANGE)) + 
  
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
    x = paste0("SUA by increasing ", SUA_ORDER), y = "Species change")  +
  
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
png(sprintf('output/figures/SUA_percent/SUA_BAR_PLOT_SUM_100_%s_%s_%s_%s_%s.png', 2070, SUA_ORDER, SUAs, KOP_ZONE, SUA_SPP),      
    10, 8, units = 'in', res = 500)

## 2070
ggplot(SUA.PLOT.70.M,  aes(x = reorder(SUA, eval(parse(text = SUA_ORDER))), fill = AREA_CHANGE)) + 
  
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
    x = paste0("SUA by increasing ", SUA_ORDER), y = "Species change")  +
  
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
## OUTSTANDING PLOT TASKS:
#########################################################################################################################