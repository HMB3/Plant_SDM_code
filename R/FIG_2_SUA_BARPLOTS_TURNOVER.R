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
## Read in Linda's data. Could calculate 
## https://cran.r-project.org/web/packages/codyn/vignettes/Temporal_Diversity_Indices.html
SUA_TURNOVER = read.csv("./data/ANALYSIS/SUA_TURNOVER.csv", stringsAsFactors = FALSE)
View(SUA_TURNOVER)
View(SUA.PLOT.30.M)


#########################################################################################################################
## How to calculate turnover from Linda's data?
## 


#########################################################################################################################
## Use ggplot to plot the percentage of species inside an LGA, which is being lost or gained
## How to calculate the gain/loss? Could do :
## Final - Now / Now *100. OR
## Gain % = gain/(stable + lost)   * 100
## Lost % = lost/(stable + lost)   * 100


## Make all the figures in the MS colorblind safe. The "Paried" color scheme is aparently safe.
## No. 1 = light blue, 4 = green, and 8 = dark orange 
SUA.plot.cols = brewer.pal(12, "Paired")


#########################################################################################################################
## Create PNG output for all SUAs for 2030, ordered by mean annual temperature
png(sprintf('output/figures/FIG_2/SUA_BAR_PLOT_SUM_100_%s_%s_%s_%s.png', 2030, SUA_ORDER, SUAs, SUA_SPP),      
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
  
  scale_fill_manual(values = rev(colorRampPalette(c(SUA.plot.cols[8], SUA.plot.cols[4]))(2))) +
  
  labs(#title = "Predicted native species gain/loss (number) within SUAs to 2030", 
    x = "", #paste0("SUA by increasing ", SUA_ORDER),
    y = "") + # "Species change")  +
  
  ## Format axes
  theme(#axis.title.x     = element_text(face = "bold", colour = "black", size = 15),
    #axis.text.x      = element_text(angle = 90, vjust = 0.5, size = 8),
    #axis.title.y     = element_text(face = "bold", colour = "black", size = 15),
    axis.text.y      = element_text(vjust = 0.5, size = 30),
    #title            = element_text(face = "bold", colour = "black", size = 15),
    #legend.title     = element_text(face = "bold", colour = "black", size = 12),
    #legend.text      = element_text(face = "bold", size = 12),
    axis.title.x     = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),  
    legend.position  = "none",
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    panel.border     = element_rect(colour = "black", fill = NA, size = 2)) #+

#guides(fill = guide_legend(title = "Current vs. 2030"))

dev.off()


#########################################################################################################################
## Create PNG output for all SUAs for 2070, ordered by mean annual temperature
png(sprintf('output/figures/FIG_2/SUA_BAR_PLOT_SUM_100_%s_%s_%s_%s.png', 2070, SUA_ORDER, SUAs, SUA_SPP),      
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
  
  scale_fill_manual(values = rev(colorRampPalette(c(SUA.plot.cols[8], SUA.plot.cols[4]))(2))) +
  
  labs(#title = "Predicted native species gain/loss (number) within SUAs to 2030", 
    x = "", #paste0("SUA by increasing ", SUA_ORDER),
    y = "") + # "Species change")  +
  
  ## Format axes
  theme(#axis.title.x     = element_text(face = "bold", colour = "black", size = 15),
    #axis.text.x      = element_text(angle = 90, vjust = 0.5, size = 8),
    #axis.title.y     = element_text(face = "bold", colour = "black", size = 15),
    axis.text.y      = element_text(vjust = 0.5, size = 30),
    #title            = element_text(face = "bold", colour = "black", size = 15),
    #legend.title     = element_text(face = "bold", colour = "black", size = 12),
    #legend.text      = element_text(face = "bold", size = 12),
    axis.title.x     = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),  
    legend.position  = "none",
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    panel.border     = element_rect(colour = "black", fill = NA, size = 2)) #+

#guides(fill = guide_legend(title = "Current vs. 2030"))

dev.off()





## How can we calcualte turnover with this dataset? 
## No. new species / baseline % (E.G.)


#########################################################################################################################
## Create PNG output for all SUAs for 2030, ordered by mean annual temperature
png(sprintf('output/figures/FIG_2/SUA_BAR_PLOT_SUM_100_%s_%s_%s_%s.png', 2030, SUA_ORDER, SUAs, SUA_SPP),      
    10, 8, units = 'in', res = 500)

## 2030
ggplot(SUA_TURNOVER,  aes(x = reorder(SUA, eval(parse(text = SUA_ORDER))), fill = AREA_CHANGE)) + 
  
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
  
  scale_fill_manual(values = rev(colorRampPalette(c(SUA.plot.cols[8], SUA.plot.cols[4]))(2))) +
  
  labs(#title = "Predicted native species gain/loss (number) within SUAs to 2030", 
    x = "", #paste0("SUA by increasing ", SUA_ORDER),
    y = "") + # "Species change")  +
  
  ## Format axes
  theme(#axis.title.x     = element_text(face = "bold", colour = "black", size = 15),
    #axis.text.x      = element_text(angle = 90, vjust = 0.5, size = 8),
    #axis.title.y     = element_text(face = "bold", colour = "black", size = 15),
    axis.text.y      = element_text(vjust = 0.5, size = 30),
    #title            = element_text(face = "bold", colour = "black", size = 15),
    #legend.title     = element_text(face = "bold", colour = "black", size = 12),
    #legend.text      = element_text(face = "bold", size = 12),
    axis.title.x     = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),  
    legend.position  = "none",
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    panel.border     = element_rect(colour = "black", fill = NA, size = 2)) #+

#guides(fill = guide_legend(title = "Current vs. 2030"))

dev.off()




#########################################################################################################################
## OUTSTANDING PLOT TASKS:
#########################################################################################################################