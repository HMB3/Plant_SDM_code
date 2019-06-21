#########################################################################################################################
## Create PNG output for all SUAs for 2030, ordered by mean annual PET.......
png(sprintf('output/figures/SUA_percent/ALL_SUA_BAR_PLOT_SUM_100_%s_%s.png', 2030, 'PET'),      
    10, 8, units = 'in', res = 500)

## 2030
ggplot(SUA.PLOT.30.M,  aes(x = reorder(SUA, CURRENT_PET), fill = AREA_CHANGE)) + 
  
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
    x = "SUA by increasing PET", y = "Species change")  +
  
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
## Create PNG output for all SUAs for 2070, ordered by PET
png(sprintf('output/figures/SUA_percent/ALL_SUA_BAR_PLOT_SUM_100_%s_%s.png', 2070, 'PET'),      
    10, 8, units = 'in', res = 500)

## 2070
ggplot(SUA.PLOT.70.M,  aes(x = reorder(SUA, CURRENT_PET), fill = AREA_CHANGE)) + 
  
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
    x = "SUA by increasing PET", y = "Species change")  +
  
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
## Create PNG output for all SUAs for 2070, ordered by area
png(sprintf('output/figures/SUA_percent/ALL_SUA_BAR_PLOT_SUM_100_%s_%s.png', 2070, 'area'),      
    10, 8, units = 'in', res = 500)

## 2070
ggplot(SUA.PLOT.70.M,  aes(x = reorder(SUA, AREASQKM16), fill = AREA_CHANGE)) + 
  
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
    x = "SUA by increasing Area", y = "Species change")  +
  
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
## Try plotting a subset of just the mapped SUAs
SUA.PLOT.30.M$SUA = gsub("Canberra - Queanbeyan", "Canberra", SUA.PLOT.30.M$SUA)
SUA.PLOT.50.M$SUA = gsub("Canberra - Queanbeyan", "Canberra", SUA.PLOT.50.M$SUA)
SUA.PLOT.70.M$SUA = gsub("Canberra - Queanbeyan", "Canberra", SUA.PLOT.70.M$SUA)

SUA.PLOT.MAJOR.30 = SUA.PLOT.30.M[SUA.PLOT.30.M$SUA %in% MAP_SUA, ]
SUA.PLOT.MAJOR.50 = SUA.PLOT.50.M[SUA.PLOT.50.M$SUA %in% MAP_SUA, ] 
SUA.PLOT.MAJOR.70 = SUA.PLOT.70.M[SUA.PLOT.70.M$SUA %in% MAP_SUA, ] 
unique(SUA.PLOT.MAJOR.30$SUA);unique(SUA.PLOT.MAJOR.70$SUA)


## Create a subset of just major cities
SUA.30.MJ.LOSS        = subset(SUA.PLOT.MAJOR.30, AREA_CHANGE %in% c("LOSS"))
SUA.30.MJ.GAIN        = subset(SUA.PLOT.MAJOR.30, AREA_CHANGE %in% c("GAIN"))
SUA.30.MJ.STABLE      = subset(SUA.PLOT.MAJOR.30, AREA_CHANGE %in% c("STABLE"))

SUA.50.MJ.LOSS        = subset(SUA.PLOT.MAJOR.50, AREA_CHANGE %in% c("LOSS"))
SUA.50.MJ.GAIN        = subset(SUA.PLOT.MAJOR.50, AREA_CHANGE %in% c("GAIN"))
SUA.50.MJ.STABLE      = subset(SUA.PLOT.MAJOR.50, AREA_CHANGE %in% c("STABLE"))

SUA.70.MJ.LOSS        = subset(SUA.PLOT.MAJOR.70, AREA_CHANGE %in% c("LOSS"))
SUA.70.MJ.GAIN        = subset(SUA.PLOT.MAJOR.70, AREA_CHANGE %in% c("GAIN"))
SUA.70.MJ.STABLE      = subset(SUA.PLOT.MAJOR.70, AREA_CHANGE %in% c("STABLE"))


## Change % cacl to include gains when plotting SUA_RECORDS >0


#########################################################################################################################
## Create PNG output for Major captuials, 2030
png(sprintf('output/figures/SUA_percent/CAPITAL_SUA_BAR_PLOT_SUM_100_%s_%s.png', 2030, 'PET'),      
    10, 8, units = 'in', res = 500)

## 2030
ggplot(SUA.PLOT.MAJOR.30,  aes(x = reorder(SUA, CURRENT_PET), fill = AREA_CHANGE)) + 
  
  ## Plot species being lost from each SUA
  geom_bar(data        = subset(SUA.PLOT.MAJOR.30, AREA_CHANGE %in% c("LOSS")),
           
           ## Create % lost: the species lost, divided by those that were both lost + those always there
           aes(y = -(SPECIES_COUNT/(SUA.30.MJ.LOSS$SPECIES_COUNT + SUA.30.MJ.STABLE$SPECIES_COUNT + SUA.30.MJ.GAIN$SPECIES_COUNT))),
           position = "stack", stat = "identity") +
  
  ## Plot species being gained inside each SUA
  geom_bar(data  = subset(SUA.PLOT.MAJOR.30, AREA_CHANGE %in% c("GAIN")), 
           
           ## Create % gained: the species gained, divided by those that were both lost + those always there
           aes(y = (SPECIES_COUNT/(SUA.30.MJ.LOSS$SPECIES_COUNT + SUA.30.MJ.STABLE$SPECIES_COUNT + SUA.30.MJ.GAIN$SPECIES_COUNT))),
           position = "stack", stat = "identity") +
  
  scale_y_continuous(labels = scales::percent) +
  
  scale_fill_manual(values = rev(colorRampPalette(c('brown1', 'seagreen3'))(2))) +
  
  labs(x = "SUA by increasing PET", y = "Species change")  +
  
  ## Format axes
  theme(axis.title.x     = element_text(face = "bold", colour = "black", size = 15),
        axis.text.x      = element_text(angle = 45, vjust = 0.5, size = 12),
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
## Create PNG output for Major captuials, 2070
png(sprintf('output/figures/SUA_percent/CAPITAL_SUA_BAR_PLOT_SUM_100_%s_%s.png', 2070, 'PET'),      
    10, 8, units = 'in', res = 500)

## 2070
ggplot(SUA.PLOT.MAJOR.70,  aes(x = reorder(SUA, CURRENT_PET), fill = AREA_CHANGE)) + 
  
  ## Plot species being lost from each SUA
  geom_bar(data        = subset(SUA.PLOT.MAJOR.70, AREA_CHANGE %in% c("LOSS")),
           
           ## Create % lost: the species lost, divided by those that were both lost + those always there
           aes(y = -(SPECIES_COUNT/(SUA.70.MJ.LOSS$SPECIES_COUNT + SUA.70.MJ.STABLE$SPECIES_COUNT + SUA.70.MJ.GAIN$SPECIES_COUNT))),
           position = "stack", stat = "identity") +
  
  ## Plot species being gained inside each SUA
  geom_bar(data  = subset(SUA.PLOT.MAJOR.70, AREA_CHANGE %in% c("GAIN")), 
           
           ## Create % gained: the species gained, divided by those that were both lost + those always there
           aes(y = (SPECIES_COUNT/(SUA.70.MJ.LOSS$SPECIES_COUNT + SUA.70.MJ.STABLE$SPECIES_COUNT + SUA.70.MJ.GAIN$SPECIES_COUNT))),
           position = "stack", stat = "identity") +
  
  scale_y_continuous(labels = scales::percent) +
  
  scale_fill_manual(values = rev(colorRampPalette(c('brown1', 'seagreen3'))(2))) +
  
  labs(x = "SUA by increasing PET", y = "Species change")  +
  
  ## Format axes
  theme(axis.title.x     = element_text(face = "bold", colour = "black", size = 15),
        axis.text.x      = element_text(angle = 45, vjust = 0.5, size = 12),
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
## Create PNG output for Major captuials, 2070
png(sprintf('output/figures/SUA_percent/CAPITAL_SUA_BAR_PLOT_SUM_100_%s_%s.png', 2070, 'area'),      
    10, 8, units = 'in', res = 500)

## 2070
ggplot(SUA.PLOT.MAJOR.70,  aes(x = reorder(SUA, AREASQKM16), fill = AREA_CHANGE)) + 
  
  ## Plot species being lost from each SUA
  geom_bar(data        = subset(SUA.PLOT.MAJOR.70, AREA_CHANGE %in% c("LOSS")),
           
           ## Create % lost: the species lost, divided by those that were both lost + those always there
           aes(y = -(SPECIES_COUNT/(SUA.70.MJ.LOSS$SPECIES_COUNT + SUA.70.MJ.STABLE$SPECIES_COUNT + SUA.70.MJ.GAIN$SPECIES_COUNT))),
           position = "stack", stat = "identity") +
  
  ## Plot species being gained inside each SUA
  geom_bar(data  = subset(SUA.PLOT.MAJOR.70, AREA_CHANGE %in% c("GAIN")), 
           
           ## Create % gained: the species gained, divided by those that were both lost + those always there
           aes(y = (SPECIES_COUNT/(SUA.70.MJ.LOSS$SPECIES_COUNT + SUA.70.MJ.STABLE$SPECIES_COUNT + SUA.70.MJ.GAIN$SPECIES_COUNT))),
           position = "stack", stat = "identity") +
  
  scale_y_continuous(labels = scales::percent) +
  
  scale_fill_manual(values = rev(colorRampPalette(c('brown1', 'seagreen3'))(2))) +
  
  labs(x = "SUA by increasing Area", y = "Species change")  +
  
  ## Format axes
  theme(axis.title.x     = element_text(face  = "bold", colour = "black", size = 15),
        axis.text.x      = element_text(angle = 45, vjust = 0.5, size = 12),
        axis.title.y     = element_text(face  = "bold", colour = "black", size = 15),
        axis.text.y      = element_text(vjust = 0.5, size = 12),
        title            = element_text(face  = "bold", colour = "black", size = 15),
        legend.title     = element_text(face  = "bold", colour = "black", size = 12),
        legend.text      = element_text(face  = "bold", size = 12),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border     = element_rect(colour = "black", fill = NA, size = 2)) +
  
  guides(fill = guide_legend(title = "Current vs. 2070"))

dev.off()