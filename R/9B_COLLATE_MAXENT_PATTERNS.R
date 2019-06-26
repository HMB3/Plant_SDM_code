#########################################################################################################################
############################################# SUMARISE MAXENT RESULTS ################################################### 
#########################################################################################################################



#########################################################################################################################
## 1). COMBINE SDM RESULTS WITH HIA CONTEXT DATA, & VISUALISE THE RELATIONSHIPS BETWEEN VARIABLES
#########################################################################################################################


#########################################################################################################################
## So we've got a range of quantitative and qualitative measures for each species (e.g. maxent output, number of records,
## model score). Before embarking on a machine learning odessy, can we see some patterns in the 450-odd species that have
## already been rated by hand? Simple categorical visualations in R should be sufficient to enough to get at the patterns.


## Do the relationships between model performance variables vary between models rated as good, fair and poor?


#########################################################################################################################
## Grouped boxplot :: 
## Are the total number of maxent training records related to the model performance?
MXT.PLOT = na.omit(MAXENT.CONTEXT)


## Are the number of background points related to model performance?
ggplot(MXT.PLOT, aes(x = Number_background_points, y = Max_tss, fill = MAXENT_RATING)) + 
  geom_boxplot()


#########################################################################################################################
## Simple correlations :: do the relationships between variables differ by maxent rating?
#########################################################################################################################


## Create a scatterplot maxrix for all the quantiatitve variables
## Need to minise these variables, which make the most sense?
MXT.COR = dplyr::select(MXT.PLOT, Max_tss, Logistic_threshold, Aus_records, KOP_count,
                        Number_var, Perm_imp, KOP_count, MAXENT_RATING)


## Try using ggpairs to summarise all the data :: need to tweak the settings, to create more informative plot
## Also, set the colours to good, fair etc
message('Creating plot of maxent results for ', length(MXT.PLOT$searchTaxon), ' species')
png(sprintf("./output/maxent/maxent_results_patterns_%s.png", save_run),
    16, 10, units = 'in', res = 500)

## Use the ggpairs default to create a plot of all the points
ggpairs(MXT.COR, aes(colour = MAXENT_RATING, alpha = 0.4)) +
  
  ## Use the classic theme
  theme_classic() +
  
  ## Change the axes sizes, etc.
  theme(axis.title.x     = element_text(colour = "black", size = 25),
        axis.text.x      = element_text(size = 12),
        
        axis.title.y     = element_text(colour = "black", size = 25),
        axis.text.y      = element_text(size = 12),
        
        panel.background = element_blank(),
        panel.border     = element_rect(colour = "black", fill = NA, size = 1.2),
        plot.title       = element_text(size   = 20, face = "bold"),
        legend.text      = element_text(size   = 10),
        legend.title     = element_text(size   = 10),
        legend.key.size  = unit(1.5, "cm")) +
  
  ## And title
  ggtitle(paste0("Relationships between maxent indicators for ", 
                 length(unique(MXT.PLOT$searchTaxon)), 
                 ' horticultural species'))

dev.off()


## Try a multiple regression, informed by the visualisations :: this doesn't account for maxent rating
## Maybe try running GAMs on each subset - good, fair, poor, or using 'rating' as factor
MXT.GD = subset(MXT.PLOT, MAXENT_RATING == "GOOD")
MXT.FA = subset(MXT.PLOT, MAXENT_RATING == "FAIR")
MXT.PR = subset(MXT.PLOT, MAXENT_RATING == "POOR")


## All ratings
summary(gam(Max_tss ~ s(Number_background_points, k = 3) +
              s(KOP_count, k = 3) +
              s(Maxent_records, k = 3) +
              s(Aus_records, k = 3) +
              s(Logistic_threshold, k = 3) +
              s(AOO, k = 3),
            data = MXT.PLOT))[24]


## Poor models
summary(gam(Max_tss ~ s(Number_background_points, k = 3) +
              s(KOP_count, k = 3) +
              s(Maxent_records, k = 3) +
              s(Aus_records, k = 3) +
              s(Logistic_threshold, k = 3) +
              s(AOO, k = 3),
            data = MXT.PR))[24]


#########################################################################################################################
## Consider how else to represent the relationships?


#########################################################################################################################
##################################################### TBC ############################################################### 
#########################################################################################################################