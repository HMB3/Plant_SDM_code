#########################################################################################################################
## 7). PLOT HISTOGRAMS AND BAR CHARTS FOR EACH SPECIES AT 1KM
#########################################################################################################################


# ##############################################################################################
# ## Plot histograms of temperature and rainfall
# ## species = spp.geo[5]
# for (species in spp.geo) {
#   
#   ## Subset the spatial dataframe into records for each species
#   SP.DF     <- NICHE.1KM.84[NICHE.1KM.84$searchTaxon %in% species , ] 
#   
#   TMP.GLO   <- subset(GLOB.NICHE,   searchTaxon == species)[c("searchTaxon", "Annual_mean_temp_q95_q05", 
#                                                               "Annual_mean_temp_q05", "Annual_mean_temp_q95")]
#   
#   TMP.AUS   <- subset(AUS.NICHE,    searchTaxon == species)[c("searchTaxon", "Annual_mean_temp_q95_q05", 
#                                                               "Annual_mean_temp_q05", "Annual_mean_temp_q95")]
#   
#   #############################################################
#   ## Now, build a df of the temperature vectors 
#   if(nrow(TMP.GLO) > 0){
#     TMP.GLO$RANGE = "GLOBAL"
#   } else {
#     message("No global data for ", species)
#   }
#   
#   if(nrow(TMP.AUS) > 0){
#     TMP.AUS$RANGE = "AUS"
#   } else {
#     message("No Australian data for ", species, " don't plot the range")
#   }
#   
#   TMP.RANGE <- rbind(TMP.GLO, TMP.AUS)
#   names(TMP.RANGE)[2] = c("Temperature_range")
#   
#   ## Subset DF into records for each species
#   DF     <- subset(COMBO.SUA.POA, searchTaxon == species)
#   DF.OCC <- subset(COMBO.SUA.POA, searchTaxon == species & SOURCE != "INVENTORY")
#   DF.INV <- subset(COMBO.SUA.POA, searchTaxon == species & SOURCE == "INVENTORY")
#   
#   #############################################################
#   ## Plot occurrence points by source for the world
#   message('Writing global occ sources for ', species)
#   png(sprintf("./data/ANALYSIS/SPECIES_RANGES/%s_%s", species, "1km_occ_points_source.png"),
#       16, 10, units = 'in', res = 500)
#   
#   plot(LAND.WGS, main = paste0("Global points for ", species),
#        lwd = 0.01, asp = 1, col = 'grey', bg = 'sky blue')
#   
#   points(SP.DF,
#          pch = ".", cex = 3.3, cex.lab = 3, cex.main = 4, cex.axis = 2,
#          xlab = "", ylab = "", asp = 1,
#          col = factor(SP.DF$SOURCE))
#   
#   dev.off()
#   
#   #############################################################
#   ## Plot temperature barchart
#   message('Writing global temp histograms for ', species)
#   png(sprintf("./data/ANALYSIS/SPECIES_RANGES/%s_%s", species, "temp_barchart_1km_records.png"),
#       16, 10, units = 'in', res = 500)
#   
#   ## Use the 'SOURCE' column to create a histogram for each source.
#   ## Back to here..............................................................................
#   max.temp = max(TMP.RANGE$Annual_mean_temp_q95)+5
#   min.temp = min(TMP.RANGE$Annual_mean_temp_q05)-5
#   
#   temp.bar = 
#     ggplot(TMP.RANGE, aes(y = Temperature_range, x = RANGE, fill = RANGE)) +
#     
#     # scale_y_discrete(limits = c(min.temp,
#     #                             max.temp)) +
#     
#     geom_bar(stat = "identity", position = "identity", width = 0.1) +
#     coord_flip() +
#     
#     ## Add some median lines : overall, ALA and GBIF
#     ## This will only work if we plot the full range of temperatures on the x-axis
#     geom_vline(aes(xintercept = POA.SYD$Annual_mean_temp),
#                col = 'blue', size = 1) +
#     geom_vline(aes(xintercept = POA.SYD$Annual_mean_temp_50),
#                col = 'red', size = 1) +
#     geom_vline(aes(xintercept = POA.SYD$Annual_mean_temp_70),
#                col = 'green', size = 1) +
#     
#     
#     ggtitle(paste0("Worldclim temperature ranges for ", species)) +
#     
#     ## Add themes
#     theme(axis.title.x     = element_text(colour = "black", size = 35),
#           axis.text.x      = element_text(size = 25),
#           
#           axis.title.y     = element_text(colour = "black", size = 35),
#           axis.text.y      = element_text(size = 25),
#           
#           panel.background = element_blank(),
#           panel.border     = element_rect(colour = "black", fill = NA, size = 3),
#           plot.title       = element_text(size   = 40, face = "bold"),
#           legend.text      = element_text(size   = 20),
#           legend.title     = element_text(size   = 20),
#           legend.key.size  = unit(1.5, "cm"))
#   
#   ## Print the plot and close the device 
#   print(temp.hist + ggtitle(paste0("Worldclim temp niches for ", species)))
#   dev.off()




#############################################################
## Plot temperature histograms
# message('Writing global temp histograms for ', species)
# png(sprintf("./data/ANALYSIS/SPECIES_RANGES/%s_%s", species, "temp_niche_histograms_1km_records.png"),
#     16, 10, units = 'in', res = 500)
# 
# ## Use the 'SOURCE' column to create a histogram for each source.
# temp.hist = ggplot(DF, aes(x = Annual_mean_temp, group = SOURCE, fill = SOURCE)) +
#   
#   geom_histogram(position = "identity", alpha = 0.5, binwidth = 0.25,
#                  aes(y =..density..))  + 
#   geom_density(col = 4, alpha = 0.5) +
#   
#   ## Add some median lines : overall, ALA and GBIF
#   geom_vline(aes(xintercept = median(DF.INV$Annual_mean_temp)),
#              col = 'blue', size = 1) +
#   geom_vline(aes(xintercept = median(DF.OCC$Annual_mean_temp)),
#              col = 'red', size = 1) +
#   
#   ggtitle(paste0("Worldclim temp niches for ", species)) +
#   
#   ## Add themes
#   theme(axis.title.x     = element_text(colour = "black", size = 35),
#         axis.text.x      = element_text(size = 25),
#         
#         axis.title.y     = element_text(colour = "black", size = 35),
#         axis.text.y      = element_text(size = 25),
#         
#         panel.background = element_blank(),
#         panel.border     = element_rect(colour = "black", fill = NA, size = 3),
#         plot.title       = element_text(size   = 40, face = "bold"),
#         legend.text      = element_text(size   = 20),
#         legend.title     = element_text(size   = 20),
#         legend.key.size  = unit(1.5, "cm"))
# 
# ## Print the plot and close the device 
# print(temp.hist + ggtitle(paste0("Worldclim temp niches for ", species)))
# dev.off()
# 
# #############################################################
# ## Plot rainfall histograms
# message('Writing global rain histograms for ', species)
# png(sprintf("./data/ANALYSIS/SPECIES_RANGES/%s_%s", species, "rain_niche_histograms_1km_records.png"),
#     16, 10, units = 'in', res = 500)
# 
# ## Use the 'SOURCE' column to create a histogram for each source.
# rain.hist = ggplot(DF, aes(x = Annual_precip, group = SOURCE, fill = SOURCE)) +
#   
#   geom_histogram(position = "identity", alpha = 0.5, binwidth = 15,
#                  aes(y =..density..))  + 
#   geom_density(col = 4, alpha = 0.5) +
#   
#   ## Add some median lines : overall, ALA and GBIF
#   geom_vline(aes(xintercept = median(DF.INV$Annual_precip)),
#              col = 'blue', size = 1) +
#   geom_vline(aes(xintercept = median(DF.OCC$Annual_precip)),
#              col = 'red', size = 1) +
#   
#   ggtitle(paste0("Worldclim rain niches for ", species)) +
#   
#   ## Add themes
#   theme(axis.title.x     = element_text(colour = "black", size = 35),
#         axis.text.x      = element_text(size = 25),
#         
#         axis.title.y     = element_text(colour = "black", size = 35),
#         axis.text.y      = element_text(size = 25),
#         
#         panel.background = element_blank(),
#         panel.border     = element_rect(colour = "black", fill = NA, size = 3),
#         plot.title       = element_text(size   = 40, face = "bold"),
#         legend.text      = element_text(size   = 20),
#         legend.title     = element_text(size   = 20),
#         legend.key.size  = unit(1.5, "cm"))
# 
# ## Print the plot and close the device 
# print(rain.hist)
# dev.off()

#}