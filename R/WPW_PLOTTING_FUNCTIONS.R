## create simple maps of all GBIF records for selected taxa, in Australia and overseas
density_GBIF_records = function (DF, taxa.list, env.1, env.col.1, env.units.1,
                                   env.2, env.col.2, env.units.2) {
  
  ###############################
  ## for all the taxa in the list
  for (taxa.n in taxa.list) {
    
    #####################################################
    ## Subset histograms for global occurences of taxa.n
    env.1.all = DF[ which(DF$searchTaxon == taxa.n), ][[env.var.1]]
    min.1     = min(env.1, na.rm = TRUE)
    max.1     = max(env.1, na.rm = TRUE)
    
    env.1.ala = DF[ which(DF$searchTaxon == taxa.n 
                          & DF$SOURCE != "INVENTORY"), ][[env.var.1]]
    env.1.inv = DF[ which(DF$searchTaxon == taxa.n 
                          & DF$SOURCE == "INVENTORY"), ][[env.var.1]]
    
    #####################################################
    ## Subset histograms for global occurences of taxa.n
    env.2.all = DF[ which(DF$searchTaxon == taxa.n), ][[env.var.2]]
    min.2     = min(env.1, na.rm = TRUE)
    max.2     = max(env.1, na.rm = TRUE)
    
    env.2.ala = DF[ which(DF$searchTaxon == taxa.n 
                          & DF$SOURCE != "INVENTORY"), ][[env.var.2]]
    env.2.inv = DF[ which(DF$searchTaxon == taxa.n 
                          & DF$SOURCE == "INVENTORY"), ][[env.var.2]]
    
    #################################
    ## 2). Save env.1 histogram to file
    message("Saving histogram for ", taxa.n)
    CairoPNG(width  = 16180, height = 12000,
             file   = paste("./data/ANALYSIS/CLEAN_GBIF/", 
                            taxa.n, "_", env.var.1, "_all_records_histo.png", sep = ""),
             canvas = "white", bg = "white", units = "px", dpi = 600)
    
    ## Create layout
    nf <- layout(mat = matrix(c(1,6), 6,1, byrow = TRUE), height = c(1,6))
    par(mar = c(7.5, 5.5, 4, 2),
        oma = c(4, 4, 4, 4),
        mgp = c(6, 3, 0),
        font.lab = 2, lwd = 2)
    
    #################################
    ## print the boxplot
    boxplot(env.1.all, horizontal = TRUE,  outline = TRUE, 
            ylim  = c(min.1, max.1), 
            frame = FALSE, col = env.col.1, axes = TRUE)
    
    ## and print the
    hist(env.1.all, xlim = c(min.1, max.1),
         breaks = 50, border = NA, col = env.col.1, main = "",
         xlab = paste0(taxa.n, " ", env.var.1, " ", "(", env.units.1, ")", sep = ""), ylab = "",
         cex.lab = 5, cex.axis = 4)
    
    #################################
    ## print the boxplot
    boxplot(env.1.ala, horizontal = TRUE,  outline = TRUE, 
            ylim  = c(min.1, max.1), 
            frame = FALSE, col = env.col.1, axes = TRUE)
    
    ## and print the
    hist(env.1.ala, xlim = c(min.1, max.1),
         breaks = 50, border = NA, col = env.col.1, main = "",
         xlab = paste0(taxa.n, " ", env.var.1, " ", "(", env.units.1, ")", sep = ""), ylab = "",
         cex.lab = 5, cex.axis = 4)
    
    #################################
    ## print the boxplot
    boxplot(env.1.inv, horizontal = TRUE,  outline = TRUE, 
            ylim  = c(min.1, max.1), 
            frame = FALSE, col = env.col.1, axes = TRUE)
    
    ## and print the
    hist(env.1.inv, xlim = c(min.1, max.1),
         breaks = 20, border = NA, col = env.col.1, main = "",
         xlab = paste0(taxa.n, " ", env.var.1, " ", "(", env.units.1, ")", sep = ""), ylab = "",
         cex.lab = 5, cex.axis = 4)
    
    ## finsh the device
    dev.off()
    
    #################################
    ## 3). Save env.2 histogram to file
    CairoPNG(width  = 16180, height = 12000,
             file   = paste("./output/figures/Traits_v_glasshouse/",
                            taxa.n, "_", env.var.2, "_world_GBIF_histo.png", sep = ""),
             canvas = "white", bg = "white", units = "px", dpi = 600)
    
    ## create layout
    nf <- layout(mat = matrix(c(1,2), 2,1, byrow = TRUE), height = c(1,3))
    par(mar = c(7.5, 5.5, 4, 2),
        oma = c(4, 4, 4, 4),
        mgp = c(6, 3, 0),
        font.lab = 2, lwd = 2)
    
    ## print the boxplot
    boxplot(env.2, horizontal = TRUE,  outline = TRUE, 
            ylim  = c(min.2, max.2), 
            frame = FALSE, col = env.col.2, axes = TRUE)
    
    ## and print the
    hist(env.2,
         #xlim = c(min.2, max.2),
         breaks = 50, border = NA,
         col = env.col.2, main = "",
         xlab = paste0(taxa.n, " ", env.var.2, " ", "(", env.units.2, ")", sep = ""),
         ylab = "", cex.lab = 5, cex.axis = 4)
    
    ## finsh the device
    dev.off()
    
  }
  
}
                                   




## create simple maps of all GBIF records for selected taxa, in Australia and overseas
histogram_GBIF_records = function (DF, taxa.list, env.1, env.col.1, env.units.1,
                                   env.2, env.col.2, env.units.2) {
  
  ###############################
  ## for all the taxa in the list
  for (taxa.n in taxa.list) {
    
    #####################################################
    ## Subset histograms for global occurences of taxa.n
    env.1.all = DF[ which(DF$searchTaxon == taxa.n), ][[env.var.1]]
    min.1     = min(env.1, na.rm = TRUE)
    max.1     = max(env.1, na.rm = TRUE)
    
    env.1.ala = DF[ which(DF$searchTaxon == taxa.n 
                          & DF$SOURCE != "INVENTORY"), ][[env.var.1]]
    env.1.inv = DF[ which(DF$searchTaxon == taxa.n 
                          & DF$SOURCE == "INVENTORY"), ][[env.var.1]]
    
    #####################################################
    ## Subset histograms for global occurences of taxa.n
    env.2.all = DF[ which(DF$searchTaxon == taxa.n), ][[env.var.2]]
    min.2     = min(env.1, na.rm = TRUE)
    max.2     = max(env.1, na.rm = TRUE)
    
    env.2.ala = DF[ which(DF$searchTaxon == taxa.n 
                          & DF$SOURCE != "INVENTORY"), ][[env.var.2]]
    env.2.inv = DF[ which(DF$searchTaxon == taxa.n 
                          & DF$SOURCE == "INVENTORY"), ][[env.var.2]]
    
    #################################
    ## 2). Save env.1 histogram to file
    message("Saving histogram for ", taxa.n)
    CairoPNG(width  = 16180, height = 12000,
             file   = paste("./data/ANALYSIS/CLEAN_GBIF/", 
                            taxa.n, "_", env.var.1, "_all_records_histo.png", sep = ""),
             canvas = "white", bg = "white", units = "px", dpi = 600)
    
    ## Create layout
    nf <- layout(mat = matrix(c(1,6), 6,1, byrow = TRUE), height = c(1,6))
    par(mar = c(7.5, 5.5, 4, 2),
        oma = c(4, 4, 4, 4),
        mgp = c(6, 3, 0),
        font.lab = 2, lwd = 2)
    
    #################################
    ## print the boxplot
    boxplot(env.1.all, horizontal = TRUE,  outline = TRUE, 
            ylim  = c(min.1, max.1), 
            frame = FALSE, col = env.col.1, axes = TRUE)
    
    ## and print the
    hist(env.1.all, xlim = c(min.1, max.1),
         breaks = 50, border = NA, col = env.col.1, main = "",
         xlab = paste0(taxa.n, " ", env.var.1, " ", "(", env.units.1, ")", sep = ""), ylab = "",
         cex.lab = 5, cex.axis = 4)
    
    #################################
    ## print the boxplot
    boxplot(env.1.ala, horizontal = TRUE,  outline = TRUE, 
            ylim  = c(min.1, max.1), 
            frame = FALSE, col = env.col.1, axes = TRUE)
    
    ## and print the
    hist(env.1.ala, xlim = c(min.1, max.1),
         breaks = 50, border = NA, col = env.col.1, main = "",
         xlab = paste0(taxa.n, " ", env.var.1, " ", "(", env.units.1, ")", sep = ""), ylab = "",
         cex.lab = 5, cex.axis = 4)
    
    #################################
    ## print the boxplot
    boxplot(env.1.inv, horizontal = TRUE,  outline = TRUE, 
            ylim  = c(min.1, max.1), 
            frame = FALSE, col = env.col.1, axes = TRUE)
    
    ## and print the
    hist(env.1.inv, xlim = c(min.1, max.1),
         breaks = 20, border = NA, col = env.col.1, main = "",
         xlab = paste0(taxa.n, " ", env.var.1, " ", "(", env.units.1, ")", sep = ""), ylab = "",
         cex.lab = 5, cex.axis = 4)
    
    ## finsh the device
    dev.off()
    
    #################################
    ## 3). Save env.2 histogram to file
    CairoPNG(width  = 16180, height = 12000,
             file   = paste("./output/figures/Traits_v_glasshouse/",
                            taxa.n, "_", env.var.2, "_world_GBIF_histo.png", sep = ""),
             canvas = "white", bg = "white", units = "px", dpi = 600)
    
    ## create layout
    nf <- layout(mat = matrix(c(1,2), 2,1, byrow = TRUE), height = c(1,3))
    par(mar = c(7.5, 5.5, 4, 2),
        oma = c(4, 4, 4, 4),
        mgp = c(6, 3, 0),
        font.lab = 2, lwd = 2)
    
    ## print the boxplot
    boxplot(env.2, horizontal = TRUE,  outline = TRUE, 
            ylim  = c(min.2, max.2), 
            frame = FALSE, col = env.col.2, axes = TRUE)
    
    ## and print the
    hist(env.2,
         #xlim = c(min.2, max.2),
         breaks = 50, border = NA,
         col = env.col.2, main = "",
         xlab = paste0(taxa.n, " ", env.var.2, " ", "(", env.units.2, ")", sep = ""),
         ylab = "", cex.lab = 5, cex.axis = 4)
    
    ## finsh the device
    dev.off()
    
  }
  
}


