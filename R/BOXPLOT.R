## create simple maps of all GBIF records for selected taxa, in Australia and overseas
Boxplot_GBIF_records = function (DF, taxa.list, env.1, env.col.1, env.units.1,
                                   env.2, env.col.2, env.units.2) {
  
  ###############################
  ## for all the taxa in the list
  for (taxa.n in taxa.list) {
    
    
    ## First, check if the file exists
    file  = paste("./output/figures/Traits_v_glasshouse/", 
                  taxa.n, "_", env.var.1, "_world_GBIF_histo.png", sep = "")
    
    ## If it's already downloaded, skip
    if (file.exists (file)) {

      print (paste ("Global histogram exists for", taxa.n, "skipping"))
      next

    }
    
    ## create vector for the environmental dimenson, and set min and max for plotting
    env.1 = DF[ which(DF$searchTaxon == taxa.n), ][[env.1]]
    min.1 = min(env.1, na.rm = TRUE)
    max.1 = max(env.1, na.rm = TRUE)
    
    ## create vector for the environmental dimenson, and set min and max for plotting
    env.2 = DF[ which(DF$searchTaxon == taxa.n), ][[env.2]]
    min.2 = min(env.2, na.rm = TRUE)
    max.2 = max(env.2, na.rm = TRUE)

    #####################################################
    ## 1). Plot histograms for global occurences of taxa.n
    
    #################################
    ## 2). Save env.1 histogram to file
    message("Saving boxplot for ", taxa.n)
    CairoPNG(width  = 16180, height = 12000,
             file   = paste("./output/figures/Traits_v_glasshouse/", 
                            taxa.n, "_", env.var.1, "_world_GBIF_histo.png", sep = ""),
             canvas = "white", bg = "white", units = "px", dpi = 600)
    
    ## create layout
    nf <- layout(mat = matrix(c(1,2), 2,1, byrow = TRUE), height = c(1,1))
    par(mar = c(7.5, 0.5, 0.5, 0.5),
        oma = c(4, 4, 4, 4),
        mgp = c(7, 3, 0),
        font.lab = 2, lwd = 2)
    
    ## print the boxplot
    boxplot(env.1, horizontal = TRUE,  outline = TRUE, 
            ylim  = c(min.1, max.1), 
            frame = FALSE, col = env.col.1, axes = TRUE, lwd = 6, cex.lab = 4, cex.axis = 4,
            xlab = paste0(env.var.1, " ", "(", env.units.1, ")", sep = ""))

    ## print the boxplot
    boxplot(env.2, horizontal = TRUE,  outline = TRUE, 
            ylim  = c(min.2, max.2), 
            frame = FALSE, col = env.col.2, axes = TRUE, lwd = 6, cex.lab = 4, cex.axis = 4,
            xlab  = paste0(env.var.2, " ", "(", env.units.2, ")", sep = ""))

    dev.off()
    
  }
  
}