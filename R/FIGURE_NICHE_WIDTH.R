#########################################################################################################################
############################################## NICHE WIDTH FIGURE ####################################################### 
#########################################################################################################################


#########################################################################################################################
## Read in data

 
NICHE.DATA.ALL       = readRDS("./data/base/HIA_LIST/COMBO/CLEAN_ONLY_HIA_SPP.rds")
aus                  = readRDS("F:/green_cities_sdm/data/base/CONTEXTUAL/aus_states.rds") %>%
  spTransform(ALB.CONICAL)
source('./R/HIA_LIST_MATCHING.R')


#########################################################################################################################
## 1). BOXPLOTS
#########################################################################################################################


## Subset data
Ficus.brachypoda = subset(COMBO.RASTER.CONTEXT, searchTaxon == "Ficus brachypoda")
dim(Ficus.brachypoda)


## Create spatial object
FICUS.POINTS   = SpatialPointsDataFrame(coords      = Ficus.brachypoda[c("lon", "lat")], 
                                        data        = Ficus.brachypoda[c("lon", "lat")],
                                        proj4string = ALB.CONICAL)

## Create PNG file
CairoPNG(width  = 10000, height = 10000,
         file   = "./output/figs/WPW_MILESTONE_1/NICHE_FIG/AUS_RECORDS_FICUS_BRACHYPODA.png", 
         canvas = "white", bg = "white", units = "px", dpi = 600)

## set margins
par(mfrow = c(2, 1),            ## c(nrows, ncols)
    mar   = c(5, 18, 10, 1),  ## b, l, t, r
    mgp   = c(13, 2.5, 0),
    oma   = c(1, 1, 1, 1))


## Plot just the Australian points
plot(aus, main = "Ficus brachypoda")
points(Ficus.brachypoda, col = "red", cex = .3, pch = 19)


## Finish the device
dev.off()


#####################################################################################################################
## now make boxplots of rain and temp niche for all records of Ficus brachypoda
CairoPNG(width  = 10000, height = 16180, ## 13000, 10000
         file   = "./output/figures/WPW_MILESTONE_1/NICHE_FIG/NICHE_BOXPLOT_FICUS_BRACHYPODA.png", 
         canvas = "white", bg = "white", units = "px", dpi = 600)

par(mfrow = c(2, 1),            ## c(nrows, ncols)
    mar   = c(5, 18, 10, 1),  ## b, l, t, r
    mgp   = c(13, 2.5, 0),
    oma   = c(1, 1, 1, 1))


#################################################################
## RAIN
## the Annual_precip boxplot
boxplot(Ficus.brachypoda$Annual_precip, ylab = "Annual precip (mm)",
        medcol = "black", staplelwd = 8, outpch = 8, outcex = 3, whisklwd = 8, 
        cex = 2.5, cex.axis = 5, cex.lab = 6, las = 1, boxlwd = 8,
        col = (c("blue")), cex.main = 8, col.main = "blue",
        main = paste("RE = ", format(quantile(Ficus.brachypoda$Annual_precip, 0.95) - 
                                       quantile(Ficus.brachypoda$Annual_precip, 0.05), digits = 2)))

## compute the required quantiles, as per the niche width calcs
r.qntl <- quantile(Ficus.brachypoda$Annual_precip, probs = c(0.05, 0.95))
#quantile(Ficus.brachypoda$Annual_precip, 0.95) - quantile(Ficus.brachypoda$Annual_precip, 0.05)

## add them as a rgu plot to the left hand side
rug(r.qntl, side = 2, col = "blue", lwd = 10)


box(lwd = 8)

# ## add the Annual_precipfall niche width as a legend
# legend("topright", bty = "n", cex = 7, pt.cex = 7,
#        text.col = "blue", 
#        legend = paste("RW = ", format(quantile(Ficus.brachypoda$Annual_precip, 0.95) - 
#                                          quantile(Ficus.brachypoda$Annual_precip, 0.05), digits = 2)))


#################################################################
## TEMP
## then Annual_mean_temp boxplot
boxplot(Ficus.brachypoda$Annual_mean_temp , ylab = "Annual mean temp (°C)",
        medcol = "black", staplelwd = 8, outpch = 8, outcex = 3, whisklwd = 8, 
        cex = 2.5, cex.axis = 5, cex.lab = 6, las = 1, boxlwd = 8,
        col = (c("orange")), cex.main = 8, col.main = "orange",
        main = paste("TE = ", format(quantile(Ficus.brachypoda$Annual_mean_temp, 0.95) - 
                                       quantile(Ficus.brachypoda$Annual_mean_temp, 0.05), digits = 2)))


## compute the required quantiles, as per the niche width calcs
t.qntl <- quantile(Ficus.brachypoda$Annual_mean_temp, probs = c(0.05, 0.95))
quantile(Ficus.brachypoda$Annual_mean_temp, 0.95) - quantile(Ficus.brachypoda$Annual_mean_temp, 0.05)


## add them as a rgu plot to the left hand side
rug(t.qntl, side = 2, col = "orange", lwd = 10)


box(lwd = 8)


## add the Annual_mean_temp niche width as a legend
# legend("topright", bty = "n", cex = 7, pt.cex = 7,
#        text.col = "orange", 
#        legend = paste("TW = ", format(quantile(Ficus.brachypoda$Annual_mean_temp, 0.95) - 
#                                          quantile(Ficus.brachypoda$Annual_mean_temp, 0.05), digits = 2)))


## finish the device
dev.off()





#########################################################################################################################
##############################  END OF CHAPTER 4 DATA PREP  ############################################################# 
#########################################################################################################################