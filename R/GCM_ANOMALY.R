#########################################################################################################################
############################################ CALCULATE ANOMALIES AND PLOT ############################################### 
#########################################################################################################################


#########################################################################################################################
## also calcualte the temp and rain anomalies for 2050 and 2070, and map them
#calculate.anomaly.2050(scen_2050)
#calculate.anomaly.2070(scen_2070)


## Now combine the anomalies across the GCMs
## BIO1.2050.anomaly  = list.files("./data/base/worldclim/aus/0.5/bio/anomalies/BIO1",  pattern = 'bi50.tif')
anomalies = list.files("./data/base/worldclim/aus/0.5/bio/anomalies/")


## Create anomaly rasters in the global environment
setwd("F:/green_cities_sdm/data/base/worldclim/aus/0.5/bio/anomalies/")
for(i in anomalies ) { assign(unlist(strsplit(i, "[.]"))[1], raster(i)) } 
setwd("F:/green_cities_sdm")


#########################################################################################################################
## Make a panel plot showing 2050 anomalies for MAT (BIO1)
## Create a PNG
png(sprintf('./data/base/worldclim/aus/0.5/bio/anomalies/BI01_GCM_2050_anomalies.png'),
    13, 8, units = 'in', res = 500)


## Create a raster stack
print(levelplot(stack(BIO1_anomaly_ac85bi50,
                      BIO1_anomaly_cn85bi50,
                      BIO1_anomaly_gf85bi50,
                      BIO1_anomaly_hg85bi50,
                      BIO1_anomaly_mc85bi50,
                      BIO1_anomaly_no85bi50), margin = FALSE,
                
                ## Create a colour scheme using colbrewer: 100 is to make it continuos
                ## Also, make it a one-directional colour scheme
                scales      = list(draw = FALSE), 
                at = seq(1, 3.5, length = 100),
                col.regions = colorRampPalette(rev(brewer.pal(11, 'Spectral'))),
                
                ## Give each plot a name
                names.attr = c(names(BIO1_anomaly_ac85bi50), 
                               names(BIO1_anomaly_cn85bi50), 
                               names(BIO1_anomaly_gf85bi50),
                               
                               names(BIO1_anomaly_hg85bi50),
                               names(BIO1_anomaly_mc85bi50),
                               names(BIO1_anomaly_no85bi50)),
                
                ## Check with 6 panels?
                colorkey   = list(height = 0.5, width = 3), xlab = '', ylab = '',
                main       = "BIO1 GCM 2050 - current anomalies") +
        
        ## Plot the Aus shapefile 
        layer(sp.polygons(aus)))

## Finish the device
dev.off()



#########################################################################################################################
## Make a panel plot showing 2050 anomalies for Precip (BIO12)
## Create a PNG
png(sprintf('./data/base/worldclim/aus/0.5/bio/anomalies/BI012_GCM_2050_anomalies.png'),
    13, 8, units = 'in', res = 500)


## Create a raster stack
print(levelplot(stack(BIO12_anomaly_ac85bi50,
                      BIO12_anomaly_cn85bi50,
                      BIO12_anomaly_gf85bi50,
                      BIO12_anomaly_hg85bi50,
                      BIO12_anomaly_mc85bi50,
                      BIO12_anomaly_no85bi50), margin = FALSE,
                
                ## Check the colour schemes..........................................................................
                
                ## Create a colour scheme using colbrewer: 100 is to make it continuos
                ## Also, make it a one-directional colour scheme
                scales      = list(draw = FALSE), 
                at = seq(-100, 800, length = 100),
                col.regions = colorRampPalette(rev(brewer.pal(11, 'Spectral'))),
                
                ## Give each plot a name
                names.attr = c(names(BIO12_anomaly_ac85bi50), 
                               names(BIO12_anomaly_cn85bi50), 
                               names(BIO12_anomaly_gf85bi50),
                               
                               names(BIO12_anomaly_hg85bi50),
                               names(BIO12_anomaly_mc85bi50),
                               names(BIO12_anomaly_no85bi50)),
                
                ## Check with 6 panels?
                colorkey   = list(height = 0.5, width = 3), xlab = '', ylab = '',
                main       = "BIO12 GCM 2050 - current anomalies") +
        
        ## Plot the Aus shapefile
        layer(sp.polygons(aus)))

## Finish the device
dev.off()





#########################################################################################################################
## Make a panel plot showing 2070 anomalies for MAT (BIO1)
## Create a PNG
png(sprintf('./data/base/worldclim/aus/0.5/bio/anomalies/BI01_GCM_2070_anomalies.png'),
    13, 8, units = 'in', res = 700)


## Create a raster stack
print(levelplot(stack(BIO1_anomaly_ac85bi70,
                      BIO1_anomaly_cn85bi70,
                      BIO1_anomaly_gf85bi70,
                      BIO1_anomaly_hg85bi70,
                      BIO1_anomaly_mc85bi70,
                      BIO1_anomaly_no85bi70), margin = FALSE,
                
                ## Create a colour scheme using colbrewer: 100 is to make it continuos
                ## Also, make it a one-directional colour scheme?
                scales      = list(draw = FALSE), 
                at = seq(1.5, 6, length = 100),
                col.regions = colorRampPalette(rev(brewer.pal(11, 'Spectral'))),
                
                ## Check the no85bi70 anomaly..........................................................................
                
                ## Give each plot a name
                names.attr = c(names(BIO1_anomaly_ac85bi70), 
                               names(BIO1_anomaly_cn85bi70), 
                               names(BIO1_anomaly_gf85bi70),
                               
                               names(BIO1_anomaly_hg85bi70),
                               names(BIO1_anomaly_mc85bi70),
                               names(BIO1_anomaly_no85bi70)),
                
                ## Check with 6 panels?
                colorkey   = list(height = 0.5, width = 3), xlab = '', ylab = '',
                main       = "BIO1 GCM 2070 - current anomalies") +
        
        ## Plot the Aus shapefile 
        layer(sp.polygons(aus)))

## Finish the device
dev.off()



#########################################################################################################################
## Make a panel plot showing 2070 anomalies for Precip (BIO12)
## Create a PNG
png(sprintf('./data/base/worldclim/aus/0.5/bio/anomalies/BI012_GCM_2070_anomalies.png'),
    13, 8, units = 'in', res = 700)


## Create a raster stack
print(levelplot(stack(BIO12_anomaly_ac85bi70,
                      BIO12_anomaly_cn85bi70,
                      BIO12_anomaly_gf85bi70,
                      BIO12_anomaly_hg85bi70,
                      BIO12_anomaly_mc85bi70,
                      BIO12_anomaly_no85bi70), margin = FALSE,
                
                ## Create a colour scheme using colbrewer: 100 is to make it continuos
                ## Also, make it a one-directional colour scheme
                scales      = list(draw = FALSE), 
                at = seq(-100, 800, length = 100),
                col.regions = colorRampPalette(rev(brewer.pal(11, 'Spectral'))),
                
                ## Give each plot a name
                names.attr = c(names(BIO12_anomaly_ac85bi70), 
                               names(BIO12_anomaly_cn85bi70), 
                               names(BIO12_anomaly_gf85bi70),
                               
                               names(BIO12_anomaly_hg85bi70),
                               names(BIO12_anomaly_mc85bi70),
                               names(BIO12_anomaly_no85bi70)),
                
                ## Check with 6 panels?
                colorkey   = list(height = 0.5, width = 3), xlab = '', ylab = '',
                main       = "BIO12 GCM 2070 - current anomalies") +
        
        ## Plot the Aus shapefile
        layer(sp.polygons(aus)))

## Finish the device
dev.off()





#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################