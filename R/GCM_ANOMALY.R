#########################################################################################################################
############################################ CALCULATE ANOMALIES AND PLOT ############################################### 
#########################################################################################################################


## This code calculate the anomalies between the current environment and the GCMs for mean annual temp and precipitation
## 


#########################################################################################################################
## 1). CALCULATE ANOMALIES
#########################################################################################################################


#########################################################################################################################
## Calcualte the temp and rain anomalies for 2030, 2050 and 2070, and map them
calculate.anomaly.2030 = calculate.anomaly(scen_list    = scen_2030,
                                           time_slice   = 30,
                                           climate_path = "./data/base/worldclim/aus/0.5/bio")

calculate.anomaly.2050 = calculate.anomaly(scen_list    = scen_2050,
                                           time_slice   = 50,
                                           climate_path = "./data/base/worldclim/aus/0.5/bio")

calculate.anomaly.2070 = calculate.anomaly(scen_list    = scen_2070,
                                           time_slice   = 70,
                                           climate_path = "./data/base/worldclim/aus/0.5/bio")


## Now combine the anomalies across the GCMs
anomaly.list = list.files("./data/base/worldclim/aus/0.5/bio/anomalies/", pattern = '.tif')


## Create anomaly rasters in the global environment. Clunky but it works
setwd("F:/green_cities_sdm/data/base/worldclim/aus/0.5/bio/anomalies/")
for(i in anomaly.list ) { assign(unlist(strsplit(i, "[.]"))[1], raster(i)) } 
setwd("H:/green_cities_sdm")


## Check the extent is the same :: the level plot function doesn't like different extents
extent(BIO1_anomaly_cc85bi30)
extent(BIO12_anomaly_ac85bi30)

extent(BIO1_anomaly_cc85bi50)
extent(BIO12_anomaly_ac85bi50)

extent(BIO1_anomaly_cc85bi70)
extent(BIO12_anomaly_ac85bi70)


## create an Australia shapefile
aus <- ne_states(country = 'Australia') %>% 
  subset(!grepl('Island', name))





#########################################################################################################################
## 2). CALCULATE STANDARD DEVIATION IN VALUES ACROSS ALL GCMS
#########################################################################################################################


#########################################################################################################################
## Create lists of rasters
Bio1.30 = stack(list.files("./data/base/worldclim/aus/0.5/bio/", 
                          pattern    = '301.tif$', 
                          recursive  = TRUE,
                          full.names = TRUE))

Bio1.50 = stack(list.files("./data/base/worldclim/aus/0.5/bio/", 
                          pattern    = '501.tif$', 
                          recursive  = TRUE,
                          full.names = TRUE))

Bio1.70 = stack(list.files("./data/base/worldclim/aus/0.5/bio/", 
                           pattern    = '501.tif$', 
                           recursive  = TRUE,
                           full.names = TRUE))


## Also just calcualte the standard deviation across all gcms for a few variables
Bio1.30 = stack(list.files("./data/base/worldclim/aus/0.5/bio/", 
                           pattern    = '301.tif$', 
                           recursive  = TRUE,
                           full.names = TRUE))

Bio1.50 = stack(list.files("./data/base/worldclim/aus/0.5/bio/", 
                           pattern    = '501.tif$', 
                           recursive  = TRUE,
                           full.names = TRUE))

Bio1.70 = stack(list.files("./data/base/worldclim/aus/0.5/bio/", 
                           pattern    = '501.tif$', 
                           recursive  = TRUE,
                           full.names = TRUE))


## Calcualte the standard deviation of a raster stack
bio1.30.sd = calc(Bio1.30, sd)
bio1.50.sd = calc(Bio1.50, sd)
bio1.70.sd = calc(Bio1.70, sd)





#########################################################################################################################
#########################################################################################################################
## Make a panel plot showing 2030 anomalies for MAT (BIO1)
## Create a PNG
png(sprintf('./data/base/worldclim/aus/0.5/bio/anomalies/BI01_GCM_2030_anomalies.png'),
    13, 8, units = 'in', res = 300)

## Create a raster stack
print(levelplot(stack(BIO1_anomaly_ac85bi30,
                      BIO1_anomaly_cc85bi30,
                      BIO1_anomaly_gf85bi30,
                      BIO1_anomaly_hg85bi30,
                      BIO1_anomaly_mc85bi30,
                      BIO1_anomaly_no85bi30, quick = TRUE), margin = FALSE,
                
                ## Create a colour scheme using colbrewer: 100 is to make it continuos
                ## Also, make it a one-directional colour scheme
                scales      = list(draw = FALSE), 
                at = seq(1, 3.5, length = 100),
                col.regions = colorRampPalette(rev(brewer.pal(11, 'Spectral'))),
                
                ## Give each plot a name
                names.attr = c(names(BIO1_anomaly_ac85bi30), 
                               names(BIO1_anomaly_cc85bi30), 
                               names(BIO1_anomaly_gf85bi30),
                               
                               names(BIO1_anomaly_hg85bi30),
                               names(BIO1_anomaly_mc85bi30),
                               names(BIO1_anomaly_no85bi30)),
                
                ## Check with 6 panels?
                colorkey   = list(height = 0.5, width = 3), xlab = '', ylab = '',
                main       = "BIO1 GCM 2030 - current anomalies") +
        
        ## Plot the Aus shapefile 
        layer(sp.polygons(aus)))

## Finish the device
dev.off()


#########################################################################################################################
## Make a panel plot showing 2030 anomalies for Precip (BIO12)
## Create a PNG
png(sprintf('./data/base/worldclim/aus/0.5/bio/anomalies/BI012_GCM_2030_anomalies.png'),
    13, 8, units = 'in', res = 300)


## Create a raster stack
print(levelplot(stack(BIO12_anomaly_ac85bi30,
                      BIO12_anomaly_cc85bi30,
                      BIO12_anomaly_gf85bi30,
                      BIO12_anomaly_hg85bi30,
                      BIO12_anomaly_mc85bi30,
                      BIO12_anomaly_no85bi30, quick = TRUE), margin = FALSE,
                
                ## Make the colour scale for temperature uni-directional..................................................
                
                ## Create a colour scheme using colbrewer: 100 is to make it continuos
                ## Also, make it a one-directional colour scheme
                scales      = list(draw = FALSE), 
                at = seq(-100, 800, length = 100),
                col.regions = colorRampPalette(rev(brewer.pal(11, 'Spectral'))),
                
                ## Give each plot a name
                names.attr = c(names(BIO12_anomaly_ac85bi30), 
                               names(BIO12_anomaly_cc85bi30), 
                               names(BIO12_anomaly_gf85bi30),
                               
                               names(BIO12_anomaly_hg85bi30),
                               names(BIO12_anomaly_mc85bi30),
                               names(BIO12_anomaly_no85bi30)),
                
                ## Check with 6 panels?
                colorkey   = list(height = 0.5, width = 3), xlab = '', ylab = '',
                main       = "BIO12 GCM 2030 - current anomalies") +
        
        ## Plot the Aus shapefile
        layer(sp.polygons(aus)))

## Finish the device
dev.off()





#########################################################################################################################
#########################################################################################################################
## Make a panel plot showing 2050 anomalies for MAT (BIO1)
## Create a PNG
png(sprintf('./data/base/worldclim/aus/0.5/bio/anomalies/BI01_GCM_2050_anomalies.png'),
    13, 8, units = 'in', res = 500)

## Create a raster stack
print(levelplot(stack(BIO1_anomaly_ac85bi50,
                      BIO1_anomaly_cc85bi50,
                      BIO1_anomaly_gf85bi50,
                      BIO1_anomaly_hg85bi50,
                      BIO1_anomaly_mc85bi50,
                      BIO1_anomaly_no85bi50, quick = TRUE), margin = FALSE,
                
                ## Create a colour scheme using colbrewer: 100 is to make it continuos
                ## Also, make it a one-directional colour scheme
                scales      = list(draw = FALSE), 
                at = seq(1, 3.5, length = 100),
                col.regions = colorRampPalette(rev(brewer.pal(11, 'Spectral'))),
                
                ## Give each plot a name
                names.attr = c(names(BIO1_anomaly_ac85bi50), 
                               names(BIO1_anomaly_cc85bi50), 
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
                      BIO12_anomaly_cc85bi50,
                      BIO12_anomaly_gf85bi50,
                      BIO12_anomaly_hg85bi50,
                      BIO12_anomaly_mc85bi50,
                      BIO12_anomaly_no85bi50), margin = FALSE,
                
                ## Make the colour scale for temperature uni-directional..................................................
                
                ## Create a colour scheme using colbrewer: 100 is to make it continuos
                ## Also, make it a one-directional colour scheme
                scales      = list(draw = FALSE), 
                at = seq(-100, 800, length = 100),
                col.regions = colorRampPalette(rev(brewer.pal(11, 'Spectral'))),
                
                ## Give each plot a name
                names.attr = c(names(BIO12_anomaly_ac85bi50), 
                               names(BIO12_anomaly_cc85bi50), 
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
                      BIO1_anomaly_cc85bi70,
                      BIO1_anomaly_gf85bi70,
                      BIO1_anomaly_hg85bi70,
                      BIO1_anomaly_mc85bi70,
                      BIO1_anomaly_no85bi70), margin = FALSE,
                
                ## Create a colour scheme using colbrewer: 100 is to make it continuos
                ## Also, make it a one-directional colour scheme?
                scales      = list(draw = FALSE), 
                at = seq(1.5, 6, length = 100),
                col.regions = colorRampPalette(rev(brewer.pal(11, 'Spectral'))),
                
                ## Check the no85bi70 anomaly...........................................................................
                
                ## Give each plot a name
                names.attr = c(names(BIO1_anomaly_ac85bi70), 
                               names(BIO1_anomaly_cc85bi70), 
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
                      BIO12_anomaly_cc85bi70,
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
                               names(BIO12_anomaly_cc85bi70), 
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