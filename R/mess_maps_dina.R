library(raster)
library(dismo)
library(gdalUtils)
library(RColorBrewer)
library(rasterVis)
library(latticeExtra)
library(rnaturalearth)
library(magrittr)
library(devtools)
source_gist('26e8091f082f2b3dd279', filename = 'polygonizer.R')
source_gist('c6a1cb61b8b6616143538950e6ec34aa', filename = 'hatch.R')
###

setwd("C:/data")
aus <- ne_states(country='Australia') %>% 
  subset(!grepl('Island', name)) %>% 
  spTransform(CRS('+init=epsg:3577'))


# First, make sure that all bioclim raster files for variables 1-9 have leading
# zeroes. This renames files on your computer... make sure that's what you ## 
# want to do.
ff <- list.files('C:/data', 'bio_\\d\\.', full=TRUE, recursive=TRUE) 
file.rename(ff, sub('_(\\d)\\.', '_0\\1.', ff))


## Read in the files for B. aquilonis
swd <- readRDS('C:/data/max_result_3rd_chapter/Bactrocera_jarvisi_04_07_11_12_15/swd.rds')
model <- readRDS('C:/data/max_result_3rd_chapter/Bactrocera_jarvisi_04_07_11_12_15/maxent_fitted.rds')
vars <- names(model$me_full@presence)


## Read in HS predictions
hs_current <- raster('C:/data/max_result_3rd_chapter/masked_continuous/Bactrocera_jarvisi_04_07_11_12_15_current_logistic.tif')
hs_current_binary <- raster('C:/data/max_result_3rd_chapter/masked_binary/Bactrocera_jarvisi_04_07_11_12_15_current_logistic_binary_10pctTraining_masked.tif')


## Read in climate data for current
s_current <- stack(sprintf('C:/data/worldclim/0.5/current/%s.bil', vars))


## Read in climate data for future
s_cccma_canesm2_rcp8_5_2030    <- stack(sprintf('C:/data/Future_climates/cccma_canesm2_rcp8_5_2030/%s.bil',   vars))
s_cccma_canesm2_rcp8_5_2050    <- stack(sprintf('C:/data/Future_climates/cccma_canesm2_rcp8_5_2050/%s.bil',   vars))
s_cccma_canesm2_rcp8_5_2070    <- stack(sprintf('C:/data/Future_climates/cccma_canesm2_rcp8_5_2070/%s.bil',   vars))
s_csiro_access1_0_rcp8_5_2030  <- stack(sprintf('C:/data/Future_climates/csiro_access1_0_rcp8_5_2030/%s.bil', vars))
s_csiro_access1_0_rcp8_5_2050  <- stack(sprintf('C:/data/Future_climates/csiro_access1_0_rcp8_5_2050/%s.bil', vars))
s_csiro_access1_0_rcp8_5_2070  <- stack(sprintf('C:/data/Future_climates/csiro_access1_0_rcp8_5_2070/%s.bil', vars))
s_gfdl_esm2m_rcp8_5_2030       <- stack(sprintf('C:/data/Future_climates/gfdl_esm2m_rcp8_5_2030/%s.bil', vars))
s_gfdl_esm2m_rcp8_5_2050       <- stack(sprintf('C:/data/Future_climates/gfdl_esm2m_rcp8_5_2050/%s.bil', vars))
s_gfdl_esm2m_rcp8_5_2070       <- stack(sprintf('C:/data/Future_climates/gfdl_esm2m_rcp8_5_2070/%s.bil', vars))
s_miroc_miroc5_rcp8_5_2030     <- stack(sprintf('C:/data/Future_climates/miroc_miroc5_rcp8_5_2030/%s.bil', vars))
s_miroc_miroc5_rcp8_5_2050     <- stack(sprintf('C:/data/Future_climates/miroc_miroc5_rcp8_5_2050/%s.bil', vars))
s_miroc_miroc5_rcp8_5_2070     <- stack(sprintf('C:/data/Future_climates/miroc_miroc5_rcp8_5_2070/%s.bil', vars))
s_mohc_hadgem2_cc_rcp8_5_2030  <- stack(sprintf('C:/data/Future_climates/mohc_hadgem2_cc_rcp8_5_2030/%s.bil', vars))
s_mohc_hadgem2_cc_rcp8_5_2050  <- stack(sprintf('C:/data/Future_climates/mohc_hadgem2_cc_rcp8_5_2050/%s.bil', vars))
s_mohc_hadgem2_cc_rcp8_5_2070  <- stack(sprintf('C:/data/Future_climates/mohc_hadgem2_cc_rcp8_5_2070/%s.bil', vars))
s_ncc_noresm1_m_rcp8_5_2030    <- stack(sprintf('C:/data/Future_climates/ncc_noresm1_m_rcp8_5_2030/%s.bil', vars))
s_ncc_noresm1_m_rcp8_5_2050    <- stack(sprintf('C:/data/Future_climates/ncc_noresm1_m_rcp8_5_2050/%s.bil', vars))
s_ncc_noresm1_m_rcp8_5_2070    <- stack(sprintf('C:/data/Future_climates/ncc_noresm1_m_rcp8_5_2070/%s.bil', vars))


# .. repeat for each scenario
# Create MESS maps for B.jarvisi - repeat for each scenario
mess_current <- mess(s_current, swd)            # ~7.3 minutes
novel_current <- mess_current < 0                ## "0"-not clear
novel_current[novel_current==0] <- NA            ## "NA"- not clear
novel_current <- mask(novel_current, hs_current) ## just to mask the grid to exclude PNG etc.
writeRaster(novel_current,'mess/bacjar_current_mess.tif', datatype='INT2S')

mess_cccma_canesm2_rcp8_5_2030 <- mess(s_cccma_canesm2_rcp8_5_2030, swd)
novel_cccma_canesm2_rcp8_5_2030 <- mess_cccma_canesm2_rcp8_5_2030 < 0
novel_cccma_canesm2_rcp8_5_2030[novel_cccma_canesm2_rcp8_5_2030==0] <- NA
novel_cccma_canesm2_rcp8_5_2030 <- mask(
  novel_cccma_canesm2_rcp8_5_2030, hs_current) # you can leave this as hs_current.. it's just to mask the grid to exclude PNG etc.
writeRaster(novel_cccma_canesm2_rcp8_5_2030,
            'mess/bacjar_cccma_canesm2_rcp8_5_2030_mess.tif', 
            datatype='INT2S')


### for 2050
mess_cccma_canesm2_rcp8_5_2050 <- mess(s_cccma_canesm2_rcp8_5_2050, swd)
novel_cccma_canesm2_rcp8_5_2050 <- mess_cccma_canesm2_rcp8_5_2050 < 0
novel_cccma_canesm2_rcp8_5_2050[novel_cccma_canesm2_rcp8_5_2050==0] <- NA
novel_cccma_canesm2_rcp8_5_2050 <- mask(
  novel_cccma_canesm2_rcp8_5_2050, hs_current) # you can leave this as hs_current.. it's just to mask the grid to exclude PNG etc.
writeRaster(novel_cccma_canesm2_rcp8_5_2050,
            'mess/bacjar_cccma_canesm2_rcp8_5_2050_mess.tif', 
            datatype='INT2S')

### for 2070
mess_cccma_canesm2_rcp8_5_2070 <- mess(s_cccma_canesm2_rcp8_5_2070, swd)
novel_cccma_canesm2_rcp8_5_2070 <- mess_cccma_canesm2_rcp8_5_2070 < 0
novel_cccma_canesm2_rcp8_5_2070[novel_cccma_canesm2_rcp8_5_2070==0] <- NA
novel_cccma_canesm2_rcp8_5_2070 <- mask(
  novel_cccma_canesm2_rcp8_5_2070, hs_current) # you can leave this as hs_current.. it's just to mask the grid to exclude PNG etc.
writeRaster(novel_cccma_canesm2_rcp8_5_2070,
            'mess/bacjar_cccma_canesm2_rcp8_5_2070_mess.tif', 
            datatype='INT2S')

#### for CSIRO scenarios:

mess_csiro_access1_0_rcp8_5_2030 <- mess(s_csiro_access1_0_rcp8_5_2030, swd)
novel_csiro_access1_0_rcp8_5_2030 <- mess_csiro_access1_0_rcp8_5_2030 < 0
novel_csiro_access1_0_rcp8_5_2030[novel_csiro_access1_0_rcp8_5_2030==0] <- NA
novel_csiro_access1_0_rcp8_5_2030 <- mask(
  novel_csiro_access1_0_rcp8_5_2030, hs_current) # you can leave this as hs_current.. it's just to mask the grid to exclude PNG etc.
writeRaster(novel_csiro_access1_0_rcp8_5_2030,
            'mess/bacjar_csiro_access1_0_rcp8_5_2030_mess.tif', 
            datatype='INT2S')

mess_csiro_access1_0_rcp8_5_2050 <- mess(s_csiro_access1_0_rcp8_5_2050, swd)
novel_csiro_access1_0_rcp8_5_2050 <- mess_csiro_access1_0_rcp8_5_2050 < 0
novel_csiro_access1_0_rcp8_5_2050[novel_csiro_access1_0_rcp8_5_2050==0] <- NA
novel_csiro_access1_0_rcp8_5_2050 <- mask(
  novel_csiro_access1_0_rcp8_5_2050, hs_current) # you can leave this as hs_current.. it's just to mask the grid to exclude PNG etc.
writeRaster(novel_csiro_access1_0_rcp8_5_2050,
            'mess/bacjar_csiro_access1_0_rcp8_5_2050_mess.tif', 
            datatype='INT2S')

mess_csiro_access1_0_rcp8_5_2070 <- mess(s_csiro_access1_0_rcp8_5_2070, swd)
novel_csiro_access1_0_rcp8_5_2070 <- mess_csiro_access1_0_rcp8_5_2070 < 0
novel_csiro_access1_0_rcp8_5_2070[novel_csiro_access1_0_rcp8_5_2070==0] <- NA
novel_csiro_access1_0_rcp8_5_2070 <- mask(
  novel_csiro_access1_0_rcp8_5_2070, hs_current) # you can leave this as hs_current.. it's just to mask the grid to exclude PNG etc.
writeRaster(novel_csiro_access1_0_rcp8_5_2070,
            'mess/bacjar_csiro_access1_0_rcp8_5_2070_mess.tif', 
            datatype='INT2S')

###For GFDL
mess_gfdl_esm2m_rcp8_5_2030 <- mess(s_gfdl_esm2m_rcp8_5_2030, swd)
novel_gfdl_esm2m_rcp8_5_2030 <- mess_gfdl_esm2m_rcp8_5_2030 < 0
novel_gfdl_esm2m_rcp8_5_2030[novel_gfdl_esm2m_rcp8_5_2030==0] <- NA
novel_gfdl_esm2m_rcp8_5_2030 <- mask(
  novel_gfdl_esm2m_rcp8_5_2030, hs_current) # you can leave this as hs_current.. it's just to mask the grid to exclude PNG etc.
writeRaster(novel_gfdl_esm2m_rcp8_5_2030,
            'mess/bacjar_gfdl_esm2m_rcp8_5_2030_mess.tif', 
            datatype='INT2S')

mess_gfdl_esm2m_rcp8_5_2050 <- mess(s_gfdl_esm2m_rcp8_5_2050, swd)
novel_gfdl_esm2m_rcp8_5_2050 <- mess_gfdl_esm2m_rcp8_5_2050 < 0
novel_gfdl_esm2m_rcp8_5_2050[novel_gfdl_esm2m_rcp8_5_2050==0] <- NA
novel_gfdl_esm2m_rcp8_5_2050 <- mask(
  novel_gfdl_esm2m_rcp8_5_2050, hs_current) # you can leave this as hs_current.. it's just to mask the grid to exclude PNG etc.
writeRaster(novel_gfdl_esm2m_rcp8_5_2050,
            'mess/bacjar_gfdl_esm2m_rcp8_5_2050_mess.tif', 
            datatype='INT2S')

mess_gfdl_esm2m_rcp8_5_2070 <- mess(s_gfdl_esm2m_rcp8_5_2070, swd)
novel_gfdl_esm2m_rcp8_5_2070 <- mess_gfdl_esm2m_rcp8_5_2070 < 0
novel_gfdl_esm2m_rcp8_5_2070[novel_gfdl_esm2m_rcp8_5_2070==0] <- NA
novel_gfdl_esm2m_rcp8_5_2070 <- mask(
  novel_gfdl_esm2m_rcp8_5_2070, hs_current) # you can leave this as hs_current.. it's just to mask the grid to exclude PNG etc.
writeRaster(novel_gfdl_esm2m_rcp8_5_2070,
            'mess/bacjar_gfdl_esm2m_rcp8_5_2070_mess.tif', 
            datatype='INT2S')

###for MIROC
mess_miroc_miroc5_rcp8_5_2030 <- mess(s_miroc_miroc5_rcp8_5_2030, swd)
novel_miroc_miroc5_rcp8_5_2030 <- mess_miroc_miroc5_rcp8_5_2030 < 0
novel_miroc_miroc5_rcp8_5_2030[novel_miroc_miroc5_rcp8_5_2030==0] <- NA
novel_miroc_miroc5_rcp8_5_2030 <- mask(
  novel_miroc_miroc5_rcp8_5_2030, hs_current) # you can leave this as hs_current.. it's just to mask the grid to exclude PNG etc.
writeRaster(novel_miroc_miroc5_rcp8_5_2030,
            'mess/bacjar_miroc_miroc5_rcp8_5_2030.tif', 
            datatype='INT2S')

mess_miroc_miroc5_rcp8_5_2050 <- mess(s_miroc_miroc5_rcp8_5_2050, swd)
novel_miroc_miroc5_rcp8_5_2050 <- mess_miroc_miroc5_rcp8_5_2050 < 0
novel_miroc_miroc5_rcp8_5_2050[novel_miroc_miroc5_rcp8_5_2050==0] <- NA
novel_miroc_miroc5_rcp8_5_2050 <- mask(
  novel_miroc_miroc5_rcp8_5_2050, hs_current) # you can leave this as hs_current.. it's just to mask the grid to exclude PNG etc.
writeRaster(novel_miroc_miroc5_rcp8_5_2050,
            'mess/bacjar_miroc_miroc5_rcp8_5_2050_mess.tif', 
            datatype='INT2S')

mess_miroc_miroc5_rcp8_5_2070 <- mess(s_miroc_miroc5_rcp8_5_2070, swd)
novel_miroc_miroc5_rcp8_5_2070 <- mess_miroc_miroc5_rcp8_5_2070 < 0
novel_miroc_miroc5_rcp8_5_2070[novel_miroc_miroc5_rcp8_5_2070==0] <- NA
novel_miroc_miroc5_rcp8_5_2070 <- mask(
  novel_miroc_miroc5_rcp8_5_2070, hs_current) # you can leave this as hs_current.. it's just to mask the grid to exclude PNG etc.
writeRaster(novel_miroc_miroc5_rcp8_5_2070,
            'mess/bacjar_miroc_miroc5_rcp8_5_2070_mess.tif', 
            datatype='INT2S')
###MOHC
mess_mohc_hadgem2_cc_rcp8_5_2030 <- mess(s_mohc_hadgem2_cc_rcp8_5_2030, swd)
novel_mohc_hadgem2_cc_rcp8_5_2030 <- mess_mohc_hadgem2_cc_rcp8_5_2030 < 0
novel_mohc_hadgem2_cc_rcp8_5_2030[novel_mohc_hadgem2_cc_rcp8_5_2030==0] <- NA
novel_mohc_hadgem2_cc_rcp8_5_2030 <- mask(
  novel_mohc_hadgem2_cc_rcp8_5_2030, hs_current) # you can leave this as hs_current.. it's just to mask the grid to exclude PNG etc.
writeRaster(novel_mohc_hadgem2_cc_rcp8_5_2030,
            'mess/bacjar_mohc_hadgem2_cc_rcp8_5_2030.tif', 
            datatype='INT2S')

mess_mohc_hadgem2_cc_rcp8_5_2050 <- mess(s_mohc_hadgem2_cc_rcp8_5_2050, swd)
novel_mohc_hadgem2_cc_rcp8_5_2050 <- mess_mohc_hadgem2_cc_rcp8_5_2050 < 0
novel_mohc_hadgem2_cc_rcp8_5_2050[novel_mohc_hadgem2_cc_rcp8_5_2050==0] <- NA
novel_mohc_hadgem2_cc_rcp8_5_2050 <- mask(
  novel_mohc_hadgem2_cc_rcp8_5_2050, hs_current) # you can leave this as hs_current.. it's just to mask the grid to exclude PNG etc.
writeRaster(novel_mohc_hadgem2_cc_rcp8_5_2050,
            'mess/bacjar_mohc_hadgem2_cc_rcp8_5_2050_mess.tif', 
            datatype='INT2S')

mess_mohc_hadgem2_cc_rcp8_5_2070 <- mess(s_mohc_hadgem2_cc_rcp8_5_2070, swd)
novel_mohc_hadgem2_cc_rcp8_5_2070 <- mess_mohc_hadgem2_cc_rcp8_5_2070 < 0
novel_mohc_hadgem2_cc_rcp8_5_2070[novel_mohc_hadgem2_cc_rcp8_5_2070==0] <- NA
novel_mohc_hadgem2_cc_rcp8_5_2070 <- mask(
  novel_mohc_hadgem2_cc_rcp8_5_2070, hs_current) # you can leave this as hs_current.. it's just to mask the grid to exclude PNG etc.
writeRaster(novel_mohc_hadgem2_cc_rcp8_5_2070,
            'mess/bacjar_mohc_hadgem2_cc_rcp8_5_2070_mess.tif', 
            datatype='INT2S')

###NCC
mess_ncc_noresm1_m_rcp8_5_2030 <- mess(s_ncc_noresm1_m_rcp8_5_2030, swd)
novel_ncc_noresm1_m_rcp8_5_2030 <- mess_ncc_noresm1_m_rcp8_5_2030 < 0
novel_ncc_noresm1_m_rcp8_5_2030[novel_ncc_noresm1_m_rcp8_5_2030==0] <- NA
novel_ncc_noresm1_m_rcp8_5_2030 <- mask(
  novel_ncc_noresm1_m_rcp8_5_2030, hs_current) # you can leave this as hs_current.. it's just to mask the grid to exclude PNG etc.
writeRaster(novel_ncc_noresm1_m_rcp8_5_2030,
            'mess/bacjar_ncc_noresm1_m_rcp8_5_2030.tif', 
            datatype='INT2S')

mess_ncc_noresm1_m_rcp8_5_2050 <- mess(s_ncc_noresm1_m_rcp8_5_2050, swd)
novel_ncc_noresm1_m_rcp8_5_2050 <- mess_ncc_noresm1_m_rcp8_5_2050 < 0
novel_ncc_noresm1_m_rcp8_5_2050[novel_ncc_noresm1_m_rcp8_5_2050==0] <- NA
novel_ncc_noresm1_m_rcp8_5_2050 <- mask(
  novel_ncc_noresm1_m_rcp8_5_2050, hs_current) # you can leave this as hs_current.. it's just to mask the grid to exclude PNG etc.
writeRaster(novel_ncc_noresm1_m_rcp8_5_2050,
            'mess/bacjar_ncc_noresm1_m_rcp8_5_2050_mess.tif', 
            datatype='INT2S')

mess_ncc_noresm1_m_rcp8_5_2070 <- mess(s_ncc_noresm1_m_rcp8_5_2070, swd)
novel_ncc_noresm1_m_rcp8_5_2070 <- mess_ncc_noresm1_m_rcp8_5_2070 < 0
novel_ncc_noresm1_m_rcp8_5_2070[novel_ncc_noresm1_m_rcp8_5_2070==0] <- NA
novel_ncc_noresm1_m_rcp8_5_2070 <- mask(
  novel_ncc_noresm1_m_rcp8_5_2070, hs_current) # you can leave this as hs_current.. it's just to mask the grid to exclude PNG etc.
writeRaster(novel_ncc_noresm1_m_rcp8_5_2070,
            'mess/bacjar_ncc_noresm1_m_rcp8_5_2070_mess.tif', 
            datatype='INT2S')


# Convert to polygon ... this may not work on your computer... it requires GDAL
# and Python to be set up in a specific way
poly <- polygonizer(novel_current)

# Show novel climates as hatched polygon
#plot(hs_current_masked, col=colorRampPalette(rev(brewer.pal(11, 'Spectral')))(100), zlim=c(0, 1))
#plot(poly, add=TRUE, density=40, border=NA)
poly_hatch <- hatch(poly, 60)

levelplot(hs_current, at=seq(0, 1, len=101), 
          col.regions=colorRampPalette(rev(brewer.pal(11, 'Spectral')))(100),
          margin=FALSE, scales=list(draw=FALSE), colorkey=list(height=0.5)) +
  layer(sp.polygons(aus)) +
  layer(sp.polygons(poly_hatch))


levelplot(hs_current_binary, at=seq(0, 1, len=101), 
          col.regions=c('white', brewer.pal(11, 'Spectral')[1]),
          margin=FALSE, scales=list(draw=FALSE), colorkey=FALSE) +
  layer(sp.polygons(aus, lwd=2)) +
  layer(sp.polygons(poly_hatch, lwd=1.5))

levelplot(hs_current_binary, at=seq(0, 1, len=101), 
          col.regions=c('white', brewer.pal(11, 'Spectral')[1]),
          margin=FALSE, scales=list(draw=FALSE), colorkey=FALSE) +
  layer(sp.polygons(aus, lwd=2)) +
  layer(sp.polygons(poly, fill='#00000050', col='transparent'), 
        data=list(poly=poly))
