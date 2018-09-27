#########################################################################################################################
################################################## CREATE MESS MAPS ##################################################### 
#########################################################################################################################


#########################################################################################################################
## This code creates MESS maps for problematic species


#########################################################################################################################
## 1). CREATE THRESHOLD RASTERS
#########################################################################################################################


## Read in table of species checks 
## load("SDM_COMBINE.RData")
SDM.SPAT.ALL = readRDS(paste0('data/base/HIA_LIST/COMBO/SDM_SPAT_ALL_', save_run, '.rds'))
length(intersect(unique(SDM.SPAT.ALL$searchTaxon), GBIF.spp))

  
MAXENT.CHECK = read.csv("./output/maxent/MAXENT_SUA_CHECK_OLD_ALA_300_SPAT_1209_2016.csv", stringsAsFactors = FALSE)
MESS.SPP     = subset(MAXENT.CHECK, Action == "MESS")[["searchTaxon"]]
THRESH.SPP   = subset(MAXENT.CHECK, check.map == 1)[["searchTaxon"]]
map.spp      = gsub("_", " ", map_spp)
mess_spp     = intersect(MESS.SPP, map.spp)
mess_spp     = gsub(" ", "_", mess_spp)


## Create a corresponding list of thresholds for these species
MESS.THRESH  =  as.numeric(MAXENT.CHECK[MAXENT.CHECK$searchTaxon %in% MESS.SPP , ]
                           [["X10.percentile.training.presence.Logistic.threshold"]]) 


#########################################################################################################################
## Create a list of the current rasters. Just run this on the final SUA table. 
threshold.list = list.files(maxent_path, pattern = '_current_suit_above_', full.names = TRUE, recursive = TRUE)
current.list   = list.files(maxent_path, pattern = '_current.tif$',        full.names = TRUE, recursive = TRUE)
length(threshold.list);length(current.list)


#########################################################################################################################
## Can you loop like this?
for (i in 1:length(map.spp)) {
  
  for (i in 1:length(threshold.list)) {
    
    for (i in 1:length(current.list)) {
      
      ## Create points for each species
      spp.points <- SDM.SPAT.ALL[SDM.SPAT.ALL$searchTaxon == map.spp[i], ] %>%
        spTransform(CRS('+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'))
      
      ## Create threshold raster
      spp.thresh = raster(threshold.list[i])

      ## Create current raster
      spp.current = raster(current.list[i])
      
      ## Now create the empty panel just before plotting
      empty_ras <- init(spp.thresh, function(x) NA)
      save.spp = gsub(' ', '_', map.spp[i])
      
      png(sprintf('%s/%s/full/%s_%s.png', maxent_path, save.spp, save.spp, "current_suit"),      
          11, 4, units = 'in', res = 300)
      
      message('writing map for ', 'map.spp[i]')
      ## Need an empty frame
      print(levelplot(stack(empty_ras,
                            spp.current, 
                            spp.thresh, 
                            quick = TRUE), margin = FALSE,
                      
                      ## Create a colour scheme using colbrewer: 100 is to make it continuos
                      ## Also, make it a one-directional colour scheme
                      scales      = list(draw = FALSE), 
                      at = seq(0, 1, length = 100),
                      col.regions = colorRampPalette(rev(brewer.pal(11, 'Spectral'))),
                      
                      ## Give each plot a name: the third panel is the GCM
                      names.attr = c('Australian records', 'Current', 'Current threshold'),
                      colorkey   = list(height = 0.5, width = 3), xlab = '', ylab = '',
                      main       = list(gsub('_current_suit_above_', '_', names(spp.thresh)), font = 4, cex = 2)) +
              
              ## Plot the Aus shapefile with the occurrence points for reference
              ## Can the points be made more legible for both poorly and well recorded species?
              layer(sp.polygons(aus_albers), data = list(aus_albers = aus_albers)) +
              layer(sp.points(spp.points, pch = 19, cex = 0.15, 
                              col = c('red', 'transparent', 'transparent')[panel.number()]), data = list(spp.points = spp.points)))
      dev.off()
      
    }
    
  }
  
}





#########################################################################################################################
## 2). CREATE MESS MAPS FOR SELECTED TAXA
#########################################################################################################################


## Do we need mess maps for every scenario, or just for one?


## What steps do we want the MESS map function to produce for each species?

## Read in the files for each species

## Read in HS predictions

## Read in current climate - that is aus.grids.current

## Read in future climate  - that is the stack "s" in my mapping function , "x" for scenario 

## Create MESS maps for B.jarvisi - repeat for each scenarion, i.e. for x



#########################################################################################################################
## Run function to create MESS maps over a list of species and current suitability thresholds


## The double handling of names is causing the looping errors
## Do the subsetting of the rater frame outside the loop
grid_names           = sdm.predictors
current_grids        = aus.grids.current
names(current_grids) = grid_names  ## Error in `names<-`(`*tmp*`, value = grid_names) : incorrect number of layer names
current_grids        = subset(current_grids, intersect(names(current_grids), sdm.select))
names(current_grids)


species        = mess_spp[1]
threshold      = MESS.THRESH[1]
maxent_path    = maxent_path 
current_grids  = current_grids


##  layer(sp.polygons(aus_albers), data = list(aus_albers = aus_albers))
create_maxent_mess(species_list   = mess_spp,
                   threshold_list = MESS.THRESH,
                   maxent_path    = maxent_path, 
                   current_grids  = current_grids)



