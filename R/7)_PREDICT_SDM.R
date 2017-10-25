#########################################################################################################################
######################################### PREDICT MAXENT TO FUTURE CLIMATES ############################################# 
#########################################################################################################################




## This code takes a table of all species occurrences (rows) and environmental values "columns", and runs a maxent model 
## for each species. 



env.grids.current = stack(
  file.path('//sci-7910/F/data/worldclim/aus/0.5/bio/current',
            sprintf('bio_%02d.tif', 1:19)))

env.grids.future = stack(
  sprintf('//sci-7910/F/data/worldclim/aus/0.5/bio/2050/%s/%s%s.tif',
          scen, scen, 1:19))

for(i in 1:11) {
  message(i)
  env.grids.current[[i]] <- env.grids.current[[i]]/10
  env.grids.future[[i]] <- env.grids.future[[i]]/10  
}

names(env.grids.current) <- names(env.grids.future) <- c(
  'Annual_mean_temp', 'Mean_diurnal_range',
  'Isothermality', 'Temp_seasonality', 'Max_temp_warm_month',
  'Min_temp_cold_month', 'Temp_annual_range', 'Mean_temp_wet_qu',
  'Mean_temp_dry_qu', 'Mean_temp_warm_qu', 'Mean_temp_cold_qu', 'Annual_precip',
  'Precip_wet_month', 'Precip_dry_month', 'Precip_seasonality', 'Precip_wet_qu',
  'Precip_dry_qu', 'Precip_warm_qu', 'Precip_col_qu')

aus <- ne_states(country='Australia') %>% 
  subset(!grepl('Island', name))

species_list <- basename(list.dirs('F:/green_cities_sdm/output/maxent/baseline', recursive=FALSE))
scenario_list <- basename(list.dirs('//sci-7910/F/data/worldclim/aus/0.5/bio/2050', recursive=FALSE))
lapply(species_list, function(species) {
  message('Doing ', species)
  lapply(scenario_list, function(scen)) {
    message('  Doing ', scen)
    m <- readRDS(sprintf('F:/green_cities_sdm/output/maxent/baseline/%s/maxent_fitted.rds', species))
    occ <- readRDS(sprintf('F:/green_cities_sdm/output/maxent/baseline/%s/occ.rds', species)) %>% 
      spTransform(CRS('+init=epsg:4326'))
    
    pred.current <- rmaxent::project(m$me_full, env.grids.current[[colnames(m$me_full@presence)]])
    pred.future <- rmaxent::project(m$me_full, env.grids.future[[colnames(m$me_full@presence)]])
    writeRaster(pred.current, sprintf('F:/green_cities_sdm/output/maxent/baseline/%s/full/%s_current.tif', 
                                      species, species))
    writeRaster(pred.future, sprintf('F:/green_cities_sdm/output/maxent/baseline/%s/full/%s_%s.tif', 
                                     species, species, scen))
    empty <- init(pred.future$prediction_logistic, function(x) NA)
    png(sprintf('F:/green_cities_sdm/output/maxent/baseline/%s/full/%s.png', species, species), 
        11, 4, units='in', res=300)
    levelplot(stack(empty,
                    pred.current$prediction_logistic,
                    pred.future$prediction_logistic), margin=FALSE, 
              scales=list(draw=FALSE), at=seq(0, 1, length=100),
              col.regions=colorRampPalette(rev(brewer.pal(11, 'Spectral'))),
              names.attr=c('Occurrence', 'Current', 'CSIRO Access 1.0, 2050, RCP8.5'),
              colorkey=list(height=0.6, width=3), xlab='', ylab='',
              main=list(gsub('_', ' ', species), font=4, cex=2)) +
      layer(sp.polygons(aus)) +
      layer(sp.points(occ, pch=20, cex=0.5, 
                      col=c('red', 'transparent', 'transparent')[panel.number()]))
    dev.off()
  }
})

