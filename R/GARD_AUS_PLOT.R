#########################################################################################################################
#########################################################################################################################
## 

#########################################################################################################################
## Jacaranda_mimosifoliacurrent
png(sprintf('%s/%s/full/%s_%s.png', maxent_path, 'Jacaranda_mimosifolia', 'Jacaranda_mimosifolia', "current_garden_plot"),
    11, 4, units = 'in', res = 300)

print(levelplot(stack(Jacaranda_mimosifolia_current,
                      quick = TRUE), margin = FALSE,
                
                ## Create a colour scheme using colbrewer: 100 is to make it continuos
                ## Also, make it a one-directional colour scheme
                scales      = list(draw = FALSE),
                at = seq(0, 1, length = 100),
                col.regions = colorRampPalette(rev(brewer.pal(11, 'Spectral'))),

                ## Give each plot a name: the third panel is the GCM
                names.attr = c('Current climatic suitability'),#, paste0('Suitability in 20', time_slice)),#, 'Gain/loss raster'),
                colorkey   = list(height = 0.5, width = 3), xlab = '', ylab = '',
                main       = c('Jacaranda mimosifoliacurrent suitability')) +
        
        ## Plot the Aus shapefile with the occurrence points for reference
        ## Can the points be made more legible for both poorly and well recorded species?
        latticeExtra::layer(sp.polygons(aus), data = list(aus = aus)))

dev.off()


#########################################################################################################################
## Jacaranda_mimosifolia2030
png(sprintf('%s/%s/full/%s_%s.png', maxent_path, 'Jacaranda_mimosifolia', 'Jacaranda_mimosifolia', "2030_garden_plot"),
    11, 4, units = 'in', res = 300)

print(levelplot(stack(Jacaranda_mimosifolia_2030,
                      quick = TRUE), margin = FALSE,
                
                ## Create a colour scheme using colbrewer: 100 is to make it continuos
                ## Also, make it a one-directional colour scheme
                scales      = list(draw = FALSE),
                at = seq(0, 6, length = 6),
                col.regions = colorRampPalette(rev(brewer.pal(11, 'Spectral'))),
                
                
                ## Give each plot a name: the third panel is the GCM
                names.attr = c('Current climatic suitability'),#, paste0('Suitability in 20', time_slice)),#, 'Gain/loss raster'),
                colorkey   = list(height = 0.5, width = 3), xlab = '', ylab = '',
                main       = c('Jacaranda mimosifolia2030 suitability - 6GCMs')) +
        
        ## Plot the Aus shapefile with the occurrence points for reference
        ## Can the points be made more legible for both poorly and well recorded species?
        latticeExtra::layer(sp.polygons(aus), data = list(aus = aus)))


dev.off()



#########################################################################################################################
## Jacaranda_mimosifolia2070
png(sprintf('%s/%s/full/%s_%s.png', maxent_path, 'Jacaranda_mimosifolia', 'Jacaranda_mimosifolia', "2070_current_garden_plot"),
    11, 4, units = 'in', res = 300)

print(levelplot(stack(Jacaranda_mimosifolia_2070,
                      quick = TRUE), margin = FALSE,
                
                ## Create a colour scheme using colbrewer: 100 is to make it continuos
                ## Also, make it a one-directional colour scheme
                scales      = list(draw = FALSE),
                at = seq(0, 6, length = 6),
                col.regions = colorRampPalette(rev(brewer.pal(11, 'Spectral'))),
                
                
                ## Give each plot a name: the third panel is the GCM
                names.attr = c('Current climatic suitability'),#, paste0('Suitability in 20', time_slice)),#, 'Gain/loss raster'),
                colorkey   = list(height = 0.5, width = 3), xlab = '', ylab = '',
                main       = c('Jacaranda mimosifolia2070 suitability - 6GCMs')) +
        
        ## Plot the Aus shapefile with the occurrence points for reference
        ## Can the points be made more legible for both poorly and well recorded species?
        latticeExtra::layer(sp.polygons(aus), data = list(aus = aus)))

dev.off()






