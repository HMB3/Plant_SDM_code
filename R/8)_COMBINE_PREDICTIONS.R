#########################################################################################################################
################################################ COMBINE MAXENT PREDICTIONS ############################################# 
#########################################################################################################################


## The aim of this code is to combine multiple maxent predictions into a single suitability raster for each species
## The ensemble.raster package creates a weighted average, and a consensus layer.

## Lots of different definitions of a consensus. Could start by getting a layer which is has 1 where all cells meet a
## threshold, 
## e.g. above 0.5, above 0.7, above 0.9
## e.g. below 0.5, below 0.3, etc


#########################################################################################################################
## Load packages
source('./R/HIA_LIST_MATCHING.R')
load("./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT.RData")


## Require packages
p <- c('ff',    'things',         'raster',    'dismo',        'sp',           'latticeExtra', 'data.table', 
       'rgdal', 'rgeos',          'gdalUtils', 'rmaxent',      'readr',        'dplyr',        'tidyr',
       'readr', 'rnaturalearth',  'rasterVis', 'RColorBrewer', 'latticeExtra', 'parallel')
sapply(p, require, character.only = TRUE)





#########################################################################################################################
## 1). CREATE LISTS OF GCMs AND SPECIES FOR MODEL RUNS 
#########################################################################################################################


#########################################################################################################################
## create a list of GCM scenarios (below is from CSIRO): 

## Eight of the 40 CMIP5 models assessed in this project have been selected for use in provision of application-ready data. 
## This facilitates efficient exploration of climate projections for Australia.
## https://www.climatechangeinaustralia.gov.au/en/support-and-guidance/faqs/eight-climate-models-data/
# scen = c("ip85bi50", "mc85bi50", "mg85bi50", "mi85bi50", "mp85bi50", 
#          "mr85bi50", "no85bi50", "ac85bi50", "bc85bi50", "cc85bi50", 
#          "cn85bi50", "gf85bi50", "gs85bi50", "hd85bi50", "he85bi50",
#          "hg85bi50", "in85bi50")


## Create a lookup table of GCMs
h <- read_html('http://www.worldclim.org/cmip5_30s') # Also https://www.climatechangeinaustralia.gov.au/en/support-and-guidance/faqs/eight-climate-models-data/
gcms <- h %>% 
  html_node('table') %>% 
  html_table(header=TRUE) %>% 
  filter(rcp85 != '')

id.50 <- h %>% 
  html_nodes(xpath = "//a[text()='bi']/@href") %>% 
  as.character %>% 
  grep('85bi50', ., value=TRUE) %>% 
  basename %>% 
  sub('\\.zip.', '', .)

id.70 <- h %>% 
  html_nodes(xpath = "//a[text()='bi']/@href") %>% 
  as.character %>% 
  grep('85bi70', ., value=TRUE) %>% 
  basename %>% 
  sub('\\.zip.', '', .)

gcms.50 <- cbind(gcms, id.50)
gcms.50$GCM = sub(" \\(#\\)", "", gcms$GCM)  ## sub replaces first instance in a string, gsub = global

gcms.70 <- cbind(gcms, id.70)
gcms.70$GCM = sub(" \\(#\\)", "", gcms$GCM) 


## Now create the scenario lists across which to loop
gcms.50 ; gcms.70


## Just get the 8 models picked by CSIRO for Australia, for 2050 and 2070
## Also will be a better way to get at this...
scen_2050 = c("mc85bi50", "no85bi50", "ac85bi50", "cn85bi50", "gf85bi50", "hg85bi50")
scen_2070 = c("mc85bi70", "no85bi70", "ac85bi70", "cn85bi70", "gf85bi70", "hg85bi70")

#########################################################################################################################

## All species on the growers list, and also a test list which includes Renee's species
species_list  <- basename(list.dirs('F:/green_cities_sdm/output/maxent/STD_VAR_ALL',   recursive = FALSE))
species_rev   = sort(species_list, decreasing = TRUE)


spp.all  <- unique(COMBO.NICHE.CONTEXT$searchTaxon)
str(spp.all)                 ## 6782


## The trial species
test.spp = sort(unique(c(renee.full$Species, "Betula pendula", "Fraxinus excelsior", "Quercus robur", "Fagus sylvatica")))
test.spp 


## Combine with 45 from the main list
HIA.SAMPLE = head(COMBO.NICHE.CONTEXT, 53)[, c("searchTaxon")]
test.spp   = sort(unique(c(test.spp, HIA.SAMPLE)))


## Now make the test species directory names
test_spp = gsub(" ", "_", test.spp)
test_spp = intersect(test_spp, species_list)
test_rev = sort(test_spp, decreasing = TRUE)


## Also create a list of thresholds:



#########################################################################################################################
## 2). CREATE AN AVERAGE AND A CONSENSUS SUITABILITY RASTER FOR EACH SPECIES
#########################################################################################################################


#########################################################################################################################
## First, iterate over each species, then each scenario EG: x = scen_2050[1];species = test_rev[1]
env.grids.2050 = lapply(test_rev, function(species) { ## lapply(scen_2050, function(x)
  
  lapply(scen_2050, function(x) {
    
    ## Create a raster stack for each 2050 GCM - also an empty raster for the final plot
    # sprintf('F:/green_cities_sdm/output/maxent/STD_VAR_ALL/%s/full/%s_%s.tif', 
    #         species, species, x)
    
    ########################################################################################################################
    ## First read in each prediction
    m <-   sprintf('F:/green_cities_sdm/output/maxent/STD_VAR_ALL/%s/full/%s_%s.tif', 
                   species, species, x)
    
    ########################################################################################################################
    ## Take the average
    #print(m)
    # mean <- overlay(r2000, r2001, r2002, r2003, r2004, r2005, fun=meanIgnoringZeroes)
    
    
    ########################################################################################################################
    ## If the future raster doesn't exist, create it 
    
    
    ########################################################################################################################
    ## Use the levelplot function to make a multipanel output: average, threshold 1, threshold 2
    png(sprintf('F:/green_cities_sdm/output/maxent/STD_VAR_ALL/%s/full/%s_%s.png', species, species, x),      
        11, 4, units = 'in', res = 300)
    
    ## Need an empty frame
    print(levelplot(stack(pred.average,
                          pred.thresh.05,
                          pred.thresh.05), margin = FALSE,
                    
                    ## Create a colour scheme using colbrewer: 100 is to make it continuos
                    ## Also, make it a one-directional colour scheme
                    scales      = list(draw = FALSE), 
                    at = seq(0, 1, length = 100),
                    col.regions = colorRampPalette(rev(brewer.pal(11, 'Spectral'))),
                    
                    ## Give each plot a name
                    names.attr = c('Occurrence', 'Current', sprintf('%s, 2050, RCP8.5', scen_name)),
                    colorkey   = list(height = 0.5, width = 3), xlab = '', ylab = '',
                    main       = list(gsub('_', ' ', species), font = 4, cex = 2)))
            
            ## Plot the Aus shapefile with the occurrence points for reference
            ## Can the points be made more legible for both poorly and well recorded species?
            # layer(sp.polygons(aus)) +
            # layer(sp.points(occ, pch = 20, cex = 0.4, 
            #                 col = c('red', 'transparent', 'transparent')[panel.number()]), data = list(occ = occ)))
    
    ##
    dev.off()
    
  })
  
  
})





