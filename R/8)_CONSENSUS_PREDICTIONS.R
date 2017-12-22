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
load("./data/base/HIA_LIST/COMBO/COMBO_NICHE_CONTEXT.RData")
source('./R/HIA_LIST_MATCHING.R')


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
species_list  <- basename(list.dirs('./output/maxent/STD_VAR_ALL',   recursive = FALSE))
species_rev   = sort(species_list, decreasing = TRUE)


## Now make the test species directory names
test_spp = gsub(" ", "_", test.spp)
#test_spp = intersect(test_spp, species_list)
#test_rev = sort(test_spp, decreasing = TRUE)





#########################################################################################################################
## 2). CREATE LIST OF DIRECTORIES
#########################################################################################################################


#########################################################################################################################
## Just list the test species
SDM.RESULTS.DIR <- test_spp[c(1:length(test_spp))] %>%
  
  ## pipe the list into lapply
  lapply(function(species) {
    
    ## create the character string
    m <-   sprintf('./output/maxent/STD_VAR_ALL/%s/full/', species)
    m 
    
  }) %>%
  
  ## Bind the list together
  c()


## Check the output
str(SDM.RESULTS.DIR)





#########################################################################################################################
## 3). CREATE AN AVERAGE AND A CONSENSUS SUITABILITY RASTER FOR EACH SPECIES, FOR 2050 
#########################################################################################################################


#########################################################################################################################
## Iterate over each directory: 
ensemble.2050 = lapply(SDM.RESULTS.DIR, function(DIR) { 
  
  ## And each species - although we don't want all possible combinations
  lapply(test_spp, function(species) {
    
    ## First, as a workaround for the lapply combination, check if the file combination is correct
    f_current <- paste0(DIR, species, "_current.tif")
    
    ## Then only run the raster calculations if the file path exists
    if(file.exists(f_current)) {
      
      ###################################################################################################################
      ## Create a list of the rasters in each directory, then take the mean. How long does the mean calculation take?
      ## Tidy this up with %, etc
      raster.list = list.files(as.character(DIR), pattern = "bi50.tif", full.names = TRUE)
      suit        = stack(raster.list)
      suit.list   = unstack(suit)
      mean.suit   = mean(suit)   ## plot(mean)
      
      ## Write the mean to file
      f_mean = sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_2050_suitability_mean.tif', 
                       species, species)
      
      
      if(!file.exists(f_mean)) {
        
        writeRaster(mean.suit, sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_2050_suitability_mean.tif', 
                                       species, species), overwrite = TRUE)
        
      } else {
        
        message(species, '2050 mean suitability skipped - already exists')   ## 
        
      }
      
      ###################################################################################################################
      ## Then create rasters that meet habitat suitability criteria thresholds
      #suits = list()
      for (thresh in c(0.5, 0.7, 0.9)) {
        
        ## Check if the combined suitability raster exists
        f_suit <- sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_%s%s.tif',
                          species, species, "2050_suitability_consensus_greater_", thresh)
        
        ## If it exists, create the suitability rasters
        if(!file.exists(f_suit)) {
          
          ## First create a simple function to threshold each of the rasters in raster.list
          ## So how does this change to add them up? First
          scen_greater = function (x) {x > thresh}

          ## Then apply the function to the GCM list for each species
          ## The Init function initializes a raster object with values
          #suit_ras_greater    = reduce(suit.list, thresh_greater_fun, .init = suit.list[[1]] > thresh)
          #suit_ras_less       = reduce(suit.list, thresh_less_fun,    .init = suit.list[[1]] < thresh)
          suit_ras1_greater  = scen_greater(suit.list[[1]])   ## do this better...
          suit_ras2_greater  = scen_greater(suit.list[[2]])
          suit_ras3_greater  = scen_greater(suit.list[[3]])
          suit_ras4_greater  = scen_greater(suit.list[[4]])
          suit_ras5_greater  = scen_greater(suit.list[[5]])
          suit_ras6_greater  = scen_greater(suit.list[[6]])
          
          ## Then sum them up: could use a function or magrittr to compress this..
          suit_ras_consensus_sum   =  Reduce("+", list(suit_ras1_greater, suit_ras2_greater, suit_ras3_greater,
                                                       suit_ras4_greater, suit_ras5_greater, suit_ras6_greater))
          ## plot(suit_ras_consensus_sum )
          
          ## Write the raster for each species and threshold inside the loop. But how to access the rasters for plotting?
          message('Writing ', species, '2050 suitability > ', thresh) 
          writeRaster(suit_ras_consensus_sum, sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_%s%s.tif',
                                                      species, species, "2050_suitability_consensus_greater_", thresh), overwrite = TRUE)
          
        } else {
          
          message(species, '2050 suitability consensus > ', thresh, ' skipped - already exists')   ## 
          
        }
        
        # message('Writing ', species, ' suitability < ', thresh) 
        # writeRaster(suit_ras_less, sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_%s%s.tif',
        #                                    species, species, "suitability_less_thresh", thresh), overwrite = TRUE)
        
      }
      
      ########################################################################################################################
      ## Use the levelplot function to make a multipanel output: average, threshold 1, threshold 2
      # png(sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_suitability_combo.tif', 
      #             species, species),
      #     11, 4, units = 'in', res = 300)
      # 
      # ## Need an empty frame
      # print(levelplot(stack(pred.average,
      #                       pred.thresh.05,
      #                       pred.thresh.05), margin = FALSE,
      #                 
      #                 ## Create a colour scheme using colbrewer: 100 is to make it continuos
      #                 ## Also, make it a one-directional colour scheme
      #                 scales      = list(draw = FALSE), 
      #                 at = seq(0, 1, length = 100),
      #                 col.regions = colorRampPalette(rev(brewer.pal(11, 'Spectral'))),
      #                 
      #                 ## Give each plot a name
      #                 names.attr = c('GCM average', 'GCM > 0.7', 'GCM > 0.9'),
      #                 colorkey   = list(height = 0.5, width = 3), xlab = '', ylab = '',
      #                 main       = list(gsub('_', ' ', species), font = 4, cex = 2)))
      
      ## Plot the Aus shapefile with the occurrence points for reference
      ## Can the points be made more legible for both poorly and well recorded species?
      # layer(sp.polygons(aus)) +
      # layer(sp.points(occ, pch = 20, cex = 0.4, 
      #                 col = c('red', 'transparent', 'transparent')[panel.number()]), data = list(occ = occ)))
      
      ##
      #dev.off()
      
    } else {
      
      message(species, ' ', ' skipped - incorrect directory')   ## not needed with a proper loop
      
    }
    
  })
  
})





#########################################################################################################################
## 4). CREATE AN AVERAGE AND A CONSENSUS SUITABILITY RASTER FOR EACH SPECIES, FOR 2070 
#########################################################################################################################


#########################################################################################################################
## Iterate over each directory: 
ensemble.2070 = lapply(SDM.RESULTS.DIR, function(DIR) { 
  
  lapply(test_spp, function(species) {
    
    ## First, as a workaround for the lapply combination, check if the file combination is correct
    f_current <- paste0(DIR, species, "_current.tif")
    
    ## Then only run the raster calculations if the file path exists
    if(file.exists(f_current)) {
      
      ###################################################################################################################
      ## Create a list of the rasters in each directory, then take the mean. How long does the mean calculation take?
      ## Tidy this up with %, etc
      raster.list = list.files(as.character(DIR), pattern = "bi70.tif", full.names = TRUE)
      suit        = stack(raster.list)
      suit.list   = unstack(suit)
      mean.suit   = mean(suit)   ## plot(mean)
      
      ## Write the mean to file
      f_mean = sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_2070_suitability_mean.tif', 
                       species, species)
      
      
      if(!file.exists(f_mean)) {
        
        writeRaster(mean.suit, sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_2070_suitability_mean.tif', 
                                       species, species), overwrite = TRUE)
        
      } else {
        
        message(species, '2070 mean suitability skipped - already exists')   ## not needed with a proper loop
        
      }
      
      ###################################################################################################################
      ## Then create rasters that meet habitat suitability criteria thresholds
      #suits = list()
      for (thresh in c(0.5, 0.7, 0.9)) {
        
        ## Check if the combined suitability raster exists
        f_suit <- sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_%s%s.tif',
                          species, species, "2070_suitability_greater_thresh", thresh)
        
        ## If it exists, create the suitability rasters
        if(!file.exists(f_suit)) {
          
          ## First create a simple function to threshold each of the rasters in raster.list
          thresh_greater_fun  = function (x1, x2) {x1 & x2 > thresh}
          thresh_less_fun     = function (x1, x2) {x1 & x2 < thresh}
          
          ## Then apply the function to the GCM list for each species
          ## The Init function initializes a raster object with values
          #Error in .local(x, ...) : not a valid subset is due to the fact that you are providing a logical index vector
          suit_ras_greater    = reduce(suit.list, thresh_greater_fun, .init = suit.list[[1]] > thresh)
          #suit_ras_less       = reduce(suit.list, thresh_less_fun,    .init = suit.list[[1]] < thresh)
          
          ## Write the raster for each species and threshold inside the loop. But how to access the rasters for plotting?
          message('Writing ', species, '2070 suitability > ', thresh) 
          writeRaster(suit_ras_greater, sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_%s%s.tif',
                                                species, species, "2070_suitability_greater_thresh", thresh), overwrite = TRUE)
          
        } else {
          
          message(species, '2070 suitability > ', thresh, ' skipped - already exists')   ## 
          
        }
        
        # message('Writing ', species, ' suitability < ', thresh) 
        # writeRaster(suit_ras_less, sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_%s%s.tif',
        #                                    species, species, "suitability_less_thresh", thresh), overwrite = TRUE)
        
      }
      
      ########################################################################################################################
      ## Use the levelplot function to make a multipanel output: average, threshold 1, threshold 2
      # png(sprintf('./output/maxent/STD_VAR_ALL/%s/full/%s_suitability_combo.tif', 
      #             species, species),
      #     11, 4, units = 'in', res = 300)
      # 
      # ## Need an empty frame
      # print(levelplot(stack(pred.average,
      #                       pred.thresh.05,
      #                       pred.thresh.05), margin = FALSE,
      #                 
      #                 ## Create a colour scheme using colbrewer: 100 is to make it continuos
      #                 ## Also, make it a one-directional colour scheme
      #                 scales      = list(draw = FALSE), 
      #                 at = seq(0, 1, length = 100),
      #                 col.regions = colorRampPalette(rev(brewer.pal(11, 'Spectral'))),
      #                 
      #                 ## Give each plot a name
      #                 names.attr = c('GCM average', 'GCM > 0.7', 'GCM > 0.9'),
      #                 colorkey   = list(height = 0.5, width = 3), xlab = '', ylab = '',
      #                 main       = list(gsub('_', ' ', species), font = 4, cex = 2)))
      
      ## Plot the Aus shapefile with the occurrence points for reference
      ## Can the points be made more legible for both poorly and well recorded species?
      # layer(sp.polygons(aus)) +
      # layer(sp.points(occ, pch = 20, cex = 0.4, 
      #                 col = c('red', 'transparent', 'transparent')[panel.number()]), data = list(occ = occ)))
      
      ##
      #dev.off()
      
    } else {
      
      message(species, ' ', ' skipped - incorrect directory')   ## not needed with a proper loop
      
    }
    
  })
  
})





#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################