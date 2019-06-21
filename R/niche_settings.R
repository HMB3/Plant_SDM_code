#########################################################################################################################
## SUA SPP
GBIF.spp      = native.good.models ## your list of species
GBIF.spp.rev  = sort(GBIF.spp, decreasing = TRUE)                ## the list reversed - only needed for a big list

save_run      = "SUA_ANALYSIS_NATIVE_GOOD"                       ## a variable to append the run name to the output files
map_spp_list  = gsub(" ", "_", GBIF.spp)                         ## species list with "_" for mapping
map_spp_rev   = sort(map_spp_list, decreasing = TRUE)            ## reversed, so we can run two at once

GBIF_path     = "./data/base/HIA_LIST/GBIF/OCC_SEARCH/"          ## The path where GBIF data is stored
ALA_path      = "./data/base/HIA_LIST/ALA/TREES_TEST/"           ## The path where ALA data is stored  place

maxent_path   = './output/maxent/SUA_TREES_ANALYSIS/'            ## The directory where files are saved               
maxent_dir    = 'output/maxent/SUA_TREES_ANALYSIS'               ## Another version of the path that John's coding needs to run a loop
save_data     = 'TRUE'                                           ## Arguments for saving the intermediary output - i.e. niches
read_data     = 'FALSE'                                          ## Leave these the same - saves data, but doesn't read back in
save_path     = 'data/base/HIA_LIST/COMBO'

#########################################################################################################################
## TRAIT SPP
GBIF.spp      = trait.spp ## your list of species
GBIF.spp.rev  = sort(GBIF.spp, decreasing = TRUE)                ## the list reversed - only needed for a big list

save_run      = "TRAIT_SPP"                                  ## a variable to append the run name to the output files
map_spp_list  = gsub(" ", "_", GBIF.spp)                         ## species list with "_" for mapping
map_spp_rev   = sort(map_spp_list, decreasing = TRUE)            ## reversed, so we can run two at once

GBIF_path     = "./data/base/HIA_LIST/GBIF/OCC_SEARCH/"          ## The path where GBIF data is stored
ALA_path      = "./data/base/HIA_LIST/ALA/TREES_TEST/"           ## The path where ALA data is stored  place

maxent_path   = './output/maxent/SUA_TREES_ANALYSIS/'            ## The directory where files are saved               
maxent_dir    = 'output/maxent/SUA_TREES_ANALYSIS'               ## Another version of the path that John's coding needs to run a loop
save_data     = 'TRUE'                                           ## Arguments for saving the intermediary output - i.e. niches
read_data     = 'FALSE'                                          ## Leave these the same - saves data, but doesn't read back in
save_path     = 'data/base/HIA_LIST/COMBO

#########################################################################################################################
## HIA SPP
GBIF.spp      = unique(c(TPL.HIA, TPL.CLEAN, ALL.INV.EV))[1:500] ## your list of species
GBIF.spp.rev  = sort(GBIF.spp, decreasing = TRUE)                ## the list reversed - only needed for a big list

save_run      = "EVERGREEN_500"                                  ## a variable to append the run name to the output files
map_spp_list  = gsub(" ", "_", GBIF.spp)                         ## species list with "_" for mapping
map_spp_rev   = sort(map_spp_list, decreasing = TRUE)            ## reversed, so we can run two at once

GBIF_path     = "./data/base/HIA_LIST/GBIF/OCC_SEARCH/"          ## The path where GBIF data is stored
ALA_path      = "./data/base/HIA_LIST/ALA/TREES_TEST/"           ## The path where ALA data is stored  place

maxent_path   = './output/maxent/SUA_TREES_ANALYSIS/'            ## The directory where files are saved               
maxent_dir    = 'output/maxent/SUA_TREES_ANALYSIS'               ## Another version of the path that John's coding needs to run a loop
save_data     = 'TRUE'                                           ## Arguments for saving the intermediary output - i.e. niches
read_data     = 'FALSE'                                          ## Leave these the same - saves data, but doesn't read back in
save_path     = 'data/base/HIA_LIST/COMBO'