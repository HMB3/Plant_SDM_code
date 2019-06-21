#########################################################################################################################
################################# SUMMARISE SPEICES GAINS BY SUA ######################################################## 
#########################################################################################################################



##
location     = map_spp[1:2]
new.location = "HOLLOWS"



file.copy(from      = location, 
          to        = new.location, 
          overwrite = FALSE, 
          recursive = TRUE, 
          copy.mode = TRUE)


## Full list of files
t = list.files(maxent_path, full.names = TRUE) 

Match.SN = GBIF.TRIM.TAXO  %>%
  mutate(Match.SN.ST = 
           str_detect(searchTaxon, t)) %>%  ## scientificName, New_binomial
  
  select(one_of(c("scientificName",
                  "searchTaxon",
                  "Taxonomic.status",
                  "New.Taxonomic.status",
                  "New.Genus",
                  "New.Species",
                  "country",
                  "Match.SN.ST")))

