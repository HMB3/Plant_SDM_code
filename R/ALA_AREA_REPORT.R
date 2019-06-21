#########################################################################################################################
########################################  CREATE ALA SPECIES LISTS ###################################################### 
#########################################################################################################################


## This code takes the plant lists from the ALA area report and creates a unique list
SYD.ALA = read.csv("./data/base/HIA_LIST/ALA/Syd_60km_species_list.csv",  stringsAsFactors = FALSE)
BRI.ALA = read.csv("./data/base/HIA_LIST/ALA/Bri_30km_species_list.csv",  stringsAsFactors = FALSE)
MEL.ALA = read.csv("./data/base/HIA_LIST/ALA/Mel_30km_species_list.csv",  stringsAsFactors = FALSE)
HOB.ALA = read.csv("./data/base/HIA_LIST/ALA/Hob_15km_species_list.csv",  stringsAsFactors = FALSE)
CAN.ALA = read.csv("./data/base/HIA_LIST/ALA/Can_15km_species_list.csv",  stringsAsFactors = FALSE)
ADE.ALA = read.csv("./data/base/HIA_LIST/ALA/Ade_40km_species_list.csv",  stringsAsFactors = FALSE)
DAR.ALA = read.csv("./data/base/HIA_LIST/ALA/Dar_20km_species_list.csv",  stringsAsFactors = FALSE)
PER.ALA = read.csv("./data/base/HIA_LIST/ALA/Per_30km_species_list.csv",  stringsAsFactors = FALSE)


## What are the columns?
colnames(SYD.ALA)
identical(colnames(SYD.ALA), 
          colnames(PER.ALA))


## Add a column for the city
SYD.ALA$CITY = "SYD";BRI.ALA$CITY = "BRI";MEL.ALA$CITY = "MEL";ADE.ALA$CITY = "ADE";
HOB.ALA$CITY = "HOB";CAN.ALA$CITY = "CAN";DAR.ALA$CITY = "DAR";PER.ALA$CITY = "PER"


## Merge the data together
CITY.ALA = bind_rows(SYD.ALA, BRI.ALA, MEL.ALA, HOB.ALA, CAN.ALA, ADE.ALA, DAR.ALA, PER.ALA)
