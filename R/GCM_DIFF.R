#########################################################################################################################
################################################### OPTIONS FOR GCM CALC ################################################ 
#########################################################################################################################




## Then subtract the binary current layer from the binary future layer. So need to mask the no-data from this overlay calc.  
# binary_future_minus_current     = overlay(combo_suit_binary,
#                                           current_suit_thresh,
#                                           fun = function(r1, r2) {return (r1 - r2)})

#########################################################################################################################
## Subtract the binary current layer from the discrete future layer 
discrete_future_minus_current = overlay(combo_suit_band,
                                        current_suit_thresh,
                                        fun = function(r1, r2) {return (r1 - r2)})                                           

#########################################################################################################################
## Subtract the binary current layer from the integer future layer 
integer_future_minus_current  = overlay(combo_suit_thresh,
                                        current_suit_thresh,
                                        fun = function(r1, r2) {return (r1 - r2)})

## Plot the difference between future and current layers...
plot(current_suit_thresh, main = gsub('_', ' ', (sprintf('%s current Max_train_sensit > %s', species, thresh))))
plot(combo_suit_thresh,   main = gsub('_', ' ', (sprintf('%s future Max_train_sensit > %s',  species, thresh))))


## Plot the integer difference
plot(integer_future_minus_current,
     main = gsub('_', ' ', (sprintf('%s future - current  Max_train_sensit > %s', species, thresh))))

## Plot the band difference
plot(discrete_future_minus_current,
     main = gsub('_', ' ', (sprintf('%s future - discrete  Max_train_sensit > %s', species, thresh))))


#########################################################################################################################
## Now, mask out the zero values?
# plus    = overlay(combo_suit_thresh,
#                   current_suit_thresh,
#                   fun = function(r1, r2) {return (r1 + r2)})
# 
# mask <- calc(plus, fun = rc)
# mask[mask < 0] <- NA
# 
# pol <- rasterToPolygons(mask, fun = function(x){x>0})
# 
# 
# 
# test   =  polygonizer(plus, outshape = NULL, pypath = "F:/green_cities_sdm/R")
# crs(b) <- crs(r)
# crop   <- crop(r, b)


#########################################################################################################################
## Create difference layers for all GCMs, then sum them up. 
## For each species:
GCMS = (1:6)
total = stack()
  
for (GCM in GCMS) { 
  
  GCM_minus_current = overlay(eval(parse(text = sprintf('suit_ras%s_thresh', GCM))),
                              current_suit_thresh,
                              fun = function(r1, r2) {return (r1 - r2)})
  
  plot(GCM_minus_current,
       main = gsub('_', ' ', (sprintf('%s  - current thresh > %s', names(suit.list[[GCM]]), thresh))))
  
}


## Current – GCM1 
Acacia_implexa_1_future_minus_current = overlay(suit_ras1_thresh,
                                                current_suit_thresh,
                                                fun = function(r1, r2) {return (r1 - r2)})

names(Acacia_implexa_1_future_minus_current) = names(suit.list[[1]])


## Current – GCM2 
Acacia_implexa_2_future_minus_current = overlay(suit_ras2_thresh,
                                                current_suit_thresh,
                                                fun = function(r1, r2) {return (r1 - r2)})

names(Acacia_implexa_2_future_minus_current) = names(suit.list[[2]])


## Current – GCM3 
Acacia_implexa_3_future_minus_current = overlay(suit_ras3_thresh,
                                                current_suit_thresh,
                                                fun = function(r1, r2) {return (r1 - r2)})

names(Acacia_implexa_3_future_minus_current) = names(suit.list[[3]])


## Current – GCM4
Acacia_implexa_4_future_minus_current = overlay(suit_ras4_thresh,
                                                current_suit_thresh,
                                                fun = function(r1, r2) {return (r1 - r2)})

names(Acacia_implexa_4_future_minus_current) = names(suit.list[[4]])


## Current – GCM5
Acacia_implexa_5_future_minus_current = overlay(suit_ras5_thresh,
                                                current_suit_thresh,
                                                fun = function(r1, r2) {return (r1 - r2)})

names(Acacia_implexa_5_future_minus_current) = names(suit.list[[5]])


## Current – GCM6
Acacia_implexa_6_future_minus_current = overlay(suit_ras6_thresh,
                                                current_suit_thresh,
                                                fun = function(r1, r2) {return (r1 - r2)})

names(Acacia_implexa_6_future_minus_current) = names(suit.list[[6]])


##
GCM_diffs   =  Reduce("+", list(Acacia_implexa_1_future_minus_current,
                                Acacia_implexa_2_future_minus_current,
                                Acacia_implexa_3_future_minus_current,
                                Acacia_implexa_4_future_minus_current,
                                Acacia_implexa_5_future_minus_current,
                                Acacia_implexa_6_future_minus_current))

plot(GCM_diffs)


#########################################################################################################################
## Create difference layers for one combined GCM
mean_suit_thresh  = thresh_greater(mean.suit)
plot(mean_suit_thresh)


##
plot(overlay(mean_suit_thresh,
             current_suit_thresh,
             fun = function(r1, r2) {return (r1 - r2)}))   



#########################################################################################################################
## When we subtract the current binary layer (1, 0) from the future binary layer, we have: 

## F0 - C0 =  0 (no data in either layer)
## F0 - C1 = -1 (LOSS across all GCMs)
## F1 - C0 =  1 (GAIN across all GCMs) 
## F1 - C1 =  0 (NO CHANGE: also no data before the overlay...)


## However, this changes when we subtract the current binary layer (taking values 1, 0) from the future integer layer 
## (taking values from 0 to 6). Here, any negative value (-1) is a loss: the species was predicted to be present in that cell 
## based on current conditions, but predicted to be absent under future conditions (i.e. the future layer did not meet the suitabiity
## threshold in that location at that time.

## However, positive values (1-6) can either be a gain, or, no change. The difference needs to be accounted for, because otherwise
## we don't know if the species was predicted to occur. Can we fix this by masking the overlay to just cells with data?
## Effectively this means excluding 0 values from the overlay.

## Not all these possibilities will occur, but they are : 

## F0 - C0 =  0 no data in either layer
## F0 - C0 =  0 (NO CHANGE according to all GCMs)
## F0 - C1 =  -1 (LOSS according to all GCMs)

## F1 - C1 =  0 (NO CHANGE according to one GCM: also, no data before the overlay)
## F1 - C0 =  1 (GAIN according to one GCM)

## F2 - C1 =  1 (NO CHANGE, according to two GCMs)
## F2 - C0 =  2 (GAIN according to two GCMs)

## F3 - C1 =  2 (NO CHANGE, according to three GCMs)
## F3 - C0 =  3 (GAIN, according to three GCMs)

## F4 - C1 =  3 (NO CHANGE, according to four GCMs)
## F4 - C0 =  4 (GAIN, according to four GCMs)

## F5 - C1 =  4 (NO CHANGE, according to five GCMs)
## F5 - C0 =  5 (GAIN according to five GCMs)

## F6 - C1 =  5 (NO CHANGE, according to six GCMs?)
## F6 - 0  =  6 (GAIN, according to six GCMs?)


#########################################################################################################################
## Do the different ways of calculating difference change the combinations?
unique(combo_suit_band)
unique(current_suit_thresh)
unique(discrete_future_minus_current)

## F0 - C0 =   0 (no data in either layer) 
## F0 - C0 =   0 (no change according to all GCMs) 
## F0 - C1 =  -1 (LOSS according to all GCMs) 


## F3 - C0 =   3 (GAIN according to <3 GCMs) 
## F3 - C1 =   2 (NO CHANGE according to <3 GCMs) 


## F4 - C0 =   4 (GAIN according to >4 GCMs) 
## F4 - C1 =   3 (NO CHANGE, according to >4 GCMs) 


unique(Acacia_implexa_1_future_minus_current)
unique(GCM_diffs)
## -6 = sum(-1 + -1 + -1 + -1 + -1 + -1)

## -6 = sum(-1 + -1 + -1 + -1 + -1 + -1)  Loss, according to 6 GCMs
## -5 = sum(-1 + -1 + -1 + -1 + -1 +  0)  Loss, according to 5 GCMs
## -4 = sum(-1 + -1 + -1 + -1 +  0 +  0)  Loss, according to 4 GCMs
## -3 = sum(-1 + -1 + -1 +  0 +  0 +  0)  Loss, according to 3 GCMs
##  0 = sum(0)

##  6 = sum(1 + 1 + 1 + 1 + 1 + 1)  Loss, according to 6 GCMs
##  5 = sum(1 + 1 + 1 + 1 + 1 + 0)  Loss, according to 5 GCMs
##  4 = sum(1 + 1 + 1 + 1 + 0 + 0)  Loss, according to 4 GCMs
##  3 = sum(1 + 1 + 1 + 0 + 0 + 0)  Loss, according to 3 GCMs




