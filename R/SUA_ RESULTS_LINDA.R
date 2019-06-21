#########################################################################################################################
###################### PREDICT MAXENT TO FUTURE CLIMATES AND SUMARISE THE RESULTS ####################################### 
#########################################################################################################################

  
  
##This code summarises Maxent output to be incorporated into paper 1. A draft of the paper is currently in Google Docs (https://docs.google.com/document/d/1BZUX_OlA9tbATMC_6s6R3zw9UdyZErnMgkUUw6zAD8c/edit?userstoinvite=hugh.burley@gmail.com&ts=5b272325)


# Summarise data. There should be 176 species and 82 SUAs
length(unique(my.dat$SPECIES))
length(unique(my.dat$SUA))   



## Next, subset data for 2030, 2050 and 2070. Note, data for the baseline period will be within each of these.
dat.2030 <- subset(my.dat, PERIOD == "30")
dim(dat.2030)

dat.2050 <- subset(my.dat, PERIOD == "50")
dim(dat.2050)

dat.2070 <- subset(my.dat, PERIOD == "70")
dim(dat.2070)


## The above should have the same dimensions. This means that there are records for each species/SUA for the three future time periods.
## If any are missing, this is the code to find out which ones:
missing2030 <- setdiff(test$`unique(new.dat$SPECIES)`,s30$`unique(dat.2030$SPECIES)`)


## ANALYSIS
## Section A: How many species have at least one pixel with suitable conditions in an SUA under the baseline conditions?
  sua.species <- filter(my.dat, PERIOD == "30", CURRENT_SUITABLE > "0") %>% 
  group_by(SUA) %>% 
  summarise(MAT     = mean(CURRENT_MAT),
            MAP     = mean(CURRENT_MAP),
            MAXT    = mean(CURRENT_MAXT),
            PET     = mean(CURRENT_PET),
            AI      = mean(CURRENT_AI),
            n_sua   = length(CURRENT_SUITABLE)) 

summary(sua.species)
sd(sua.species$n_sua)



## Section B: How many species have at least one pixel with suitable conditions in an SUA in 2030?
sua.species.30 <- filter(my.dat, PERIOD == "30", FUTURE_SUITABLE > "0") %>% 
  group_by(SUA) %>% 
  summarise(n_sua.30 = length(FUTURE_SUITABLE)) 

summary(sua.species.30)
sd(sua.species.30$n_sua.30)


## 2070
sua.species.70 <- filter(my.dat, PERIOD == "70", FUTURE_SUITABLE > "0") %>% 
  group_by(SUA) %>% 
  summarise(n_sua.70 = length(FUTURE_SUITABLE)) 

summary(sua.species.70)
sd(sua.species.70$n_sua.70)


Section C: How many species are 'gained'? That is, there is no suitable habitat in the baseline, but suitable habitat DOES exist in 2030 or 2070?
  ```{r}
sua.species.gain.30 <- filter(my.dat, CURRENT_SUITABLE == "0", PERIOD == "30", FUTURE_SUITABLE > "0") %>% 
  group_by(SUA) %>% 
  summarise(n_sua.gain.30 = length(FUTURE_SUITABLE)) 

summary(sua.species.gain.30)
sd(sua.species.gain.30$n_sua.gain.30)


sua.species.gain.70 <- filter(my.dat, CURRENT_SUITABLE == "0", PERIOD == "70", FUTURE_SUITABLE > "0") %>% 
  group_by(SUA) %>% 
  summarise(n_sua.gain.70 = length(FUTURE_SUITABLE)) 

summary(sua.species.gain.70)
sd(sua.species.gain.70$n_sua.gain.70)
```

Section D: How many species are 'lost'? That is, there is suitable habitat in the baseline but NO suitable habitat exists in 2030 or 2070?
  ```{r}
sua.species.loss.30 <- filter(my.dat, CURRENT_SUITABLE > "0", PERIOD == "30", FUTURE_SUITABLE == "0") %>% 
  group_by(SUA) %>% 
  summarise(n_sua.loss.30 = length(FUTURE_SUITABLE)) 

summary(sua.species.loss.30)
sd(sua.species.loss.30$n_sua.loss.30)


sua.species.loss.70 <- filter(my.dat, CURRENT_SUITABLE > "0", PERIOD == "70", FUTURE_SUITABLE == "0") %>% 
  group_by(SUA) %>% 
  summarise(n_sua.loss.70 = length(FUTURE_SUITABLE)) 

summary(sua.species.loss.70)
sd(sua.species.loss.70$n_sua.loss.70)
```

Section E: Next, bind all data
```{r}
## Clunky method to join tables, because sua.species.gain.30 has 81 SUAs not 82.
sua.species.2 <- merge(sua.species, sua.species.gain.30, all = TRUE)
all.sua.species <- cbind(sua.species[,1:7], 
                         sua.species.30[,2], 
                         sua.species.loss.30[,2], 
                         sua.species.2[,8], 
                         sua.species.70[,2],
                         sua.species.loss.70[,2], 
                         sua.species.gain.70[,2])

#all.sua.species$change.30 <- (all.sua.species$n_sua.30-all.sua.species$n_sua)/all.sua.species$n_sua
#all.sua.species$change.70 <- (all.sua.species$n_sua.70-all.sua.species$n_sua)/all.sua.species$n_sua

summary(all.sua.species)

sd(all.sua.species$n_sua)
sd(all.sua.species$n_sua.30)
sd(all.sua.species$n_sua.loss.30)
sd(all.sua.species$n_sua.gain.30)
sd(all.sua.species$n_sua.70)
sd(all.sua.species$n_sua.loss.70)
sd(all.sua.species$n_sua.gain.70)
```


Section F: Calculate the number of SUAs a species currently has suitable climate in. 
```{r}
species.records.sua <- filter(dat.2030, CURRENT_SUITABLE > "0") %>%
  group_by(SPECIES) %>%
  summarise(n = length(CURRENT_SUITABLE))

summary(species.records.sua)
sd(species.records.sua$n)
```
##### The above code shows that each species currently has suitable climate in 36 SUAs (± 16 SUAs). 


Section G: Calculate the area (km2) of suitable habitat per species, across all 82 SUAs. Do this for both 2030 and 2070.
```{r}
my.dat.2030 <-subset(my.dat, PERIOD == "30") %>%
  group_by(SPECIES) %>%
  summarise(n           = sum(CURRENT_SUITABLE),
            n.30        = sum(FUTURE_SUITABLE),
            n.30.stable = sum(STABLE))

summary(my.dat.2030)
sd(my.dat.2030$n)
sd(my.dat.2030$n.30)


my.dat.2070 <-subset(my.dat, PERIOD == "70") %>%
  group_by(SPECIES) %>%
  summarise(n           = sum(CURRENT_SUITABLE),
            n.70        = sum(FUTURE_SUITABLE),
            n.70.stable = sum(STABLE))

summary(my.dat.2070)
sd(my.dat.2070$n.70)
```

Section H: Calculate the percent change in area of suitable habitat per species, across all 82 SUAs. Do this for both 2030 and 2070. This will be put in column "n.30.per" or "n.70.per". In the code below, n = area in baseline, n.30 = area in 2030, n.30.per is this as a proportion.
```{r}
my.dat.2030$n.30     <- my.dat.2030$n.30 + 1
my.dat.2030$n        <- my.dat.2030$n    + 1
my.dat.2030$n.30.per <- my.dat.2030$n.30/my.dat.2030$n
summary(my.dat.2030, na.rm=TRUE)
sd(my.dat.2030$n.30.per)

my.dat.2070$n.70     <- my.dat.2070$n.70 + 1
my.dat.2070$n        <- my.dat.2070$n    + 1
my.dat.2070$n.70.per <- my.dat.2070$n.70/my.dat.2070$n
summary(my.dat.2070)
sd(my.dat.2070$n.70.per)

```






Leftover.
```{r}
#my.dat.2030$n <- my.dat.2030$n+1
#my.dat.2030$n.2030 <- my.dat.2030$n.2030+1
#my.dat.2030$P_CHANGE <- round((((my.dat.2030$n.2030-my.dat.2030$n)/my.dat.2030$n)*100),4)
dim(subset(my.dat.2030, P_CHANGE < -25))
dim(subset(my.dat.2030, P_CHANGE < -50))
dim(subset(my.dat.2030, P_CHANGE >= 0))
```

```{r}
my.dat.2070$n <- my.dat.2070$n+1
my.dat.2070$n.2070 <- my.dat.2070$n.2070+1
my.dat.2070$P_CHANGE <- round((((my.dat.2070$n.2070-my.dat.2070$n)/my.dat.2070$n)*100),4)
dim(subset(my.dat.2070, P_CHANGE < -25))
dim(subset(my.dat.2070, P_CHANGE < -50))
dim(subset(my.dat.2070, P_CHANGE >= 0))
```