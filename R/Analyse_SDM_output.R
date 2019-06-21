# This code is to analyse the output of species distribution models.
# It can be used to convert maps of continuous data into binary data, 
# to calculate the size of ranges based on equal area projections,
# to calculate changes in range size.

# 11 March 2015

# Make sure workspace is clean
rm(list = ls())

#dat <- c(-12,0,-2,-16,0,-33,11,-7,-29,-14,-8,-42,-13,-21,-25,4,9,-7,3,-16,-6,-2,8,3,-5,-25,-2,-3,-12)
# load library 
library(raster)

# get directories where data are located
work.dir<-"E:\\Australian_shrubs\\analysis_papers2_3\\future\\maxent_climsoil\\"
out.dir<-"E:\\Australian_shrubs\\analysis_papers2_3\\future\\threshold_climsoil\\"

# To convert continuous data to binary, a file is needed with species names and threshold values, e.g. maxentResults.csv
dat1<-read.csv("E:\\Australian_shrubs\\analysis_papers2_3\\future\\maxent_climsoil\\maxentResults.csv",header=TRUE)

# From the above file, get the names of species and threshold.
maxthresh <- dat1[grep("average",dat1$Species),c("Species","Maximum.training.sensitivity.plus.specificity.logistic.threshold")] 


for (i in 1:length(maxthresh$Species)){
  
  sp_files<-list.files(paste(work.dir),pattern=paste0(sub(" \\(average\\)","",maxthresh[i,1]),".*avg\\.asc$"),full.names = TRUE)  #Sub substitutes some, and the \\ in front of (average) means the search has to identify the brackets.
  
  s <- stack(sp_files) # Puts all files in a stack
  sthr <- s>maxthresh[i,2] # Applies true/false to the stack. If cell value > threshold, TRUE else FALSE
  
  # Write binary rasters.
  writeRaster(sthr,filename=paste0(out.dir,names(sthr),".asc"),overwrite=TRUE,NAflag=-9999,bylayer=TRUE)
  
}




###

# Function to calculate change summaries
dist_change_binary <- function(from, to) {
  require(raster)
  s <- stack(from, to)
  d <- data.frame(
    from=from,
    to=to,
    count_from = sum(values(s[[1]]), na.rm=TRUE),
    count_to = sum(values(s[[2]]), na.rm=TRUE),
    stay_good = sum(values(s[[1]] & s[[2]]), na.rm=TRUE),
    stay_bad = sum(values(!(s[[1]] | s[[2]])), na.rm=TRUE),
    gained = sum(values(!s[[1]] & s[[2]]), na.rm=TRUE),
    lost = sum(values(s[[1]] & !s[[2]]), na.rm=TRUE))
  within(d, {
    prop_from_lost <- lost/count_from
    prop_to_gained <- gained/count_to
    prop_change <- (count_to - count_from)/count_from
    prop_stay_good <- stay_good/count_from
  })
}



ff <- list.files("E:\\Australian_shrubs\\analysis_papers2_3\\future\\threshold_climsoil\\", 
                 full.names=TRUE)
ff_list <- split(ff, sub(".*\\\\threshold_climsoil\\\\([^_]*_[^_]*).*", '\\1', ff))


change_summaries <- lapply(ff_list, function(x) {
  sapply(2:length(x), function(y) {
    from <- x[1]
    to <- x[y]
    dist_change_binary(from, to)
  })
})
change_summaries_mat <- do.call(rbind, change_summaries)
write.csv(change_summaries_mat, 'E:\\Australian_shrubs\\analysis_papers2_3\\future\\threshold_climsoil\\change_summaries.csv', row.names=FALSE)
