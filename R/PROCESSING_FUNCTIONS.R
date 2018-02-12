#########################################################################################################################
########################################### FUNCTIONS FOR MAXENT POST-PROCESSING ######################################## 
#########################################################################################################################


#########################################################################################################################
## TSS FUNCTION
#########################################################################################################################


## Calculate the mean and sd of maximum test TSS from a Maxent model
maxtss <- function(x) {
  
  ## x: a file path to a directory containing Maxent outputs
  library(raster)
  library(rmaxent) ## devtools::install_github('johnbaums/rmaxent')
  
  a <- read.csv(file.path(x, '/full/absence')) #names(a)
  f <- function(m, preds) {
    
    pred_p <- subset(read.csv(preds), Test.or.train == 'test')$Logistic.prediction
    pred_a <- c(rmaxent::project(m, a, quiet = TRUE)$prediction_logistic,
                subset(read.csv(preds), Test.or.train == 'train')$Logistic.prediction)
    
    pred <- c(pred_p, pred_a)  
    obs  <- rep(1:0, c(length(pred_p), length(pred_a)))
    thr  <- seq(0, 1, len = 100)
    
    sensspec  <- t(sapply(thr, function(x) {
      pred    <- factor(pred > x, c(FALSE, TRUE), 0:1)
      confuse <- prop.table(table(obs, pred), margin = 1)
      c(sens = confuse[2, 2], spec = confuse[1, 1])
      
    }))
    
    c(max_tss = max(rowSums(sensspec)) - 1,
      thr     = thr[which.max(rowSums(sensspec))])
    
  }
  
  max_tss <- mapply(f, list.files(x, '\\.lambdas$', full.names = TRUE),             ## try changing the path, and using recursive
                    list.files(x, 'samplePredictions\\.csv$', full.names = TRUE))   ## try changing the path, and using recursive
  
  out <- t(max_tss)
  rownames(out) <- basename(rownames(out))
  list(max_tss      = out, 
       max_tss_mean = mean(out[, 'max_tss']), 
       max_tss_sd   = sd(out[, 'max_tss']))
  
}





#########################################################################################################################
## OTHER GISTS
#########################################################################################################################


similarity <- function(x, ref, full = FALSE) {
  
  # x: a list, matrix, or data.frame where each column/element represents
  #    focal values of an environmental variable.
  
  # ref: a list, matrix or data.frame where each column/element represents
  #      reference values for an environmental variable (corresponding to those 
  #      given in x).
  
  # full: (logical) should similarity values be returned for all variables? If
  #       FALSE, then only the minimum similarity scores across variables will be 
  #       returned.
  
  ## double-check the conditions below
  if(is(ref, 'list')) {
    
    ref <- as.data.frame(ref)
    
  } 
  
  if(is.null(dim(ref))) {
    
    rng <- as.data.frame(range(ref, na.rm=TRUE))
    
  } else {
    
    rng <- as.data.frame(apply(ref, 2, range, na.rm=TRUE))
    
  }
  
  pct_less <- mapply(function(x, ref) {
    
    findInterval(x, sort(ref))/length(ref)
    
  }, x, ref, SIMPLIFY=FALSE)
  
  sim <- mapply(function(f, rng, p) {
    
    ifelse(f==0, (p-rng[1])/diff(rng)*100,
           ifelse(f > 0 & f <= 0.5, f*200,
                  ifelse(f > 0.5 & f < 1, (1-f)*200,
                         (rng[2]-p)/diff(rng)*100)))
  }, pct_less, rng, x)
  
  min_sim <- if(is.matrix(sim)) apply(sim, 1, min) else(min(sim))
  
  if(isTRUE(full)) {
    
    list(similarity = sim, similarity_min=min_sim)
    
  } else min_sim
  
}






#########################################################################################################################
#####################################################  TBC ############################################################## 
#########################################################################################################################



