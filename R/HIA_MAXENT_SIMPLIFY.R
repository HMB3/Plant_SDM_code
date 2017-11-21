## This should replace John's code
HIA_SIMPLIFY = function (occ, bg, path, species_column = "searchTaxon", response_curves = FALSE, 
                         logistic_format = TRUE, type = "PI", cor_thr, pct_thr, k_thr, 
                         quiet = TRUE) 
  
{
  if (missing(path)) {
    
    save <- FALSE
    path <- tempdir()
    
  }
  
  else save <- TRUE
  occ_by_species <- split(occ, occ[[species_column]])
  bg_by_species  <- split(bg, bg[[species_column]])
  
  
  ## 
  # str(occ_by_species)
  # str(bg_by_species)
  
  ## This code breaks on my data, because the background points are taken from points that are not the species 
  ## background <- subset(SDM.DATA.ALL, searchTaxon != x)
  
  if (!identical(sort(names(occ_by_species)), sort(names(bg_by_species)))) {    ##
    print(paste0("The same set of species names must exist in occ and bg"))
  }
  
  ## 
  type <- switch(type, PI = "permutation.importance", PC = "contribution", 
                 stop("type must be either \"PI\" or \"PC\".", call. = FALSE))
  
  ## Maxent args?
  args <- c("threshold=false", "hinge=false")
  if (isTRUE(response_curves)) 
    args <- c(args, "responsecurves=TRUE")
  
  if (isTRUE(logistic_format)) 
    args <- c(args, "outputformat=logistic")
  
  
  ## Explain
  lapply(names(occ_by_species), function(name) {
    if (!quiet) 
      message("\n\nDoing ", name)
    
    ## The problem is in here: conflict between these lines and the main background argument:
    ## background <- subset(SDM.DATA.ALL, searchTaxon != x)
    ## Not sure why this dataframe needs to be matched?
    name_ <- gsub(" ", "_", name)
    swd <- rbind(occ_by_species[[name]], bg_by_species[[name]])
    #swd <- swd[, -match(species_column, names(swd))]
    ## str(swd)
    
    ##
    if (ncol(swd) < k_thr) 
      stop("Initial number of variables < k_thr")
    
    ## what does round do?
    round(cor(as.data.frame(swd)[2:20], use = "pairwise"), 2)
    swd.cor = as.data.frame(swd)[2:20]
    pa <- rep(1:0, c(nrow(occ_by_species[[name]]), nrow(bg_by_species[[name]])))  ## problem here
    ok <- as.character(usdm::vifcor(swd.cor, maxobservations = nrow(swd), 
                                    th = cor_thr)@results$Variables)              ## ok
    
    ##
    swd_uncor             <- swd@data[, ok]
    swd_uncor$searchTaxon <- swd@data$searchTaxon                                 ## 
    d <- file.path(path, name_, "full")
    m <- dismo::maxent(swd_uncor, pa, args = args, path = d)  ## problem here
    
    ##
    if (isTRUE(save)) 
      saveRDS(m, file.path(d, "model.rds"))
    
    pct <- m@results[grep(type, rownames(m@results)), ]
    pct <- sort(pct[pct > 0])
    
    names(pct) <- sub(paste0("\\.", type), "", names(pct))
    
    
    
    if (min(pct) >= pct_thr || length(pct) <= k_thr) {
      if (isTRUE(save)) {
        d_out <- file.path(path, name_, "final")
        dir.create(d_out)
        file.copy(list.files(d, full.names = TRUE), 
                  d_out, recursive = TRUE)
        
      }
      
      return(m)
      
    }
    
    ## Explain
    while (min(pct) < pct_thr && length(pct) > k_thr) {
      
      message("Dropping ", names(pct)[1])
      swd_uncor <- swd_uncor[, -match(names(pct)[1], names(swd_uncor))]
      tmp <- tempfile()
      
      if (!quiet) 
        message(sprintf("%s variables: %s", ncol(swd_uncor), 
                        paste0(colnames(swd_uncor), collapse = ", ")))
      m <- dismo::maxent(swd_uncor, pa, args = args, path = tmp)
      pct <- m@results[grep(type, rownames(m@results)), 
                       ]
      pct <- sort(pct)
      names(pct) <- sub(paste0("\\.", type), "", names(pct))
      
    }
    
    ## Explain
    if (isTRUE(save)) {
      
      d_out <- file.path(path, name_, "final")
      file.copy(tmp, file.path(path, name_), recursive = TRUE)
      file.rename(file.path(path, name_, basename(tmp)), 
                  d_out)
      saveRDS(m, file.path(path, name_, "final/model.rds"))
      
    }
    
    ## return the set of predictors which are less correlated
    return(m)
    
  })
  
}