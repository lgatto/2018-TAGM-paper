#' Calculate Cross-validation Metrics 
#' 
#' 
#' 
#' 
#'@params
#'
#'
#'@return
#'
#'
#'

bayesgmmgibbsCval <- function(object,
                         times = 100,
                         test.size = 0.2,
                         iterations = 1000L,
                         burnin = 100L,
                         thin = 5L,
                         BPPARAM = bpparam()
){
  
  
  marker.data <- markerMSnSet(object)
  X <- pRoloc:::subsetAsDataFrame(marker.data, "markers", train = TRUE)
  K <- length(getMarkerClasses(object))
  
  .testPartitions <- .cmMatrices <- vector("list", length = times)
  f1score <- matrix(0, times)
  .f1Matrices <- matrix(0, times)
  cmlist <- vector("list", length = times)
  quadloss <- vector("list", length = times) 
  .size <- ceiling(table(fData(marker.data)$markers) * test.size) #get sizes
  .size <- .size[unique(fData(marker.data)$markers)] #strata needs size to be ordered as they appear in data
  test.idx <- vector("list", length = times)
  
  
  for(i in seq_along(cmlist)){
    test.idx[[i]] <- strata(X, "markers", size = .size, method = "srswor")$ID_unit #get strata indexes
  } 
  
  .res <- bplapply(test.idx, gibbsMetric,
                             object = object,
                             fcol = "markers",
                             iterations = iterations,
                             burnin = burnin,
                             thin = thin,
                             mu0 = NULL,
                             lambda0 = 0.01,
                             nu0 = NULL,
                             S0 = NULL,
                             beta0 = NULL,
                             u = 2,
                             v = 10,
                             BPPARAM = BPPARAM
  )
  
  for(i in seq_along(cmlist)){
    cmlist[[i]] <- .res[[i]]$conf
    quadloss[[i]] <- .res[[i]]$quadloss
  }
    
    
  return(list(cmlist = cmlist, quadloss = quadloss))
}   
