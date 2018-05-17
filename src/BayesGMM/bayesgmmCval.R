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

bayesgmmCval <- function(object,
                         times = 10,
                         test.size = 0.2
                         ){
  
  
  marker.data <- markerMSnSet(object)
  X <- pRoloc:::subsetAsDataFrame(marker.data, "markers", train = TRUE)
  K <- length(getMarkerClasses(object))
  
  .testPartitions <- .cmMatrices <- vector("list", length = times)
  f1score <- matrix(0, times)
  .f1Matrices <- matrix(0, times)
  cmlist <- vector("list", length = times)
  quadloss <- vector("list", length = times)
  
  
  for(i in seq_along(cmlist)){
    
    .size <- ceiling(table(fData(marker.data)$markers) * test.size) #get sizes
    .size <- .size[unique(fData(marker.data)$markers)] #strata needs size to be ordered as they appear in data
    
    test.idx <- strata(X, "markers", size = .size, method = "srswor")$ID_unit #get strata indexes
    
    .test1   <- MSnSet(exprs(marker.data)[test.idx, ], fData(marker.data[test.idx, ]), pData(marker.data)) ## 'unseen' test set
    .train1  <- MSnSet(exprs(marker.data)[-test.idx, ], fData(marker.data[-test.idx, ]), pData(marker.data)) ## 'seen' training set
    
    test.markers <- fData(.test1)$markers #save true marker labels
    mydata <- combine(.test1, .train1) #create new combined MSnset
    levels(fData(mydata)$markers) <- c(levels(fData(mydata)$markers), "unknown")
    fData(mydata)[rownames(.test1), "markers"] <- "unknown" #hide marker labels
    
    params <- bayesgmmOptimisation(object = mydata, iterations = 100) #optimsiation
    
    if (abs(params$loglike[100]-params$loglike[99])>0.5){
      print("EM algorithm has not converged")
    }
    
    res <- bayesgmmPredict(object = mydata, optimRes = params)
    
    data <- factor(fData(res$object[rownames(.test1),])$predicted.allocation, levels = getMarkerClasses(mydata))
    reference <- factor(test.markers, levels = getMarkerClasses(mydata))
    
    cmlist[[i]] <- conf <- confusionMatrix(data = data, reference = reference)$table

    allocmatrix <- matrix(0,length(test.idx), length(getMarkerClasses(mydata)))
    
    #create allocation matrix
    for(j in seq_along(test.idx)){
      allocmatrix[j, as.numeric(factor(test.markers), seq(1,length(unique(test.markers))))[j]] <- 1
    }

    quadloss[[i]] <- sum((allocmatrix - res$orgalloc[rownames(.test1),])^2) #compute quad loss
    
    
  }
   
  return(list(cmlist = cmlist, quadloss = quadloss))
}   