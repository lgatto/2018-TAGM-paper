##' Cross-validation functions for 2018-TAGM paper
##' 
##' 
##' 
##' @title Functions to performed Cross-validation on SVM algorithms
##'  
##' 
##' @svmparams An object of class svmparms
##' @object An object of class MSnSet
##' 
##' @return list of quadratics losses

svmquadloss <- function(svmparams, object){

  svmquadloss <- vector("list", length = 100) ## quadratic loss to save
  
  for (i in seq_along(svmparams@f1Matrices)) {
    
    index <- svmparams@testPartitions[[i]] ## get test partitions
    markerSet <- markerMSnSet(object) ## create markerset
    trueclasses <- fData(markerSet)$markers[index] ## save true markers
    
    levels(fData(markerSet)$markers) <- c(getMarkerClasses(object), "unknown") ## change levels
    fData(markerSet)$markers[index] <- "unknown" ## hide markers
    
    sigma <- as.numeric(svmparams@results[i, 2]) ## get sigma
    cost <- as.numeric(svmparams@results[i, 3]) ## get cost
    
    res <- svmClassification(markerSet, fcol= "markers", sigma = sigma, cost = cost, scores = "all") ## svm classification
    scoresmatrix <- fData(res)[index, "svm.all.scores"] ## get scores as matrix
    colnames(scoresmatrix) <- unique(trueclasses) ## change names
    allocmatrix <- matrix(0, length(index), length(unique(trueclasses)))
    
    ## create allocation matrix
    for (j in seq_along(index)) {
      allocmatrix[j, as.numeric(factor(trueclasses), seq(1,length(unique(trueclasses))))[j]] <- 1 
    }
    
    svmquadloss[[i]] <- sum((scoresmatrix[,getMarkerClasses(markerSet)] - allocmatrix)^2) ## compute quadloss
  }

return(svmquadloss = svmquadloss)
}
##' 
##' @title Functions to performed Cross-validation on KNN algorithm
##'  
##' 
##' @param knnparams An object of class knnparms
##' @param object An object of class MSnSet
##' 
##' @return list of quadratics losses


knnquadloss <- function(knnparams, object){
  
  knnquadloss <- vector("list", length = 100) #quadratic loss to save
  
    for(iLoss in seq_along(knnquadloss)){
      
      markerSet <- markerMSnSet(object)
      
      ##get index
      index <- knnparams@testPartitions[[iLoss]]
      trueclasses <- fData(markerSet)$markers[index]
      
      levels(fData(markerSet)$markers) <- c(getMarkerClasses(object), "unknown") ## change levels
      fData(markerSet)$markers[index] <- "unknown" #hide markers
      
      nbr <- knnparams@results[iLoss, 2]
      train <- exprs(markerSet)[-index, ]
      test <- exprs(markerSet)[index, ]
      
      K <- length(getMarkerClasses(object))
      
      res <- FNN::knn(train, test, cl = fData(markerSet)$markers[-index], k = nbr, prob = TRUE) ## knn classification
      nn.index <- attr(res, "nn.index")
      distMat <- attr(res, "nn.dist") ## get atrributes
      
      scoresmatrix <- matrix(0, length(index), length(getMarkerClasses(object)))
      colnames(scoresmatrix) <- getMarkerClasses(object) ## change names
      
      irates <- tabulate(fData(markerSet)$markers[-index])/sum(tabulate(fData(markerSet)$markers[-index]))
      
      
      for(i in seq_along(index)){
        scoresmatrix[i, ] <- (0.5 * irates * K)/(nbr + 0.5 * K) ## compute non-parametric probability
        
        for(j in seq(1:K)[table(fData(markerSet)$markers[-index][nn.index[i,]])!=0]){
          count <- as.numeric(table(fData(markerSet)$markers[-index][nn.index[i, ]])[j])
          scoresmatrix[i,j] <-  (count + 0.5 * irates[j] * K)/(nbr + 0.5 * K)
        }
      }
      scoresmatrix <- scoresmatrix/rowSums(scoresmatrix) ## normalise
    
      #create allocation matrix
      allocmatrix <- matrix(0, length(index), length(unique(trueclasses)))
      for(j in seq_along(index)){
        allocmatrix[j, as.numeric(factor(trueclasses), seq(1,length(unique(trueclasses))))[j]] <- 1 
      }
      
      knnquadloss[[iLoss]] <- sum((scoresmatrix[,getMarkerClasses(markerSet)] - allocmatrix)^2) #compute quad loss
    }
  
  return(knnquadloss = knnquadloss)
}


##' Calculate Cross-validation Metrics using TAGM-MAP
##' 
##' 
##' 
##' 
##'@param object  An object of class MSnset
##'@param times  Number of times to run cross-validation
##'@param test.size The test size split of the data. Default is 0.2.
##'
##'@return
##'

tagmMapCval <- function(object,
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
  
  
  for (i in seq_along(cmlist)) {
    
    .size <- ceiling(table(fData(marker.data)$markers) * test.size) ## get sizes
    .size <- .size[unique(fData(marker.data)$markers)] ## strata needs size to be ordered as they appear in data
    
    test.idx <- strata(X, "markers", size = .size, method = "srswor")$ID_unit ## get strata indexes
    
    .test1   <- MSnSet(exprs(marker.data)[test.idx, ],
                       fData(marker.data[test.idx, ]), pData(marker.data)) ## 'unseen' test set
    .train1  <- MSnSet(exprs(marker.data)[-test.idx, ],
                       fData(marker.data[-test.idx, ]), pData(marker.data)) ## 'seen' training set
    
    test.markers <- fData(.test1)$markers ## save true marker labels
    mydata <- combine(.test1, .train1) ## create new combined MSnset
    levels(fData(mydata)$markers) <- c(levels(fData(mydata)$markers), "unknown")
    fData(mydata)[rownames(.test1), "markers"] <- "unknown" ## hide marker labels
    
    params <- pRoloc:::tagmMapTrain(object = mydata, numIter = 100) #optimsiation
    
    if (abs(params@posteriors$logposterior[100]-params@posteriors$logposterior[99])>0.5){
      print("EM algorithm has not converged")
    }
    
    res <- tagmMapPredict(object = mydata, params = params, probJoint = TRUE)
    
    data <- factor(fData(res)$tagm.map.allocation[seq_len(nrow(.test1))], levels = getMarkerClasses(mydata))
    reference <- factor(test.markers, levels = getMarkerClasses(mydata))
    
    cmlist[[i]] <- conf <- confusionMatrix(data = data, reference = reference)$table
    
    allocmatrix <- matrix(0, length(test.idx), length(getMarkerClasses(mydata)))
    
    #create allocation matrix
    for (j in seq_along(test.idx)) {
      allocmatrix[j, as.numeric(factor(test.markers), seq(1,length(unique(test.markers))))[j]] <- 1
    }
    quadloss[[i]] <- sum((allocmatrix - fData(res)$tagm.map.joint[rownames(.test1),])^2) #compute quad loss
    
  }
  
  return(list(cmlist = cmlist, quadloss = quadloss))
}

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

tagmMcmcCval <- function(test.idx,
                        object,
                        fcol = "markers",
                        numIter = 1000L,
                        burnin = 100L,
                        thin = 5L,
                        mu0 = NULL,
                        lambda0 = 0.01,
                        nu0 = NULL,
                        S0 = NULL,
                        beta0 = NULL,
                        u = 2,
                        v = 10
){
  
  suppressPackageStartupMessages({
    require(pRoloc)
    require(MSnbase)
  })
  
  
  marker.data <- markerMSnSet(object)
  X <- pRoloc:::subsetAsDataFrame(marker.data, "markers", train = TRUE)
  
  .test1   <- MSnSet(exprs(marker.data)[test.idx, ], fData(marker.data[test.idx, ]), pData(marker.data)) ## 'unseen' test set
  .train1  <- MSnSet(exprs(marker.data)[-test.idx, ], fData(marker.data[-test.idx, ]), pData(marker.data)) ## 'seen' training set
  
  test.markers <- fData(.test1)$markers #save true marker labels
  mydata <- combine(.test1, .train1) #create new combined MSnset
  levels(fData(mydata)$markers) <- c(levels(fData(mydata)$markers), "unknown")
  fData(mydata)[rownames(.test1), "markers"] <- "unknown" #hide marker labels
  
  params <- tagmMcmcTrain(object = mydata, numIter = numIter, burnin = burnin, thin = thin) ## optimsiation
  
  N <- nrow(X[test.idx, ])
  K <- length(getMarkerClasses(mydata))
  
  params <- tagmMcmcProcess(params)
  mydata <- tagmMcmcPredict(object = mydata, MCMCParams = params, probJoint = TRUE)
  
  data <- factor(fData(mydata[rownames(.test1),])$tagm.allocation, levels = getMarkerClasses(mydata))
  reference <- factor(test.markers, levels = getMarkerClasses(mydata))
  
  conf <- confusionMatrix(data = data, reference = reference)$table
  
  quad <- vector("numeric", length = length(test.idx))
  allocmatrix <- matrix(0, length(test.idx), K)
  
  for (j in seq_along(test.idx)) {
    allocmatrix[j, as.numeric(factor(test.markers), seq(1, K))[j]] <- 1 
  }
  quadloss <- sum(rowSums((allocmatrix - fData(mydata[rownames(.test1),])$tagm.mcmc.joint)^2))
  
  return(list(conf = conf, quadloss = quadloss))
} 


#' Parallel function for running cross-validation MCMC 
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

tagmMcmcCvalpar <- function(object,
                              times = 100,
                              test.size = 0.2,
                              numIter = 1000L,
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
  
  .res <- bplapply(test.idx, tagmMcmcCval,
                   object = object,
                   fcol = "markers",
                   numIter = numIter,
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





