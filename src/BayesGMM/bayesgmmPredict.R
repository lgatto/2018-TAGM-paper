#' BayesGMM mixture model for LOPIT dataset prediction routine
#' 
#' 
#' @param object An instance of class MSnset
#' @param optimRes Output produced by bayesgmmOptimsation
#' 
#' 
#' 
#' @return 
#' 
#' 
#' 

bayesgmmPredict <- function(object,
                            optimRes
                            ){
  
  eps <- optimRes$eps
  mu <- optimRes$mu
  sigma <- optimRes$sigma
  weights <- optimRes$weights
  
  
  X <- exprs(unknownMSnSet(object))
  K <- length(getMarkerClasses(object))
  
  a <- matrix(0, nrow = nrow(X), ncol = K)
  b <- matrix(0, nrow = nrow(X), ncol = K)
  predictionprob <- matrix(0, nrow = nrow(X), ncol =  K)
  organellealloc <- matrix(0, nrow = nrow(X), ncol = 2)
  
  M <- colMeans(exprs(object))
  V <- cov(exprs(object))/2
  
  for(j in 1:K){
    a[, j] <- log( weights[j] ) + log( 1 - eps) + mvtnorm::dmvnorm(X, mean = mu[j, ], sigma = sigma[j, , ], log = T)
    b[, j] <- log( weights[j] ) + log(eps) + mvtnorm::dmvt(X, delta = M, sigma = V, df = 4, log = T)
  }
  
  #correct for underflow by adding constant
  ab <- cbind(a,b)
  c <- apply(ab, 1, max)
  ab <- ab - c                   #add constant
  ab <- exp(ab)/rowSums(exp(ab)) #normlise
  a <- ab[, 1:K]
  b <- ab[, (K + 1):(2 * K)]
  predictionprob <- a + b
  
  organellealloc[, 1] <- getMarkerClasses(object)[apply(a, 1, which.max)]
  proballoc <- apply(a, 1, which.max)
  
  for(i in 1:nrow(X)){
    organellealloc[i, 2] <- as.numeric(a[i, proballoc[i]])
  }
  rownames(a) <- rownames(unknownMSnSet(object))
  rownames(b) <- rownames(unknownMSnSet(object))
  rownames(predictionprob) <- rownames(unknownMSnSet(object))
  rownames(organellealloc) <- rownames(unknownMSnSet(object))
  
  pred <- c(organellealloc[, 1], as.character(fData(markerMSnSet(object))[,"markers"]))
  prob <- c(organellealloc[, 2], rep(1,length(fData(markerMSnSet(object))$markers)))
  
  names(prob) <- c(rownames(unknownMSnSet(object)), rownames(fData(markerMSnSet(object))))
  names(pred) <- c(rownames(unknownMSnSet(object)), rownames(fData(markerMSnSet(object))))
  
  fData(object)$predicted.allocation <- pred[rownames(fData(object))] 
  fData(object)$predicted.probability <- prob[rownames(fData(object))]
  
  return(list(orgalloc = a, globalalloc = b, predictionprob = predictionprob, finalassignementprob = organellealloc, object = object))
}