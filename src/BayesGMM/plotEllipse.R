#' A function to plot probabiltiy ellipses on marker data
#' 
#' 
#' 
#' 
#' @param object An instance of class MSnset
#' 
#' 
#' @return A PCA plot of the marker data with probability ellipises
#' 

plotEllipse <- function(object,
                     components = c(1,2),
                     mu0 = NULL,
                     lambda0 = 0.01,
                     nu0 = NULL,
                     S0 = NULL,
                     method = "MAP",
                     ...){
  
  #make marker MSnSet
  markersubset <- markerMSnSet(object)
  marker.names <- getMarkerClasses(object)
  mydata <- exprs(markersubset)
  K <- length(marker.names)
  
  #create dataset
  D <- ncol(mydata)
  N <- nrow(mydata)
  
  if (is.null(nu0)) {
    nu0 <- D + 2
  }
  if (is.null(S0)) { 
    S0 <- diag( colSums(( mydata - mean( mydata)) ^ 2) / N)/( K ^ (1/D))
  }
  if(is.null(mu0)){
    mu0 <- colMeans( mydata)
  }
  
  
  #create storage for posterior parameters
  mk <- matrix(0, nrow = K, ncol = D)
  lambdak <- matrix(0, nrow = K, ncol = 1)
  nuk <- matrix(0, nrow = K, ncol = 1)
  sk <- array(0, dim = c(K, D, D))
  
  #create storage for cluster parameters
  muk <- matrix(0, nrow = K, ncol = D)
  sigmak <- array(0, dim = c(K, D, D))
  xk <- matrix(0, nrow = K, ncol = D)
  
  #update prior with training data
  nk <- tabulate(fData(markersubset)[, "markers"])
  for(j in 1:K){
    xk[j, ] <- colMeans(mydata[fData(markersubset)[, "markers"] == getMarkerClasses(markersubset)[j], ])
  }
  lambdak <- lambda0 + nk
  nuk <- nu0 + nk
  mk <- (nk * xk + lambda0 * mu0) / lambdak
  
  for(j in 1:K){
    sk[j, , ] <- S0 + t(mydata[fData(markersubset)[,"markers"] == getMarkerClasses(markersubset)[j], ]) %*% 
      mydata[fData(markersubset)[, "markers"] == getMarkerClasses(markersubset)[j],] +
      lambda0 * mu0 %*% t(mu0) - lambdak[j] * mk[j, ] %*% t(mk[j, ]) 
  }
  #initial posterior mode
  muk <- mk
  if(method == "MAP"){
  for(j in 1:K){
    sigmak[j, , ] <- sk[j, , ] / (nuk[j] + D + 1)
  }
  }else if(method == "MCMC"){
  for(j in 1:K){
    sigmak[j, , ] <- sk[j, , ] / (nuk[j] - D - 1)
  }
  }
  
  
  #scale data correctly and format
  sigmaROT <- array(0, c(K, D, D))
  par(mfrow = c(1,1))
  pca <- prcomp(mydata, center = TRUE, scale = FALSE)
  plot2D(markersubset, dims = components, methargs = list(scale = FALSE, center= TRUE),...)
  meanROT <- scale(muk, center = pca$center, scale = pca$scale) %*% pca$rotation
  points(meanROT[, components], pch=19, col='black')
  for(j in 1:K){
    sigmaROT[j, , ] <- scale(sigmak[j, , ], scale=pca$scale, center = FALSE)
  }
  
  #ellipses from Mixtools package
  for(j in 1:K){
    ellipse(meanROT[j,components],sigmaROT[j,components,components], alpha = .1, npoints=250, newplot = FALSE)
    ellipse(meanROT[j,components],sigmaROT[j,components,components], alpha = .01, npoints=250, newplot = FALSE)
    ellipse(meanROT[j,components],sigmaROT[j,components,components], alpha = .05, npoints=250, newplot = FALSE)
  }
  
}  