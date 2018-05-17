#' BayesGMM mixture model for LOPIT datasets Optimisation routine using EM
#' 
#' 
#' @param Object An instance of class MSnset
#' @param fcol The feature metadata containing marker defintions
#' @param iterations The number of iterations of the algorithm. Default is \code{100}
#' @param mu0 The prior mean. Default is colmeans of expression data.
#' @param lambda0 The prior shrinkage. Default is 0.01
#' @param nu0 The prior degreed of freedom. Default is \code(ncol(exprs(object))+2)
#' @param beta0 The prior Dirichlet distribution concentration. Default is 1 for each class
#' @param u The prior shape parameter for Beta(u, v). Default is 2
#' @param v The prior shape parameter for Beta(u, v). Default is 999.
#' 
#' 
#' @return Maximum a posteriori (MAP) estimates of the mean, variance, component weights, structure weights
#' and the log likelihood.
#' 


bayesgmmOptimisation <- function(object,
                                 fcol = "markers",
                                 iterations = 100,
                                 mu0 = NULL,
                                 lambda0 = 0.01,
                                 nu0 = NULL,
                                 S0 = NULL,
                                 beta0 = NULL,
                                 u = 2,
                                 v = 10
                                 ) {
  
  #get expression marker data
  markersubset <- markerMSnSet(object)
  mydata <- exprs(markersubset)
  X <- exprs(unknownMSnSet(object))
  
  #get data dize
  N <- nrow(mydata)
  D <- ncol(mydata)
  K <- length(getMarkerClasses(object))
  
  if (is.null(nu0)) {
    nu0 <- D + 2
  }
  if (is.null(S0)) { 
   S0 <- diag( colSums(( mydata - mean( mydata)) ^ 2) / N)/( K ^ (1/D))
  }
  if(is.null(mu0)){
    mu0 <- colMeans( mydata)
  }
  if(is.null(beta0)){
    beta0 <- rep(1, K)
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
    xk[j, ] <- colSums(mydata[fData(markersubset)[, "markers"] == getMarkerClasses(markersubset)[j], ])/nk[j]
  }
  lambdak <- lambda0 + nk
  nuk <- nu0 + nk
  mk <- (nk * xk + lambda0 * mu0) / lambdak
  
  for(j in 1:K){
    sk[j, , ] <- S0 + t(mydata[fData(markersubset)[,"markers"] == getMarkerClasses(markersubset)[j], ]) %*% 
      mydata[fData(markersubset)[, "markers"] == getMarkerClasses(markersubset)[j],] +
      lambda0 * mu0 %*% t(mu0) - lambdak[j] * mk[j, ] %*% t(mk[j, ]) 
  }
  betak <- beta0 + nk
  
  #initial posterior mode
  muk <- mk
  for(j in 1:K){
    sigmak[j, , ] <- sk[j, , ] / (nuk[j] + D + 1)
  }
  #initial cluster probabilty weights
  pik <- (betak - 1) / (sum(betak) - K)
  
  #global parameters
  M <- colMeans(exprs(object))
  V <- cov(exprs(object))/2
  eps <- (u - 1) / (u + v - 2)
   
  #storage for Estep
  a <- matrix(0, nrow = nrow(X), ncol = K)
  b <- matrix(0, nrow = nrow(X), ncol = K)
  w <- matrix(0, nrow = nrow(X), ncol = K)
  #storage for Mstep
  xbar <- matrix(0, nrow = nrow(X), ncol = K)
  lambda <- matrix(0, K)
  nu <- matrix(0, K)
  m <- matrix(0, K, D)
  S <- array(0, c(K, D, D))
  loglike <- vector(mode = "numeric", length = iterations)
  
  for (t in 1:iterations){
    
    #E-Step
    for(k in 1:K){
      a[, k] <- log( pik[k] ) + log( 1 - eps) + mvtnorm::dmvnorm(X, mean = muk[k, ], sigma = sigmak[k, , ], log = T)
      b[, k] <- log( pik[k] ) + log(eps) + mvtnorm::dmvt(X, delta = M, sigma = V, df = 4, log = T)
    }
    
    #correct for underflow by adding constant
    ab <- cbind(a,b)
    c <- apply(ab, 1, max)
    ab <- ab - c                   #add constant
    ab <- exp(ab)/rowSums(exp(ab)) #normlise
    a <- ab[, 1:K]
    b <- ab[, (K + 1):(2 * K)]
    w <- a + b
    r <- colSums(w)
    
    #M-Step
    
    #structure weights
    eps <- (u + sum(b) - 1) / ( (sum(a) + sum(b)) + (u + v) - 2)
    
    xbar <- apply(a, 2, function(x){colSums( x * X )})
    xbar[, colSums(xbar)!=0] <- t(t(xbar[, colSums(xbar)!=0])/colSums(a)[colSums(xbar)!=0]) 
    
    #component weights
    pik <- (r + betak - 1)/(nrow(X) + sum(betak) - K)
    
    #component parameters
    lambda <- lambdak + colSums(a)
    nu <- nuk + colSums(a)
    m <- (colSums(a) * t(xbar) + lambdak * mk)/lambda

    TS <- array(0, c(K, D, D))
    for(j in 1:K){
      for(i in 1:nrow(X)){
        TS[j, , ] <- TS[j, , ] + a[i, j] * (X[i, ] - xbar[, j]) %*% t((X[i, ] - xbar[, j]))
      }
    }

    vv <- (lambdak * colSums(a))/ lambda
    for(j in 1:K){
      S[j, , ] <- sk[j, , ] + vv[j] * (xbar[, j] - mk[j,]) %*% t((xbar[, j] - mk[j,])) + TS[j, , ] 
      sigmak[j, , ] <- S[j, , ]/(nu[j] + D + 2) 
    }
    muk <- m
    
    #compute log-likelihood
    for (j in 1:K){
      loglike[t] <- loglike[t] + sum( a[, j] * mvtnorm::dmvnorm(X, mean = muk[j, ], sigma = sigmak[j, , ], log = T)) +
        sum( w[,j] * log( pik[j]) )
    }
    loglike[t] <- loglike[t] + sum(a) * log(1 - eps) + sum(b) * log(eps) + 
      sum(rowSums(b) *  mvtnorm::dmvt(X, delta = M, sigma = V, df = 4, log = T))
    
  }
  
  return(list(mu = muk, sigma = sigmak, epsilon = eps, weights = pik, loglike = loglike))
}
                                  
                                  
                                  