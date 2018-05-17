#' BayesGMM mixture model for LOPIT datasets Optimisation routine using Gibbs sampling
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
#' @return An MCMC samples of the allocation probabilities and final parameters values.
#' 


bayesgmmOptimisationGibbs <- function(object,
                                 fcol = "markers",
                                 iterations = 1000L,
                                 burnin = 100L,
                                 thin = 5L,
                                 mu0 = NULL,
                                 lambda0 = 0.01,
                                 nu0 = NULL,
                                 S0 = NULL,
                                 beta0 = NULL,
                                 u = 2,
                                 v = 10
) {
  
  if(burnin >= iterations){
    stop("Burnin is larger than iterations you will not retain any samples")
  }
  
  retained <- seq.int(burnin + 1L, iterations , thin)
  
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
  
  #global parameters
  M <- colMeans(exprs(object))
  V <- cov(exprs(object))/2
  
  #storage
  allocprob <- array(0, c(nrow(X), K, iterations))
  allocstruc <- matrix(0, nrow = nrow(X), ncol = iterations)
  allocstrucprob <- array(0, c(nrow(X), iterations, 2))
  alloc <- array(0, c(nrow(X), iterations))
  
  #initially assigned all unlabelled points to clusters greedily initially
  for(j in 1:K){
    allocprob[, j, 1] <- dmvtCpp(X, mu = mk[j, ], sigma = (1 + lambdak[j]) * sk[j,,] / (lambdak[j] * (nuk[j] - D + 1)),
                     df = nuk[j] - D + 1, log = T, ncores_ = 1, isChol_ = F)
  }

  alloc[, 1] <- apply(X = allocprob[,,1], 1, FUN = which.max) 
  
  #initially assign all components to global component i.e. allocstruc = 0
  allocstruc[, 1] <- 0
  
  #initial allocation statistics
  tau1 <- sum(allocstruc[, 1] == 1) + N
  tau2 <- sum(allocstruc[, 1] == 0)
  
  for(t in seq.int(2L, iterations)){
    
    #consider each protein in turn

    for(i in 1:nrow(X)){

      #if assigned to a cluster remove statistics
      if( allocstruc[i, t - 1] == 1){
        idx <- alloc[i, t - 1] #temporary variable
        temp <- mk[idx,] %*% t(mk[idx,])
        #sk[idx, , ] <- sk[idx, , ] - ((lambdak[idx]) / lambdak[idx] - 1) * (X[i, ] - mk[idx,]) %*% t((X[i, ] - mk[idx,]))
        mk[idx, ] <- (lambdak[idx] * mk[idx, ] - X[i, ]) / (lambdak[idx] - 1)
        lambdak[idx] <- lambdak[idx] - 1
        nuk[idx] <- nuk[idx] - 1
        nk[idx] <- nk[idx] - 1
        tau1 <- tau1 - 1
        sk[idx,,] <- sk[idx,,] - (X[i,] %*% t(X[i,])) + (lambdak[idx] + 1) * temp  - lambdak[idx] * mk[idx,] %*% t(mk[idx,])
      }else{
        if(t > 2){
        tau2 <- tau2 - 1
        }
      }
      
      #compute probability of belonging to each organelle
      #precompute terms for speed
      weight <- (nk + betak)/(sum(nk) + sum(betak) - 1) #weight
      sigmak <- ((1 + lambdak) * sk)/(lambdak * (nuk - D + 1)) #scale matrix
      degf <- nuk - D + 1 #degrees freedom
      for(j in 1:K){
        allocprob[i, j, t] <- log(weight[j]) + dmvtCpp(X[i, ,drop = FALSE], 
                                                       mu = mk[j, ], 
                                                       sigma = sigmak[j,,], 
                                                       df = degf[j], 
                                                       log = T,
                                                       ncores_ = 1, 
                                                       isChol_ = F)
      }
      
      
      #normalise with underflow correction
      c <-  max(allocprob[i, , t])
      allocprob[i, , t] <- exp(allocprob[i, , t] - c) / sum(exp(allocprob[i, , t] - c))
      
      #sample allocation
      alloc[i, t] <- sample(x = 1:K, size = 1, prob = allocprob[i, , t])
      
      #compute structure component allocation
      n <- nrow(object)
      idk <- alloc[i, t] #temporary variable
      allocstrucprob[i, t, 1] <- log((tau1 + v)/(n + u + v - 1)) + dmvtCpp(X[i, ,drop=FALSE], 
                                                                      mu = mk[idk, ], 
                                                                      sigma = sigmak[idk,,], 
                                                                      df = degf[idk], 
                                                                      log = T, ncores_ = 1, isChol_ = F)
      allocstrucprob[i, t, 2] <- log((tau2 + u)/(n + u + v - 1)) + dmvtCpp(X[i, ,drop = FALSE], mu = M, sigma = V, df = 4, log = T,
                                                                           ncores_ = 1, isChol_ = F)
      
      #normalise and sample
      allocstrucprob[i, t, ] <- exp(allocstrucprob[i, t, ])/sum(exp(allocstrucprob[i, t, ]))
      allocstruc[i, t] <- sample(x = c(1, 0), 1, prob = allocstrucprob[i, t, ]) # reversed sample so 2nd entry is prob of 0
      
      #allocation statistics
      if( allocstruc[i, t] == 1){
        idx <- alloc[i, t] #temporary variable
        temp <- mk[idx,] %*% t(mk[idx,])
        #sk[idx, , ] <- sk[idx, , ] + ((lambdak[idx] ) / lambdak[idx] + 1 ) * (X[i, ] - mk[idx,]) %*% t((X[i, ] - mk[idx,]))
        mk[idx, ] <- (lambdak[idx] * mk[idx, ] + X[i, ]) / (lambdak[idx] + 1)
        lambdak[idx] <- lambdak[idx] + 1
        nuk[idx] <- nuk[idx] + 1
        nk[idx] <- nk[idx] + 1
        tau1 <- tau1 + 1
        if( t == 2){
        tau2 <- tau2 - 1  
        } 
        sk[idx,,] <- sk[idx,,] + (X[i,] %*% t(X[i,])) + (lambdak[idx] - 1) * temp  - lambdak[idx] * mk[idx,] %*% t(mk[idx,])
      }else{
        if(t > 2){ # default on first round in phi = 0
          tau2 <- tau2 + 1
        }
      } 
    }
  }
 
  return(list(allocstrucprob = allocstrucprob[, retained, ], 
              allocstruc = allocstruc[, retained], 
              alloc = alloc[, retained], 
              allocprob = allocprob[, , retained], mk = mk, lambdak = lambdak, nuk = nuk, sk = sk))
  
}
