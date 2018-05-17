# Metroplis-Hastings for GP hyperparameters with normal prior


metropolisGP <- function(inith,
                         X,
                         tau,
                         nk,
                         D,
                         niter
                         ){
  
  
  h <- matrix(0, ncol = niter + 1, nrow = length(inith))
  h[, 1] <- inith
  ar <- 0 
  
  for(i in seq.int(niter)){
   xi <- rnorm(length(inith), mean = 0, sd = 0.1) # sample random walk steps
   oldHypers <- h[, i]
   proposedHypers <- h[, i] + xi #random walk proposal
   #compute metropolis ratio, likelihoodGPcpp return negative logliklihood
   mhratio <- -likelihoodGPcpp(X, tau, proposedHypers, nk, D) + likelihoodGPcpp(X, tau, oldHypers, nk, D) + 
     sum(dnorm(proposedHypers, mean = 0, sd = 1, log = TRUE)) - sum(dnorm(oldHypers, mean = 0, sd = 1, log = TRUE))
  
   if(mhratio > log(runif(1, 0, 1))){
    h[, i + 1] <- proposedHypers
    ar <- ar + 1
   }else{
    h[, i + 1] <- oldHypers
   }
   
  }
  ar <- ar/niter
  
  return(list(h = h, ar = ar))
  
}


