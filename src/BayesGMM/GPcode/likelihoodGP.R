

likelihoodGP <- function(Xk, tau, h, nk, D){
  
  negloglike <- likelihoodGPcpp(Xk, tau, h, nk, D)
  
  return(negloglike)
}