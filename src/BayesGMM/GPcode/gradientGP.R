
gradientGP <- function(Xk, tau, h, nk, D){
  
  grad <- gradientGPcpp(Xk, tau, h, nk, D)
  
  return(unlist(grad))
  
}