




X <- exprs(hlnorm[fData(hlnorm)$markers == getMarkerClasses(hlnorm)[1],])
idx <- c(1,11,2,12,3,13,4,14,5,15,6,16,7,17,8,18,9,19,10,20)
X <- X[,idx]
tau <- seq_along(idx)
nk <- length(X)/max(idx)
D <- 20
matplot(t(X), pch = 19, col = "black", ylim = c(-0.2, 0.5))

ar <- 0 
mh <- matrix(0, 10001, 3)
mh[1, ] <- c(-2,2, -4)

for(i in 1:10000){

  metropolisRes <- metropolisGP(mh[i,], X, tau, nk, D, ar)
  mh[i+1, ] <- c(metropolisRes$h[1:2],-4)
  ar <- metroplisRes$ar
  
}

S <- matrix(rep(1:length(tau),length(tau)), nrow = length(tau))
sigmak <- exp(2 * mean(h[5000:10001,3]))
a <- exp(mean(h[5000:10001,2]))
l <- exp(mean(h[5000:10001,1]))


covA <- a * exp( - (S - t(S))^ 2 / l)
cov <- kronecker(matrix(1,nk, nk), covA) + sigmak*diag(1, nk*D)
Kstar <- a*exp(-(matrix(rep(tau, nk*D),  nrow = D, byrow = F) - matrix(rep(tau, nk*D),nrow = D, byrow = T))^2/l)
Kstarstar <- rep(a+sigmak, length(tau))


#Y <- matrix(X, nrow = 20 , byrow = TRUE)
#matplot(t(Y), pch = 19, col = "black", ylim = c(min(X)-1,max(X)+1))

M <- Kstar %*% solve(cov) %*% as.vector(t(X))
V <- sqrt(diag(diag(Kstarstar, length(tau)) - Kstar %*% solve(cov) %*% t(Kstar)))
points(seq_along(tau), M, col = "red", pch = 19, cex = 1.2, type = "b")
arrows(seq_along(tau), M-1.96*V, seq_along(tau), M+2*V, length=0.1, angle=90, code=3, col = "red", lwd = 2)

gridx <- seq(-7, 3.2, 0.05)
gridy <- seq(-8,3, 0.05)
likelisurf <- matrix(0, length(gridx), length(gridy))

for(i in 1:length(gridx)){
  for(j in 1:length(gridy)){
   likelisurf[i,j] <- likelihoodGPcpp(X, tau, h = c(gridx[i], gridy[j], -4), nk, D) - dmvnorm(c(gridx[i],gridy[j], -4), mean = rep(0,3), sigma = diag(1,3), log = TRUE)
  }
}

likelisurf[is.na(likelisurf)] <- 0.1
cols = rev(colorRampPalette(c('yellow','green','lightblue','darkblue'))(40))
filled.contour(x = gridx, y = gridy, z = log(likelisurf - min(likelisurf)), col = cols, nlevels = 10, frame.plot = FALSE)
points(unique(h)[,1],unique(h)[,2], col = "red", pch = 19, cex =1)
lines(unique(h)[,1], unique(h)[,2], col = "red")

require(graphics)
require(grDevices)



for(i in 1:35){
png(filename = paste0(i, c("hmc.png")))
filled.contour(x = gridx, y = gridy, z = log(likelisurf - min(likelisurf)), col = cols, nlevels = 10, frame.plot = FALSE)
points(unique(h)[1:i,1],unique(h)[1:i,2], col = "red", pch = 19, cex =1)
lines(unique(h)[1:i,1], unique(h)[1:i,2], col = "red")
dev.off()
}

for(i in 1:300){
  png(filename = paste0(i, c("mhmcmc.png")))
  filled.contour(x = gridx, y = gridy, z = log(likelisurf - min(likelisurf)), col = cols, nlevels = 10, frame.plot = FALSE)
  points(h[30*seq(1:300)[1:i],1],h[30*seq(1:300)[1:i],2], col = "red", pch = 19)
  lines(h[30*seq(1:300)[1:i],1],h[30*seq(1:300)[1:i],2], col = "red")
  dev.off()
}


ar <- 0 
h <- matrix(0, 301, 3)
h[1, ] <- c(-2,2, -4)

for(i in 1:300){
  
  if(i == 1){
    hamiltonianRes <- hamiltonianGP(h[i,], X, tau, nk, D, ar)
   h[i+1, ] <- c(hamiltonianRes$h[1:2],-4)
  ar <- hamiltonianRes$ar
  }else{
    hamiltonianRes <- hamiltonianGP(h[i,], X, tau, nk, D, ar, hamiltonianRes$residualp)
    h[i+1, ] <- c(hamiltonianRes$h[1:2],-4)
    ar <- hamiltonianRes$ar  
  }

}

ar <- 0 
h <- matrix(0, 301, 3)
h[1, ] <- c(0.5,-4, -6)

for(i in 1:300){
  
  if(i == 1){
  hamiltonianRes <- hamiltonianGP(h[i,], X, tau, nk, D, ar)
  h[i+1, ] <- c(0.5,hamiltonianRes$h[2:3])
  ar <- hamiltonianRes$ar 
  }else{
    hamiltonianRes <- hamiltonianGP(h[i,], X, tau, nk, D, ar, hamiltonianRes$residualp)
    h[i+1, ] <- c(0.5,hamiltonianRes$h[2:3])
    ar <- hamiltonianRes$ar 
  }
}

gridx <- seq(-7, 3.2, 0.05)
gridy <- seq(-8,3, 0.05)
likelisurffixedl <- matrix(0, length(gridx), length(gridy))

for(i in 1:length(gridx)){
  for(j in 1:length(gridy)){
    likelisurffixedl[i,j] <- likelihoodGPcpp(X, tau, h = c(0.5, gridx[i], gridy[j]), nk, D) - dmvnorm(c(0.5, gridx[i], gridy[j]), mean = rep(0,3), sigma = diag(1,3), log = TRUE)
  }
}

likelisurffixedl[is.na(likelisurf)] <- 0.1
cols = rev(colorRampPalette(c('yellow','green','lightblue','darkblue'))(40))
filled.contour(x = gridx, y = gridy, z = log(likelisurffixedl - min(likelisurffixedl)), col = cols, nlevels = 20, frame.plot = FALSE)
points(h[,2], h[,3], col = "red", pch = 19)
lines(h[,2], h[,3], col = "red")


gridx <- seq(-7, 3.2, 0.05)
gridy <- seq(-8,3, 0.05)
likelisurffixeda <- matrix(0, length(gridx), length(gridy))

for(i in 1:length(gridx)){
  for(j in 1:length(gridy)){
    likelisurffixeda[i,j] <- likelihoodGPcpp(X, tau, h = c(gridx[i], -2.5 , gridy[j]), nk, D) - dmvnorm(c(0.5, gridx[i], gridy[j]), mean = rep(0,3), sigma = diag(1,3), log = TRUE)
  }
}
likelisurffixeda[is.na(likelisurf)] <- 0.1
cols = rev(colorRampPalette(c('yellow','green','lightblue','darkblue'))(40))
filled.contour(x = gridx, y = gridy, z = log(likelisurffixeda - min(likelisurffixeda)), col = cols, nlevels = 20, frame.plot = FALSE)

ar <- 0 
ha <- matrix(0, 301, 3)
ha[1, ] <- c(-1,-2.5, -4)

for(i in 1:300){
  
  if(i == 1){
    hamiltonianRes <- hamiltonianGP(ha[i,], X, tau, nk, D, ar)
    ha[i+1, ] <- c(hamiltonianRes$h[1], -2.5 ,hamiltonianRes$h[3])
    ar <- hamiltonianRes$ar 
  }else{
    hamiltonianRes <- hamiltonianGP(ha[i,], X, tau, nk, D, ar, hamiltonianRes$residualp)
    ha[i+1, ] <- c(hamiltonianRes$h[1], -2.5, hamiltonianRes$h[3])
    ar <- hamiltonianRes$ar 
  }
}
filled.contour(x = gridx, y = gridy, z = log(likelisurffixeda - min(likelisurffixeda)), 
               col = cols, nlevels = 20, frame.plot = FALSE, plot.axes = {points(ha[,1], ha[,3], col = "red", pch = 19)})




ar <- 0 
hl <- matrix(0, 301, 3)
hl[1, ] <- c(0.5,-4, -6)

for(i in 1:300){
  
  if(i == 1){
    hamiltonianRes <- hamiltonianGP(hl[i,], X, tau, nk, D, ar)
    hl[i+1, ] <- c(0.5,hamiltonianRes$h[2:3])
    ar <- hamiltonianRes$ar 
  }else{
    hamiltonianRes <- hamiltonianGP(hl[i,], X, tau, nk, D, ar, hamiltonianRes$residualp)
    hl[i+1, ] <- c(0.5,hamiltonianRes$h[2:3])
    ar <- hamiltonianRes$ar 
  }
}
filled.contour(x = gridx, y = gridy, z = log(likelisurffixedl - min(likelisurffixedl)),
               col = cols, nlevels = 20, frame.plot = FALSE, plot.axes = {points(hl[,2], hl[,3], col = "red", pch = 19)})




ar <- 0 
hs <- matrix(0, 301, 3)
hs[1, ] <- c(-2,2, bestparams$par[3])

for(i in 1:300){
  
  if(i == 1){
    hamiltonianRes <- hamiltonianGP(hs[i,], X, tau, nk, D, ar)
    hs[i+1, ] <- c(hamiltonianRes$h[1:2], bestparams$par[3])
    ar <- hamiltonianRes$ar
  }else{
    hamiltonianRes <- hamiltonianGP(hs[i,], X, tau, nk, D, ar, hamiltonianRes$residualp)
    hs[i+1, ] <- c(hamiltonianRes$h[1:2], bestparams$par[3])
    ar <- hamiltonianRes$ar  
  }
  
}

filled.contour(x = gridx, y = gridy, z = log(likelisurf - min(likelisurf)),
               col = cols, nlevels = 10, frame.plot = FALSE, plot.axes = {points(hs[,1], hs[,2], col = "red", pch = 19)})


ar <- 0 
hpar <- matrix(0, 301, 3)
hpar[1, ] <- bestparams$par

hamiltonianRes <- hamiltonianGP(c(1,-2, -0.6),tau = tau, X = Xk, nk = 20, D = 5, niter = 1000, stepsize = c(0.03,0.04), steps = 35)

metropolisRes <- metropolisGP(inith = c(1,-2, -0.6), X = Xk, tau = tau,nk = 20,D = 5,niter = 20000)


res <- apply(initz, 1,function(z){lbfgs(likelihoodGP,
                                        gradientGP,
                                        vars = z,
                                        invisible=1,
                                        epsilon = 1e-8,
                                        Xk = as.vector(t(X)),
                                        tau = seq.int(20),
                                        nk = 43,
                                        D = 20)})

bestparams <- res[[which.min(lapply(res, function(x){max(x$value)}))]]

hamiltonianRes <- hamiltonianGP(bestparams$par, tau = seq.int(20), X = as.vector(t(X)), nk = 43, D = 20, niter = 500, stepsize = c(0.01,0.02), steps = 20)

metropolisRes <- metropolisGP(inith = bestparams$par, tau = seq.int(20), X = as.vector(t(X)), nk = 43, D = 20, niter = 20000)


hlnorm <- normalise(hyperLOPIT2015, method = "quantiles")
hlnorm <- normalise(hyperLOPIT2015, method = "center.mean")
X <- exprs(hlnorm[fData(hlnorm)$markers == getMarkerClasses(hlnorm)[4],])
idx <- c(1,11,2,12,3,13,4,14,5,15,6,16,7,17,8,18,9,19,10,20)
X <- X[,idx]
matplot(t(X), pch = 19, col = "black", ylim = c(-0.2, 0.5))

S <- matrix(rep(1:length(tau),length(tau)), nrow = length(tau))
sigmak <- exp(2 * hpar[3])
a <- exp(2 * hpar[2])
l <- exp(hpar[1])

covA <- a * exp( - (S - t(S))^ 2 / l)
cov <- kronecker(matrix(1,nk, nk), covA) + sigmak*diag(1, nk*D)
Kstar <- a*exp(-(matrix(rep(tau, nk*D),  nrow = D, byrow = F) - matrix(rep(tau, nk*D),nrow = D, byrow = T))^2/l)
Kstarstar <- rep(a+sigmak, length(tau))


#Y <- matrix(X, nrow = 20 , byrow = TRUE)
#matplot(t(Y), pch = 19, col = "black", ylim = c(min(X)-1,max(X)+1))

M <- Kstar %*% solve(cov) %*% as.vector(t(X))
V <- sqrt(diag(diag(Kstarstar, length(tau)) - Kstar %*% solve(cov) %*% t(Kstar)))
points(seq_along(tau), M, col = "red", pch = 19, cex = 1.2, type = "b")
arrows(seq_along(tau), M-1.96*V, seq_along(tau), M+2*V, length=0.1, angle=90, code=3, col = "red", lwd = 2)

S <- matrix(rep(1:length(tau),length(tau)), nrow = length(tau))
hampar <- rowMeans(hamiltonianRes$h)
sigmak <- exp(2 * hampar[3])
a <- exp(2 * hampar[2])
l <- exp(hampar[1])

metpar <- rowMeans(metropolisRes$h)

gridx <- seq(-2, 2, 0.05)
gridy <- seq(-3,3, 0.05)
likelisurf <- matrix(0, length(gridx), length(gridy))

for(i in 1:length(gridx)){
  for(j in 1:length(gridy)){
    likelisurf[i,j] <- likelihoodGPcpp(as.vector(t(X)), tau, h = c(gridx[i], gridy[j], -3.6), nk, D) - dmvnorm(c(gridx[i],gridy[j], -3.6), mean = rep(0,3), sigma = diag(1,3), log = TRUE)
  }
}

gridx <- seq(-4, 0, 0.05)
gridy <- seq(-5,1, 0.05)
likelisurffixedl <- matrix(0, length(gridx), length(gridy))
for(i in 1:length(gridx)){
  for(j in 1:length(gridy)){
    likelisurffixedl[i,j] <- likelihoodGPcpp(as.vector(t(X)), tau, h = c(0.7, gridx[i], gridy[j]), nk, D) - dmvnorm(c(0.5, gridx[i], gridy[j]), mean = rep(0,3), sigma = diag(1,3), log = TRUE)
  }
}

gridx <- seq(-2, 5, 0.05)
gridy <- seq(-5,-1, 0.05)
likelisurffixeda <- matrix(0, length(gridx), length(gridy))
for(i in 1:length(gridx)){
  for(j in 1:length(gridy)){
    likelisurffixeda[i,j] <- likelihoodGPcpp(as.vector(t(X)), tau, h = c(gridx[i], -2.2 , gridy[j]), nk, D) - dmvnorm(c(0.5, gridx[i], gridy[j]), mean = rep(0,3), sigma = diag(1,3), log = TRUE)
  }
}


cols = rev(colorRampPalette(c('yellow','green','lightblue','darkblue'))(40))
filled.contour(x = gridx, y = gridy, z = log(likelisurf - min(likelisurf)),
               col = cols, nlevels = 20, frame.plot = FALSE,
plot.axes = {points(metropolisRes$h[1,-(1:18000)],metropolisRes$h[2,-(1:18000)] , col = "red", pch = 19)})
filled.contour(x = gridx, y = gridy, z = log(likelisurf - min(likelisurf)),
               col = cols, nlevels = 20, frame.plot = FALSE,
               plot.axes = {points(unique(hamiltonianRes$h[1,]),unique(hamiltonianRes$h[2,]) , col = "red", pch = 19)})
filled.contour(x = gridx, y = gridy, z = log( likelisurffixedl - min(likelisurffixedl)),
               col = cols, nlevels = 10, frame.plot = FALSE,
               plot.axes = {points(metropolisRes$h[2,-(1:18000)],metropolisRes$h[3,-(1:18000)] , col = "red", pch = 19)})
filled.contour(x = gridx, y = gridy, z = log( likelisurffixedl - min( likelisurffixedl)),
               col = cols, nlevels = 10, frame.plot = FALSE,
               plot.axes = {points(unique(hamiltonianRes$h[2,]),unique(hamiltonianRes$h[3,]) , col = "red", pch = 19)})




filled.contour(x = gridx, y = gridy, z = log(likelisurffixeda - min(likelisurffixeda)),
               col = cols, nlevels = 10, frame.plot = FALSE,
               plot.axes = {points(metropolisRes$h[1,-(1:18000)],metropolisRes$h[3,-(1:18000)] , col = "red", pch = 19)})
filled.contour(x = gridx, y = gridy, z = log(likelisurffixeda - min(likelisurffixeda)),
               col = cols, nlevels = 10, frame.plot = FALSE,
               plot.axes = {points(unique(hamiltonianRes$h[1,]),unique(hamiltonianRes$h[3,]) , col = "red", pch = 19, type = "b")})



