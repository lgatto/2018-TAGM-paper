

Orgparamsmetropolis <- vector(mode = "list", length(getMarkerClasses(hlnorm)))

for(j in seq.int(length(getMarkerClasses(hlnorm)))){
  
  D <- 20
  
  exprs <- as.vector(t(exprs(hlnorm[fData(hlnorm)$markers == getMarkerClasses(hlnorm)[j],idx])))
  
metropolisRes <- metropolisGP(inith = getallpars[,j], tau = seq.int(20), X = exprs, nk = length(exprs)/D, D = D, niter = 20000)

Orgparamsmetropolis[[j]] <- metropolisRes$h[, 15000 + 20*seq(1:250)]

}

save(Orgparamsmetropolis, file = "metropolisparams.rda")

methypers <- sapply(Orgparamsmetropolis, rowMeans)
colnames(methypers) <- getMarkerClasses(hlnorm)
rownames(methypers) <- c("Memory", "Amplitude", "Noise")

methyperconf <- array(0, c(14,3,2))

for(i in 1:14){
  for(j in 1:3){
    methyperconf[i,j,] <- quantile(Orgparamsmetropolis[[i]][j,], probs = c(0.025, 0.975))
  }
}

par(mfrow = c(4,4))
for(j in 1:14){
  Orgdata <- t(exprs(hlnorm[fData(hlnorm)$markers == getMarkerClasses(hlnorm)[j],idx]))
  matplot(Orgdata, col = getStockcol()[j],pch = 19, type = "b", lty = 1, lwd = 1.5, main = paste(getMarkerClasses(hlnorm)[j]),
          xlab = "Fraction", ylab = "", cex.main = 2)
  
  nk <- table(fData(hlnorm)$markers)[j]
  S <- matrix(rep(1:length(tau),length(tau)), nrow = length(tau))
  params <- methypers[,j]
  sigmak <- exp(2 * params[3])
  a <- exp(2 * params[2])
  l <- exp(params[1])
  
  
  covA <- a * exp( - (S - t(S))^ 2 / l)
  R <- diag(1,D) + (nk * covA)/sigmak;
  trenchres <- trenchDetcpp(R[1,])
  Z <- trenchInvcpp(trenchres$v)
  invcov <- diag(1, nk*D)/sigmak - kronecker(matrix(1, nk, nk), Z %*% covA)/sigmak^2
  Kstar <- a*exp(-(matrix(rep(tau, nk*D),  nrow = D, byrow = F) - matrix(rep(tau, nk*D),nrow = D, byrow = T))^2/l)
  Kstarstar <- rep(a+sigmak, length(tau))
  M <- Kstar %*% invcov %*% as.vector(Orgdata)
  V <- sqrt(diag(diag(Kstarstar, length(tau)) - Kstar %*% invcov %*% t(Kstar)))
  points(seq_along(tau), M, col = "black", pch = 19, cex = 1.3, type = "b", lwd = 5, lty = 1)
  arrows(seq_along(tau), M-1.96*V, seq_along(tau), M+1.96*V, length=0.1, angle=90, code=3, col = "black", lwd = 3)
}

Orgparamsmetropolis <- vector(mode = "list", length(getMarkerClasses(hlnorm)))

for(j in seq.int(length(getMarkerClasses(hlnorm)))){
  
  D <- 20
  
  exprs <- as.vector(t(exprs(hlnorm[fData(hlnorm)$markers == getMarkerClasses(hlnorm)[j],idx])))
  
  metropolisRes <- metropolisGP(inith = getallpars[,j], tau = seq.int(20), X = exprs, nk = length(exprs)/D, D = D, niter = 20000)
  
  Orgparamsmetropolis[[j]] <- metropolisRes$h[, 15000 + 20*seq(1:250)]
  
}

require(mvtnorm)
hlnorm <- normalise(hyperLOPIT2015, method = "quantiles")
hlnorm <- normalise(hyperLOPIT2015, method = "center.mean")
idx <- c(1,11,2,12,3,13,4,14,5,15,6,16,7,17,8,18,9,19,10,20)
exprs <- as.vector(t(exprs(hlnorm[fData(hlnorm)$markers == getMarkerClasses(hlnorm)[3],idx])))
nk <- length(exprs)/20

system.time(hamiltonianRes <- hamiltonianGP(c(0,-2,-4), tau = seq.int(20), X = exprs, nk = nk, D = 20, niter = 500, stepsize = c(0.13,0.13), steps = 10))

system.time(MetroplisRes <- metropolisGP(c(0,-2,-4), tau = seq.int(20), X = exprs, nk = nk, D = 20, niter = 50000))

aY <- exprs(unknownMSnSet(hlnorm))



plot(unique(hamiltonianRes$h[3,-(1:20)]), type = "b", col = "blue", cex = 0.1)
effectiveSize(mcmc(unique(MetroplisRes$h[1,])))
effectiveSize(mcmc(unique(MetroplisRes$h[2,])))
effectiveSize(mcmc(unique(MetroplisRes$h[3,])))
effectiveSize(mcmc(unique(hamiltonianRes$h[1,])))
effectiveSize(mcmc(unique(hamiltonianRes$h[2,])))
effectiveSize(mcmc(unique(hamiltonianRes$h[3,])))
