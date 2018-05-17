hlnorm <- normalise(hyperLOPIT2015[,1:10], method = "quantiles")
hlnorm <- normalise(hyperLOPIT2015[,1:10], method = "center.mean")

init <- seq(-3,3, 2)
initz <- matrix(0,length(init),3)

for(i in seq_along(init)){
  initz[i,] <- init[sample.int(length(init), size = 3, replace = T)]
}

idx <- seq.int(10)
tau <- seq.int(10)
nk <- length(X)/max(idx)
D <- 10


Orgparams <- vector(mode = "list", length(getMarkerClasses(hlnorm)))

system.time(for(j in seq.int(length(getMarkerClasses(hlnorm)))){
  exprs <- as.vector(t(exprs(hlnorm[fData(hlnorm)$markers == getMarkerClasses(hlnorm)[j],idx])))
  
  res <- apply(initz, 1,function(z){lbfgs(likelihoodGP,
                                          gradientGP,
                                          vars = z,
                                          invisible=1,
                                          epsilon = 1e-8,
                                          Xk = exprs,
                                          tau =  seq.int(D),
                                          nk = length(exprs)/D,
                                          D = 10)})
  Orgparams[[j]] <- res[[which.min(lapply(res, function(x){max(x$value)}))]]
  
})


#save(Orgparams, file = "lbfgsOrg.rda")

getallpars <- sapply(Orgparams, function(x) {x$par})
colnames(getallpars) <- getMarkerClasses(hlnorm)
rownames(getallpars) <- c("Memory", "Amplitude", "Noise")

par(mfrow = c(4,4))
for(j in 1:14){
  Orgdata <- t(exprs(hlnorm[fData(hlnorm)$markers == getMarkerClasses(hlnorm)[j],idx]))
  matplot(Orgdata, col = getStockcol()[j],pch = 19, type = "b", lty = 1, lwd = 1.5, main = paste(getMarkerClasses(hlnorm)[j]),
          xlab = "Fraction", ylab = "", cex.main = 2, ylim = c(min(Orgdata) - 0.05, max(Orgdata) + 0.05))
  
  nk <- table(fData(hlnorm)$markers)[j]
  S <- matrix(rep(1:length(tau),length(tau)), nrow = length(tau))
  params <- Orgparams[[j]]$par
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


andyGP <- mixGP(mydata = hyperLOPIT2015[,1:10], hypLearn = "HMC", numIter = 20000, hypIter = 50)
