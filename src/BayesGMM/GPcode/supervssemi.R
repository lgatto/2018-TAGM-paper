
#Tan LBFGS
#random grid sampling for starting values
initialvalues <- seq(-3,3, 2)
init <- matrix(0, length(initialvalues), 3)
for(i in seq_along(initialvalues)){
  init[i,] <- initialvalues[sample.int(length(initialvalues), size = 3, replace = T)]
}
tanparams <- vector(mode = "list", length(getMarkerClasses(tan2009r1)))
for(j in seq.int(length(getMarkerClasses(tan2009r1)))){
  exprs <- as.vector(t(exprs(tan2009r1[fData(tan2009r1)$markers == getMarkerClasses(tan2009r1)[j],])))
  D <- ncol(tan2009r1)
  
  res <- apply(init, 1,function(z){lbfgs(likelihoodGP,
                                          gradientGP,
                                          vars = z,
                                          invisible=1,
                                          epsilon = 1e-8,
                                          Xk = exprs,
                                          tau =  seq.int(D),
                                          nk = length(exprs)/D,
                                          D = D)})
  
  tanparams[[j]] <- res[[which.min(lapply(res, function(x){max(x$value)}))]]
  
}

getallpars <- sapply(tanparams, function(x) {x$par})


tanparamsmetropolis <- vector(mode = "list", length(getMarkerClasses(tan2009r1)))

for(j in seq.int(length(getMarkerClasses(tan2009r1)))){
  
  
  
  exprs <- as.vector(t(exprs(tan2009r1[fData(tan2009r1)$markers == getMarkerClasses(tan2009r1)[j],])))
  
  metropolisRes <- metropolisGP(inith = getallpars[,j], tau = seq.int(D), X = exprs, nk = length(exprs)/D, D = D, niter = 20000)
  
  tanparamsmetropolis[[j]] <- metropolisRes$h[,10*seq(1:2000)]
  
}
tanGP <- mixGP(tan2009r1, hypLearn = "MH", numIter = 20000, hypIter = 5)

par(mfrow= c(3,4))
for(i in 1:11){
hist(tanparamsmetropolis[[i]][3,], breaks = 25, xlim = c(-4.5,-1.5), ylim = c(0, 450), col = rgb(0,0,1,0.5),
     main = paste0(getMarkerClasses(tan2009r1)[i]), xlab = "log-noise")
abline(v=tanparams[[i]]$par[3], col="red", lwd = 2)


hist(sapply(tanGP$hypers, function(x)x[i,3])[2500 + 5 *seq(1:3500)], add = T, col = rgb(1,0,0,0.5), breaks = 25)


#yfit <- rnorm(2000,mean=0,sd=1) 
#hist(yfit, col="red", lwd=2, add = T, breaks = 30)
}
plot.new()
legend("top", legend = c("supervised", "semi-supervised"), fill = c(rgb(0,0,1,0.5), rgb(1,0,0,0.5)), cex = 1)
legend("bottom", legend = c("L-BFGS"), col= "red", lty = 1 )



tan <- tan2009r1
plotDist(object[fData(object)$predicted.allocation == "ER", ],
         pcol = "grey", xlab= "", lwd = 1, type = "l") 
title("Labelled data for ER along with proteins predicted to the ER")
matlines(t(exprs(object[fData(object)$markers == "ER" , ])),lty = 1, col = getStockcol()[1], lwd = 1, type="l" )


legend("topright", c("ER"),
       lty = 1, col = c(getStockcol()[1]), bty = "n")

tau <- seq.int(1:D)
nk <- sum(fData(object)$predicted.allocation == "ER")
S <- matrix(rep(1:length(tau),length(tau)), nrow = length(tau))
params <- rowMeans(sapply(tanGP$hypers, function(x)x[1,]))
sigmak <- exp(2 * params[3])
a <- exp(2 * params[2])
l <- exp(params[1])


covA <- a * exp( - (S - t(S))^ 2 / l)
R <- diag(1,D) + (nk * covA)/sigmak;
trenchres <- trenchDetcpp(R[1,])
Z <- trenchInvcpp(trenchres$v)
invcov <- diag(1, nk*D)/sigmak - kronecker(matrix(1, nk, nk), Z %*% covA)/sigmak^2
Kstar <- a*exp(-(matrix(rep(tau, nk*D),  nrow = D, byrow = F) - matrix(rep(tau, nk*D),nrow = D, byrow = T))^2/l)
Kstarstar <- rep(a+sigmak, length(D))
M <- Kstar %*% invcov %*% as.vector(t(exprs(object[fData(object)$predicted.allocation == "ER", ])))
V <- sqrt(diag(diag(Kstarstar, length(tau)) - Kstar %*% invcov %*% t(Kstar)))
points(seq_along(tau), M, col = "black", pch = 19, cex = 1.3, type = "b", lwd = 5, lty = 1)
arrows(seq_along(tau), M-1.96*V, seq_along(tau), M+1.96*V, length=0.1, angle=90, code=3, col = "black", lwd = 3)


plotDist(object[fData(object)$predicted.allocation == "Nucleus", ],
         pcol = "grey", xlab= "", lwd = 1, type = "l") 
title("Labelled data for Nucleus along with proteins predicted to the Nucleus")
matlines(t(exprs(object[fData(object)$markers == "Nucleus" , ])),lty = 1, col = getStockcol()[6], lwd = 1, type="l" )


legend("topright", c("Nucleus"),
       lty = 1, col = c(getStockcol()[6]), bty = "n")

tau <- seq.int(1:D)
nk <- sum(fData(object)$predicted.allocation == "Nucleus")
S <- matrix(rep(1:length(tau),length(tau)), nrow = length(tau))
params <- rowMeans(sapply(tanGP$hypers, function(x)x[1,]))
sigmak <- exp(2 * params[3])
a <- exp(2 * params[2])
l <- exp(params[1])


covA <- a * exp( - (S - t(S))^ 2 / l)
R <- diag(1,D) + (nk * covA)/sigmak;
trenchres <- trenchDetcpp(R[1,])
Z <- trenchInvcpp(trenchres$v)
invcov <- diag(1, nk*D)/sigmak - kronecker(matrix(1, nk, nk), Z %*% covA)/sigmak^2
Kstar <- a*exp(-(matrix(rep(tau, nk*D),  nrow = D, byrow = F) - matrix(rep(tau, nk*D),nrow = D, byrow = T))^2/l)
Kstarstar <- rep(a+sigmak, length(D))
M <- Kstar %*% invcov %*% as.vector(t(exprs(object[fData(object)$predicted.allocation == "Nucleus", ])))
V <- sqrt(diag(diag(Kstarstar, length(tau)) - Kstar %*% invcov %*% t(Kstar)))
points(seq_along(tau), M, col = "black", pch = 19, cex = 1.3, type = "b", lwd = 5, lty = 1)
arrows(seq_along(tau), M-1.96*V, seq_along(tau), M+1.96*V, length=0.1, angle=90, code=3, col = "black", lwd = 3)


#assessing convergence Hamiltonian Monte-Carlo
tanGP1 <- mixGP(tan2009r1, hypLearn = "HMC", numIter = 20000, hypIter = 50, seed = 1)
tanGP2 <- mixGP(tan2009r1, hypLearn = "HMC", numIter = 20000, hypIter = 50, seed = 2)

mcmc1 <- mcmc(sapply(tanGP1$hypers, function(x)x[1,3])[50 *seq(1:400)])
mcmc2 <- mcmc(sapply(tanGP2$hypers, function(x)x[1,3])[50 *seq(1:400)])
mcmc1and2 <- mcmc.list(mcmc1,mcmc2)
res <- gelman.diag(mcmc1and2, autoburnin = F)

plot(sapply(tanGP1$hypers, function(x)x[1,3])[1000 + 50 *seq(1:400)], cex = 0.1, col = "blue", type = "b", xlab = "sample", ylab ="Noise", main = "assessing convergence using parallel chains")
points(sapply(tanGP2$hypers, function(x)x[1,3])[1000+ 50 *seq(1:400)], cex = 0.1, col = "red", type = "b")

plot(sapply(tanGP1$hypers, function(x)x[2,2])[200 + 50*seq(1:400)], cex = 0.1, col = "blue", type = "b", xlab = "sample", ylab ="Length-Scale", main = "assessing convergence using parallel chains")
points(sapply(tanGP2$hypers, function(x)x[2,2])[200 + 50* seq(1:400)], cex = 0.1, col = "red", type = "b")

mcmc1 <- mcmc(sapply(tanGP1$hypers, function(x)x[2,2])[50 *seq(1:400)])
mcmc2 <- mcmc(sapply(tanGP2$hypers, function(x)x[2,2])[50 *seq(1:400)])
mcmc1and2 <- mcmc.list(mcmc1,mcmc2)
res <- gelman.diag(mcmc1and2, autoburnin = F)

save(tanGP1, file = "tanGPHMC1.rda")
save(tanGP2, file = "tanGPHMC2.rda")

#assesing metropolis hastings convergence
tanGP3 <- mixGP(tan2009r1, hypLearn = "MH", numIter = 20000, hypIter = 10, seed = 1)
tanGP4 <- mixGP(tan2009r1, hypLearn = "MH", numIter = 20000, hypIter = 10, seed = 2)

mcmc1 <- mcmc(sapply(tanGP3$hypers, function(x)x[1,3])[5000 + 10 *seq(1:1500)])
mcmc2 <- mcmc(sapply(tanGP4$hypers, function(x)x[1,3])[5000 + 10 *seq(1:1500)])
mcmc1and2 <- mcmc.list(mcmc1,mcmc2)
res <- gelman.diag(mcmc1and2, autoburnin = F)

plot(sapply(tanGP3$hypers, function(x)x[1,3])[5000 + 10 *seq(1:1500)], cex = 0.1, col = "blue", type = "b", xlab = "sample", ylab ="Noise", main = "assessing convergence using parallel chains")
points(sapply(tanGP4$hypers, function(x)x[1,3])[5000 + 10 *seq(1:1500)], cex = 0.1, col = "red", type = "b")

plot(sapply(tanGP3$hypers, function(x)x[2,2])[5000 + 20 *seq(1:750)], cex = 0.1, col = "blue", type = "b", xlab = "sample", ylab ="Length-Scale", 
     main = "assessing convergence using parallel chains", ylim =c(-3,-0.5))
points(sapply(tanGP4$hypers, function(x)x[2,2])[5000 + 20 *seq(1:750)], cex = 0.1, col = "red", type = "b")

mcmc3 <- mcmc(sapply(tanGP3$hypers, function(x)x[2,2])[5000 + 20 *seq(1:750)])
mcmc4 <- mcmc(sapply(tanGP4$hypers, function(x)x[2,2])[5000 + 20 *seq(1:750)])
mcmc3and4 <- mcmc.list(mcmc3,mcmc4)
res <- gelman.diag(mcmc1and2, autoburnin = F)
