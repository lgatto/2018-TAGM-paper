#supervised mixture of gaussian process regression models


sprvsdmixGP <- function(mydata,
                        hypLearn = "LBFGS"){
  


hlnorm <- normalise(mydata, method = "quantiles")
hlnorm <- normalise(mydata, method = "center.mean")
D <- ncol(exprs(hlnorm))
alpha <- 1
numIter <- 2
u <- 2
v <- 50
M <- colMeans(exprs(hlnorm))
V <- cov(exprs(hlnorm))/2
K <- length(getMarkerClasses(hlnorm))

componenthypers <- vector(mode = "list", length(getMarkerClasses(hlnorm)))

#random grid sampling for starting values
initialvalues <- seq(-3,3, 2)
init <- matrix(0, length(initialvalues), 3)
for(i in seq_along(initialvalues)){
  init[i,] <- initialvalues[sample.int(length(initialvalues), size = 3, replace = T)]
}

idx <- seq.int(1:D)

#LBFGS routine to get hypers
for(j in seq.int(length(getMarkerClasses(hlnorm)))){
  
  exprs <- as.vector(t(exprs(hlnorm[fData(hlnorm)$markers == getMarkerClasses(hlnorm)[j],idx])))
  
  res <- apply(init, 1,function(z){lbfgs(likelihoodGP,
                                          gradientGP,
                                          vars = z,
                                          invisible=1,
                                          epsilon = 1e-8,
                                          Xk = exprs,
                                          tau =  seq.int(D),
                                          nk = length(exprs)/D,
                                          D = D)})
  componenthypers[[j]] <- res[[which.min(lapply(res, function(x){max(x$value)}))]]$par
  
}

#seperate data into known and unknown
unknownData <- unknownMSnSet(hlnorm)
exprsUnknown <- exprs(unknownData)
exprsKnown <- exprs(markerMSnSet(hlnorm))
allocKnown <- seq.int(1:K)[fData(markerMSnSet(hlnorm))$markers]
numProtein <- nrow(unknownData)

#some storage
alloc <- matrix(0, nrow = numProtein, ncol = numIter)
allocOut <- matrix(0, nrow = numProtein, ncol = numIter)
outlierprob <- matrix(0, nrow = numProtein, ncol = numIter)
allocOutprob <- array(0, c(nrow(unknownData), numIter, 2))
allocprob <- array(0, c(nrow(unknownData), numIter, K))
logMarginalNum <- vector(mode = "numeric", length = K)
logMarginalDen <- vector(mode = "numeric", length = K)

#random allocation of unknown Proteins
alloc[, 1] <- sample.int(n = K, size = numProtein, replace = TRUE)

#number of proteins allocated to each component
nk <- table(getMarkers(hlnorm, verbose = FALSE))[getMarkerClasses(hlnorm)]

#number of proteins allocated to each component with unknowns
nk <- nk

#precalculate component marginal likelihood here
for(j in seq.int(K)){
 componentsdata <- c(as.vector(t(exprsKnown[allocKnown==j, ])))  
 logMarginalNum[j] <- -likelihoodGPcpp(componentsdata, seq.int(D), componenthypers[[j]], nk[j], D) 
}
logMarginalDen <- logMarginalNum


#number initial allocated to outlier component
allocOut[, 1] <- 0
tau1 <- sum(allocOut[, 1] == 1) + nrow(exprsKnown)
tau2 <- sum(allocOut[, 1] == 0)

for(t in 2:numIter){

 alloctemp <- alloc[, t - 1]
 outlier <- allocOut[, t - 1]
  
 for(i in seq.int(numProtein)){
  
  priorPredictive <- (nk + alpha/K)/(sum(nk) + alpha - 1)
  #if outlier not remove from cluster
  if(outlier[i] != 0){
   priorPredictive[alloctemp[i]] <- (nk[alloctemp[i]] - 1 + alpha/K)/(sum(nk[alloctemp[i]]) + alpha - 1) 
  }
  
  Y <- exprsUnknown[-i, ] #remove ith protein from dataset 
  if(outlier[i] == 0 ){
    for(j in seq.int(K)){
     componentswithProtein <- c(as.vector(t(exprsKnown[allocKnown == j, ] )), 
                                 as.vector(t(Y[(alloctemp[-i] == j)&(outlier[-i] != 0),])), exprsUnknown[i, ])
     componentswithoutProtein <- c(as.vector(t(exprsKnown[allocKnown == j, ] )), 
                                    as.vector(t(Y[(alloctemp[-i] == j)&(outlier[-i] != 0),])))
      
     logMarginalNum[j] <- -likelihoodGPcpp(componentswithProtein, seq.int(D), componenthypers[[j]], nk[j] + 1, D)
     logMarginalDen[j] <- -likelihoodGPcpp(componentswithoutProtein, seq.int(D), componenthypers[[j]], nk[j], D)   
    }
  }else if (outlier[i]!=0){
    for(j in seq.int(K)){
      componentswithProtein <- c(as.vector(t(exprsKnown[allocKnown == j, ] )), 
                                 as.vector(t(Y[(alloctemp[-i] == j)&(outlier[-i] != 0),])), exprsUnknown[i, ])
      componentswithoutProtein <- c(as.vector(t(exprsKnown[allocKnown == j, ] )), 
                                    as.vector(t(Y[(alloctemp[-i] == j)&(outlier[-i] != 0),])))
      if(alloctemp[i] == j){
        logMarginalNum[j] <- -likelihoodGPcpp(componentswithProtein, seq.int(D), componenthypers[[j]], nk[j], D)
        logMarginalDen[j] <- -likelihoodGPcpp(componentswithoutProtein, seq.int(D), componenthypers[[j]], nk[j] - 1, D)
      }else{
        logMarginalNum[j] <- -likelihoodGPcpp(componentswithProtein, seq.int(D), componenthypers[[j]], nk[j] + 1, D)
        logMarginalDen[j] <- -likelihoodGPcpp(componentswithoutProtein, seq.int(D), componenthypers[[j]], nk[j], D)
      }
    }
  }
  logratioML <- logMarginalNum  - logMarginalDen
  print(logratioML)
  unnormliasedlogallocprob <- log(priorPredictive)  + logratioML
  unnormliasedlogallocprob <- unnormliasedlogallocprob - max(unnormliasedlogallocprob) #add consant for numerical stability
  allocprob[i, t, ] <- exp(unnormliasedlogallocprob)/sum(exp(unnormliasedlogallocprob))
  
  alloctemp[i] <- sample.int(n = K, size = 1, prob = allocprob[i, t, ])

  #compute allocation to outlier component, allocation to component already computed  
  n <- numProtein
  
  tau <- c(tau1, tau2) #temp variable for allocation
  if(allocOut[i, t - 1] == 0){
    tau[2] <- tau[2] - 1
  }else{
     tau[1] <- tau[1] - 1
  }
  priorOutlier <- (tau + c(v,u))/(sum(nk) +  u + v - 1)
  allocOutprob[i, t, 1] <- log(priorOutlier[1]) + logratioML[alloctemp[i]]
  allocOutprob[i, t, 2] <- log(priorOutlier[2]) + dmvtCpp(exprsUnknown[i,, drop =FALSE], mu = M, sigma = V, df = 4, log = T,
                                                                       ncores_ = 1, isChol_ = F)
  print(allocOutprob[i, t, ])
  #still need to check dmtCpp and normalise if needed
  allocOutprob[i, t, ] <- allocOutprob[i,t,] - max(allocOutprob[i,t, ])
  #normalise and sample
  allocOutprob[i, t, ] <- exp(allocOutprob[i, t, ])/sum(exp(allocOutprob[i, t, ]))
  outlier[i] <- sample(x = c(1, 0), 1, prob = allocOutprob[i, t, ]) # reversed sample so 2nd entry is prob of 0
 
  #reallocate
  if(outlier[i] - allocOut[i, t - 1] == 1 ){
   tau1 <- tau1 + 1
   tau2 <- tau2 - 1
  } else if (outlier[i] - allocOut[i, t - 1] == -1){
   tau1 <- tau1 - 1
   tau2 <- tau2 + 1
  }  
  
  #reallocated component allocations
  if((outlier[i]==0)&(allocOut[i, t - 1] == 0)){
    #no reallocations needed
  }else if((outlier[i]==0)&(allocOut[i, t - 1]==1)){
     nk[alloc[i, t - 1]] <- nk[alloc[i, t - 1]] - 1
  }else if((outlier[i] == 1)&(allocOut[i, t - 1]==1)){
    if(alloctemp[i] - alloc[i, t-1] != 0 ){
     nk[alloctemp[i]] <- nk[alloctemp[i]] + 1
     nk[alloc[i, t - 1]] <- nk[alloc[i, t - 1]] - 1
    }# if reallocated to same component do nothing
  }else if ((outlier[i] == 1)&(allocOut[i, t - 1]==0)){
     nk[alloctemp[i]] <- nk[alloctemp[i]] + 1
  }
  
 }
 
 #update hypers
  #if(hypLearn != "LBFGS"){
    #do something here
  #}
 
  alloc[,t] <- alloctemp
  allocOut[,t] <- outlier

}  

return(list(allocOutprob = allocOutprob, allocOut = allocOut, allocprob = allocprob, alloc = alloc))


}