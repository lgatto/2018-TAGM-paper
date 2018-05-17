
knnquadloss <- function(knnparams, object){

knnquadloss <- vector("list", length = 100) #quadratic loss to save

for(iLoss in seq_along(knnquadloss)){
 
 markerhl <- markerMSnSet(object)

 index <- knnparams@testPartitions[[iLoss]]
 trueclasses <- fData(markerhl)$markers[index]

 levels(fData(markerhl)$markers) <- c(getMarkerClasses(object), "unknown") #change levels
 fData(markerhl)$markers[index] <- "unknown" #hide markers

 nbr <- knnparams@results[iLoss,2]
 train <- exprs(markerhl)[-index,]
 test <- exprs(markerhl)[index,]

 K <- length(getMarkerClasses(object))

 res <- FNN::knn(train, test, cl = fData(markerhl)$markers[-index], k = nbr, prob = TRUE) #knn classification
 nn.index <- attr(res, "nn.index")
 distMat <- attr(res, "nn.dist") #get atrributes

 scoresmatrix <- matrix(0, length(index), length(getMarkerClasses(object)))
 colnames(scoresmatrix) <- getMarkerClasses(object) #change names

 irates <- tabulate(fData(markerhl)$markers[-index])/sum(tabulate(fData(markerhl)$markers[-index]))


  for(i in seq_along(index)){
    scoresmatrix[i, ] <- (0.5 * irates * K)/(nbr + 0.5 * K)
    
    for(j in seq(1:K)[table(fData(markerhl)$markers[-index][nn.index[i,]])!=0]){
     count <- as.numeric(table(fData(markerhl)$markers[-index][nn.index[i,]])[j])
     scoresmatrix[i,j] <-  (count + 0.5 * irates[j] * K)/(nbr + 0.5 * K)
     
   }
  }

  scoresmatrix <- scoresmatrix/rowSums(scoresmatrix)
  
  allocmatrix <- matrix(0, length(index), length(unique(trueclasses)))
  
  #create allocation matrix
  for(j in seq_along(index)){
    allocmatrix[j, as.numeric(factor(trueclasses), seq(1,length(unique(trueclasses))))[j]] <- 1 
  }
  
  knnquadloss[[iLoss]] <- sum((scoresmatrix[,getMarkerClasses(markerhl)] - allocmatrix)^2) #compute quad loss
}

return(knnquadloss = knnquadloss)
}
