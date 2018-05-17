
svmquadloss <- function(svmparams, object){
  
svmquadloss <- vector("list", length = 100) #quadratic loss to save

for(i in seq_along(svmparams@f1Matrices)){
  
  index <- svmparams@testPartitions[[i]] #gettestpartitions
  
  
  markerhl <- markerMSnSet(object) #create markerset
  trueclasses <- fData(markerhl)$markers[index] #save true markers
  
  levels(fData(markerhl)$markers) <- c(getMarkerClasses(object), "unknown") #change levels
  fData(markerhl)$markers[index] <- "unknown" #hide markers
  
  sigma <- as.numeric(svmparams@results[i, 2]) #get sigma
  cost <- as.numeric(svmparams@results[i, 3]) #getcost
  
  res <- svmClassification(markerhl, fcol= "markers", sigma = sigma, cost = cost, scores = "all" ) #svm classification
  scoresmatrix <- fData(res)[index, "svm.all.scores"] #get scores as matrix
  colnames(scoresmatrix) <- unique(trueclasses) #change names
  
  allocmatrix <- matrix(0, length(index), length(unique(trueclasses)))
  
  #create allocation matrix
  for(j in seq_along(index)){
    allocmatrix[j, as.numeric(factor(trueclasses), seq(1,length(unique(trueclasses))))[j]] <- 1 
  }
  
  svmquadloss[[i]] <- sum((scoresmatrix[,getMarkerClasses(markerhl)] - allocmatrix)^2) #compute quad loss
}

  return(svmquadloss = svmquadloss)
}




