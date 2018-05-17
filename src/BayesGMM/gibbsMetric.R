#' Calculate Cross-validation Metrics 
#' 
#' 
#' 
#' 
#'@params
#'
#'
#'@return
#'
#'
#'

gibbsMetric <- function(test.idx,
                        object,
                        fcol = "markers",
                        iterations = 1000L,
                        burnin = 100L,
                        thin = 5L,
                        mu0 = NULL,
                        lambda0 = 0.01,
                        nu0 = NULL,
                        S0 = NULL,
                        beta0 = NULL,
                        u = 2,
                        v = 10
){
  
  suppressPackageStartupMessages({
    require(pRoloc)
    require(MSnbase)
  })
  
  
  marker.data <- markerMSnSet(object)
  X <- pRoloc:::subsetAsDataFrame(marker.data, "markers", train = TRUE)
  
  .test1   <- MSnSet(exprs(marker.data)[test.idx, ], fData(marker.data[test.idx, ]), pData(marker.data)) ## 'unseen' test set
  .train1  <- MSnSet(exprs(marker.data)[-test.idx, ], fData(marker.data[-test.idx, ]), pData(marker.data)) ## 'seen' training set
  
  test.markers <- fData(.test1)$markers #save true marker labels
  mydata <- combine(.test1, .train1) #create new combined MSnset
  levels(fData(mydata)$markers) <- c(levels(fData(mydata)$markers), "unknown")
  fData(mydata)[rownames(.test1), "markers"] <- "unknown" #hide marker labels
  
  params <- cmp_Gibbs(object = mydata, iterations = iterations, burnin = burnin) #optimsiation
  
  N <- nrow(X[test.idx, ])
  K <- length(getMarkerClasses(mydata))
  
  membershipProb <- matrix(0, N, K)
  
  for(i in seq_len(N)){
    membershipProb[i, ] <- rowMeans(params$allocprob[i, , ])
  }
  
  organellealloc <- matrix(0, nrow = N, ncol = 2)
  organellealloc[, 1] <- getMarkerClasses(mydata)[apply(membershipProb, 1, which.max)]
  proballoc <- apply(membershipProb, 1, which.max)
  
  for(i in seq_len(N)){
    organellealloc[i, 2] <- as.numeric(membershipProb[i, proballoc[i]])
  }
  rownames(membershipProb) <- rownames(unknownMSnSet(mydata))
  rownames(organellealloc) <- rownames(unknownMSnSet(mydata))
  
  pred <- c(organellealloc[, 1], as.character(fData(markerMSnSet(mydata))[,"markers"]))
  prob <- c(organellealloc[, 2], rep(1,length(fData(markerMSnSet(mydata))$markers)))
  
  names(prob) <- c(rownames(unknownMSnSet(mydata)), rownames(fData(markerMSnSet(mydata))))
  names(pred) <- c(rownames(unknownMSnSet(mydata)), rownames(fData(markerMSnSet(mydata))))
  
  fData(mydata)$predicted.allocation <- pred[rownames(fData(mydata))] 
  fData(mydata)$predicted.probability <- prob[rownames(fData(mydata))]
  
  data <- factor(fData(mydata[rownames(.test1),])$predicted.allocation, levels = getMarkerClasses(mydata))
  reference <- factor(test.markers, levels = getMarkerClasses(mydata))
  
  conf <- confusionMatrix(data = data, reference = reference)$table
  
  quad <- vector("numeric", length = length(test.idx))
  allocmatrix <- matrix(0, length(test.idx), K)
  
  for(j in seq_along(test.idx)){
    allocmatrix[j, as.numeric(factor(test.markers), seq(1,K))[j]] <- 1 
  }
  quadloss <- sum(rowSums((allocmatrix - membershipProb[rownames(.test1),])^2))
  
return(list(conf = conf, quadloss = quadloss))
} 

  
