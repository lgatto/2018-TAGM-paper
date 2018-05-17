
membershipProb <- matrix(0,N,K)

for(i in 1:N){
   membershipProb[i, ] <- rowMeans(res2$allocprob[i,,-(1:100)])
}

X <- exprs(unknownMSnSet(object))
organellealloc <- matrix(0, nrow = nrow(X), ncol = 2)
organellealloc[, 1] <- getMarkerClasses(object)[apply(membershipProb, 1, which.max)]
proballoc <- apply(membershipProb, 1, which.max)

for(i in 1:nrow(X)){
  organellealloc[i, 2] <- as.numeric(membershipProb[i, proballoc[i]])
}
rownames(membershipProb) <- rownames(unknownMSnSet(object))
rownames(organellealloc) <- rownames(unknownMSnSet(object))

pred <- c(organellealloc[, 1], as.character(fData(markerMSnSet(object))[,"markers"]))
prob <- c(organellealloc[, 2], rep(1,length(fData(markerMSnSet(object))$markers)))

names(prob) <- c(rownames(unknownMSnSet(object)), rownames(fData(markerMSnSet(object))))
names(pred) <- c(rownames(unknownMSnSet(object)), rownames(fData(markerMSnSet(object))))

fData(object)$predicted.allocation <- pred[rownames(fData(object))] 
fData(object)$predicted.probability <- prob[rownames(fData(object))]

nml <- rowMeans(res2$allocstrucprob[,-(1:100),1])<0.995
names(nml) <- rownames(unknownMSnSet(object))

ptsze <- as.numeric(fData(object)$predicted.probability)^2
names(ptsze) <- rownames(fData(object))
ptsze[names(nml)[nml==T]] <- 0.01
ptsze[names(uncertain)[uncertain >10^{-4}]] <- 0.01

setStockcol(paste0(getStockcol(), 90))
.pca <- plot2D(object, fcol="predicted.allocation", dims = c(1,2), method = "t-SNE", cex = ptsze, main="Prediction with pointer size scaled with probability of membership" )
addLegend(hl, cex=0.6, where="topleft")
setStockcol(getLisacol())

probquantilesAll<-array(0, c(N,K,2))
for(i in 1:N){
  for(j in 1:K){
    probquantilesAll[i,j,] <- quantile(res2$allocprob[i,j,-(1:100)], probs=c(0.025,0.975))
  }
}

#calculate quantile differences, mean quantile distance and per organelle quantile distance.
quantDiff<-probquantilesAll[,,2]-probquantilesAll[,,1]
mean(quantDiff)
boxplot(quantDiff)

uncertain <- rowSums(quantDiff)/1.5
uncertain <- c(uncertain, rep(0,length(fData(markerMSnSet(object))$markers)))
names(uncertain) <- c(rownames(unknownMSnSet(object)), rownames(fData(markerMSnSet(object))))
fData(object)$uncertain <- uncertain[rownames(fData(object))] 
setStockcol(paste0(getStockcol(), 90))
ptze <- fData(object)$uncertain 
plot2D(object, fcol = "predicted.allocation", dims = c(1,2), cex = ptze, main = "Visualising whole proteome uncertainty", mirrorY = FALSE)
