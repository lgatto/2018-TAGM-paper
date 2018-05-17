
prob <- t(apply(testtan$allocprob[,100:1000,],1,colMeans))
proboutlier
membershipProb <- prob

object <- tan2009r1
rownames(prob) <- rownames(unknownMSnSet(object))


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

ptsze <- exp(as.numeric(fData(object)$predicted.probability)) - 1

setStockcol(paste0(getStockcol(), 90))
plot2D(object, fcol="predicted.allocation", cex = ptsze, main="Prediction with pointer size scaled with probability of membership" )
addLegend(object, cex=0.6, where="bottomleft")
setStockcol(getLisacol())


hlo <- dunkley2006[,1:4]
plotDist(hlo[fData(hlo)$markers == "Golgi", ],
         pcol = getStockcol()[3], xlab= "", lwd = 10, type = "l")
matlines(t(exprs(hlo[fData(object)$predicted.allocation == "Golgi" , ])),lty = 1, col = "grey", lwd = 1, type="b", alpha = 0.1)

prob <- t(apply(testandy$allocprob[,-1,],1,colMeans))
proboutlier <- t(apply(testandy$allocOutprob[,-1,],1,colMeans))[,1]
membershipProb <- proboutlier * prob

object <- hyperLOPIT2015
rownames(prob) <- rownames(unknownMSnSet(object))


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

ptsze <- exp(as.numeric(fData(object)$predicted.probability)) - 1

setStockcol(paste0(getStockcol(), 90))
plot2D(object, fcol="predicted.allocation", cex = ptsze, main="Prediction with pointer size scaled with probability of membership" )
addLegend(object, cex=0.6, where="bottomleft")
setStockcol(getLisacol())