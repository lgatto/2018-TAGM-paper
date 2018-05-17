
andyrun1 <- andy2015mcmc[[1]]
andyrun2 <- andy2015mcmc[[2]]
andyrun3 <- andy2015mcmc[[3]]
andyrun4 <- andy2015mcmc[[4]]
andyrun5 <- andy2015mcmc[[5]]
andyrun6 <- andy2015mcmc[[6]]
andyrun7 <- andy2015mcmc[[7]]
andyrun8 <- andy2015mcmc[[8]]
andyrun9 <- andy2015mcmc[[9]]

rm(andy2015mcmc)

par(mfrow = c(2,3))
plot(colSums(andyrun1$allocstruc), type = "b", cex = 0.1, col = "blue", ylim = c(3200,3450))
plot(colSums(andyrun3$allocstruc), type = "b", cex = 0.1, col = "blue", ylim = c(3200,3450))
plot(colSums(andyrun4$allocstruc), type = "b", cex = 0.1, col = "blue", ylim = c(3200,3450))
plot(colSums(andyrun5$allocstruc), type = "b", cex = 0.1, col = "blue", ylim = c(3200,3450))
plot(colSums(andyrun8$allocstruc), type = "b", cex = 0.1, col = "blue", ylim = c(3200,3450))
plot(colSums(andyrun9$allocstruc), type = "b", cex = 0.1, col = "blue", ylim = c(3200,3450))

andydiag <- list(andymcmc1 = andyrun1$allocstruc,
                 andymcmc2 = andyrun3$allocstruc,
                 andymcmc3 = andyrun4$allocstruc,
                 andymcmc4 = andyrun5$allocstruc,
                 andymcmc5 = andyrun8$allocstruc,
                 andymcmc6 = andyrun9$allocstruc)

save(andydiag, file = "andydiagtest.rda")


require(coda)
listmcmc <- list(mcmc(colSums(andyrun1$allocstruc)[501:1100]),
                 mcmc(colSums(andyrun3$allocstruc)[501:1100]),
                 mcmc(colSums(andyrun4$allocstruc)[501:1100]),
                 mcmc(colSums(andyrun8$allocstruc)[501:1100]),
                 mcmc(colSums(andyrun9$allocstruc)[501:1100]))
listmcmc <- as.mcmc.list(listmcmc)
gd <- gelman.diag(x = listmcmc, autoburnin = F)
gd$psrf

dim(andyrun2$allocprob)

andycmb <- abind(andyrun1$allocprob[,,501:1100], andyrun3$allocprob[,,501:1100], andyrun4$allocprob[,,501:1100]
  ,andyrun8$allocprob[,,501:1100], andyrun9$allocprob[,,501:1100])

andycmbstruc <- abind(andyrun1$allocstrucprob[,501:1100,1], andyrun3$allocstrucprob[,501:1100,1], andyrun4$allocstrucprob[,501:1100,1]
                      , andyrun8$allocstrucprob[,501:1100,1], andyrun9$allocstrucprob[,501:1100,1])

andycmbstrucalloc <- cbind(andyrun1$allocstruc[,501:1100],andyrun3$allocstruc[,501:1100], andyrun4$allocstruc[,501:1100]
                           ,andyrun8$allocstruc[,501:1100], andyrun9$allocstruc[,501:1100])



membershipProb <- t(apply(andycmb[,], 1, FUN = rowMeans))
data("hyperLOPIT2015")
object <- hyperLOPIT2015

rownames(andycmb) <- rownames(unknownMSnSet(object))
  
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

outlier.allocation <- c(rowMeans(andycmbstrucalloc), rep(1,length(fData(markerMSnSet(object))$markers)))
names(outlier.allocation) <- c(rownames(unknownMSnSet(object)), rownames(fData(markerMSnSet(object))))

fData(object)$outlier.allocation <- outlier.allocation[rownames(fData(object))]


ptsze <- exp(as.numeric(fData(object)$predicted.probability)) - 1
names(ptsze) <- rownames(fData(object))
ptsze[as.numeric(fData(object)$outlier.allocation) <0.99999] <- 0.01

par(mfrow = c(1,1))
setStockcol(paste0(getStockcol(), 90))
.pca <- plot2D(object, fcol="predicted.allocation", cex = ptsze, main="Prediction with pointer size scaled with probability of membership" )
addLegend(object, cex=0.6, where="bottomleft")
setStockcol(getLisacol())

N <- nrow(X)
K <- length(getMarkerClasses(object))
probquantilesAll<-array(0, c(N, K, 2))
for(i in 1:N){
  for(j in 1:K){
    probquantilesAll[i,j,] <- quantile(andycmb[i,j,], probs=c(0.025,0.975))
  }
}

#calculate quantile differences, mean quantile distance and per organelle quantile distance.
quantDiff <- probquantilesAll[, , 2] - probquantilesAll[, , 1]
mean(quantDiff)
boxplot(quantDiff)

shannon <- - apply(andycmb * log(andycmb), c(1,3), sum)
shannon[is.na(shannon)] <- 0
meanshannon <- rowMeans(shannon)
meanshannon <- c(meanshannon, rep(0,length(fData(markerMSnSet(object))$markers)))
names(meanshannon) <- c(rownames(unknownMSnSet(object)), rownames(fData(markerMSnSet(object))))
 
fData(andyMCMChl)$meanshannon <- meanshannon[rownames(andyMCMChl)]

plot(as.numeric(fData(andyMCMChl)$predicted.probability), meanshannon[rownames(andyMCMChl)], col = getStockcol()[as.factor(fData(andyMCMChl)$predicted.allocation)], pch = 19)

uncertain <- rowSums(quantDiff)
uncertain <- c(uncertain, rep(0,length(fData(markerMSnSet(object))$markers)))
names(uncertain) <- c(rownames(unknownMSnSet(object)), rownames(fData(markerMSnSet(object))))
fData(object)$uncertain <- uncertain[rownames(fData(object))] 
setStockcol(paste0(getStockcol(), 90))
ptsze <- 5*fData(andyMCMChl)$meanshannon
plot2D(andyMCMChl, fcol = "predicted.allocation", cex = ptsze, main = "Visualising whole proteome uncertainty", mirrorY = FALSE)
plot(as.numeric(fData(andyMCMChl)$predicted.probability), fData(andyMCMChl)$uncertain, col = getStockcol()[as.factor(fData(andyMCMChl)$predicted.allocation)], pch = 19)


andyMCMChl <- object
setwd("C:/Users/OllyC/Desktop/bayesian-spatial-proteomics/code/data/")
save(andyMCMChl, file = "andyMCMChl.rda")

rowMeans(andycmbstrucalloc)[31]



dfG5E870 <- t(andycmb[31,,])
save(dfG5E870, file = "probdistE3.rda")

boxplot(df)
colnames(df) <- getMarkerClasses(hl)
df <- melt(df)
colnames(df) <- c("One","Organelle","Probability")

gg <- ggplot(df, aes(Organelle, Probability, width = (Probability))) + geom_violin(aes(fill = Organelle), scale = "width")
gg <- gg + scale_fill_manual(values = getStockcol()[1:14]) + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x = element_blank())
gg <- gg + ylab("Membership Probability") + ggtitle(paste0("Distribution of Subcellular Membership for Protein ", rownames(andycmb)[31]))
gg <- gg + theme(legend.position="none")
gg

prots <- c(rownames(andycmb)[33])
foi13s <- FeaturesOfInterest(description = "",
                             fnames = prots,
                             object = hl)

foi13s
setStockcol(paste0(getStockcol(), 90))
ptsze <- exp(as.numeric(fData(object)$predicted.probability)) - 1
names(ptsze) <- rownames(fData(object))
ptsze[as.numeric(fData(object)$outlier.allocation) <0.99995] <- 0.01
plot2D(andyMCMChl, fcol = "predicted.allocation", cex = ptsze, dims = c(1,2), main = "PCA plot with Protein Q9WUA2 indicated")
setStockcol(getLisacol())
addLegend(hl, cex = .75)
highlightOnPlot(hl, foi13s, cex = 1.5, col = "yellow", pch = 19)
highlightOnPlot(hl, foi13s, labels = TRUE, cex = 2, col = "black", pch = 19, pos = 4)






dfQ924C1 <- t(andycmb["Q924C1",,])
save(dfQ924C1, file = "probdistQ924C1.rda")
boxplot(df2)
colnames(df2) <- getMarkerClasses(hl)
df2 <- melt(df2)
colnames(df2) <- c("One","Organelle","Probability")

gg2 <- ggplot(df2, aes(Organelle, Probability, width = (Probability))) + geom_violin(aes(fill = Organelle), scale = "width")
gg2 <- gg2 + scale_fill_manual(values = getStockcol()[1:14]) + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x = element_blank())
gg2 <- gg2 + ylab("Membership Probability") + ggtitle(paste0("Distribution of Subcellular Membership for Protein Q924C1" ))
gg2 <- gg2 + theme(legend.position="none")
gg2   

setStockcol()
gg <- ggplot(df2[df2$Organelle == c("Cytosol","Nucleus - Non-chromatin"),], aes(Probability, fill = Organelle)) + geom_density(aes(fill = Organelle), alpha = 0.1) + xlim(c(0,1))
gg <- gg + scale_fill_manual(values = getStockcol()[c(4,11)]) + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x = element_blank())
gg <- gg + ylab("Density") + xlab("Membership Probability") + ggtitle("Density plot of Subcellular Membership for Protein Q924C1")
gg <- gg + theme_bw()
gg


prots <- "Q924C1"
foi13s <- FeaturesOfInterest(description = "",
                             fnames = prots,
                             object = hl)

foi13s
setStockcol(paste0(getStockcol(), 90))
ptsze <- exp(as.numeric(fData(object)$predicted.probability)) - 1
names(ptsze) <- rownames(fData(object))
ptsze[as.numeric(fData(object)$outlier.allocation) <0.99995] <- 0.01
plot2D(andyMCMChl, fcol = "predicted.allocation", cex = ptsze, dims = c(1,2), main = "PCA plot with Protein Q924C1 indicated")
setStockcol(getLisacol())
addLegend(hl, cex = .75)
highlightOnPlot(hl, foi13s, cex = 1.5, col = "black", pch = 19)
highlightOnPlot(hl, foi13s, labels = TRUE, cex = 2, col = "black", pch = 19, pos = 4)






dfQ9WUA2 <- t(andycmb["Q9WUA2",,])
save(dfQ9WUA2, file = "probdistQ9WUA2.rda")
boxplot(df3)
colnames(df3) <- getMarkerClasses(hl)
df3 <- melt(df3)
colnames(df3) <- c("One","Organelle","Probability")

gg3 <- ggplot(df3, aes(Organelle, Probability, width = (Probability))) + geom_violin(aes(fill = Organelle), scale = "width")
gg3 <- gg3 + scale_fill_manual(values = getStockcol()[1:14]) + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x = element_blank())
gg3 <- gg3 + ylab("Membership Probability") + ggtitle(paste0("Distribution of Subcellular Membership for Protein Q9WUA2" ))
gg3 <- gg3 + theme(legend.position="none")
gg3   


prots <- "Q9WUA2"
foi13s <- FeaturesOfInterest(description = "",
                             fnames = prots,
                             object = hl)

foi13s
setStockcol(paste0(getStockcol(), 90))
ptsze <- exp(as.numeric(fData(object)$predicted.probability)) - 1
names(ptsze) <- rownames(fData(object))
ptsze[as.numeric(fData(object)$outlier.allocation) <0.99995] <- 0.01
plot2D(andyMCMChl, fcol = "predicted.allocation", cex = ptsze, dims = c(1,2), main = "PCA plot with Protein Q9WUA2 indicated")
setStockcol(getLisacol())
addLegend(hl, cex = .75)
highlightOnPlot(hl, foi13s, cex = 1.5, col = "black", pch = 19)
highlightOnPlot(hl, foi13s, labels = TRUE, cex = 2, col = "black", pch = 19, pos = 4)



dfQ8VDR9 <- t(andycmb["Q8VDR9",,])
save(dfQ8VDR9, file = "probdistQ8VDR9.rda")
boxplot(df4)
colnames(df4) <- getMarkerClasses(hl)
df4 <- melt(df4)
colnames(df4) <- c("One","Organelle","Probability")

gg4 <- ggplot(df4, aes(Organelle, Probability, width = (Probability))) + geom_violin(aes(fill = Organelle), scale = "width")
gg4 <- gg4 + scale_fill_manual(values = getStockcol()[1:14]) + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x = element_blank())
gg4 <- gg4 + ylab("Membership Probability") + ggtitle(paste0("Distribution of Subcellular Membership for Protein Q8VDR9" ))
gg4 <- gg4 + theme(legend.position="none")
gg4   

prots <- "P97287"
foi13s <- FeaturesOfInterest(description = "",
                             fnames = prots,
                             object = hl)

foi13s
setStockcol(paste0(getStockcol(), 90))
ptsze <- exp(as.numeric(fData(object)$predicted.probability)) - 1
names(ptsze) <- rownames(fData(object))
ptsze[as.numeric(fData(object)$outlier.allocation) <0.99995] <- 0.01
plot2D(andyMCMChl, fcol = "predicted.allocation", cex = ptsze, dims = c(1,2), main = "PCA plot with Protein Q8VDR9 indicated")
setStockcol(getLisacol())
addLegend(hl, cex = .75)
highlightOnPlot(hl, foi13s, cex = 1.5, col = "black", pch = 19)
highlightOnPlot(hl, foi13s, labels = TRUE, cex = 2, col = "black", pch = 19, pos = 4)


uncertain[uncertain > 0]

data("hyperLOPIT2015goCC")
hlgoCC <- hyperLOPIT2015goCC

numCC <- rowSums(exprs(hlgoCC))

plot(uncertain[uncertain > 0.1], numCC[names(uncertain[uncertain > 0.1])])

source("https://bioconductor.org/biocLite.R")
require("clusterProfiler")
nmlgo <- rowMeans(andycmbstruc)<0.95
names(nmlgo) <- rownames(unknownMSnSet(object))

ego <- enrichGO(gene          = rownames(fData(andyMCMChl))[fData(andyMCMChl)$outlier.allocation < 0.95],
                universe      = rownames(exprs(object)),
                OrgDb         = org.Mm.eg.db,
                ont           = "CC",
                keyType       = "UNIPROT",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

egobp <- enrichGO(gene        = rownames(fData(andyMCMChl))[fData(andyMCMChl)$outlier.allocation < 0.95],
                universe      = rownames(exprs(andyMCMChl)),
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                keyType       = "UNIPROT",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

egoMF <- enrichGO(gene          = rownames(fData(andyMCMChl))[fData(andyMCMChl)$outlier.allocation < 0.95],
                  universe      = rownames(exprs(object)),
                  OrgDb         = org.Mm.eg.db,
                  ont           = "MF",
                  keyType       = "UNIPROT",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)


DOSE::dotplot(ego, showCategory=10, title = "Cellular Compartment Enrichment", font.size = 16)
DOSE::dotplot(egobp, showCategory=10, title = "Biological Process Enrichment", font.size = 16)
DOSE::dotplot(egoMF, showCategory = 10, title = "Molecular Function Enrichment", font.size = 16)

unkhl <- unknownMSnSet(andyMCMChl)

X <- fData(andyMCMChl)[fData(andyMCMChl)$predicted.probability>0.99,"svm.classification"]

Y <- fData(andyMCMChl)$predicted.allocation[fData(andyMCMChl)$predicted.probability>0.99]

for(j in seq_along(getMarkerClasses(andyMCMChl))){
  Y[Y==j]<-getMarkerClasses(hyperLOPIT2015)[j]
}

cont <- table(X,Y)

mycol <- colorRampPalette(c("white", "red"))(10)


library(lattice)
levelplot(cont)
cont <- cont/rowSums(cont)
levelplot(cont, xlab="SVM", ylab="BayesGMM (MCMC)", main="Agreement of TAGM at 0.95 threshold compared with SVM", scales = list(x = list(draw = FALSE)), 
          cuts = 50, pretty = T, col.regions = mycol, at = seq(0, 10)/10)

X2 <- fData(andyMCMChl)[fData(andyMCMChl)[,"final.assignment"]!="unknown","final.assignment"]
Y2 <- fData(andyMCMChl)$predicted.allocation[fData(andyMCMChl)[,"final.assignment"]!="unknown"]

for(j in seq_along(getMarkerClasses(andyMCMChl))){
  Y2[Y2==j] <- getMarkerClasses(hyperLOPIT2015)[j]
}
df3<-table(X2,Y2)
levelplot(df3)
df4 <- df3/rowSums(df3)
levelplot(df4, xlab="SVM", ylab="BGMM (MCMC)", main="Final assignments from SVM compared with TAGM", scales = list(x = list(draw = FALSE)),
          cuts = 50, pretty = T, col.regions = mycol, at = seq(0,10)/10)



