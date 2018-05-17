params2 <- bayesgmmOptimisation(object = hl, iterations = 100)
res <- bayesgmmPredict(object = hl, optimRes = params2)
metrics <- bayesgmmCval(object = hl, times = 100)

f1score<-matrix(0, 100)

for(i in 1:100){
  
  conf <- metrics$cmlist[[i]]
  
  .p <- MLInterfaces:::precision(conf)
  .r <- MLInterfaces:::recall(conf)
  
  f1 <- MLInterfaces:::.macroF1(.p, .r, naAs0 = TRUE) ## macro F1 score for .time's iteration
  
  f1score[i] <- f1  ##store macro F1 score
  
  
}

boxplot(f1score, xlab="Macro F1 score BGMM", ylim=c(0,1), main="Mouse Stem Cell Data")

#Per organelle F1 scores

f1scoreperO<-matrix(0,100,dim(metrics$cmlist[[1]])[1])

for(i in 1:100){
  
  conf<-metrics$cmlist[[i]]
  
  f1perO <- MLInterfaces:::F1(conf, naAs0 = TRUE) ## macro F1 score for .time's iteration
  
  f1scoreperO[i,]<-f1perO  ##store macro F1 score
  
  
}




setStockcol(paste0(getStockcol(), 90))
ptsze <- exp(as.numeric(fData(res$object)$predicted.probability))-0.99
.pca <- plot2D(res$object, fcol="predicted.allocation", dims = c(1,2), method = "t-SNE", cex = ptsze, main="Prediction with pointer size scaled with probability of membership" )
addLegend(hl, cex=0.6, where="topleft")
setStockcol(getLisacol())


rownames(res$globalalloc)[rowSums(res$globalalloc)>0.5]
source("https://bioconductor.org/biocLite.R")
biocLite("clusterProfiler")
library(clusterProfiler)
source("https://bioconductor.org/biocLite.R")
biocLite("org.Mm.eg.db")
library(org.Mm.eg.db)

ego <- enrichGO(gene          = names(nml)[nml==T],
                universe      = rownames(object),
                OrgDb         = org.Mm.eg.db,
                keyType = "UNIPROT",
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
dotplot(ego, showCategory = 20)

egomf <- enrichGO(gene          = names(nml)[nml==T],
                universe      = rownames(object),
                OrgDb         = org.Mm.eg.db,
                keyType = "UNIPROT",
                ont           = "MF",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

dotplot(egomf, showCategory = 20)

egobp <- enrichGO(gene          = names(nml)[nml==T],
                  universe      = rownames(object),
                  OrgDb         = org.Mm.eg.db,
                  keyType = "UNIPROT",
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)

dotplot(egobp, showCategory = 20)

  