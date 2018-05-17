source('C:/Users/OllyC/Desktop/bayesian-spatial-proteomics/code/BayesGMM/bayesgmmCval.R')
source('C:/Users/OllyC/Desktop/bayesian-spatial-proteomics/code/BayesGMM/bayesgmmOptimisation.R')
source('C:/Users/OllyC/Desktop/bayesian-spatial-proteomics/code/BayesGMM/bayesgmmPredict.R')

set.seed(1)

p <- lapply(hirstconcvalmcmc$cmlist, function(x) MLInterfaces:::precision(x))
r <- lapply(hirstconcvalmcmc$cmlist, function(x) MLInterfaces:::recall(x))
hirstconF1 <- vector("numeric", length = 100)
for(i in 1:100){
  hirstconF1[i] <- MLInterfaces:::.macroF1(p[[i]], r[[i]], naAs0 = TRUE)
}
mean(hirstconF1)


p <- lapply(hirstc2cvalmcmc$cmlist, function(x) MLInterfaces:::precision(x))
r <- lapply(hirstc2cvalmcmc$cmlist, function(x) MLInterfaces:::recall(x))
hirstc2F1 <- vector("numeric", length = 100)
for(i in 1:100){
  hirstc2F1[i] <- MLInterfaces:::.macroF1(p[[i]], r[[i]], naAs0 = TRUE)
}
mean(hirstc2F1)

p <- lapply(hirstc6cvalmcmc$cmlist, function(x) MLInterfaces:::precision(x))
r <- lapply(hirstc6cvalmcmc$cmlist, function(x) MLInterfaces:::recall(x))
hirstc6F1 <- vector("numeric", length = 100)
for(i in 1:100){
  hirstc6F1[i] <- MLInterfaces:::.macroF1(p[[i]], r[[i]], naAs0 = TRUE)
}
mean(hirstc6F1)

p <- lapply(beltranM24cvalmcmc$cmlist, function(x) MLInterfaces:::precision(x))
r <- lapply(beltranM24cvalmcmc$cmlist, function(x) MLInterfaces:::recall(x))
beltranM24F1 <- vector("numeric", length = 100)
for(i in 1:100){
  beltranM24F1[i] <- MLInterfaces:::.macroF1(p[[i]], r[[i]], naAs0 = TRUE)
}
mean(beltranM24F1)

p <- lapply(beltranM48cvalmcmc$cmlist, function(x) MLInterfaces:::precision(x))
r <- lapply(beltranM48cvalmcmc$cmlist, function(x) MLInterfaces:::recall(x))
beltranM48F1 <- vector("numeric", length = 100)
for(i in 1:100){
  beltranM48F1[i] <- MLInterfaces:::.macroF1(p[[i]], r[[i]], naAs0 = TRUE)
}
mean(beltranM48F1)

p <- lapply(beltranM72cvalmcmc$cmlist, function(x) MLInterfaces:::precision(x))
r <- lapply(beltranM72cvalmcmc$cmlist, function(x) MLInterfaces:::recall(x))
beltranM72F1 <- vector("numeric", length = 100)
for(i in 1:100){
  beltranM72F1[i] <- MLInterfaces:::.macroF1(p[[i]], r[[i]], naAs0 = TRUE)
}
mean(beltranM72F1)

p <- lapply(beltranM96cvalmcmc$cmlist, function(x) MLInterfaces:::precision(x))
r <- lapply(beltranM96cvalmcmc$cmlist, function(x) MLInterfaces:::recall(x))
beltranM96F1 <- vector("numeric", length = 100)
for(i in 1:100){
  beltranM96F1[i] <- MLInterfaces:::.macroF1(p[[i]], r[[i]], naAs0 = TRUE)
}
mean(beltranM96F1)

p <- lapply(beltranM120cvalmcmc$cmlist, function(x) MLInterfaces:::precision(x))
r <- lapply(beltranM120cvalmcmc$cmlist, function(x) MLInterfaces:::recall(x))
beltranM120F1 <- vector("numeric", length = 100)
for(i in 1:100){
  beltranM120F1[i] <- MLInterfaces:::.macroF1(p[[i]], r[[i]], naAs0 = TRUE)
}
mean(beltranM120F1)

p <- lapply(beltranH24cvalmcmc$cmlist, function(x) MLInterfaces:::precision(x))
r <- lapply(beltranH24cvalmcmc$cmlist, function(x) MLInterfaces:::recall(x))
beltranH24F1 <- vector("numeric", length = 100)
for(i in 1:100){
  beltranH24F1[i] <- MLInterfaces:::.macroF1(p[[i]], r[[i]], naAs0 = TRUE)
}
mean(beltranH24F1)

p <- lapply(beltranH48cvalmcmc$cmlist, function(x) MLInterfaces:::precision(x))
r <- lapply(beltranH48cvalmcmc$cmlist, function(x) MLInterfaces:::recall(x))
beltranH48F1 <- vector("numeric", length = 100)
for(i in 1:100){
  beltranH48F1[i] <- MLInterfaces:::.macroF1(p[[i]], r[[i]], naAs0 = TRUE)
}
mean(beltranH48F1)

p <- lapply(beltranH72cvalmcmc$cmlist, function(x) MLInterfaces:::precision(x))
r <- lapply(beltranH72cvalmcmc$cmlist, function(x) MLInterfaces:::recall(x))
beltranH72F1 <- vector("numeric", length = 100)
for(i in 1:100){
  beltranH72F1[i] <- MLInterfaces:::.macroF1(p[[i]], r[[i]], naAs0 = TRUE)
}
mean(beltranH72F1)

p <- lapply(beltranH96cvalmcmc$cmlist, function(x) MLInterfaces:::precision(x))
r <- lapply(beltranH96cvalmcmc$cmlist, function(x) MLInterfaces:::recall(x))
beltranh96F1 <- vector("numeric", length = 100)
for(i in 1:100){
  beltranh96F1[i] <- MLInterfaces:::.macroF1(p[[i]], r[[i]], naAs0 = TRUE)
}
mean(beltranh96F1)

p <- lapply(beltranH120cvalmcmc$cmlist, function(x) MLInterfaces:::precision(x))
r <- lapply(beltranH120cvalmcmc$cmlist, function(x) MLInterfaces:::recall(x))
beltranH120F1 <- vector("numeric", length = 100)
for(i in 1:100){
  beltranH120F1[i] <- MLInterfaces:::.macroF1(p[[i]], r[[i]], naAs0 = TRUE)
}
mean(beltranH120F1)


p <- lapply(groencvalmcmc$cmlist, function(x) MLInterfaces:::precision(x))
r <- lapply(groencvalmcmc$cmlist, function(x) MLInterfaces:::recall(x))
groenF1 <- vector("numeric", length = 100)
for(i in 1:100){
  groenF1[i] <- MLInterfaces:::.macroF1(p[[i]], r[[i]], naAs0 = TRUE)
}
mean(groenF1)


p <- lapply(E14TG2aRcvalmcmc$cmlist, function(x) MLInterfaces:::precision(x))
r <- lapply(E14TG2aRcvalmcmc$cmlist, function(x) MLInterfaces:::recall(x))
E14TG2aRF1 <- vector("numeric", length = 100)
for(i in 1:100){
  E14TG2aRF1[i] <- MLInterfaces:::.macroF1(p[[i]], r[[i]], naAs0 = TRUE)
}
mean(E14TG2aRF1)


w <- table(getMarkers(hirstcon, verbose = F))
w <- 1/w[names(w) != "unknown"]
hirstconsvmparams <- svmOptimisation(hirstcon, fcol = "markers", times = 100, xval = 5, class.weights = w)
hirstconknnparams <- knnOptimisation(hirstcon, fcol = "markers", times = 100, xval = 5)
hirstconMAPparams <- bayesgmmCval(object = hirstcon, times = 100)
hirstconsvmquadloss <- svmquadloss(svmparam = hirstconsvmparams, object = hirstcon)
hirstconknnquadloss <- knnquadloss(knnparams = hirstconknnparams, object = hirstcon)

mean(hirstconsvmparams@results[,1])
mean(hirstconknnparams@results[,1])

phirstconMAPparams <- lapply(hirstconMAPparams$cmlist, function(x) MLInterfaces:::precision(x))
rhirstconMAPparams <- lapply(hirstconMAPparams$cmlist, function(x) MLInterfaces:::recall(x))
hirstconMAPF1 <- vector("numeric", length = 100)
for(i in 1:100){
  hirstconMAPF1[i] <- MLInterfaces:::.macroF1(phirstconMAPparams[[i]], rhirstconMAPparams[[i]], naAs0 = TRUE)
}
mean(hirstconMAPF1)

save(hirstconsvmparams, file = "hirstconsvm.rda")
save(hirstconknnparams, file = "hirstconknn.rda")
save(hirstconMAPparams, file = "hirstconMAP.rda")
save(hirstconsvmquadloss, file = "hirstconsvmquadloss.rda")
save(hirstconknnquadloss, file = "hirstconknnquadloss.rda")


w <- table(getMarkers(hirstc2, verbose = F))
w <- 1/w[names(w) != "unknown"]
hirstc2svmparams <- svmOptimisation(hirstc2, fcol = "markers", times = 100, xval = 5, class.weights = w)
hirstc2knnparams <- knnOptimisation(hirstc2, fcol = "markers", times = 100, xval = 5)
hirstc2MAPparams <- bayesgmmCval(object = hirstc2, times = 100)
hirstc2svmquadloss <- svmquadloss(svmparam = hirstc2svmparams, object = hirstc2)
hirstc2knnquadloss <- knnquadloss(knnparams = hirstc2knnparams, object = hirstc2)

mean(hirstc2svmparams@results[,1])
mean(hirstc2knnparams@results[,1])

phirstc2MAPparams <- lapply(hirstc2MAPparams$cmlist, function(x) MLInterfaces:::precision(x))
rhirstc2MAPparams <- lapply(hirstc2MAPparams$cmlist, function(x) MLInterfaces:::recall(x))
hirstc2MAPF1 <- vector("numeric", length = 100)
for(i in 1:100){
  hirstc2MAPF1[i] <- MLInterfaces:::.macroF1(phirstc2MAPparams[[i]], rhirstc2MAPparams[[i]], naAs0 = TRUE)
}
mean(hirstc2MAPF1)

save(hirstc2svmparams, file = "hirstc2svm.rda")
save(hirstc2knnparams, file = "hirstc2knn.rda")
save(hirstc2MAPparams, file = "hirstc2MAP.rda")
save(hirstc2svmquadloss, file = "hirstc2svmquadloss.rda")
save(hirstc2knnquadloss, file = "hirstc2knnquadloss.rda")


w <- table(getMarkers(hirstc6, verbose = F))
w <- 1/w[names(w) != "unknown"]
hirstc6svmparams <- svmOptimisation(hirstc6, fcol = "markers", times = 100, xval = 5, class.weights = w)
hirstc6knnparams <- knnOptimisation(hirstc2, fcol = "markers", times = 100, xval = 5)
hirstc6MAPparams <- bayesgmmCval(object = hirstc6, times = 100)
hirstc6svmquadloss <- svmquadloss(svmparam = hirstc6svmparams, object = hirstc6)
hirstc6knnquadloss <- knnquadloss(knnparams = hirstc6knnparams, object = hirstc6)

mean(hirstc6svmparams@results[,1])
mean(hirstc6knnparams@results[,1])

phirstc6MAPparams <- lapply(hirstc6MAPparams$cmlist, function(x) MLInterfaces:::precision(x))
rhirstc6MAPparams <- lapply(hirstc6MAPparams$cmlist, function(x) MLInterfaces:::recall(x))
hirstc6MAPF1 <- vector("numeric", length = 100)
for(i in 1:100){
  hirstc6MAPF1[i] <- MLInterfaces:::.macroF1(phirstc6MAPparams[[i]], rhirstc6MAPparams[[i]], naAs0 = TRUE)
}
mean(hirstc6MAPF1)

save(hirstc6svmparams, file = "hirstc6svm.rda")
save(hirstc6knnparams, file = "hirstc6knn.rda")
save(hirstc6MAPparams, file = "hirstc6MAP.rda")
save(hirstc6svmquadloss, file = "hirstc6svmquadloss.rda")
save(hirstc6knnquadloss, file = "hirstc6knnquadloss.rda")


w <- table(getMarkers(beltran2016MOCK24, verbose = F))
w <- 1/w[names(w) != "unknown"]
beltranM24svmparams <- svmOptimisation(beltran2016MOCK24, fcol = "markers", times = 100, xval = 5, class.weights = w)
beltranM24knnparams <- knnOptimisation(beltran2016MOCK24, fcol = "markers", times = 100, xval = 5)
beltranM24MAPparams <- bayesgmmCval(object = beltran2016MOCK24, times = 100)
beltranM24svmquadloss <- svmquadloss(svmparam = beltranM24svmparams, object = beltran2016MOCK24)
beltranM24knnquadloss <- knnquadloss(knnparams = beltranM24knnparams, object = beltran2016MOCK24)

mean(beltranM24svmparams@results[,1])
mean(beltranM24knnparams@results[,1])

pbeltranM24MAPparams <- lapply(beltranM24MAPparams$cmlist, function(x) MLInterfaces:::precision(x))
rbeltranM24MAPparams <- lapply(beltranM24MAPparams$cmlist, function(x) MLInterfaces:::recall(x))
beltranM24MAPF1 <- vector("numeric", length = 100)
for(i in 1:100){
  beltranM24MAPF1[i] <- MLInterfaces:::.macroF1(pbeltranM24MAPparams[[i]], rbeltranM24MAPparams[[i]], naAs0 = TRUE)
}
mean(beltranM24MAPF1)

save(beltranM24svmparams, file = "beltranM24nsvm.rda")
save(beltranM24knnparams, file = "beltranM24knn.rda")
save(beltranM24MAPparams, file = "beltranM24MAP.rda")
save(beltranM24svmquadloss, file = "beltranM24svmquadloss.rda")
save(beltranM24knnquadloss, file = "beltranM24knnquadloss.rda")

data("beltran2016MOCK48")
w <- table(getMarkers(beltran2016MOCK48, verbose = F))
w <- 1/w[names(w) != "unknown"]
beltranM48svmparams <- svmOptimisation(beltran2016MOCK48, fcol = "markers", times = 100, xval = 5, class.weights = w)
beltranM48knnparams <- knnOptimisation(beltran2016MOCK48, fcol = "markers", times = 100, xval = 5)
beltranM48MAPparams <- bayesgmmCval(object = beltran2016MOCK48, times = 100)
beltranM48svmquadloss <- svmquadloss(svmparam = beltranM48svmparams, object = beltran2016MOCK48)
beltranM48knnquadloss <- knnquadloss(knnparams = beltranM48knnparams, object = beltran2016MOCK48)


mean(beltranM48svmparams@results[,1])
mean(beltranM48knnparams@results[,1])

pbeltranM48MAPparams <- lapply(beltranM48MAPparams$cmlist, function(x) MLInterfaces:::precision(x))
rbeltranM48MAPparams <- lapply(beltranM48MAPparams$cmlist, function(x) MLInterfaces:::recall(x))
beltranM48MAPF1 <- vector("numeric", length = 100)
for(i in 1:100){
  beltranM48MAPF1[i] <- MLInterfaces:::.macroF1(pbeltranM48MAPparams[[i]], rbeltranM48MAPparams[[i]], naAs0 = TRUE)
}
mean(beltranM48MAPF1)

save(beltranM48svmparams, file = "beltranM48nsvm.rda")
save(beltranM48knnparams, file = "beltranM48knn.rda")
save(beltranM48MAPparams, file = "beltran48MAP.rda")
save(beltranM48svmquadloss, file = "beltranM48svmquadloss.rda")
save(beltranM48knnquadloss, file = "beltranM48knnquadloss.rda")

data("beltran2016MOCK72")
w <- table(getMarkers(beltran2016MOCK72, verbose = F))
w <- 1/w[names(w) != "unknown"]
beltranM72svmparams <- svmOptimisation(beltran2016MOCK72, fcol = "markers", times = 100, xval = 5, class.weights = w)
beltranM72knnparams <- knnOptimisation(beltran2016MOCK72, fcol = "markers", times = 100, xval = 5)
beltranM72MAPparams <- bayesgmmCval(object = beltran2016MOCK72, times = 100)
beltranM72svmquadloss <- svmquadloss(svmparam = beltranM72svmparams, object = beltran2016MOCK72)
beltranM72knnquadloss <- knnquadloss(knnparams = beltranM72knnparams, object = beltran2016MOCK72)


mean(beltranM72svmparams@results[,1])
mean(beltranM72knnparams@results[,1])

pbeltranM72MAPparams <- lapply(beltranM72MAPparams$cmlist, function(x) MLInterfaces:::precision(x))
rbeltranM72MAPparams <- lapply(beltranM72MAPparams$cmlist, function(x) MLInterfaces:::recall(x))
beltranM72MAPF1 <- vector("numeric", length = 100)
for(i in 1:100){
  beltranM72MAPF1[i] <- MLInterfaces:::.macroF1(pbeltranM72MAPparams[[i]], rbeltranM72MAPparams[[i]], naAs0 = TRUE)
}
mean(beltranM72MAPF1)



save(beltranM72svmparams, file = "beltranM72nsvm.rda")
save(beltranM72knnparams, file = "beltranM72knn.rda")
save(beltranM72MAPparams, file = "beltranM72MAP.rda")
save(beltranM72svmquadloss, file = "beltranM72svmquadloss.rda")
save(beltranM72knnquadloss, file = "beltranM72knnquadloss.rda")

data("beltran2016MOCK96")
w <- table(getMarkers(beltran2016MOCK96, verbose = F))
w <- 1/w[names(w) != "unknown"]
beltranM96svmparams <- svmOptimisation(beltran2016MOCK96, fcol = "markers", times = 100, xval = 5, class.weights = w)
beltranM96knnparams <- knnOptimisation(beltran2016MOCK96, fcol = "markers", times = 100, xval = 5)
beltranM96MAPparams <- bayesgmmCval(object = beltran2016MOCK96, times = 100)
beltranM96svmquadloss <- svmquadloss(svmparam = beltranM96svmparams, object = beltran2016MOCK96)
beltranM96knnquadloss <- knnquadloss(knnparams = beltranM96knnparams, object = beltran2016MOCK96)

mean(beltranM96svmparams@results[,1])
mean(beltranM96knnparams@results[,1])

pbeltranM96MAPparams <- lapply(beltranM96MAPparams$cmlist, function(x) MLInterfaces:::precision(x))
rbeltranM96MAPparams <- lapply(beltranM96MAPparams$cmlist, function(x) MLInterfaces:::recall(x))
beltranM96MAPF1 <- vector("numeric", length = 100)
for(i in 1:100){
  beltranM96MAPF1[i] <- MLInterfaces:::.macroF1(pbeltranM96MAPparams[[i]], rbeltranM96MAPparams[[i]], naAs0 = TRUE)
}
mean(beltranM96MAPF1)

save(beltranM96svmparams, file = "beltranM96nsvm.rda")
save(beltranM96knnparams, file = "beltranM96knn.rda")
save(beltranM96MAPparams, file = "beltranM96MAP.rda")
save(beltranM96svmquadloss, file = "beltranM96svmquadloss.rda")
save(beltranM96knnquadloss, file = "beltranM96knnquadloss.rda")


data("beltran2016MOCK120")
w <- table(getMarkers(beltran2016MOCK120, verbose = F))
w <- 1/w[names(w) != "unknown"]
beltranM120svmparams <- svmOptimisation(beltran2016MOCK120, fcol = "markers", times = 100, xval = 5, class.weights = w)
beltranM120knnparams <- knnOptimisation(beltran2016MOCK120, fcol = "markers", times = 100, xval = 5)
beltranM120MAPparams <- bayesgmmCval(object = beltran2016MOCK120, times = 100)
beltranM120svmquadloss <- svmquadloss(svmparam = beltranM120svmparams, object = beltran2016MOCK120)
beltranM120knnquadloss <- knnquadloss(knnparams = beltranM120knnparams, object = beltran2016MOCK120)

mean(beltranM120svmparams@results[,1])
mean(beltranM120knnparams@results[,1])

pbeltranM120MAPparams <- lapply(beltranM120MAPparams$cmlist, function(x) MLInterfaces:::precision(x))
rbeltranM120MAPparams <- lapply(beltranM120MAPparams$cmlist, function(x) MLInterfaces:::recall(x))
beltranM120MAPF1 <- vector("numeric", length = 100)
for(i in 1:100){
  beltranM120MAPF1[i] <- MLInterfaces:::.macroF1(pbeltranM120MAPparams[[i]], rbeltranM120MAPparams[[i]], naAs0 = TRUE)
}
mean(beltranM120MAPF1)

save(beltranM120svmparams, file = "beltranM120nsvm.rda")
save(beltranM120knnparams, file = "beltranM120knn.rda")
save(beltranM120MAPparams, file = "beltranM120MAP.rda")
save(beltranM120svmquadloss, file = "beltranM120svmquadloss.rda")
save(beltranM120knnquadloss, file = "beltranM120knnquadloss.rda")

data("beltran2016HCMV24")
w <- table(getMarkers(beltran2016HCMV24, verbose = F))
w <- 1/w[names(w) != "unknown"]
beltranH24svmparams <- svmOptimisation(beltran2016HCMV24, fcol = "markers", times = 100, xval = 5, class.weights = w)
beltranH24knnparams <- knnOptimisation(beltran2016HCMV24, fcol = "markers", times = 100, xval = 5)
beltranH24MAPparams <- bayesgmmCval(object = beltran2016HCMV24, times = 100)
beltranH24svmquadloss <- svmquadloss(svmparam = beltranH24svmparams, object = beltran2016HCMV24)
beltranH24knnquadloss <- knnquadloss(knnparams = beltranH24knnparams, object = beltran2016HCMV24)

mean(beltranH24svmparams@results[,1])
mean(beltranH24knnparams@results[,1])
pbeltranH24MAPparams <- lapply(beltranH24MAPparams$cmlist, function(x) MLInterfaces:::precision(x))
rbeltranH24MAPparams <- lapply(beltranH24MAPparams$cmlist, function(x) MLInterfaces:::recall(x))
beltranH24MAPF1 <- vector("numeric", length = 100)
for(i in 1:100){
  beltranH24MAPF1[i] <- MLInterfaces:::.macroF1(pbeltranH24MAPparams[[i]], rbeltranH24MAPparams[[i]], naAs0 = TRUE)
}
mean(beltranH24MAPF1)



save(beltranH24svmparams, file = "beltranH24nsvm.rda")
save(beltranH24knnparams, file = "beltranH24knn.rda")
save(beltranH24MAPparams, file = "beltranH24MAP.rda")
save(beltranH24svmquadloss, file = "beltranH24svmquadloss.rda")
save(beltranH24knnquadloss, file = "beltranH24knnquadloss.rda")

data("beltran2016HCMV48")
w <- table(getMarkers(beltran2016HCMV48, verbose = F))
w <- 1/w[names(w) != "unknown"]
beltranH48svmparams <- svmOptimisation(beltran2016HCMV48, fcol = "markers", times = 100, xval = 5, class.weights = w)
beltranH48knnparams <- knnOptimisation(beltran2016HCMV48, fcol = "markers", times = 100, xval = 5)
beltranH48MAPparams <- bayesgmmCval(object = beltran2016HCMV48, times = 100)
beltranH48svmquadloss <- svmquadloss(svmparam = beltranH48svmparams, object = beltran2016HCMV48)
beltranH48knnquadloss <- knnquadloss(knnparams = beltranH48knnparams, object = beltran2016HCMV48)


mean(beltranH48svmparams@results[,1])
mean(beltranH48knnparams@results[,1])
pbeltranH48MAPparams <- lapply(beltranH48MAPparams$cmlist, function(x) MLInterfaces:::precision(x))
rbeltranH48MAPparams <- lapply(beltranH48MAPparams$cmlist, function(x) MLInterfaces:::recall(x))
beltranH48MAPF1 <- vector("numeric", length = 100)
for(i in 1:100){
  beltranH48MAPF1[i] <- MLInterfaces:::.macroF1(pbeltranH48MAPparams[[i]], rbeltranH48MAPparams[[i]], naAs0 = TRUE)
}
mean(beltranH48MAPF1)

save(beltranH48svmparams, file = "beltranH48nsvm.rda")
save(beltranH48knnparams, file = "beltranH48knn.rda")
save(beltranH48MAPparams, file = "beltranH48MAP.rda")
save(beltranH48svmquadloss, file = "beltranH48svmquadloss.rda")
save(beltranH48knnquadloss, file = "beltranH48knnquadloss.rda")

data("beltran2016HCMV72")
w <- table(getMarkers(beltran2016HCMV72, verbose = F))
w <- 1/w[names(w) != "unknown"]
beltranH72svmparams <- svmOptimisation(beltran2016HCMV72, fcol = "markers", times = 100, xval = 5, class.weights = w)
beltranH72knnparams <- knnOptimisation(beltran2016HCMV72, fcol = "markers", times = 100, xval = 5)
beltranH72MAPparams <- bayesgmmCval(object = beltran2016HCMV72, times = 100)
beltranH72svmquadloss <- svmquadloss(svmparam = beltranH72svmparams, object = beltran2016HCMV72)
beltranH72knnquadloss <- knnquadloss(knnparams = beltranH72knnparams, object = beltran2016HCMV72)

mean(beltranH72svmparams@results[,1])
mean(beltranH72knnparams@results[,1])
pbeltranH72MAPparams <- lapply(beltranH72MAPparams$cmlist, function(x) MLInterfaces:::precision(x))
rbeltranH72MAPparams <- lapply(beltranH72MAPparams$cmlist, function(x) MLInterfaces:::recall(x))
beltranH72MAPF1 <- vector("numeric", length = 100)
for(i in 1:100){
  beltranH72MAPF1[i] <- MLInterfaces:::.macroF1(pbeltranH72MAPparams[[i]], rbeltranH72MAPparams[[i]], naAs0 = TRUE)
}
mean(beltranH72MAPF1)

save(beltranH72svmparams, file = "beltranH72nsvm.rda")
save(beltranH72knnparams, file = "beltranH72knn.rda")
save(beltranH72MAPparams, file = "beltranH72MAP.rda")
save(beltranH72svmquadloss, file = "beltranH72svmquadloss.rda")
save(beltranH72knnquadloss, file = "beltranH72knnquadloss.rda")

data("beltran2016HCMV96")
w <- table(getMarkers(beltran2016HCMV96, verbose = F))
w <- 1/w[names(w) != "unknown"]
beltranH96svmparams <- svmOptimisation(beltran2016HCMV96, fcol = "markers", times = 100, xval = 5, class.weights = w)
beltranH96knnparams <- knnOptimisation(beltran2016HCMV96, fcol = "markers", times = 100, xval = 5)
beltranH96MAPparams <- bayesgmmCval(object = beltran2016HCMV96, times = 100)
beltranH96svmquadloss <- svmquadloss(svmparam = beltranH96svmparams, object = beltran2016HCMV96)
beltranH96knnquadloss <- knnquadloss(knnparams = beltranH96knnparams, object = beltran2016HCMV96)


mean(beltranH96svmparams@results[,1])
mean(beltranH96knnparams@results[,1])
pbeltranH96MAPparams <- lapply(beltranH96MAPparams$cmlist, function(x) MLInterfaces:::precision(x))
rbeltranH96MAPparams <- lapply(beltranH96MAPparams$cmlist, function(x) MLInterfaces:::recall(x))
beltranH96MAPF1 <- vector("numeric", length = 100)
for(i in 1:100){
  beltranH96MAPF1[i] <- MLInterfaces:::.macroF1(pbeltranH96MAPparams[[i]], rbeltranH96MAPparams[[i]], naAs0 = TRUE)
}
mean(beltranH96MAPF1)


save(beltranH96svmparams, file = "beltranH96nsvm.rda")
save(beltranH96knnparams, file = "beltranH96knn.rda")
save(beltranH96MAPparams, file = "beltranH96MAP.rda")
save(beltranH96svmquadloss, file = "beltranH96svmquadloss.rda")
save(beltranH96knnquadloss, file = "beltranH96knnquadloss.rda")

data("beltran2016HCMV120")
w <- table(getMarkers(beltran2016HCMV120, verbose = F))
w <- 1/w[names(w) != "unknown"]
beltranH120svmparams <- svmOptimisation(beltran2016HCMV120, fcol = "markers", times = 100, xval = 5, class.weights = w)
beltranH120knnparams <- knnOptimisation(beltran2016HCMV120, fcol = "markers", times = 100, xval = 5)
beltranH120MAPparams <- bayesgmmCval(object = beltran2016HCMV120, times = 100)
beltranH120svmquadloss <- svmquadloss(svmparam = beltranH120svmparams, object = beltran2016HCMV120)
beltranH120knnquadloss <- knnquadloss(knnparams = beltranH120knnparams, object = beltran2016HCMV120)

mean(beltranH120svmparams@results[,1])
mean(beltranH120knnparams@results[,1])

pbeltranH120MAPparams <- lapply(beltranH120MAPparams$cmlist, function(x) MLInterfaces:::precision(x))
rbeltranH120MAPparams <- lapply(beltranH120MAPparams$cmlist, function(x) MLInterfaces:::recall(x))
beltranH120MAPF1 <- vector("numeric", length = 100)
for(i in 1:100){
  beltranH120MAPF1[i] <- MLInterfaces:::.macroF1(pbeltranH120MAPparams[[i]], rbeltranH120MAPparams[[i]], naAs0 = TRUE)
}
mean(beltranH120MAPF1)

save(beltranH120svmparams, file = "beltranH120nsvm.rda")
save(beltranH120knnparams, file = "beltranH120knn.rda")
save(beltranH120MAPparams, file = "beltranH120MAP.rda")
save(beltranH120svmquadloss, file = "beltranH120svmquadloss.rda")
save(beltranH120knnquadloss, file = "beltranH120knnquadloss.rda")

data("E14TG2aR")
w <- table(getMarkers(E14TG2aR, verbose = F))
w <- 1/w[names(w) != "unknown"]
E14TG2aRsvmparams <- svmOptimisation(E14TG2aR, fcol = "markers", times = 100, xval = 5, class.weights = w)
E14TG2aRknnparams <- knnOptimisation(E14TG2aR, fcol = "markers", times = 100, xval = 5)
E14TG2aRMAPparams <- bayesgmmCval(object = E14TG2aR, times = 100)
E14TG2aRsvmquadloss <- svmquadloss(svmparam = E14TG2aRsvmparams, object = E14TG2aR)
E14TG2aRknnquadloss <- knnquadloss(knnparams = E14TG2aRknnparams, object = E14TG2aR)

mean(E14TG2aRsvmparams@results[,1])
mean(E14TG2aRknnparams@results[,1])
pE14TG2aRMAPparams <- lapply(E14TG2aRMAPparams$cmlist, function(x) MLInterfaces:::precision(x))
rE14TG2aRMAPparams <- lapply(E14TG2aRMAPparams$cmlist, function(x) MLInterfaces:::recall(x))
E14TG2aRMAPF1 <- vector("numeric", length = 100)
for(i in 1:100){
  E14TG2aRMAPF1[i] <- MLInterfaces:::.macroF1(pE14TG2aRMAPparams[[i]], rE14TG2aRMAPparams[[i]], naAs0 = TRUE)
}
mean(E14TG2aRMAPF1)


save(E14TG2aRsvmparams, file = "E14TG2aRsvm.rda")
save(E14TG2aRknnparams, file = "E14TG2aRknn.rda")
save(E14TG2aRMAPparams, file = "E14TG2aRMAP.rda")
save(E14TG2aRsvmquadloss, file = "E14TG2aRsvmquadloss.rda")
save(E14TG2aRknnquadloss, file = "E14TG2aRknnquadloss.rda")

hirstcondfmetrics <- NA

hirstcondfmetrics <- as.data.frame(cbind(hirstconsvmparams@results[,1], hirstconknnparams@results[,1], hirstconMAPF1, hirstconF1))
colnames(hirstcondfmetrics) <- c("SVM", "KNN", "MAP", "MCMC")
hirstcondfmetrics <- melt(hirstcondfmetrics)
hirstcondfmetrics$set <- factor(rep("HeLa Wild (Hirst et al. 2018)", 400 ))

hirstc2dfmetrics <- NA

hirstc2dfmetrics <- as.data.frame(cbind(hirstc2svmparams@results[,1], hirstc2knnparams@results[,1], hirstc2MAPF1, hirstc2F1))
colnames(hirstc2dfmetrics) <- c("SVM", "KNN", "MAP", "MCMC")
hirstc2dfmetrics <- melt(hirstc2dfmetrics)
hirstc2dfmetrics$set <- factor(rep("HeLa KO1 (Hirst et al. 2018)", 400 ))


hirstc6dfmetrics <- as.data.frame(cbind(hirstc6svmparams@results[,1], hirstc6knnparams@results[,1], hirstc6MAPF1, hirstc6F1))
colnames(hirstc6dfmetrics) <- c("SVM", "KNN", "MAP", "MCMC")
hirstc6dfmetrics <- melt(hirstc6dfmetrics)
hirstc6dfmetrics$set <- factor(rep("HeLa KO2 (Hirst et al. 2018)", 400 ))

beltranM24dfmetrics <- as.data.frame(cbind(beltranM24svmparams@results[,1], beltranM24knnparams@results[,1], beltranM24MAPF1, beltranM24F1 ))
colnames(beltranM24dfmetrics) <- c("SVM", "KNN", "MAP", "MCMC")
beltranM24dfmetrics <- melt(beltranM24dfmetrics)
beltranM24dfmetrics$set <- factor(rep("Primary Fibroblast Mock 24hpi", 400 ))

beltranM48dfmetrics <- as.data.frame(cbind(beltranM48svmparams@results[,1], beltranM48knnparams@results[,1], beltranM48MAPF1, beltranM48F1 ))
colnames(beltranM48dfmetrics) <- c("SVM", "KNN", "MAP", "MCMC")
beltranM48dfmetrics <- melt(beltranM48dfmetrics)
beltranM48dfmetrics$set <- factor(rep("Primary Fibroblast Mock 48hpi", 400 ))

beltranM72dfmetrics <- as.data.frame(cbind(beltranM72svmparams@results[,1], beltranM72knnparams@results[,1], beltranM72MAPF1, beltranM72F1 ))
colnames(beltranM72dfmetrics) <- c("SVM", "KNN", "MAP", "MCMC")
beltranM72dfmetrics <- melt(beltranM72dfmetrics)
beltranM72dfmetrics$set <- factor(rep("Primary Fibroblast Mock 72hpi", 400 ))

beltranM96dfmetrics <- as.data.frame(cbind(beltranM96svmparams@results[,1], beltranM96knnparams@results[,1], beltranM96MAPF1, beltranM96F1 ))
colnames(beltranM96dfmetrics) <- c("SVM", "KNN", "MAP", "MCMC")
beltranM96dfmetrics <- melt(beltranM96dfmetrics)
beltranM96dfmetrics$set <- factor(rep("Primary Fibroblast Mock 96hpi", 400 ))

beltranM120dfmetrics <- as.data.frame(cbind(beltranM120svmparams@results[,1], beltranM120knnparams@results[,1], beltranM120MAPF1, beltranM120F1 ))
colnames(beltranM120dfmetrics) <- c("SVM", "KNN", "MAP", "MCMC")
beltranM120dfmetrics <- melt(beltranM120dfmetrics)
beltranM120dfmetrics$set <- factor(rep("Primary Fibroblast Mock 120hpi", 400 ))

beltranH24dfmetrics <- as.data.frame(cbind(beltranH24svmparams@results[,1], beltranH24knnparams@results[,1], beltranH24MAPF1, beltranH24F1 ))
colnames(beltranH24dfmetrics) <- c("SVM", "KNN", "MAP", "MCMC")
beltranH24dfmetrics <- melt(beltranH24dfmetrics)
beltranH24dfmetrics$set <- factor(rep("Primary Fibroblast HCMV 24hpi", 400 ))

beltranH48dfmetrics <- as.data.frame(cbind(beltranH48svmparams@results[,1], beltranH48knnparams@results[,1], beltranH48MAPF1, beltranH48F1 ))
colnames(beltranH48dfmetrics) <- c("SVM", "KNN", "MAP", "MCMC")
beltranH48dfmetrics <- melt(beltranH48dfmetrics)
beltranH48dfmetrics$set <- factor(rep("Primary Fibroblast HCMV 48hpi", 400 ))

beltranH72dfmetrics <- as.data.frame(cbind(beltranH72svmparams@results[,1], beltranH72knnparams@results[,1], beltranH72MAPF1, beltranH72F1 ))
colnames(beltranH72dfmetrics) <- c("SVM", "KNN", "MAP", "MCMC")
beltranH72dfmetrics <- melt(beltranH72dfmetrics)
beltranH72dfmetrics$set <- factor(rep("Primary Fibroblast HCMV 72hpi", 400 ))

beltranH96dfmetrics <- as.data.frame(cbind(beltranH96svmparams@results[,1], beltranH96knnparams@results[,1], beltranH96MAPF1, beltranh96F1 ))
colnames(beltranH96dfmetrics) <- c("SVM", "KNN", "MAP", "MCMC")
beltranH96dfmetrics <- melt(beltranH96dfmetrics)
beltranH96dfmetrics$set <- factor(rep("Primary Fibroblast HCMV 96hpi", 400 ))

beltranH120dfmetrics <- NA

beltranH120dfmetrics <- as.data.frame(cbind(beltranH120svmparams@results[,1], beltranH120knnparams@results[,1], beltranH120MAPF1, beltranH120F1 ))
colnames(beltranH120dfmetrics) <- c("SVM", "KNN", "MAP", "MCMC")
beltranH120dfmetrics <- melt(beltranH120dfmetrics)
beltranH120dfmetrics$set <- factor(rep("Primary Fibroblast HCMV 120hpi", 400 ))

E14TG2aRdfmetrics <- as.data.frame(cbind(E14TG2aRsvmparams@results[,1], E14TG2aRknnparams@results[,1], E14TG2aRMAPF1, E14TG2aRF1))
colnames(E14TG2aRdfmetrics) <- c("SVM", "KNN", "MAP", "MCMC")
E14TG2aRdfmetrics <- melt(E14TG2aRdfmetrics)
E14TG2aRdfmetrics$set <- factor(rep("E14TG2aR (Breckels et al. 2016)", 400 ))

dfmetrics <- rbind(hirstcondfmetrics, hirstc2dfmetrics, hirstc6dfmetrics,
                   beltranM24dfmetrics, beltranM48dfmetrics,beltranM72dfmetrics,
                   beltranM96dfmetrics, beltranM120dfmetrics, beltranH24dfmetrics,
                   beltranH48dfmetrics, beltranH72dfmetrics, beltranH96dfmetrics,
                   beltranH120dfmetrics, E14TG2aRdfmetrics)

additionalF1compare <- dfmetrics

save(additionalF1compare, file = "additionalF1compare.rda")

gg <- ggplot(data = dfmetrics, aes(x=variable, y=value)) +
  geom_boxplot(aes(fill=variable)) +
  scale_y_continuous(limits = c(0, 1)) +
  ggtitle("Boxplot of Macro F1 scores") + 
  theme_bw() + labs(y = "Macro F1-Score") +  labs(x  = "Classifier") + facet_grid(. ~ set) + facet_wrap(~set) + scale_x_discrete(labels = element_blank()) + scale_fill_discrete(name = "Classifier")
gg

hirstcondfmetricsquad <- NA

hirstcondfmetricsquad <- as.data.frame(cbind(unlist(hirstconsvmquadloss), unlist(hirstconknnquadloss), 
                                             unlist(hirstconMAPparams$quadloss), unlist(hirstconcvalmcmc$quadloss)))
colnames(hirstcondfmetricsquad) <- c("SVM", "KNN", "MAP", "MCMC")
hirstcondfmetricsquad <- melt(hirstcondfmetricsquad)
hirstcondfmetricsquad$set <- factor(rep("HeLa Wild (Hirst et al. 2018)", 400 ))

hirstc2dfmetricsquad <- NA

hirstc2dfmetricsquad <- as.data.frame(cbind(unlist(hirstc2svmquadloss), unlist(hirstc2knnquadloss), 
                                             unlist(hirstc2MAPparams$quadloss), unlist(hirstc2cvalmcmc$quadloss)))
colnames(hirstc2dfmetricsquad) <- c("SVM", "KNN", "MAP", "MCMC")
hirstc2dfmetricsquad <- melt(hirstc2dfmetricsquad)
hirstc2dfmetricsquad$set <- factor(rep("HeLa KO1 (Hirst et al. 2018)", 400 ))


hirstc6dfmetricsquad <- as.data.frame(cbind(unlist(hirstc6svmquadloss), unlist(hirstc6knnquadloss), 
                                            unlist(hirstc6MAPparams$quadloss), unlist(hirstc6cvalmcmc$quadloss)))
colnames(hirstc6dfmetricsquad) <- c("SVM", "KNN", "MAP", "MCMC")
hirstc6dfmetricsquad <- melt(hirstc6dfmetricsquad)
hirstc6dfmetricsquad$set <- factor(rep("HeLa KO2 (Hirst et al. 2018)", 400 ))

beltranM24dfmetricsquad <- NA

beltranM24dfmetricsquad <- as.data.frame(cbind(unlist(beltranM24svmquadloss), unlist(beltranM24knnquadloss),
                                         unlist(beltranM24MAPparams$quadloss), unlist(beltranM24cvalmcmc$quadloss)))
colnames(beltranM24dfmetricsquad) <- c("SVM", "KNN", "MAP", "MCMC")
beltranM24dfmetricsquad <- melt(beltranM24dfmetricsquad)
beltranM24dfmetricsquad$set <- factor(rep("Primary Fibroblast Mock 24hpi", 400 ))

beltranM48dfmetricsquad <- NA

beltranM48dfmetricsquad <- as.data.frame(cbind(unlist(beltranM48svmquadloss), unlist(beltranM48knnquadloss),
                                     unlist(beltranM48MAPparams$quadloss), unlist(beltranM48cvalmcmc$quadloss)))
colnames(beltranM48dfmetricsquad) <- c("SVM", "KNN", "MAP", "MCMC")
beltranM48dfmetricsquad <- melt(beltranM48dfmetricsquad)
beltranM48dfmetricsquad$set <- factor(rep("Primary Fibroblast Mock 48hpi", 400 ))

beltranM72dfmetricsquad <- as.data.frame(cbind(unlist(beltranM72svmquadloss), unlist(beltranM72knnquadloss),
                                     unlist(beltranM72MAPparams$quadloss), unlist(beltranM72cvalmcmc$quadloss)))
colnames(beltranM72dfmetricsquad) <- c("SVM", "KNN", "MAP", "MCMC")
beltranM72dfmetricsquad <- melt(beltranM72dfmetricsquad)
beltranM72dfmetricsquad$set <- factor(rep("Primary Fibroblast Mock 72hpi", 400 ))

beltranM96dfmetricsquad <- as.data.frame(cbind(unlist(beltranM96svmquadloss), unlist(beltranM96knnquadloss),
                                     unlist(beltranM96MAPparams$quadloss), unlist(beltranM96cvalmcmc$quadloss)))
colnames(beltranM96dfmetricsquad) <- c("SVM", "KNN", "MAP", "MCMC")
beltranM96dfmetricsquad <- melt(beltranM96dfmetricsquad)
beltranM96dfmetricsquad$set <- factor(rep("Primary Fibroblast Mock 96hpi", 400 ))

beltranM120dfmetricsquad <-  as.data.frame(cbind(unlist(beltranM120svmquadloss), unlist(beltranM120knnquadloss),
                                       unlist(beltranM120MAPparams$quadloss), unlist(beltranM120cvalmcmc$quadloss)))
colnames(beltranM120dfmetricsquad) <- c("SVM", "KNN", "MAP", "MCMC")
beltranM120dfmetricsquad <- melt(beltranM120dfmetricsquad)
beltranM120dfmetricsquad$set <- factor(rep("Primary Fibroblast Mock 120hpi", 400 ))

beltranH24dfmetricsquad <- as.data.frame(cbind(unlist(beltranH24svmquadloss), unlist(beltranH24knnquadloss),
                                         unlist(beltranH24MAPparams$quadloss), unlist(beltranH24cvalmcmc$quadloss)))
colnames(beltranH24dfmetricsquad) <- c("SVM", "KNN", "MAP", "MCMC")
beltranH24dfmetricsquad <- melt(beltranH24dfmetricsquad)
beltranH24dfmetricsquad$set <- factor(rep("Primary Fibroblast HCMV 24hpi", 400))

beltranH48dfmetricsquad <- as.data.frame(cbind(unlist(beltranH48svmquadloss), unlist(beltranH48knnquadloss),
                                     unlist(beltranH48MAPparams$quadloss), unlist(beltranH48cvalmcmc$quadloss)))
colnames(beltranH48dfmetricsquad) <- c("SVM", "KNN", "MAP", "MCMC")
beltranH48dfmetricsquad <- melt(beltranH48dfmetricsquad)
beltranH48dfmetricsquad$set <- factor(rep("Primary Fibroblast HCMV 48hpi", 400 ))

beltranH72dfmetricsquad <- as.data.frame(cbind(unlist(beltranH72svmquadloss), unlist(beltranH72knnquadloss),
                                     unlist(beltranH72MAPparams$quadloss), unlist(beltranH72cvalmcmc$quadloss)))
colnames(beltranH72dfmetricsquad) <- c("SVM", "KNN", "MAP", "MCMC")
beltranH72dfmetricsquad <- melt(beltranH72dfmetricsquad)
beltranH72dfmetricsquad$set <- factor(rep("Primary Fibroblast HCMV 72hpi", 400 ))

beltranH96dfmetricsquad <- as.data.frame(cbind(unlist(beltranH96svmquadloss), unlist(beltranH96knnquadloss),
                                         unlist(beltranH96MAPparams$quadloss), unlist(beltranH96cvalmcmc$quadloss)))
colnames(beltranH96dfmetricsquad) <- c("SVM", "KNN", "MAP", "MCMC")
beltranH96dfmetricsquad <- melt(beltranH96dfmetricsquad)
beltranH96dfmetricsquad$set <- factor(rep("Primary Fibroblast HCMV 96hpi", 400 ))

beltranH120dfmetricsquad <- NA

beltranH120dfmetricsquad <- as.data.frame(cbind(unlist(beltranH120svmquadloss), unlist(beltranH120knnquadloss),
                                      unlist(beltranH120MAPparams$quadloss), unlist(beltranH120cvalmcmc$quadloss)))
colnames(beltranH120dfmetricsquad) <- c("SVM", "KNN", "MAP", "MCMC")
beltranH120dfmetricsquad <- melt(beltranH120dfmetricsquad)
beltranH120dfmetricsquad$set <- factor(rep("Primary Fibroblast HCMV 120hpi", 400 ))

E14TG2aRdfmetricsquad <- as.data.frame(cbind(unlist(E14TG2aRsvmquadloss), unlist(E14TG2aRknnquadloss),
                                             unlist(E14TG2aRMAPparams$quadloss), unlist(E14TG2aRcvalmcmc$quadloss)))
colnames(E14TG2aRdfmetricsquad) <- c("SVM", "KNN", "MAP", "MCMC")
E14TG2aRdfmetricsquad <- melt(E14TG2aRdfmetricsquad)
E14TG2aRdfmetricsquad$set <- factor(rep("E14TG2aR (Breckels et al. 2016)", 400 ))



dfmetrics2 <- rbind(hirstcondfmetricsquad, hirstc2dfmetricsquad, hirstc6dfmetricsquad,
                   beltranM24dfmetricsquad, beltranM48dfmetricsquad,beltranM72dfmetricsquad,
                   beltranM96dfmetricsquad, beltranM120dfmetricsquad, beltranH24dfmetricsquad,
                   beltranH48dfmetricsquad, beltranH72dfmetricsquad, beltranH96dfmetricsquad,
                   beltranH120dfmetricsquad, E14TG2aRdfmetricsquad)

additionalquadlosscompare <- dfmetrics2

save(additionalquadlosscompare, file = "additionalquadlosscompare.rda")


gg2 <- ggplot(data = dfmetrics2, aes(x=variable, y=value)) +
  geom_boxplot(aes(fill=variable)) +
  scale_y_continuous() +
  ggtitle("Boxplot of Quadratic Losses") +
  theme_bw() + labs(y = "Quadratic Loss")  +  labs(x  = "Classifier") + facet_grid(. ~ set) + facet_wrap(~set) + scale_x_discrete(labels = element_blank()) + scale_fill_discrete(name = "Classifier")
gg2

gg <- ggplot(data = dfmetrics[1:3200,], aes(x=variable, y=value)) +
  geom_boxplot(aes(fill=variable)) +
  scale_y_continuous(limits = c(0, 1)) +
  ggtitle("Boxplot of Macro F1 scores") + 
  theme_bw() + labs(y = "Macro F1-Score") +  labs(x  = "Classifier") + facet_grid(. ~ set) + facet_wrap(~ set, ncol = 4) + scale_x_discrete(labels = element_blank()) + scale_fill_discrete(name = "Classifier")
gg

gg <- ggplot(data = dfmetrics[3201:7600,], aes(x=variable, y=value)) +
  geom_boxplot(aes(fill=variable)) +
  scale_y_continuous(limits = c(0, 1)) +
  ggtitle("Boxplot of Macro F1 scores") + 
  theme_bw() + labs(y = "Macro F1-Score") +  labs(x  = "Classifier") + facet_grid(. ~ set) + facet_wrap(~ set, ncol = 4) + scale_x_discrete(labels = element_blank()) + scale_fill_discrete(name = "Classifier")
gg

gg2 <- ggplot(data = dfmetrics2[1:3200,], aes(x=variable, y=value)) +
  geom_boxplot(aes(fill=variable)) +
  scale_y_continuous() +
  ggtitle("Boxplot of Quadratic Losses") +
  theme_bw() + labs(y = "Quadratic Loss")  +  labs(x  = "Classifier") + facet_grid(. ~ set) + facet_wrap(~set, ncol = 4) + scale_x_discrete(labels = element_blank()) + scale_fill_discrete(name = "Classifier")
gg2

gg2 <- ggplot(data = dfmetrics2[3201:7600,], aes(x=variable, y=value)) +
  geom_boxplot(aes(fill=variable)) +
  scale_y_continuous() +
  ggtitle("Boxplot of Quadratic Losses") +
  theme_bw() + labs(y = "Quadratic Loss")  +  labs(x  = "Classifier") + facet_grid(. ~ set) + facet_wrap(~set, ncol = 4) + scale_x_discrete(labels = element_blank()) + scale_fill_discrete(name = "Classifier")
gg2

