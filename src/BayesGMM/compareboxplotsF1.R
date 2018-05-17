setwd("C:/Users/OllyC/Desktop/bayesian-spatial-proteomics/code/data/")
list.files()
files <- list("hallknn.rda","hallknnqloss.rda","hallMAPcv.rda", "hallsvm.rda", "hallsvmqloss.rda","itzhakknn.rda",
"itzhakknnqloss.rda", "itzhakMAPcv.rda", "itzhaksvm.rda", "itzhaksvmqloss.rda", "tanknn.rda", "tanknnqloss.rda",
"tanMAPcv.rda", "tansvm.rda", "tansvmqloss.rda", "u2osknn.rda", "u2osknnqloss.rda", "u2osMAPcv.rda", "u2ossvm.rda" 
,"u2ossvmqloss.rda", "andy2015cvalmcmc-01.Rda", "andyknn.rda", "andyknnqloss.rda", "andyMAPcv.rda", "andysvm.rda",
"andysvmqloss.rda", "u2oscvalmcmc-01.Rda", "hallcvalmcmc-01.Rda", "tancvalmcmc-01.Rda","itzhakcvalmcmc-01.Rda")
for(i in seq_along(files)){
 load(files[[i]]) 
}


setwd("C:/Users/OllyC/Desktop/bayesian-spatial-proteomics/code/BayesGMM/")



p <- lapply(andyRes$cmlist, function(x) MLInterfaces:::precision(x))
r <- lapply(andyRes$cmlist, function(x) MLInterfaces:::recall(x))
andyMAPf1 <- vector("numeric", length = 100)
for(i in 1:100){
  andyMAPf1[i] <- MLInterfaces:::.macroF1(p[[i]], r[[i]], naAs0 = TRUE)
}

p <- lapply(andy2015cvalmcmc$cmlist, function(x) MLInterfaces:::precision(x))
r <- lapply(andy2015cvalmcmc$cmlist, function(x) MLInterfaces:::recall(x))
andyMCMCf1 <- vector("numeric", length = 100)
for(i in 1:100){
  andyMCMCf1[i] <- MLInterfaces:::.macroF1(p[[i]], r[[i]], naAs0 = TRUE)
}



p <- lapply(itzhakres$cmlist, function(x) MLInterfaces:::precision(x))
r <- lapply(itzhakres$cmlist, function(x) MLInterfaces:::recall(x))
itzhakMAPf1 <- vector("numeric", length = 100)
for(i in 1:100){
  itzhakMAPf1[i] <- MLInterfaces:::.macroF1(p[[i]], r[[i]], naAs0 = TRUE)
}

p <- lapply(itzhakcvalmcmc$cmlist, function(x) MLInterfaces:::precision(x))
r <- lapply(itzhakcvalmcmc$cmlist, function(x) MLInterfaces:::recall(x))
itzhakMCMCf1 <- vector("numeric", length = 100)
for(i in 1:100){
  itzhakMCMCf1[i] <- MLInterfaces:::.macroF1(p[[i]], r[[i]], naAs0 = TRUE)
}


p <- lapply(u2osres$cmlist, function(x) MLInterfaces:::precision(x))
r <- lapply(u2osres$cmlist, function(x) MLInterfaces:::recall(x))
u2osMAPf1 <- vector("numeric", length = 100)
for(i in 1:100){
  u2osMAPf1[i] <- MLInterfaces:::.macroF1(p[[i]], r[[i]], naAs0 = TRUE)
}

p <- lapply(u2oscvalmcmc$cmlist, function(x) MLInterfaces:::precision(x))
r <- lapply(u2oscvalmcmc$cmlist, function(x) MLInterfaces:::recall(x))
u2osMCMCf1 <- vector("numeric", length = 100)
for(i in 1:100){
  u2osMCMCf1[i] <- MLInterfaces:::.macroF1(p[[i]], r[[i]], naAs0 = TRUE)
}



p <- lapply(tanres$cmlist, function(x) MLInterfaces:::precision(x))
r <- lapply(tanres$cmlist, function(x) MLInterfaces:::recall(x))
tanMAPf1 <- vector("numeric", length = 100)
for(i in 1:100){
  tanMAPf1[i] <- MLInterfaces:::.macroF1(p[[i]], r[[i]], naAs0 = TRUE)
}

p <- lapply(tancvalmcmc$cmlist, function(x) MLInterfaces:::precision(x))
r <- lapply(tancvalmcmc$cmlist, function(x) MLInterfaces:::recall(x))
tanMCMCf1 <- vector("numeric", length = 100)
for(i in 1:100){
  tanMCMCf1[i] <- MLInterfaces:::.macroF1(p[[i]], r[[i]], naAs0 = TRUE)
}


p <- lapply(hallres$cmlist, function(x) MLInterfaces:::precision(x))
r <- lapply(hallres$cmlist, function(x) MLInterfaces:::recall(x))
hallMAPf1 <- vector("numeric", length = 100)
for(i in 1:100){
  hallMAPf1[i] <- MLInterfaces:::.macroF1(p[[i]], r[[i]], naAs0 = TRUE)
}

p <- lapply(hallcvalmcmc$cmlist, function(x) MLInterfaces:::precision(x))
r <- lapply(hallcvalmcmc$cmlist, function(x) MLInterfaces:::recall(x))
hallMCMCf1 <- vector("numeric", length = 100)
for(i in 1:100){
  hallMCMCf1[i] <- MLInterfaces:::.macroF1(p[[i]], r[[i]], naAs0 = TRUE)
}




andydfmetrics <- NA

andydfmetrics <- as.data.frame(cbind(andysvmparams@results[,1], andyknnparams@results[,1], andyMAPf1, andyMCMCf1))
colnames(andydfmetrics) <- c("SVM", "KNN", "MAP", "MCMC")
andydfmetrics <- melt(andydfmetrics)
andydfmetrics$set <- factor(rep("Mouse", 400 ))

itzhakdfmetrics <- NA

itzhakdfmetrics <- as.data.frame(cbind(itzhaksvmparams@results[,1], itzhakknnparams@results[,1], itzhakMAPf1, itzhakMCMCf1))
colnames(itzhakdfmetrics) <- c("SVM", "KNN", "MAP", "MCMC")
itzhakdfmetrics <- melt(itzhakdfmetrics)
itzhakdfmetrics$set <- factor(rep("HeLa", 400 ))


u2osdfmetrics <- as.data.frame(cbind(u2ossvmparams@results[,1], u2osknnparams@results[,1], u2osMAPf1, u2osMCMCf1))
colnames(u2osdfmetrics) <- c("SVM", "KNN", "MAP", "MCMC")
u2osdfmetrics <- melt(u2osdfmetrics)
u2osdfmetrics$set <- factor(rep("U2-OS", 400 ))

tandfmetrics <- as.data.frame(cbind(tansvmparams@results[,1], tanknnparams@results[,1], tanMAPf1, tanMCMCf1 ))
colnames(tandfmetrics) <- c("SVM", "KNN", "MAP", "MCMC")
tandfmetrics <- melt(tandfmetrics)
tandfmetrics$set <- factor(rep("Drosophila", 400 ))

halldfmetrics <- as.data.frame(cbind(hallsvmparams@results[,1], hallknnparams@results[,1], hallMAPf1, hallMCMCf1))
colnames(halldfmetrics) <- c("SVM", "KNN", "MAP", "MCMC")
halldfmetrics <- melt(halldfmetrics)
halldfmetrics$set <- factor(rep("Chicken DT40", 400 ))

dfmetrics <- rbind(tandfmetrics, halldfmetrics, andydfmetrics, itzhakdfmetrics, u2osdfmetrics)

gg <- ggplot(data = dfmetrics, aes(x=variable, y=value)) +
  geom_boxplot(aes(fill=variable)) +
  scale_y_continuous(limits = c(0, 1))  + 
  ggtitle("Boxplot of Macro F1 scores") + theme_bw() + labs(y = "Macro F1-Score") + 
  facet_grid(. ~ set) + theme(strip.background = element_blank(), strip.text.x = element_blank()) + scale_x_discrete(labels = element_blank())
gg


#ttest for f1 score and quadratic losses
f1ttesttan <- pairwise.t.test(tandfmetrics$value, tandfmetrics$variable , p.adjust.method = "BH",
                           pool.sd = FALSE, paired = FALSE,
                           alternative = "two.sided")
f1ttest$p.value
require(xtable)
xtable(f1ttest$p.value, digits = 5,
       caption = "Adjusted P-values for pairwise T-tests for Macro F-1 score classifier evaluation on the drosphila dataset")

f1ttesthall <- pairwise.t.test(halldfmetrics$value, halldfmetrics$variable , p.adjust.method = "BH",
                           pool.sd = FALSE, paired = FALSE,
                           alternative = "two.sided", mu = 0, var.equal = FALSE)
f1ttest$p.value
require(xtable)
xtable(f1ttest$p.value, digits = 5,
       caption = "Adjusted P-values for pairwise T-tests for Macro F-1 score classifier evaluation on the Chicken DT40 dataset")

f1ttestandy <- pairwise.t.test(andydfmetrics$value, andydfmetrics$variable , p.adjust.method = "BH",
                           pool.sd = FALSE, paired = FALSE,
                           alternative = "less")
f1ttest$p.value
require(xtable)
xtable(f1ttest$p.value, digits = 5,
       caption = "Adjusted P-values for pairwise T-tests for Macro F-1 score classifier evaluation on the mouse dataset")

f1ttestitshak <- pairwise.t.test(itzhakdfmetrics$value, itzhakdfmetrics$variable , p.adjust.method = "BH",
                           pool.sd = F, paired = FALSE,
                           alternative = "two.sided")
f1ttest$p.value
require(xtable)
xtable(f1ttest$p.value, digits = 5,
       caption = "Adjusted P-values for pairwise T-tests for Macro F-1 score classifier evaluation on the HeLa dataset")

f1ttestu2os <- pairwise.t.test(u20sdfmetrics$value, u20sdfmetrics$variable , p.adjust.method = "BH",
                           pool.sd = F, paired = FALSE,
                           alternative = "two.sided")
f1ttest$p.value
require(xtable)
xtable(f1ttest$p.value, digits = 5,
       caption = "Adjusted P-values for pairwise T-tests for Macro F-1 score classifier evaluation on the U2-OS dataset")
         