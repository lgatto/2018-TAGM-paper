setwd("C:/Users/OllyC/Desktop/bayesian-spatial-proteomics/code/data/")
list.files()
files <- list("hallknnqloss.rda", "hallsvmqloss.rda", "hallcvalmcmc-01.Rda", "hallMAPcv.rda",
              "itzhakknnqloss.rda", "itzhakMAPcv.rda", "itzhaksvmqloss.rda", "itzhakcvalmcmc-01.Rda",
              "tanknnqloss.rda", "tanMAPcv.rda", "tansvmqloss.rda", "tancvalmcmc-01.Rda",
              "u2osknnqloss.rda", "u2osMAPcv.rda", "u2oscvalmcmc-01.Rda", "u2ossvmqloss.rda", 
              "andy2015cvalmcmc-01.Rda", "andyknnqloss.rda", "andyMAPcv.rda", "andysvmqloss.rda")
for(i in seq_along(files)){
  load(files[[i]]) 
}


setwd("C:/Users/OllyC/Desktop/bayesian-spatial-proteomics/code/BayesGMM/")





andydfmetricsquad <- as.data.frame(cbind(unlist(andysvmquadloss), unlist(andyknnquadloss), unlist(andyRes$quadloss), unlist(andy2015cvalmcmc$quadloss)))
colnames(andydfmetricsquad) <- c("SVM", "KNN", "MAP", "MCMC")
andydfmetricsquad <- melt(andydfmetricsquad)
andydfmetricsquad$set <- factor(rep("Mouse", 400 ))


itzhakdfmetricsquad <- as.data.frame(cbind(unlist(itzhaksvmquadloss), unlist(itzhakknnquadloss), unlist(itzhakres$quadloss), unlist(itzhakcvalmcmc$quadloss)))
colnames(itzhakdfmetricsquad) <- c("SVM", "KNN", "MAP", "MCMC")
itzhakdfmetricsquad <- melt(itzhakdfmetricsquad)
itzhakdfmetricsquad$set <- factor(rep("HeLa", 400 ))


u2osdfmetricsquad <- as.data.frame(cbind(unlist(u2ossvmquadloss), unlist(u2osknnquadloss), unlist(u2osres$quadloss), unlist(u2oscvalmcmc$quadloss)))
colnames(u2osdfmetricsquad) <- c("SVM", "KNN", "MAP", "MCMC")
u2osdfmetricsquad <- melt(u2osdfmetricsquad)
u2osdfmetricsquad$set <- factor(rep("U2-OS", 400 ))

tandfmetricsquad <- as.data.frame(cbind(unlist(tansvmquadloss), unlist(tanknnquadloss), unlist(tanres$quadloss), unlist(tancvalmcmc$quadloss) ))
colnames(tandfmetricsquad) <- c("SVM", "KNN", "MAP", "MCMC")
tandfmetricsquad <- melt(tandfmetricsquad)
tandfmetricsquad$set <- factor(rep("Drosophila", 400 ))

halldfmetricsquad <- as.data.frame(cbind(unlist(hallsvmquadloss), unlist(hallknnquadloss), unlist(hallres$quadloss), unlist(hallcvalmcmc$quadloss)))
colnames(halldfmetricsquad) <- c("SVM", "KNN", "MAP", "MCMC")
halldfmetricsquad <- melt(halldfmetricsquad)
halldfmetricsquad$set <- factor(rep("Chicken DT40", 400 ))

dfmetrics2 <- rbind(tandfmetricsquad, halldfmetricsquad, andydfmetricsquad, itzhakdfmetricsquad, u2osdfmetricsquad)

gg2 <- ggplot(data = dfmetrics2, aes(x=variable, y=value)) +
  geom_boxplot(aes(fill=variable)) +
  scale_y_continuous() +
  ggtitle("Boxplot of Quadratic Losses") + 
  theme_bw() + labs(y = "Quadratic Loss") +  labs(x  = "Classifier") + facet_grid(. ~ set) + scale_x_discrete(labels = element_blank()) + scale_fill_discrete(name = "Classifier")
gg2



#ttest for f1 score and quadratic losses
ttesttanquad <- pairwise.t.test(tandfmetricsquad$value, tandfmetricsquad$variable , p.adjust.method = "BH",
                              pool.sd = FALSE, paired = FALSE,
                              alternative = "two.sided")
ttesttanquad$p.value
require(xtable)
xtable(ttesttanquad$p.value, digits = 5,
       caption = "Adjusted P-values for pairwise T-tests for Quadratic Loss classifier evaluation on the drosphila dataset")

ttesthallquad <- pairwise.t.test(halldfmetricsquad$value, halldfmetricsquad$variable , p.adjust.method = "BH",
                               pool.sd = FALSE, paired = FALSE,
                               alternative = "two.sided", mu = 0, var.equal = FALSE)
ttesthallquad$p.value
require(xtable)
xtable(ttesthallquad$p.value, digits = 5,
       caption = "Adjusted P-values for pairwise T-tests for Macro F-1 score classifier evaluation on the Chicken DT40 dataset")

ttestandyquad <- pairwise.t.test(andydfmetricsquad$value, andydfmetricsquad$variable , p.adjust.method = "BH",
                               pool.sd = FALSE, paired = FALSE,
                               alternative = "less")
ttestandyquad$p.value
require(xtable)
xtable(ttestandyquad$p.value, digits = 5,
       caption = "Adjusted P-values for pairwise T-tests for Macro F-1 score classifier evaluation on the mouse dataset")

ttestitzhakquad <- pairwise.t.test(itzhakdfmetricsquad$value, itzhakdfmetricsquad$variable , p.adjust.method = "BH",
                                 pool.sd = F, paired = FALSE,
                                 alternative = "two.sided")
ttestitzhakquad$p.value
require(xtable)
xtable(ttestitzhakquad$p.value, digits = 5,
       caption = "Adjusted P-values for pairwise T-tests for Macro F-1 score classifier evaluation on the HeLa dataset")

ttestu2osquad <- pairwise.t.test(u2osdfmetricsquad$value, u2osdfmetricsquad$variable , p.adjust.method = "BH",
                               pool.sd = F, paired = FALSE,
                               alternative = "two.sided")
ttestu2osquad$p.value
require(xtable)
xtable(ttestu2osquad$p.value, digits = 5,
       caption = "Adjusted P-values for pairwise T-tests for Macro F-1 score classifier evaluation on the U2-OS dataset")

