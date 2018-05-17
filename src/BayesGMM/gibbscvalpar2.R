source("https://bioconductor.org/biocLite.R")
biocLite()
require(devtools)
require(MSnbase)
require(pRoloc)
require(pRolocdata)
require(compiler)
require(Rcpp)
require(sampling)
require(caret)
install_github("lgatto/pRolocdata")

source("bayesgmmOptimisationGibbs.R")
source("gibbsMetric.R")
sourceCpp("dmvtCpp.cpp")
source("bayesgmmgibbsCval.R")

cmp_Gibbs <- cmpfun(bayesgmmOptimisationGibbs)

data("hirst2018")
hirst <- hirst2018

fData(hirst)$markers[fData(hirst)$markers == "ER_Tubular"] <- "unknown"
fData(hirst)$markers[fData(hirst)$markers == "Nuclear pore complex"] <- "unknown"
fData(hirst)$markers[fData(hirst)$markers == "Peroxisome"] <- "unknown"

hirstcon <- hirst[,c(11,12,13,14,15,26,27,28,29,30,41,42,43,44,45)]
hirstc2 <- hirst[,c(1,2,3,4,5,16,17,18,19,20,31,32,33,34,35)]
hirstc6 <- hirst[,c(1,2,3,4,5,16,17,18,19,20,31,32,33,34,35)+5]



'hirstconcvalmcmc <- bayesgmmgibbsCval(object = hirstcon,
                                       times = 100,
                                       test.size = 0.2,
                                       iterations = 10000L,
                                       burnin = 1000L,
                                       thin = 5L,
                                       BPPARAM = bpparam())

save(hirstconcvalmcmc, file = "hirstconcvalmcmc-01.Rda")

hirstc2cvalmcmc <- bayesgmmgibbsCval(object = hirstc2,
                                      times = 100,
                                      test.size = 0.2,
                                      iterations = 10000L,
                                      burnin = 1000L,
                                      thin = 5L,
                                      BPPARAM = bpparam())

save(hirstc2cvalmcmc, file = "hirstc2cvalmcmc-01.Rda")

hirstc6cvalmcmc <- bayesgmmgibbsCval(object = hirstc6,
                                     times = 100,
                                     test.size = 0.2,
                                     iterations = 10000L,
                                     burnin = 1000L,
                                     thin = 5L,
                                     BPPARAM = bpparam())

save(hirstc6cvalmcmc, file = "hirstc6cvalmcmc-01.Rda")'



'data("beltran2016MOCK24")
beltranM24cvalmcmc <- bayesgmmgibbsCval(object = beltran2016MOCK24,
                                        times = 100,
                                        test.size = 0.2,
                                        iterations = 10000L,
                                        burnin = 100L,
                                        thin = 5L,
                                        BPPARAM = bpparam()) 

save(beltranM24cvalmcmc, file = "beltranM24cvalmcmc-01.Rda")


data("beltran2016MOCK48")
beltranM48cvalmcmc <- bayesgmmgibbsCval(object = beltran2016MOCK48,
                                        times = 100,
                                        test.size = 0.2,
                                        iterations = 10000L,
                                        burnin = 1000L,
                                        thin = 5L,
                                        BPPARAM = bpparam()) 

save(beltranM48cvalmcmc, file = "beltranM48cvalmcmc-01.Rda")

data("beltran2016MOCK72")
beltranM72cvalmcmc <- bayesgmmgibbsCval(object = beltran2016MOCK72,
                                        times = 100,
                                        test.size = 0.2,
                                        iterations = 10000L,
                                        burnin = 1000L,
                                        thin = 5L,
                                        BPPARAM = bpparam()) 

save(beltranM72cvalmcmc, file = "beltranM72cvalmcmc-01.Rda")

data("beltran2016MOCK96")
beltranM96cvalmcmc <- bayesgmmgibbsCval(object = beltran2016MOCK96,
                                        times = 100,
                                        test.size = 0.2,
                                        iterations = 10000L,
                                        burnin = 1000L,
                                        thin = 5L,
                                        BPPARAM = bpparam()) 

save(beltranM96cvalmcmc, file = "beltranM96cvalmcmc-01.Rda")

data("beltran2016MOCK120")
beltranM120cvalmcmc <- bayesgmmgibbsCval(object = beltran2016MOCK120,
                                        times = 100,
                                        test.size = 0.2,
                                        iterations = 10000L,
                                        burnin = 1000L,
                                        thin = 5L,
                                        BPPARAM = bpparam()) 

save(beltranM120cvalmcmc, file = "beltranM120cvalmcmc-01.Rda")


data("beltran2016HCMV24")
beltranH24cvalmcmc <- bayesgmmgibbsCval(object = beltran2016HCMV24,
                                         times = 100,
                                         test.size = 0.2,
                                         iterations = 10000L,
                                         burnin = 1000L,
                                         thin = 5L,
                                         BPPARAM = bpparam()) 

save(beltranH24cvalmcmc, file = "beltranH24cvalmcmc-01.Rda")

data("beltran2016HCMV48")
beltranH48cvalmcmc <- bayesgmmgibbsCval(object = beltran2016HCMV48,
                                        times = 100,
                                        test.size = 0.2,
                                        iterations = 10000L,
                                        burnin = 1000L,
                                        thin = 5L,
                                        BPPARAM = bpparam()) 

save(beltranH48cvalmcmc, file = "beltranH48cvalmcmc-01.Rda")

data("beltran2016HCMV72")
beltranH72cvalmcmc <- bayesgmmgibbsCval(object = beltran2016HCMV72,
                                        times = 100,
                                        test.size = 0.2,
                                        iterations = 10000L,
                                        burnin = 1000L,
                                        thin = 5L,
                                        BPPARAM = bpparam()) 

save(beltranH72cvalmcmc, file = "beltranH72cvalmcmc-01.Rda")

data("beltran2016HCMV96")
beltranH96cvalmcmc <- bayesgmmgibbsCval(object = beltran2016HCMV96,
                                        times = 100,
                                        test.size = 0.2,
                                        iterations = 10000L,
                                        burnin = 1000L,
                                        thin = 5L,
                                        BPPARAM = bpparam()) 

save(beltranH96cvalmcmc, file = "beltranH96cvalmcmc-01.Rda")

data("beltran2016HCMV120")
beltranH120cvalmcmc <- bayesgmmgibbsCval(object = beltran2016HCMV120,
                                        times = 100,
                                        test.size = 0.2,
                                        iterations = 10000L,
                                        burnin = 1000L,
                                        thin = 5L,
                                        BPPARAM = bpparam()) 

save(beltranH120cvalmcmc, file = "beltranH120cvalmcmc-01.Rda")'

'data("HEK293T2011")
fData(HEK293T2011)$markers <- fData(HEK293T2011)$markers.tl 

HEKcvalmcmc <- bayesgmmgibbsCval(object = HEK293T2011,
                                         times = 1,
                                         test.size = 0.2,
                                         iterations = 1000L,
                                         burnin = 100L,
                                         thin = 5L,
                                         BPPARAM = bpparam()) 

save(HEKcvalmcmc, file = "HEKcvalmcmc-01.Rda")'


data("groen2014cmb")

groencvalmcmc <- bayesgmmgibbsCval(object = groen2014cmb,
                                 times = 100,
                                 test.size = 0.2,
                                 iterations = 10000L,
                                 burnin = 1000L,
                                 thin = 5L,
                                 BPPARAM = bpparam()) 

save(groencvalmcmc, file = "groencvalmcmc-01.Rda")


data("E14TG2aR")

E14TG2aRcvalmcmc <- bayesgmmgibbsCval(object = E14TG2aR,
                                   times = 100,
                                   test.size = 0.2,
                                   iterations = 10000L,
                                   burnin = 1000L,
                                   thin = 5L,
                                   BPPARAM = bpparam()) 

save(E14TG2aRcvalmcmc, file = "E14TG2aRcvalmcmc-01.Rda")

