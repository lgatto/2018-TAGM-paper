source("https://bioconductor.org/biocLite.R")
biocLite()

require(MSnbase)
require(pRoloc)
require(pRolocdata)
require(compiler)
require(Rcpp)
require(sampling)
require(caret)

source("bayesgmmOptimisationGibbs.R")
source("gibbsMetric.R")
sourceCpp("dmvtCpp.cpp")
source("bayesgmmgibbsCval.R")

cmp_Gibbs <- cmpfun(bayesgmmOptimisationGibbs)

data("itzhak2016stcSILAC")
itzhak <- itzhak2016stcSILAC

itzhak <- filterNA(itzhak)
fData(itzhak)$markers[fData(itzhak)$markers == "ER_high curvature"] <- "unknown"
fData(itzhak)$markers[fData(itzhak)$markers == "Large Protein Complex"] <- "unknown"

itzhakcvallpcmcmc <- bayesgmmgibbsCval(object = itzhak,
                                  times = 100,
                                  test.size = 0.2,
                                  iterations = 10000L,
                                  burnin = 1000L,
                                  thin = 5L,
                                  BPPARAM = bpparam())

save(itzhakcvallpcmcmc, file = "itzhakcvallpcmcmc-01.Rda")

