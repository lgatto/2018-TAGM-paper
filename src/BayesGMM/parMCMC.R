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
sourceCpp("dmvtCpp.cpp")

cmp_Gibbs <- cmpfun(bayesgmmOptimisationGibbs)


data("hyperLOPIT2015")
hl <- hyperLOPIT2015

andy2015mcmc <- bplapply(1:4, cmp_Gibbs, object = hl,
                                      iterations = 10000L,
                                      burnin = 1000L,
                                      thin = 5L,
                                      BPPARAM = MulticoreParam(workers = 4))

save(andy2015mcmc, file = "andy2015mcmc-01.Rda")