source("https://bioconductor.org/biocLite.R")
biocLite()

require(MSnbase)
require(pRoloc)
require(pRolocdata)
install.packages("lbfgs")
require(lbfgs)
require(Rcpp)
require(RcppArmadillo)


sourceCpp("dmvtCpp.cpp")
source("gradientGP.R")
sourceCpp("gradienthyp.cpp")
source("hamiltonianGP.R")
sourceCpp("leapfrogGPcpp.cpp")
source("likelihoodGP.R")
source("metropolisGP.R")
source("mixGP.R")
sourceCpp("trenchDetcpp.cpp")

data("hyperLOPIT2015")

andyGP <- mixGP(mydata = hyperLOPIT2015[,1:10], hypLearn = "HMC", numIter = 20000, hypIter = 50)

save(andyGP, file = "andyGP.rda")