## Description of the data

### Data concerning application of methods to mouse pluripotent stem cell data (Christoferou et al 2016).

 .rda files can be found with both raw and processed data. All files are in the form andy[].rda
 
 andy2015cvalmcmc-01.rda contains lists with the confusion matricies and quadratic losses produced from cross-validation using TAGM-MCMC.
 
 andydiagtest is to reproduce convergence diagnostics for TAGM-MCMC. The data contains 6 parallel chains (andydiag$andymcmc[1-6]) for 4106 protein each with 1100 samples. The proteins have unknown locations and value is a sample  from a binary value indicates whether the protein was allocated to the outlier component or not. WARNING this file will load 200MB into working memory.
 
 andyknn(qloss).rda is cross-validation results using KNN algorithm for F1 scores (quadratic loss)
 
 andysvm(qloss).rda is cross-validation results using SVM algorithm for F1 scores (quadratic loss)
 
 andyMAPparams.rda are the MAP estimates for the TAGM model along with the loglikelihood at each iteration.
 
 andyMAPres.rda prediction results from applying TAGM-MAP have been appended to the fData in an MSnSet. Full prediction results available in matricies
 
 andyMCMChl.rda prediction results from applyhing TAGM-MAP have been appened to the fData in an MSnSet. These include outlier allocations and shannon entropies as well as prediction probabilities and classes.
 
 andyoutliers a numeric matrix of outlier allocations
 
 ## Cross-validation results
 
 The following naming convention is used [dataset][method]
 
 e.g. BeltranH24svmquadloss.rda provides cross-validation results from the Beltran et al. 2016 data with infection at 24 using SVM with quadratic loss. If the quadloss label is dropped, the file contains F1 scores.
 
 e.g hirstc6knn provides cross-validation F1 scores for the Hirst et al 2018 data using the KNN algorithm.
 
 for convenience additionalF1compare.rda contains F1 cross-validation results for the following datasets 
 
 HeLa Wild (Hirst et al. 2018)   HeLa KO1 (Hirst et al. 2018)    HeLa KO2 (Hirst et al. 2018)   Primary Fibroblast Mock 24hpi   Primary Fibroblast Mock 48hpi   Primary Fibroblast Mock 72hpi  Primary Fibroblast Mock 96hpi   Primary Fibroblast Mock 120hpi  Primary Fibroblast HCMV 24hpi  Primary Fibroblast HCMV 48hpi   Primary Fibroblast HCMV 72hpi   Primary Fibroblast HCMV 96hpi  Primary Fibroblast HCMV 120hpi  E14TG2aR (Breckels et al. 2016)
 
 The first column in additionalF1compare contains the method, the second column the F1 score, the final column the dataset
 
 additionalquadlosscompare.rda contains quadratic loss cross-validation results for the same datasets in the same format with F1 replaced with quadratic loss
 
 For the U2-OS, drosophila (Tan), HeLa (Itzhak), chicken (hall) datasets the naming includes the names of the dataset and the names
 of the method. For the SVM and KNN the quadloss results are included with the name tag qloss.
