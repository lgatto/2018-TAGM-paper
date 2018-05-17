# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec loglikeGPcpp(arma::vec Xk,
                    arma::mat Z,
                    arma::mat A,
                    double logcovDet,
                    double sigmak,
                    double nk,
                    double D){
  
  arma::vec loglike;
  arma::mat J = arma::ones(nk, nk);
  
  loglike = - dot(Xk,Xk) / (2 * sigmak) + (Xk.t() * arma::kron(J, Z * A) * Xk) / (2 * pow(sigmak, 2)) 
    -  logcovDet / 2 - (nk * D *  log(2 * M_PI)) / 2 ;

 return(loglike); 
}



