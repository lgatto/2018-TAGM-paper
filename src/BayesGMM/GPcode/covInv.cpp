# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat covInvcpp(arma::mat Z,
                    arma::mat A,
                    double sigmak,
                    double nk,
                    double D){
  
  arma::mat covInv;
  arma::mat J = arma::ones(nk, nk);
  
  covInv =  arma::eye(nk * D, nk * D) / sigmak - arma::kron(J, Z * A) / (pow(sigmak, 2)) ;
  
  return(covInv); 
}
