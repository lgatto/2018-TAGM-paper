# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


arma::vec gradienthyp(arma::vec Xk,
                      arma::mat covInv,
                      arma::mat drv,
                      double nk){
  
  arma::mat B;
  arma::mat J = arma::ones(nk, nk);
  arma::mat C;
  arma::vec gradienthyp; 
  
  B = (((Xk.t() * covInv) * arma::kron(J, drv)) * (covInv * Xk)) / 2 ;
  C = arma::accu(arma::kron(J, drv) % covInv)/2; //use trace is sum of entrywise product
  
  gradienthyp = B - C;  
  
  return(gradienthyp); 
}

arma::vec gradientNoise(arma::vec Xk,
                        arma::mat covInv,
                        double drv,
                        double nk){
  
  arma::mat B;
  arma::mat J = arma::ones(nk, nk);
  arma::mat C;
  arma::vec gradientNoise; 
  
  B = drv * ((Xk.t() * covInv) * (covInv * Xk))/2 ;
  C = drv * trace(covInv)/2; //use trace is sum of entrywise product
  
  gradientNoise = B - C;  
  
  return(gradientNoise); 
}

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

List trenchDetcpp(arma::vec c) {
  
  int N = c.size();
  int i;
  arma::vec xi(N-1);
  arma::vec z(N-1);
  arma::vec v(N);
  double beta = 1;
  double alpha;
  double l = 0;
  double logdet;
  
  xi = c.subvec(1, N-1)/c(0);
  z(0) = - xi(0);
  alpha = -xi(0);
  int K = xi.size(); 
  
  for(i = 0; i < (K - 1); i++){
    beta = (1 - pow(alpha, 2)) * beta ;
    l = l + log(beta) ;
    if(i == 0){
      alpha = - (xi(i + 1) + xi(0) * z(0)) / beta ;
      z(0) = z(0) + alpha * z(0) ;
    }else{
      alpha = - (xi(i + 1) + dot(flipud(xi.subvec(0, i)), z.subvec(0, i)) ) / beta ;
      z.subvec(0, i) = z.subvec(0, i) + alpha * flipud(z.subvec(0, i)) ;
    }
    
    z(i+1) = alpha  ;
  }
  
  beta = (1 - pow(alpha, 2)) * beta ;
  l = l + log(beta) ;
  
  logdet = l + N * log(c(0));
  
  v(N-1) = 1 / ((1 + dot(xi, z)) * c(0)) ;
  v.subvec(0,N-2) = v(N-1) * flipud(z.subvec(0, N - 2));
  
  return List::create(Rcpp::Named("logdet") = logdet,
                      Rcpp::Named("z") = z,
                      Rcpp::Named("v") = v);
}

arma::mat trenchInvcpp(arma::vec v) {
  
  int N = v.size();
  int i;
  int j;
  arma::mat C(N, N);
  arma::mat trenchInv;
  
  
  C.row(0) = flipud(v).t();
  C.col(0) = flipud(v);
  C.row(N - 1) = v.t();
  C.col(N - 1) = v;
  for(i = 1; i < floor( (N - 1) / 2 ) + 1; i++){
    for(j = 1; j < N - i; j++){
      C(i, j) = C(i - 1, j - 1) + (v(N - 1 - j) * v(N - 1 - i) - v(i - 1) * v(j - 1)) / v(N - 1) ;
      C(j, i) = C(i, j);
      C(N - 1 - i , N - 1 - j ) = C(i, j);
      C(N - 1 - j , N - 1 - i ) = C(i, j);
    }
  } 
  
  trenchInv = C;
  return(trenchInv) ; 
  
}



arma::vec gradientGPcpp(arma::vec Xk,
                        arma::vec tau,
                        arma::vec h,
                        int nk,
                        int D
){
  
  //terms are made unconstrained
  double sigmak; sigmak = exp(2 * h(2));
  double a; a = exp(2 * h(1));
  double l; l = exp(h(0));
  double negloglike;
  double logdet;
  arma::mat Z(D,D);
  arma::mat S(D,D); S = arma::zeros(D,D);
  arma::mat A(D,D);
  arma::mat R(D,D);
  arma::mat covInv;
  arma::mat drvl;
  arma::mat drva;
  double drvsigma;
  arma::vec temp;
  arma::vec v;
  arma::vec grdl;
  arma::vec grda;
  arma::vec grdsigma;
  arma::vec grad(3);
  List trenchres = List::create(Rcpp::Named("logdet"),
                                Rcpp::Named("z"),
                                Rcpp::Named("v"));
  
  S = S.each_col() + tau;
  A = a * exp(- pow(S - S.t(), 2)/l);
  R = arma::eye(D,D) + (nk * A)/sigmak;
  
  trenchres = trenchDetcpp(R.row(0).t());
  logdet = trenchres["logdet"];
  v = as<arma::vec>(trenchres["v"]);
  
  
  if(R_IsNA(logdet)){
    negloglike = R_NegInf;
    grad.row(0) = negloglike;
    return(grad);
  }
  
  Z = trenchInvcpp(v);
  covInv = covInvcpp(Z, A, sigmak, nk, D);
  
  drvl = A % (pow(S - S.t(), 2)/l);
  drva = 2 * A ;
  drvsigma = 2 * sigmak ;
  
  grdl = gradienthyp(Xk, covInv, drvl, nk);
  grda = gradienthyp(Xk, covInv, drva, nk);
  grdsigma = gradientNoise(Xk, covInv, drvsigma, nk);
  
  //return gradient of negative log like
  grad(0) = -grdl(0);
  grad(1) = -grda(0);
  grad(2) = -grdsigma(0);
  
  return(grad);
}

// [[Rcpp::export]]
List LeapfrogGPcpp(arma::vec Xk,
                        arma::vec tau,
                        arma::vec p,
                        arma::vec x,
                        arma::vec m,
                        int nk,
                        int D,
                        int L,
                        double delta){
  
  arma::vec gradprior(3);
  gradprior(0) = x(0); gradprior(1) = 2*x(1); gradprior(2) = 2*x(2);

for (int t=0; t < L; t++){
  //half step for momenta
  p = p - delta * (gradientGPcpp(Xk, tau, x , nk, D) + gradprior)/2 ;
  //full step for position
  x = x + delta * p / m  ;
  //another half step for momenta
  p = p - delta * (gradientGPcpp(Xk, tau, x , nk, D) + gradprior)/2 ;
} 


return List::create(Rcpp::Named("p") = p,
                    Rcpp::Named("x") = x
                    );
  
}  
  