#include "RcppArmadillo.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec rarefy_C(arma::vec freq, int n, int K, bool verbose) {
  arma::vec out(n); // Allocate a vector of dimension n
  double n1 =  n;
  double K1 = K;
  int R_show = 10000;
  double ldiv = lgamma(n1 + 1);
  for (int i = 0; i < n; i++) {
     double i1 = i + 1.0;
     arma::vec value = lgamma(n1 - freq + 1) - lgamma(n1 - freq - i1 + 1) - ldiv + lgamma(n1 - i1 + 1);
     out(i) = K1 - sum(exp(value));
     if((i+1)%R_show==0 and verbose == TRUE) {Rprintf("Iteration: %i \n", i+1);}
  }
  return out;
}
