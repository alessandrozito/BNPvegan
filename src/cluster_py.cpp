#include "RcppArmadillo.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double cluster_py_C(int n, double sigma, double theta){

  // Initialization
  double K_n = 1;

  // Cycle over the samples
  for(int i=1; i < n; i++) {
    K_n = K_n + R::rbinom(1, (theta + sigma*K_n)/(theta + i));
  }
  return K_n;
}
