#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
double EPPF_Dirichlet(IntegerVector counts, double alpha){

  // Initialization
  int k = counts.size();
  int n = sum(counts);

  double logLik = k*log(alpha) + lgamma(alpha) - lgamma(alpha+n);
  return logLik;
}

// [[Rcpp::export]]
double EPPF_PitmanYor(IntegerVector counts, double alpha, double sigma){

  // Initialization
  int k = counts.size();
  int n = sum(counts);

  if(sigma<0){
    double Inf = arma::datum::inf;
    return(-Inf);

  } else if (alpha <= -sigma){
    double Inf = arma::datum::inf;
    return(-Inf);

  } else {
    double prod_numerator = 0;
    for(int i=1; i < k; i++){
      prod_numerator += log(alpha + i * sigma);
    }

    double product_sigma = 0;
    for(int i=0; i < k; i++){
      product_sigma += lgamma(counts[i] - sigma) - lgamma(1 - sigma);
    }

    double logLik = prod_numerator -  lgamma(alpha + n) + lgamma(alpha + 1) + product_sigma;
    return logLik;
  }
}

///*** R
//species_counts = c(10, 30,12, 22, 10,15, 1,1,1,3,4,5)
//EPPF_Dirichlet(counts = c(10, 30,12, 22, 10,15, 1,1,1,3,4,5), alpha = 3)
//EPPF_PitmanYor(counts =K, alpha = 3, sigma = 0.5)

//compute_Loglik_PitmanYor(data = list("n" = sum(K), "k" = length(K), "K"  = K), pars = c(3, 0.5))
////sum(lgamma(K))
////*/
