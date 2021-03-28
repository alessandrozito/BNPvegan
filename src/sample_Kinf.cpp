#include "RcppArmadillo.h"
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double prob_LL3_Cpp(double n, double alpha, double sigma, double phi){
  double num = alpha * pow(phi, n);
  double prob = num/(num + pow(n, 1-sigma));
  return prob;
}

// [[Rcpp::export]]
double prob_Weibull_Cpp(double n, double phi, double lambda){
  double prob = pow(phi, pow(n, lambda));
  return prob;
}

/// Functions to compute the truncation point

// [[Rcpp::export]]
double truncationpoint_LL3_Cpp(double n, double alpha, double sigma, double phi, double tolerance){
  double diff = prob_LL3_Cpp(n, alpha, sigma, phi);
  double m = 0;
  while(tolerance < diff){
    m += 1;
    diff = prob_LL3_Cpp(n + m, alpha, sigma, phi);
  }
  return m;
}

// [[Rcpp::export]]
double truncationpoint_Weibull_Cpp(double n, double phi, double lambda, double tolerance){
  double diff = prob_Weibull_Cpp(n, phi, lambda);
  double m = 0;
  while(tolerance < diff){
    m += 1;
    diff = prob_Weibull_Cpp(n + m, phi, lambda);
  }
  return m;
}


/// Functions to sample from the posterior of Kinf
// [[Rcpp::export]]
arma::vec sample_Kinf_LL3_Cpp(int n_samples, double n, double k, double alpha, double sigma, double phi, double tolerance){

  // Step 0 - Obtain the optimal m to truncate the summation
  int m = truncationpoint_LL3_Cpp(n, alpha, sigma, phi, tolerance);

  // Step 1 - Compute the vector of probabilities
  arma::vec probs(m);
  for(int i=0; i<m; i++){
    probs(i) = prob_LL3_Cpp(n+i, alpha, sigma, phi);
  }

  arma::vec Kinf(n_samples);
  // Step 2 - Draw the random samples
  for(int s=0; s<n_samples; s++){
    double d = k;
    for(int j = 0; j < m; j++){
       d += R::rbinom(1, probs(j));
    }
    Kinf(s) = d;
  }
  return Kinf;
}

// [[Rcpp::export]]
arma::vec sample_Kinf_Weibull_Cpp(int n_samples, double n, double k, double phi, double lambda, double tolerance){

  // Step 0 - Obtain the optimal m to truncate the summation
  int m = truncationpoint_Weibull_Cpp(n, phi, lambda, tolerance);

  // Step 1 - Compute the vector of probabilities
  arma::vec probs(m);
  for(int i=0; i<m; i++){
    probs(i) = prob_Weibull_Cpp(n+i, phi, lambda);
  }

  arma::vec Kinf(n_samples);
  // Step 2 - Draw the random samples
  for(int s=0; s<n_samples; s++){
    double d = k;
    for(int j = 0; j < m; j++){
      d += R::rbinom(1, probs(j));
    }
    Kinf(s) = d;
  }
  return Kinf;
}



// [[Rcpp::export]]
double get_n_target_saturation_LL3(double n, double k, double alpha, double sigma, double phi, double Kinf, double target){
  double m = 0;
  while(k/Kinf < target){
    m = m + 1;
    k = k + prob_LL3_Cpp(n+m,alpha, sigma, phi);
  }
  return m;
}

// [[Rcpp::export]]
double get_n_target_saturation_Weibull(double n, double k, double phi, double lambda, double Kinf, double target){
  double m = 0;
  while(k/Kinf < target){
    m = m + 1;
    k = k + prob_Weibull_Cpp(n+m, phi,lambda);
  }
  return m;
}


