#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec logEPPF_DM_cpp(double sigma, double n, double K,  arma::vec &counts, arma::vec &H){
  arma::vec theta = - sigma*H;
  IntegerVector K_seq = seq(1, K-1);
  arma::vec Kall = as<arma::vec>(K_seq);

  int totH = H.size();
  double logLik_noH = sum(lgamma(counts - sigma)) - K * lgamma(1 - sigma);
  arma::vec logLik(totH);

  for(int i=0;i<totH; i++){
    logLik(i) =  sum(log(theta(i) + Kall * sigma)) - lgamma(theta(i) + n) +
      lgamma(theta(i) + 1) + logLik_noH;
  }
  return logLik;
}

// [[Rcpp::export]]
arma::vec logNegBin_cpp(arma::vec x, double K, double m, double r){
  int n = x.size();
  arma::vec logLik_NegBin(n);
  logLik_NegBin = lgamma(r+x-K) - lgamma(x-K+1) - lgamma(r) +
    r*(log(r) - log(r+m)) + (x-K)*(log(m) - log(m+r));
  //for(int i=0; i<n; i++){
  //  logLik_NegBin(i) = lgamma(r+x(i)-K) - lgamma(x(i)-K+1) - lgamma(r) + r*(log(r) -
  //    log(r+m)) + (x(i)-K)*(log(m) - log(m+r));
  //}
  return logLik_NegBin;
}

// [[Rcpp::export]]
double sampleH_cpp(arma::vec &counts, double sigma, double m, double r, double K, double n, double H_max){
  IntegerVector Hsome = seq(K, K+H_max);
  arma::vec H_seq = as<arma::vec>(Hsome);
  int Hsize = H_seq.size();

  arma::vec H_probs(Hsize);
  H_probs = logEPPF_DM_cpp(sigma, n, K,  counts, H_seq) + logNegBin_cpp(H_seq, K, m, r);
  arma::vec p = exp(H_probs - max(H_probs));
  H_probs = p/sum(p);
  bool replace = false;
  arma::vec H_out;
  H_out = RcppArmadillo:: sample(H_seq, 1, replace, H_probs);
  return H_out(0);
}

// [[Rcpp::export]]
arma::vec sample_phi_cpp(arma::vec &counts, int iter, double phi_old, double H,
                         double n, double K, double var_samples, double sds){
  arma::vec Hv(1);
  Hv(0) = H;
  NumericVector phioldv(1);
  phioldv(0) = phi_old;
  arma::vec out(2);

  // Propose the value via SCAM as in Haario et al 2005
  NumericVector phi;
  if(iter <= 10){
    phi = rnorm(1, phi_old, 5);
  } else {
    phi = rnorm(1, phi_old, 2.4 * sqrt(var_samples + 0.05));
  }

  double logratio;
  NumericVector prior_new = dnorm(phi, 0.0, sds, true);

  NumericVector prior_old = dnorm(phioldv, 0.0, sds, true);

  logratio = logEPPF_DM_cpp(-exp(phi(0)), n, K, counts, Hv)(0) -
    logEPPF_DM_cpp(-exp(phi_old), n, K, counts, Hv)(0) + prior_new(0) - prior_old(0);

  double u = runif(1)(0);
  if(log(u) < logratio){
    out(0) =phi(0);
    out(1) = 1;
  } else {
    out(0) = phi_old;
    out(1) = 0;
  }
  return out;
}

// [[Rcpp::export]]
NumericMatrix Gibbs_DirMult(arma::vec & counts, int R, int burnin,
                                double m, double r, double sds, double H_max,
                                bool verbose = true){
  NumericMatrix Gibbs_out(R+burnin, 3);
  int verbose_step = (R+burnin)/10 ;

  // Starting values
  double H = counts.size() + 1.0;
  double phi = 0.0;
  double K = counts.size() + 0.0;
  double n = sum(counts);
  arma::vec outphi(2);
  double var_samples = 0.0;
  //double i1; double tot;

  for(int i=0; i<R+burnin; i++){
    // Monitor the state of the sampler
    if(i%verbose_step ==0){
      Rprintf("Samling iteration: %i \n", i);
    }

    // Sample H
    H = sampleH_cpp(counts, -exp(phi), m, r, K, n,H_max);
    // Sample one phi (note that sigma = -exp(phi))
    if(i>9){
      //arma::vec s;
      NumericVector s;
      //s = Gibbs_out.rows(0, i-1).col(0);
      s = Gibbs_out.column(0);
      var_samples =  var(s[Range(0, i-1)]);
    }
    outphi = sample_phi_cpp(counts,i, phi, H, n,K,var_samples,sds);
    phi = outphi(0);

    // Store values
    Gibbs_out(i,0) = phi;
    Gibbs_out(i,1) = H;
    Gibbs_out(i,2) = outphi(1);
  }
  Gibbs_out.column(0) = -exp(Gibbs_out.column(0));

  return Gibbs_out(Range(burnin, burnin+R-1), _);

}






