#' Sample from the posterior distribution of the asymptotic species richness in a species discovery model
#'
#' @param object An object of class \code{sdm}.
#' @param n_samples Number of random samples to draw.
#' @param tolerance Tolerance to determine the truncation point. Default \code{= 1e-7}.
#' @param ... Additional model parameters
#'
#' @return A vector of size \code{n_samples}
#' @details The function sample from the posterior distribution of the asymptotic
#'          species richness under a species discovery model. Specifyng the tolerance determines the truncation point.
#' @export
sample_Kinf <- function(object, n_samples = 100, tolerance = 1e-7,...){
  # Step 1 - Determine the truncation point to stop the sampler
  n <- length(object$discoveries)
  k <- sum(object$discoveries)
  if(object$model == "LL3"){
    alpha <- object$par[1]
    sigma <- object$par[2]
    phi <- object$par[3]
    # Call the Cpp function to sample from Kinf under LL3
    Kinf <- c(sample_Kinf_LL3_Cpp(n_samples = n_samples, n = n, k = k, alpha = alpha,
                                  sigma = sigma, phi = phi, tolerance = tolerance))
  } else if (object$model == "Weibull"){
    phi <- object$par[1]
    lambda <- object$par[2]
    # Call the Cpp function to sample from Kinf under Weibull
    Kinf <- c(sample_Kinf_Weibull_Cpp(n_samples = n_samples, n = n, k = k,
                                      phi = phi, lambda = lambda, tolerance = tolerance))
  }
  return(Kinf)
}
