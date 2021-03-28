#' Sample-based or model-based Gini heterogeniety index
#'
#' @param object An object of class \code{numeric}, \code{ssm} or \code{sdm}
#' @param ... Additional parameters
#'
#' @export
Gini <- function(object, ...) {
  UseMethod("Gini", object)
}

#' Gini heterogeneity index for a vector pf species abundances
#' @param object An object of class \code{numeric}
#' @param ... Additional parameters
#'
#' @export
Gini.numeric <- function(object, ...) {
  freq_rel <- object / sum(object)
  out <- 1 - sum(freq_rel^2)
  out
}

#' Posterior Gini heterogeneity index for the Dirichlet process.
#' @param object An object of class \code{ssm, DP}.
#' @param ... Additional parameters
#'
#' @export
Gini.DP <- function(object, ...) {
  Poch2 <- function(x) x * (x + 1)

  alpha <- object$param[1]
  freq <- object$abundances
  n <- sum(freq)
  out <- 1 - 1 / Poch2(alpha + n) * (alpha + sum(Poch2(freq)))
  out
}

#' Posterior Gini heterogeneity index for the Pitman-Yor process
#' @param object An object of class \code{ssm, PY}
#' @param ... Additional parameters
#'
#' @export
Gini.PY <- function(object, ...) {
  Poch2 <- function(x) x * (x + 1)

  alpha <- object$param[1]
  sigma <- object$param[2]
  freq <- object$abundances
  n <- sum(freq)
  K <- length(freq)

  out <- 1 - 1 / Poch2(alpha + n) * ((1 - sigma) * (alpha + K * sigma) + sum(Poch2(freq - sigma)))
  out
}

######## Posterior samples for the Gini coefficient

#' Random sampler for the posterior distribution of the Gini heterogeneity index
#'
#' @param object An object of class \code{ssm}
#' @param ... Additional parameters
#'
#' @export
rGini <- function(object, ...) {
  UseMethod("rGini", object)
}

#' Random sampler for the posterior distribution of the Gini heterogeneity index under a Pitman-Yor process
#'
#' @param object An object of class \code{ssm, PY}
#' @param n_samples Number of posterior samples to draw
#' @param truncation_method Truncation method for the stick breaking sampler
#' @param prob_cutoff Probabilistic truncation cutoff for the stick breaking sampler. Valid only for \code{truncation_method = "dynamic"}
#' @param static_cutoff Integer cutoff for the stick breaking sampler. Valid only for \code{truncation_method = "static"}
#' @param ... Additional model parameters
#'
#' @export
rGini.PY <- function(object, n_samples = 100, truncation_method = "dynamic", prob_cutoff = 0.99, static_cutoff = 1000, ...){
  # Sample the weights from a posterior distribution of the Pitman-Yor given the abundances
  n <- sum(object$abundances)
  k <- length(object$abundances)
  alpha <- object$param[1]
  sigma <- object$param[2]

  Gini_out <- rep(0, n_samples)
  ## STEP 1 - random samples from a Dirichlet distribution
  alpha_dir <- c(object$abundances - sigma, alpha + sigma*k)
  dir_samples <- rdirichlet(size = n_samples, alpha = alpha_dir)

  # STEP 2 -random samples from a StickBreaking
  for(i in 1:n_samples){
    stick <- rStickBreaking(alpha = alpha + sigma*k, sigma = sigma,
                            truncation_method = truncation_method,
                            prob_cutoff = prob_cutoff, static_cutoff = static_cutoff)
    weights_squared <- c(dir_samples[i, 1:k], dir_samples[i, k+1]*stick)^2
    Gini_out[i] <- 1 - sum(weights_squared)
  }
  return(Gini_out)
}

#' Random sampler for the posterior distribution of the Gini heterogeneity index under a Dirichlet process
#'
#' @param object An object of class \code{ssm, DP}
#' @param n_samples Number of posterior samples to draw
#' @param truncation_method Truncation method for the stick breaking sampler
#' @param prob_cutoff Probabilistic truncation cutoff for the stick breaking sampler. Valid only for \code{truncation_method = "dynamic"}
#' @param static_cutoff Integer cutoff for the stick breaking sampler. Valid only for \code{truncation_method = "static"}
#' @param ... Additional model parameters
#'
#' @export
rGini.DP <- function(object, n_samples = 100, truncation_method = "dynamic", prob_cutoff = 0.99, static_cutoff = 1000, ...){
  # Sample the weights from a posterior distribution of the Pitman-Yor given the abundances
  n <- sum(object$abundances)
  k <- length(object$abundances)
  alpha <- object$param[1]
  sigma <- 0

  Gini_out <- rep(0, n_samples)
  ## STEP 1 - random samples from a Dirichlet distribution
  alpha_dir <- c(object$abundances - sigma, alpha + sigma*k)
  dir_samples <- rdirichlet(size = n_samples, alpha = alpha_dir)

  # STEP 2 -random samples from a StickBreaking
  for(i in 1:n_samples){
    stick <- rStickBreaking(alpha = alpha + sigma*k, sigma = sigma,
                            truncation_method = truncation_method,
                            prob_cutoff = prob_cutoff, static_cutoff = static_cutoff)
    weights_squared <- c(dir_samples[i, 1:k], dir_samples[i, k+1]*stick)^2
    Gini_out[i] <- 1 - sum(weights_squared)
  }
  return(Gini_out)
}





