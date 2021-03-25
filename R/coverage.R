################################################
# Functions to the sample coverage.
# Available only for both sample and model based types.
################################################

#' Sample-based or model-based coverage index
#' @param object An object of class \code{numeric}, \code{ssm} or \code{sdm}
#' @param ... Additional parameters
#' @export
coverage <- function(object, ...) {
  UseMethod("coverage", object)
}


#' Sample-based coverage estimator
#'
#' @param object An object of class \code{numeric}
#' @param ... Additional parameters
#'
#' @export
#' @details Compute the Touring estimator (Good, 1953) for the sample coverage.
coverage.numeric <- function(object, ...) {
  n <- sum(object)
  m1 <- sum(object == 1)
  1 - m1 / n
}

#' Coverage estimator for a Dirichlet process
#'
#' @param object An object of class \code{ssm, DP}
#' @param ... Additional parameters
#'
#' @export
#' @details Compute the posterior coverage estimator for the Dirichlet process.
coverage.DP <- function(object, ...) {
  alpha <- object$param[1]
  freq <- object$frequencies
  n <- sum(freq)
  K <- length(freq)

  n / (alpha + n)
}

#' Coverage estimator for a Pitman-Yor process
#'
#' @param object An object of class \code{ssm, PY}
#' @param ... Additional parameters
#'
#' @export
#' @details Compute the posterior coverage estimator for the Pitman-Yor process.
coverage.PY <- function(object, ...) {
  alpha <- object$param[1]
  sigma <- object$param[2]
  freq <- object$frequencies
  n <- sum(freq)
  K <- length(freq)

  (n - sigma * K) / (alpha + n)
}

#' Coverage estimator for a sequential discovery model
#'
#' @param object An object of class \code{sdm}
#' @param ... Additional parameters
#'
#' @export
#' @details Compute the posterior coverage estimator for a sequential discovery model
coverage.sdm <- function(object, ...) {
  if (object$model == "LL3") {
    cv <- 1 - prob_LL3(n = length(object$discoveries), alpha = object$par[1], sigma = object$par[2], phi = object$par[3])
    # Prediction under Weibull
  } else if (object$model == "Weibull") {
    cv <- 1 - prob_Weibull(n = length(object$discoveries), phi = object$par[1], lambda = object$par[2])
  }
  return(unname(cv))
}


#' Random sampler for the posterior distribution of a model-based coverage
#'
#' @param object An object of class \code{sdm}
#' @param ... Additional parameters
#'
#' @export
rcoverage <- function(object, ...) {
  UseMethod("rcoverage", object)
}

#' Random samples from the coverage posterior distribution for a Dirichlet process
#' @param object An object of class \code{ssm, DP}
#' @param n_samples Number of samples to draw
#' @param ... Additional parameters
#'
#' @export
rcoverage.DP <- function(object, n_samples = 1000, ...) {
  alpha <- object$param[1]
  freq <- object$frequencies
  n <- sum(freq)

  rbeta(n_samples, n, alpha)
}

#' Random samples from the coverage posterior distribution for a Pitman-Yor process
#' @param object An object of class \code{ssm, PY}
#' @param n_samples Number of samples to draw
#' @param ... Additional parameters
#'
#' @export
rcoverage.PY <- function(object, n_samples = 1000, ...) {
  alpha <- object$param[1]
  sigma <- object$param[2]
  freq <- object$frequencies
  n <- sum(freq)
  K <- length(freq)

  rbeta(n_samples, n - sigma * K, alpha + sigma * K)
}
