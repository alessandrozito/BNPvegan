################################################
# Functions to the sample coverage.
# Available only for both sample and model based types.
################################################

#' @export
coverage <- function(object, ...) {
  UseMethod("coverage", object)
}


#' Sample-based coverage estimator
#' @export
#' @details Compute the Touring estimator (Good, 1953) for the sample coverage.
coverage.numeric <- function(object, ...) {
  n <- sum(object)
  m1 <- sum(object == 1)
  1 - m1 / n
}

#' Coverage estimator for a Dirichlet process
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


#' @export
rcoverage <- function(object, ...) {
  UseMethod("rcoverage", object)
}

#' @export
rcoverage.DP <- function(object, R = 1000, ...) {
  alpha <- object$param[1]
  freq <- object$frequencies
  n <- sum(freq)

  rbeta(R, n, alpha)
}

#' @export
rcoverage.PY <- function(object, R = 1000, ...) {
  alpha <- object$param[1]
  sigma <- object$param[2]
  freq <- object$frequencies
  n <- sum(freq)
  K <- length(freq)

  rbeta(R, n - sigma * K, alpha + sigma * K)
}
