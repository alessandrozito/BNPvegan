################################################
# Functions to compute the extrapolation curves.
# Available only for model based types
################################################

#' Model base extrapolation curve and prediction
#'
#' @param object An object of class \code{ssm} or \code{sdm}
#' @param ... Additional parameters
#'
#' @export
extrapolation <- function(object, ...) {
  UseMethod("extrapolation", object)
}

#' Extrapolation function for a species discovery model.
#' @param object An object of class \code{sdm}.
#' @param m Additional number of samples to predict.
#' @param ... Additional parameters
#' @export
extrapolation.sdm <- function(object, m, ...) {
  m_max <- max(m)

  n <- length(object$discoveries)
  k <- sum(object$discoveries)
  if (object$model == "LL3") {
    extr <- k + cumsum(prob_LL3(c(n:(n + m_max - 1)), alpha = object$par[1], sigma = object$par[2], phi = object$par[3]))
  } else if (object$model == "Weibull") {
    extr <- k + cumsum(prob_Weibull(c(n:(n + m_max - 1)), phi = object$par[1], lambda = object$par[2]))
  }

  return(extr[m])
}


#' @export
extrapolation.DP <- function(object, m, ...) {
  alpha <- object$param[1]
  freq <- object$abundances
  n <- sum(freq)
  K <- length(freq)
  extrapolate_cl_py(m = m, n = n, K = K, sigma = 0, alpha = alpha)
}

#' @export
extrapolation.PY <- function(object, m, ...) {
  alpha <- object$param[1]
  sigma <- object$param[2]
  freq <- object$abundances
  n <- sum(freq)
  K <- length(freq)
  extrapolate_cl_py(m = m, n = n, K = K, sigma = sigma, alpha = alpha)
}
