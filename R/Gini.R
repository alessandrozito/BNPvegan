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
  freq <- object$frequencies
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
  freq <- object$frequencies
  n <- sum(freq)
  K <- length(freq)

  out <- 1 - 1 / Poch2(alpha + n) * ((1 - sigma) * (alpha + K * sigma) + sum(Poch2(freq - sigma)))
  out
}


