################################################
# Functions to compute the rarefaction curves. Note that rarefaction curves
# can be frequency based (ie. the average accumulation curve) or model-based
# (ie. estimated form a model - PY, DP or sdm)
################################################


#' Sample-based or model-based rarefaction curve
#'
#' @param object An object of class \code{numeric}, \code{ssm} or \code{sdm}
#' @param ... Additional parameters
#'
#' @export
rarefaction <- function(object, ...) {
  UseMethod("rarefaction", object)
}

#' Sample-based rarefaction curve for a vector of species abundances
#' @param object  An vector of abundances
#' @param verbose If TRUE, the estimation process will be printed
#' @param ... Additional parameters
#' @export
rarefaction.numeric <- function(object, verbose = TRUE, ...){
  # Filter out the 0 counts and make sure that we have integers
  object <- as.integer(object)
  freq <- object[object>0]

  # Call the C function
  c(rarefy_C(freq = freq, n =  sum(freq), K = length(freq), verbose = verbose))
}

#' Model-based rarefaction curve for a Sequential Discovery Model
#' @param object  An object of class \code{sdm}.
#' @param ... Additional parameters
#'
#' @export
rarefaction.sdm <- function(object, ...){
  predict(object)
}

#' Model-based rarefaction curve for a Dirichlet process model
#' @param object An object of class \code{ssm, DP}
#' @param ... Additional parameters
#' @export
rarefaction.DP <- function(object, ...) {
  n <- sum(object$abundances)
  expected_cl_py(1:n, sigma = 0, alpha = object$param)
}

#' Model-based rarefaction curve for a Dirichlet process model
#' @param object An object of class \code{ssm, PY}
#' @param ... Additional parameters
#' @export
rarefaction.PY <- function(object, ...) {
  n <- sum(object$abundances)
  expected_cl_py(1:n, sigma = object$param[2], alpha = object$param[1])
}





