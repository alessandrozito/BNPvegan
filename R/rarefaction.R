################################################
# Functions to compute the rarefaction curves. Note that rarefaction curves
# can be frequency based (ie. the average accumulation curve) or model-based
# (ie. estimated form a model - PY, DP or sdm)
################################################

#' @export
rarefaction <- function(x, ...) {
  UseMethod("rarefaction", x)
}

#' Sample-based rarefaction curve for a vector of species counts
#' @export
rarefaction.numeric <- function(frequencies, verbose = TRUE, ...){
  # Filter out the 0 counts and make sure that we have integers
  frequencies <- as.integer(frequencies)
  freq <- frequencies[frequencies>0]

  # Call the C function
  c(rarefy_C(freq = freq, n =  sum(freq), K = length(freq), verbose = verbose))
}

#' Model-based rarefaction curve for a Sequential Discovery Model
#' @param object  An object of class \code{\link[sdm]{sdm}}.
#' @param ... Additional parameters
#' @export
rarefaction.sdm <- function(object, ...){
  predict(object)
}

#' Model-based rarefaction curve for a Dirichlet process model
#' @param object An object of class \code{\link[sdm]{ssm}, \link[sdm]{DP}}
#' @export
rarefaction.DP <- function(object, ...) {
  n <- sum(object$frequencies)
  expected_cl_py(1:n, sigma = 0, alpha = object$param)
}

#' Model-based rarefaction curve for a Dirichlet process model
#' @param object An object of class \code{\link[sdm]{ssm}, \link[sdm]{PY}}
#' @export
rarefaction.PY <- function(object, ...) {
  n <- sum(object$frequencies)
  expected_cl_py(1:n, sigma = object$param[2], alpha = object$param[1])
}





