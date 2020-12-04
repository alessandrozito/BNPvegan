#' Fit a species sampling model
#'
#'
#' @param frequencies A \code{K}-dimensional vector of frequencies
#' @param model Model to fit. Available models are "DP" (Dirichlet Process) and "PY" (Pitman-Yor process)
#' #'
#' @return An object of class "ssm"
#'
#' @export
ssm <- function(frequencies, model) {
  if (model == "DP") {
    fit <- max_EPPF_DP(frequencies)
    out <- list(frequencies = frequencies, param = fit$par, loglik = -fit$objective)
    class(out) <- c("ssm", "DP")
    return(out)
  }

  if (model == "PY") {
    fit <- max_EPPF_PY(frequencies)
    out <- list(frequencies = frequencies, param = fit$par, loglik = -fit$objective)
    class(out) <- c("ssm", "PY")
    return(out)
  }
}


#' Predict method for the DP model
#'
#'
#' Predict method for a DP model, producing either the rarefaction curve or the extrapolation curve.
#'
#' @param object An object of class \code{\link[ssm]{ssm}}.
#' @param newdata A new data frame containing
#' @param ... Further arguments passed to or from other methods.
#'
#' @details The method...
#'
#' @export
#'
predict.DP <- function(object, newdata = NULL, ...) {
  n <- length(object$frequencies)
  expected_cl_py(1:n, sigma = 0, alpha = object$param)
}

#' Predict method for the DP model
#'
#'
#' Predict method for a DP model, producing either the rarefaction curve or the extrapolation curve.
#'
#' @param object An object of class \code{\link[ssm]{ssm}}.
#' @param newdata A new data frame containing
#' @param ... Further arguments passed to or from other methods.
#'
#' @details The method...
#'
#' @export
#'
predict.PY <- function(object, newdata = NULL, ...) {
  n <- length(object$frequencies)
  expected_cl_py(1:n, sigma = object$param[2], alpha = object$param[1])
}
