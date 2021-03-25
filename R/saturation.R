#' Saturation for a sequential discovery model
#'
#' @param object an object of class \code{sdm}.
#' @param method method to obtain the saturation. Available options are \code{"approximate"} and \code{"montecarlo"}
#' @param n_samples number of Monte Carlo samples from the posterior distribution of the saturation. Valid only for \code{method = "montecarlo"}
#' @param tolerance tolerance parameter to determine the truncation point for simulating K_inf. Valid only for \code{method = "montecarlo"}
#' @param ... Additional parameters.
#'
#' @export
saturation <- function(object, method = "approximate", n_samples = 100, tolerance = 1e-7, ...){
  if(method == "approximate"){
    # Approximate saturation is just k/EK_inf
    sat <- unname(object$saturation)
  } else {
    # Use a Monte Carlo simulation to obtain samples from the posterior random variable.
    # This is done by truncation of the sequence of discoveries
    Kinf <- sample_Kinf(object, n_samples = n_samples, tolerance = tolerance)
    sat <- sum(object$discoveries)/Kinf
  }
  return(sat)
}






