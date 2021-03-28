#' Saturation for a sequential discovery model
#'
#' @param object an object of class \code{sdm}.
#' @param method method to obtain the saturation. Available options are \code{"approximate"}, \code{"montecarlo"} and \code{"target"}
#' @param n_samples number of Monte Carlo samples from the posterior distribution of the saturation. Valid only for \code{method = "montecarlo"}
#' @param target Desired level of saturation to reach. Valid only for \code{"target"}
#' @param tolerance tolerance parameter to determine the truncation point for simulating K_inf. Valid only for \code{method = "montecarlo"}
#' @param ... Additional parameters.
#'
#' @export
saturation <- function(object, method = "approximate", n_samples = 100, target = NULL,  tolerance = 1e-7, ...){
  if(method == "approximate"){
    # Approximate saturation is just k/EK_inf
    sat <- unname(object$saturation)
    return(sat)
  } else if(method == "montecarlo"){
    # Use a Monte Carlo simulation to obtain samples from the posterior random variable.
    # This is done by truncation of the sequence of discoveries
    Kinf <- sample_Kinf(object, n_samples = n_samples, tolerance = tolerance)
    sat <- sum(object$discoveries)/Kinf
    return(sat)
  } else if (method == "target"){
    # Return the average number of additional samples to reach the target level of saturation
    if(is.null(target)){
      cat("Please specify a target value\n")
    } else if(target <0 | target >0.999){
      cat("Target value has to be between 0 and 0.999\n")
    } else {

      if(target < object$saturation){
        cat("The approximate saturation is already", round(unname(object$saturation), 4), "\n Please specify a higher target value \n")
      } else {
        # Return the number of additional samples
        n <- sum(object$abundances)
        k <- length(object$abundances)
        if(object$model == "LL3"){
          alpha <- object$par[1]
          sigma <- object$par[2]
          phi <- object$par[3]
          # Get the exact level for Kinf. Not rounded.
          Kinf <- cubature::cubintegrate(f = prob_LL3, alpha = alpha, sigma = sigma, phi = phi, lower = n, upper = Inf)$integral + k
          m <- get_n_target_saturation_LL3(n, k, alpha, sigma, phi, Kinf, target)
        } else if (object$model == "Weibull"){
          phi <- object$par[1]
          lambda <- object$par[2]
          Kinf <- cubature::cubintegrate(f = prob_Weibull, phi = phi, lambda = lambda, lower = n, upper = Inf)$integral + k
          m <- get_n_target_saturation_Weibull(n, k, phi, lambda, Kinf, target)
        }
        #cat("Expected number of additional samples to reach a target saturation equal to ", target, ": ", m, "\n")
        return(m)
      }
    }
 }
}






