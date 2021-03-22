get_m_saturation <- function(object, n, k, tolerance){
  ## Determine the number of additional samples to reach convergence of E(K_n)
  if(object$model == "LL3"){
    alpha <- object$par[1]
    sigma <- object$par[2]
    phi <- object$par[3]
    value_old <- k + prob_LL3(n, alpha, sigma, phi)
    value_new <- value_old + prob_LL3(n+1, alpha, sigma, phi)
    diff <- value_new - value_old
    m <- 0
    while(tolerance <= diff){
      m <- m + 1
      value_old <- value_new
      value_new <- value_old + prob_LL3(n+m, alpha, sigma, phi)
      diff <- value_new - value_old
    }
  } else if (object$model == "Weibull"){
    phi <- object$par[1]
    lambda <- object$par[2]
    value_old <- k + prob_Weibull(n, phi, lambda)
    value_new <- value_old + prob_Weibull(n+1, phi, lambda)
    diff <- value_new - value_old
    m <- 0
    while(tolerance <= diff){
      m <- m + 1
      value_old <- value_new
      value_new <- value_old + prob_Weibull(n+m, phi, lambda)
      diff <- value_new - value_old
    }
  }
  return(m)
}



#' Saturation for a sequential discovery model
#'
#' @param object an object of class \code{\link[sdm]{sdm}}.
#' @param method method to obtain the saturation. Available options are \code{"approximate"} and \code{"montecarlo"}
#' @param n_samples number of Monte Carlo samples from the posterior distribution of the saturation. Valid only for \code{method = "montecarlo"}
#' @param tolerance tolerance parameter to determine the truncation point for simulating K_inf. Valid only for \code{method = "montecarlo"}
#' @param ... Additional parameters
#'
#' @export
saturation <- function(object, method = "approximate", n_samples = 100, tolerance = 1e-7, ...){
  if(method == "approximate"){
    # Approximate saturation is just k/EK_inf
    sat <- unname(object$saturation)
    return(sat)
  } else {
    # Use a monte carlo simulation to obtain samples from the posterior random variable.
    # This is done by truncation of the sequence of discoveries
    k <- length(object$frequencies)
    n <- sum(object$frequencies)
    # Determining the truncation point
    m <- get_m_saturation(object=object, n=n, k=k, tolerance = tolerance)

    # Determing the probabilities up to the truncation point
    if(object$model == "LL3"){
      alpha <- object$par[1]
      sigma <- object$par[2]
      phi <- object$par[3]
      probs <- prob_LL3(n = c(n:(n+m)), alpha, sigma, phi)
    } else if(object$model == "Weibull"){
      phi <- object$par[1]
      lambda <- object$par[2]
      probs <- prob_Weibull(n = c(n:(n+m)), phi, lambda)
    }

    # Random sampling up to the truncation point
    sat <- rep(NA, n_samples)
    for(i in 1:n_samples){
      sat[i] <- k/(k + sum(rbinom(size = 1, n = length(probs), prob = probs)))
    }
  }
  return(sat)
}






