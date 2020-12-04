#' Find the Empirical Bayes estimates of Bayesian Nonparametric Species Sampling model
#'
#'
#' @param frequencies A \code{K}-dimensional vector of frequencies
#' @param model Model to fit. Available models are "DP" (Dirichlet Process) and "PY" (Pitman-Yor process)
#' #'
#' @return a list (?)
#' @export
ssm <- function(frequencies, model){

  if (model == "DP") {
    fit <- max_EPPF_DP(frequencies)
    out <- list(freq = freq, param = fit$par, loglik = - fit$objective)
    class(out) <- c("ssm","DP")
    return(out)
  }

  if (model == "PY") {
    fit <- max_EPPF_PY(frequencies)
    out <- list(freq = freq, param = fit$par, loglik = - fit$objective)
    class(out) <- c("ssm","PY")
    return(out)
  }
}



