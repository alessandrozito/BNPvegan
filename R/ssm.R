#' Find the Empirical Bayes estimates of Bayesian Nonparametric Species Sampling model
#'
#'
#' @param freq A \code{K}-dimensional vector of frequencies
#' @param model Model to fit. Available models are "DP" (Dirichlet Process) and "PY" (Pitman-Yor process)
#' #'
#' @return a list (?)
#' @export
ssm <- function(freq, model){

  if (model == "Dirichlet") {
    fit <-
    out <- list(freq = freq, param = fit$par, loglik = - fit$objective)
    class(out) <- c("ssm","DP")
    return(out)
  }

  if (model == "Pitman-Yor") {
    out <- nlminb(
      start = c(1, 0.5),
      objective = function(param) - EPPF_PitmanYor(freq = freq, alpha = param[1], sigma = param[2]),
      lower = c(-Inf, 1e-16),
      upper = c(Inf, 1 - 1e-10)
    )
    pars <- out$par
    loglik <- -out$objective
    return(list(
      "model" = "Pitman-Yor",
      "species_counts" = species_counts,
      "exchangeability" = TRUE,
      "parameters" = c("alpha" = pars[1], "sigma" = pars[2]),
      "loglik" = loglik,
      "Out_of_bound" = NULL,
      "n_resamples" = 1
    ))
  }
}



