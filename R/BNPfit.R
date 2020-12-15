#' Find the Empirical Bayes estimates of Bayesian Nonparametric Species Sampling model
#'
#'
#' @param species_counts Vector of integer counts of the observed clusters (species)
#' @param model Model to fit. Available models are "Dirichlet", "Pitman-Yor", "LL2" and "LL3".
#' @param n_resamples Number of resamples of the original sequence. Needed for Non-exchangeable models only. Default = 500
#'
#' @return a list (?)
#' @export
#'
BNPfit <- function(species_counts, model, n_resamples = 500L) {

  ################################################
  # Part 2 - Non-Exchangeable models
  ################################################
  #---------------------- LL-2 - Two parameter log-logistic of Zito et al. (2020)
  if (model == "LL2") {
    species_sequence <- extract_sequence(species_counts)
    discoveries <- extract_discoveries(species_sequence)
    pars <- matrix(NA, nrow = n_resamples, ncol = 2)

    # Initialize
    colnames(pars) <- c("alpha", "sigma")
    out_of_bound <- rep(NA, n_resamples)

    # Start with one sequence
    ll_out <- fit_LL2(discoveries = discoveries)
    pars[1, ] <- ll_out$parameters
    out_of_bound[1] <- ll_out$Out_of_bound

    # Iterate through the all the other re ordering of the sequence.
    for (r in 2:n_resamples) {
      species_sequence <- sample(species_sequence, size = length(species_sequence), replace = FALSE)
      discoveries <- extract_discoveries(species_sequence)
      ll_out <- fit_LL2(discoveries = discoveries)
      pars[r, ] <- ll_out$parameters
      out_of_bound[r] <- ll_out$Out_of_bound
    }

    # Return
    return(list(
      "model" = "LL2",
      "species_counts" = species_counts,
      "exchangeability" = FALSE,
      "parameters" = pars,
      "loglik" = "Still to add", ### TO DO: compute the loglikelihood for every reshuffle
      "Out_of_bound" = out_of_bound,
      "n_resamples" = n_resamples
    ))

    #---------------------- LL3 - Three parameter log-logistic of Zito et al. (2020)
  } else if (model == "LL3") {
    species_sequence <- extract_sequence(species_counts)
    pars <- matrix(NA, nrow = n_resamples, ncol = 3)

    # Initialize
    colnames(pars) <- c("alpha", "sigma", "phi")
    out_of_bound <- rep(NA, n_resamples)

    # Start with one sequence
    discoveries <- extract_discoveries(species_sequence)
    ll_out <- fit_LL3(discoveries = discoveries)
    pars[1, ] <- ll_out$parameters
    out_of_bound[1] <- ll_out$Out_of_bound

    # Iterate through the all the other re ordering of the sequence.
    for (r in 2:n_resamples) {
      species_sequence <- sample(species_sequence, size = length(species_sequence), replace = FALSE)
      discoveries <- extract_discoveries(species_sequence)
      ll_out <- fit_LL3(discoveries = discoveries)
      pars[r, ] <- ll_out$parameters
      out_of_bound[r] <- ll_out$Out_of_bound
    }

    # Return
    return(list(
      "model" = "LL3",
      "species_counts" = species_counts,
      "exchangeability" = FALSE,
      "parameters" = pars,
      "loglik" = "Still to add", ### TO DO: compute the loglikelihood for every reshuffle
      "Out_of_bound" = out_of_bound,
      "n_resamples" = n_resamples
    ))
  } else {
    print(paste("No available method for the specified model:", model))
    return(NA)
  }
}
