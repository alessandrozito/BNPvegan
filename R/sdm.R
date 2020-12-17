#' Fit the sequential discovery model of Zito et al. (2020+)
#'
#' @param frequencies Counts of the species
#' @param n_resamples number of accumulation curves to generate
#' @param verbose Whether to monitor the function or not
#'
#' @details This function.
#' @export
#'
sdm <- function(frequencies, n_resamples = 1000L, verbose = TRUE) {

  # Initialize an empty matrix for the parameters
  param <- matrix(NA, nrow = n_resamples, ncol = 4)
  colnames(param) <- c("alpha", "sigma", "phi", "loglik")
  # List to store the position of the discoveries in the list
  discoveries_indexes <- matrix(NA, nrow = n_resamples, ncol = length(frequencies))
  if(n_resamples <= 10){
    verbose <- FALSE
  }
  verbose_step <- round(n_resamples / 10)

  # Initialize the matrix of predictors
  n <- sum(frequencies) - 1
  X <- cbind(1, log(1:n), c(1:n))

  for (i in 1:n_resamples) {
    # Obtain the discovery sequence
    seqD <- sample_sequence(frequencies)
    d <- extract_discoveries(seqD)
    discoveries_indexes[i, ] <- which(d == 1)

    # Fit the three-parameter log-logistic
    fit <- logit_regression(y = d[-1], X = X)
    loglik <- max(fit$Convergence[, 2])
    if (fit$par[2] > 0 | fit$par[3] > 0) {
      fit <- max_logLik_LL3(d = d[-1], X = X)
      loglik <- -fit$objective
    }
    param[i, ] <- c(exp(fit$par[1]), 1 + fit$par[2], exp(fit$par[3]), loglik)


    # Monitor output
    if (verbose) {
      if (i %% verbose_step == 0) {
        cat(paste0("Number of accumulation curves fitted: ", i, " [", round(100 * i / n_resamples), "%]"), sep = "\n")
      }
    }
  }

  # Choose the accumulation curve with the highest loglikelihood
  selected_curve <- which.max(param[, 4])
  d <- rep(0, sum(frequencies))
  d[discoveries_indexes[selected_curve, ]] <- 1
  par <- param[selected_curve, -4]

  # Compute E(Kinf) and var(Kinf)
  Asymp_moments <- moments_Kinf(alpha = par[1], sigma = par[2], phi = par[3])

  # List to store the re-sampling output
  resampling_output <- list(param = param, discoveries_indexes = discoveries_indexes, selected_curve = selected_curve)

  # Return the output
  out <- list(
    frequencies = frequencies,
    n_resamples = n_resamples,
    resampling_output = resampling_output,
    discoveries = d,
    par = par,
    loglik = param[selected_curve, 4],
    Asymp_moments = Asymp_moments
  )
  class(out) <- "sdm"
  return(out)
}


#' Summary for a species discovery model
#'
#' @param object An object of class \code{\link[sdm]{sdm}}.
#' @param ... additional parameters
#'
#' @details A function to summarize the three parameter log-logistic output
#' @export
summary.sdm <- function(object, plot = TRUE, ...) {

  # Sample abundance and richness
  abundance <- sum(object$frequencies)
  richness <- length(object$frequencies)
  asy_rich <- unname(round(object$Asymp_moments, 2))
  # Summary
  # Print the summary
  cat("Model:",
    "\t Three-parameter log-logistic (LL3)",
    paste0("\t Number of resamples: ", object$n_resamples),
    "\nQuantities:",
    paste0("\t Abundance: ", abundance),
    paste0("\t Richness: ", richness),
    paste0("\t Expected species at infinity: ", asy_rich[1]),
    paste0("\t Standard deviation at infinity: ", asy_rich[2]),
    "\nParameters:",
    paste0("\t ", knitr::kable(t(c(object$par, object$loglik)), "simple")),
    sep = "\n"
  )

  # Output
  # out <- list(
  #  Abundance = abundance,
  #  Richness = richness,
  #  n_resamples = object$n_resamples,
  #  n_div = n_div,
  #  Asymp_richness = EK,
  #  Asymp_richness_summary = asymp_tab,
  #  Asymp_richness_plot = richness_plot,
  #  param_plot = param_plot,
  #  param_summary = pars_tab
  # )
  # class(out) <- union("summary", class(object))
  # return(invisible(out))
}

#' Predict sdm
#'
#' @param object object of class \code{\link[sdm]{sdm}}
#' @param ... additional values
#' @export
predict.sdm <- function(object, ...) {
  n <- length(object$discoveries) - 1
  N <- c(0:n)
  pred <- cumsum(prob_LL3(n = N, alpha = object$par[1], sigma = object$par[2], phi = object$par[3]))
  return(pred)
}

#' Print method for the summary
#'
#' @param x object of class \code{\link[summary.sdm]{summary.sdm}}
#' @param ... other parameters
#' @export
print.summary.sdm <- function(x, ...) {
  cat("Model:",
    "\t Three-parameter log-logistic (LL3)",
    "\nQuantities:",
    paste0("\t Abundance: ", x$Abundance),
    paste0("\t Richness: ", x$Richness),
    paste0("\t Number of divergent accumulation curves: ", x$n_div, " out of ", x$n_resamples, " sampled [", round(100 * x$n_div / x$n_resamples, 2), "%]"),
    "\nEst. Richness for non-divergent accumulation curves:",
    knitr::kable(round(x$Asymp_richness_summary, 2), "simple"),
    "\nParameters:",
    knitr::kable(x$param_summary, "simple"),
    sep = "\n"
  )
  invisible(x)
}


#' Compute rarefaction curves
#'
#' @param object An object of class \code{\link[sdm]{sdm}}.
#' @param n_points Number of points to plot in the accumulation curve
#' @param ... other parameters

#' @export
plot.sdm <- function(object, n_points = 100, ...) {
  # Step 1 - compute the average accumulation curve by averaging across re-samples
  N <- sum(object$frequencies)
  K_mat <- matrix(0, nrow = object$n_resamples, ncol = N)
  for (i in 1:object$n_resamples) {
    K_mat[i, object$discoveries_indexes[i, ]] <- 1
    K_mat[i, ] <- cumsum(K_mat[i, ])
  }
  avg_accumulation <- colMeans(K_mat)

  # Step 2 - compute the average rarefaction curve
  avg_rarefaction <- rowMeans(apply(object$param, 1, function(p) expected_rarefaction(N = N, alpha = p[1], sigma = p[2], phi = p[3])))
  df <- data.frame("n" = c(1:N), "accum" = avg_accumulation, "rar" = avg_rarefaction)

  if (nrow(df) > n_points) {
    seqX <- 1:nrow(df)
    seqY <- split(seqX, sort(seqX %% n_points))
    df <- df[unlist(lapply(seqY, function(a) tail(a, 1))), ]
  }

  ggplot2::ggplot(df) +
    ggplot2::geom_point(ggplot2::aes(x = n, y = accum), shape = 1) +
    ggplot2::geom_line(ggplot2::aes(x = n, y = rar), color = "red", size = 0.9) +
    ggplot2::theme_bw() +
    ggplot2::facet_wrap(~"Average rarefaction") +
    ggplot2::ylab("Number of species")
}
