#' Fit the sequential discovery model of Zito et al. (2020+)
#'
#' @param frequencies Counts of the species
#' @param n_resamples number of accumulation curves to generate
#' @param model model to fit. Options are "LL3" and "Weibull"
#' @param verbose Whether to monitor the state of the simulation or not
#'
#' @details This function runs the sequential discovery model is Zito et al. (2020+)
#' @export
#'
sdm <- function(frequencies, n_resamples = 1000L, model = "LL3", verbose = TRUE) {
  # Step 0 - filter out the frequencies equal to 0
  frequencies <- frequencies[frequencies > 0]

  # Initialize an empty matrix for the parameters
  if (model == "LL3") {
    param <- matrix(NA, nrow = n_resamples, ncol = 4)
    colnames(param) <- c("alpha", "sigma", "phi", "loglik")
    # Initialize the matrix of predictors
    n <- sum(frequencies) - 1
    X <- cbind(1, log(1:n), c(1:n))
  } else if (model == "Weibull") {
    param <- matrix(NA, nrow = n_resamples, ncol = 3)
    colnames(param) <- c("phi", "lambda", "loglik")
  }

  # List to store the position of the discoveries in the list
  discoveries_indexes <- matrix(NA, nrow = n_resamples, ncol = length(frequencies))
  if (n_resamples <= 10) {
    verbose <- FALSE
  }
  verbose_step <- round(n_resamples / 10)

  for (i in 1:n_resamples) {
    # Obtain the discovery sequence
    seqD <- sample_sequence(frequencies)
    d <- extract_discoveries(seqD)
    discoveries_indexes[i, ] <- which(d == 1)

    use_constrained_max <- FALSE

    if (model == "LL3") {
      # Fit the three-parameter log-logistic

      if (use_constrained_max == TRUE) {
        fit <- try(logit_regression(y = d[-1], X = X), silent = TRUE)
        if (inherits(fit, "try-error")) {
          use_constrained_max <- FALSE
          fit <- max_logLik_LL3(d = d[-1], X = X)
          loglik <- -fit$objective
        } else {
          loglik <- max(fit$Convergence[, 2])
          if (fit$par[2] > 0 | fit$par[3] > 0) {
            fit <- max_logLik_LL3(d = d[-1], X = X)
            loglik <- -fit$objective
          }
        }
      } else {
        fit <- max_logLik_LL3(d = d[-1], X = X)
        loglik <- -fit$objective
      }

      param[i, ] <- c(exp(fit$par[1]), 1 + fit$par[2], exp(fit$par[3]), loglik)
    } else if (model == "Weibull") {
      # Fit the Weibull model
      fit <- max_logLik_Weibull(d = d)
      loglik <- -fit$objective
      param[i, ] <- c(fit$par[1], fit$par[2], loglik)
    } else {
      cat(paste("No method exists for specified model:", model))
      return(NA)
    }

    # Monitor output
    if (verbose) {
      if (i %% verbose_step == 0) {
        cat(paste0("Number of accumulation curves fitted: ", i, " [", round(100 * i / n_resamples), "%]"), sep = "\n")
      }
    }
  }

  # Choose the accumulation curve with the highest loglikelihood
  if (model == "LL3") {
    selected_curve <- which.max(param[, 4])
    par <- param[selected_curve, -4]
    loglik <- param[selected_curve, 4]
  } else if (model == "Weibull") {
    selected_curve <- which.max(param[, 3])
    par <- param[selected_curve, -3]
    loglik <- param[selected_curve, 3]
  }
  d <- rep(0, sum(frequencies))
  d[discoveries_indexes[selected_curve, ]] <- 1

  # Compute E(Kinf) and var(Kinf)
  Asymp_moments <- moments_Kinf(par, n = length(d), k = sum(d), model = model)

  # List to store the re-sampling output
  resampling_output <- list(param = param, discoveries_indexes = discoveries_indexes, selected_curve = selected_curve)

  # Return the output
  out <- list(
    model = model,
    frequencies = frequencies,
    n_resamples = n_resamples,
    resampling_output = resampling_output,
    discoveries = d,
    par = par,
    loglik = loglik,
    Asymp_moments = Asymp_moments,
    saturation = unname(sum(d)/Asymp_moments[1])
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
summary.sdm <- function(object, ...) {

  # Sample abundance and richness
  abundance <- sum(object$frequencies)
  richness <- length(object$frequencies)
  saturation <- object$saturation
  asy_rich <- unname(round(object$Asymp_moments, 2))
  if (object$model == "LL3") {
    mod_print <- "\t Three-parameter log-logistic (LL3)"
  } else if (object$model == "Weibull") {
    mod_print <- "\t Weibull"
  }
  mod_str <-
    # Summary
    # Print the summary
    cat("Model:",
      mod_print,
      paste0("\t Number of resamples: ", object$n_resamples),
      "\nQuantities:",
      paste0("\t Abundance: ", abundance),
      paste0("\t Richness: ", richness),
      paste0("\t Expected species at infinity: ", asy_rich[1]),
      paste0("\t Standard deviation at infinity: ", asy_rich[2]),
      paste0("\t Sample saturation: ", round(saturation, 4)),
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
#' @param newdata Indexes to predict
#' @param ... additional values
#' @export
predict.sdm <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) {
    n <- c(1:length(object$discoveries)) - 1
    index_to_store <- newdata <- n + 1
  } else {
    n <- c(1:max(newdata)) - 1
    index_to_store <- newdata
  }

  if (object$model == "LL3") {
    pred <- cumsum(prob_LL3(n = n, alpha = object$par[1], sigma = object$par[2], phi = object$par[3]))
  } else if (object$model == "Weibull") {
    pred <- cumsum(prob_Weibull(n = n, phi = object$par[1], lambda = object$par[2]))
  }

  return(unname(pred)[index_to_store])
}


#' Plot rarefaction curves
#'
#' @param object An object of class \code{\link[sdm]{sdm}}.
#' @param n_points Number of points to plot in the accumulation curve
#' @param ... additional parameters

#' @export
plot.sdm <- function(object, n_points = 100, ...) {
  # Step 1 - Plot the accumulation curve with the highest likelihood
  accum <- cumsum(object$discoveries)

  # Plot the rarefaction by predicting
  rar <- predict(object)

  # Step 3 - compute the average rarefaction curve
  df <- data.frame("n" = c(1:length(accum)), "accum" = accum, "rar" = rar)

  if (nrow(df) > n_points) {
    seqX <- 1:nrow(df)
    seqY <- split(seqX, sort(seqX %% n_points))
    df <- df[unlist(lapply(seqY, function(a) tail(a, 1))), ]
  }

  ggplot2::ggplot(df) +
    ggplot2::geom_point(ggplot2::aes(x = n, y = accum), shape = 1) +
    ggplot2::geom_line(ggplot2::aes(x = n, y = rar), color = "red", size = 0.9) +
    ggplot2::theme_bw() +
    ggplot2::facet_wrap(~"Rarefaction curve") +
    ggplot2::ylab(expression(K[n]))
}

#' Extract the coefficients of the sdm output
#'
#' @param object an object of class \code{\link[sdm]{sdm}}.
#' @param ... additional parameters
#'
#' @export
#'
coef.sdm <- function(object, ...) {
  return(object$par)
}

#' Extract the loglikelihood of the sdm output
#'
#' @param object an object of class \code{\link[sdm]{sdm}}.
#' @param ... additional parameters
#'
#' @export
#'
logLik.sdm <- function(object, ...) {
  return(object$loglik)
}

#' Extract the species richness estimates of the sdm output
#'
#' @param object an object of class \code{\link[sdm]{sdm}}.
#' @param ... additional parameters
#'
#' @export
#'
asymptotic_richness <- function(object, ...) {
  return(c(object$Asymp_moments, "saturation" = object$saturation))
}
