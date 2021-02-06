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
sdm <- function(frequencies, model = "LL3", n_resamples = 1000, verbose = TRUE) {
  # Step 0 - filter out the frequencies equal to 0
  frequencies <- frequencies[frequencies > 0]

  # Average output
  d <- c(1, diff(rarefy_C(freq, sum(frequencies), length(frequencies))))


  # Initialize an empty matrix for the parameters
  if (model == "LL3") {
    # Initialize the matrix of predictors
    n <- sum(frequencies) - 1
    X <- cbind(1, log(1:n), c(1:n))
  }

  if (model == "LL3") {
    # Fit the three-parameter log-logistic
    D <- cbind(1 - d, d)
    fitglm <- glmnet::glmnet(
      x = X[, -1], y = D[-1, ], family = "binomial", upper.limits = 0, lambda = 0,
      standardize = FALSE, tresh = 1e-7
    )
    coeffs <- c(fitglm$a0, as.numeric(fitglm$beta))
    if (coeffs[3] == 0) {
      coeffs[3] <- -1e-7 # Force convergence
    }

    # Compute loglikelihood and save parameters
    loglik <- -logLik_LL3(d = d[-1], X = X, beta_0 = coeffs[1], beta_1 = coeffs[2], beta_2 = coeffs[3])
    par <- c(exp(coeffs[1]), 1 + coeffs[2], exp(coeffs[3]))

  } else if (model == "Weibull") {
    # Fit the Weibull model
    fit <- max_logLik_Weibull(d = d)
    loglik <- -fit$objective
    par <- c(fit$par[1], fit$par[2], loglik)
  } else {
    cat(paste("No method exists for specified model:", model))
    return(NA)
  }

  # Compute E(Kinf) and var(Kinf)
  asym_p <- moments_Kinf(par[-4], n = length(d), k = sum(d), model = model)

  # Return the output
  out <- list(
    model = model,
    frequencies = frequencies,
    n_resamples = n_resamples,
    discoveries = d,
    par = par,
    loglik = loglik,
    Asymp_moments = asym_p,
    saturation = sum(d) / asym_p[1]
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
      paste0("\t Expected new species to discover: ", asy_rich[1] - richness),
      paste0("\t Sample saturation: ", round(saturation, 4)),
      "\nParameters:",
      paste0("\t ", knitr::kable(t(c(object$par, object$loglik)), "simple")),
      sep = "\n"
    )
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
plot.sdm <- function(object, n_points = 100, type = "rarefaction", m = NULL, ...) {
  # Step 1 - Plot the accumulation curve with the highest likelihood
  accum <- cumsum(object$discoveries)
  # Rarefaction curve
  rar <- rarefaction(object)

  if (type == "rarefaction") {

    # Step 3 - compute the average rarefaction curve
    df <- data.frame("n" = c(1:length(accum)), "accum" = accum, "rar" = rar)

    if (nrow(df) > n_points) {
      seqX <- 1:nrow(df)
      seqY <- split(seqX, sort(seqX %% n_points))
      df <- df[unlist(lapply(seqY, function(a) tail(a, 1))), ]
    }

    p <- ggplot2::ggplot(df) +
      ggplot2::geom_point(ggplot2::aes(x = n, y = accum), shape = 1) +
      ggplot2::geom_line(ggplot2::aes(x = n, y = rar), color = "red", size = 0.9) +
      ggplot2::theme_bw() +
      ggplot2::facet_wrap(~"Rarefaction curve") +
      ggplot2::ylab(expression(K[n]))

    return(p)
  } else if (type == "extrapolation") {
    # Extrapolate up to m. if unspecified, m = n
    if (is.null(m)) {
      m <- length(accum)
    }

    ext <- extrapolation(object, m = m)

    df <- data.frame("n" = c(1:(length(rar) + length(ext))), "curve" = c(rar, ext))
    p <- ggplot2::ggplot(df) +
      ggplot2::geom_line(ggplot2::aes(x = n, y = curve), color = "red", size = 0.9) +
      ggplot2::theme_bw() +
      ggplot2::facet_wrap(~"Rarefaction and Extrapolation curve") +
      ggplot2::geom_segment(x = m, xend = m, y = 0, yend = Inf, linetype = "dashed") +
      ggplot2::ylab(expression(K[n]))
    return(p)
  }
}

#' @export
extrapolation <- function(x, ...) {
  UseMethod("extrapolation", x)
}

#' @export
#'
extrapolation.sdm <- function(object, m, ...) {
  n <- length(object$discoveries)
  k <- sum(object$discoveries)
  if (object$model == "LL3") {
    extr <- k + cumsum(prob_LL3(c(n:(n + m - 1)), alpha = object$par[1], sigma = object$par[2], phi = object$par[3]))
  } else if (object$model == "Weibull") {
    extr <- k + cumsum(prob_Weibull(c(n:(n + m - 1)), phi = object$par[1], lambda = object$par[2]))
  }

  return(extr)
}

#' @export
#'
coef.sdm <- function(object, ...) {
  return(object$par)
}

#' @export
#'
rarefaction.sdm <- function(object, ...) {
  predict(object)
}

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
