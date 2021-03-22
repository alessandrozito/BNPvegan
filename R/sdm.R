#' Sequential discovery model as in Zito et al. (2020+)
#'
#' @param frequencies Vector of frequencies of the observed species
#' @param model Model to fit. Options are "LL3" and "Weibull"
#' @param verbose if TRUE, monitor the construction of the average rarefaction curve.
#'
#' @details This function fits the sequential discovery model is Zito et al. (2020+) on the average rarefaction curve obtained from a vector of species frequencies.
#'          The models available are "LL3" and "Weibull".
#'
#' @export
#' @examples # Fit the model
#' frequencies <- fungalOTU
#' fit <- sdm(frequencies, model = "LL3")
#' summary(fit)
#'
sdm <- function(frequencies, model = "LL3", verbose = TRUE) {
  # Step 0 - filter out the frequencies equal to 0
  frequencies <- frequencies[frequencies > 0]

  # Extract the rarefaction curves
  d <- c(1, diff(rarefy_C(frequencies, sum(frequencies), length(frequencies), verbose)))

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
    names(par) <- c("alpha", "sigma", "phi")

  } else if (model == "Weibull") {
    # Fit the Weibull model
    fit <- max_logLik_Weibull(d = d[-1])
    loglik <- -fit$objective
    par <- c("phi" = fit$par[1], "lambda" = fit$par[2])
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
    discoveries = d,
    par = par,
    loglik = loglik,
    Asymp_moments = asym_p,
    saturation = c("saturation" = sum(d) / asym_p[1])
  )
  class(out) <- "sdm"
  return(out)
}


#' Summary for a species discovery model
#'
#' @param object An object of class \code{\link[sdm]{sdm}}.
#' @param ... additional parameters
#'
#' @details A function to print out a summary of a species discovery model on a given vector of frequencies.
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
      "\nQuantities:",
      paste0("\t Abundance: ", abundance),
      paste0("\t Richness: ", richness),
      paste0("\t Expected species at infinity: ", asy_rich[1]),
      paste0("\t Standard deviation at infinity: ", asy_rich[2]),
      paste0("\t Expected new species to discover: ", asy_rich[1] - richness),
      paste0("\t Sample saturation: ", round(saturation, 4)),
      "\nParameters:",
      paste0("\t ", knitr::kable(t(c(object$par, "logLik" = object$loglik)), "simple")),
      sep = "\n"
    )
}

#' Predict function for a species discovery model
#'
#' @param object object of class \code{\link[sdm]{sdm}}
#' @param newdata Indexes (or additional new samples) to predict
#' @param ... additional values
#'
#' @export
#'
predict.sdm <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) {
    n <- c(1:length(object$discoveries)) - 1
    index_to_store <- newdata <- n + 1
  } else {
    n <- c(1:max(newdata)) - 1
    index_to_store <- newdata
  }

  # Prediction under LL3
  if (object$model == "LL3") {
    pred <- cumsum(prob_LL3(n = n, alpha = object$par[1], sigma = object$par[2], phi = object$par[3]))
  # Prediction under Weibull
  } else if (object$model == "Weibull") {
    pred <- cumsum(prob_Weibull(n = n, phi = object$par[1], lambda = object$par[2]))
  }

  return(unname(pred)[index_to_store])
}


#' Plot the average rarefaction curve and the fitted one (in red)
#'
#' @param object An object of class \code{\link[sdm]{sdm}}.
#' @param n_points Number of points to plot in the accumulation curve
#' @param type Type of curve to plot. Available options are "rarefaction" and "extrapolation".
#'             In the second case, one needs to provide the additional number of samples m up to which extrapolate.
#' @param m Additional number of samples. Required for extrapolation only
#' @param ... additional parameters
#'
#' @export
plot.sdm <- function(object, n_points = 100, type = "rarefaction", m = NULL, ...) {
  # Step 1 - Sample-based rarefaction curve
  accum <- cumsum(object$discoveries)
  # Rarefaction curve
  rar <- rarefaction(object)

  if (type == "rarefaction") {

    # Step 3 - create the dataframe
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

    ext <- extrapolation(object, m = 1:m)
    cutoff <- length(rar)

    df <- data.frame("n" = c(1:(length(rar) + length(ext))), "curve" = c(rar, ext), "accum" = c(accum, rep(NA, length(ext))))
    if (nrow(df) > n_points) {
      seqX <- 1:nrow(df)
      seqY <- split(seqX, sort(seqX %% n_points))
      df <- df[unlist(lapply(seqY, function(a) tail(a, 1))), ]
    }

    p <- ggplot2::ggplot(df) +
      ggplot2::geom_line(ggplot2::aes(x = n, y = curve), color = "red", size = 0.9) +
      ggplot2::geom_point(ggplot2::aes(x = n, y = accum), shape = 1, na.rm=TRUE) +
      ggplot2::theme_bw() +
      ggplot2::facet_wrap(~"Rarefaction and extrapolation curve") +
      ggplot2::geom_segment(x = cutoff, xend = cutoff, y = 0, yend = Inf, linetype = "dashed") +
      ggplot2::ylab(expression(K[n]))
    return(p)
  }
}


#' Asymptotic species richness estimates for the species discovery model.
#'
#' @param object an object of class \code{\link[sdm]{sdm}}.
#' @param ... additional parameters
#' @export
asym_richness <- function(object, ...) {
  rich <- unname(c(object$Asymp_moments))
  rich[1] <- as.integer(rich[1])
  names(rich) <- c("Exp. species at Inf","sd at Inf")
  return(rich)
}



#' @export
coef.sdm <- function(object, ...) {
  return(object$par)
}

#' @export
logLik.sdm <- function(object, ...) {
  return(object$loglik)
}


