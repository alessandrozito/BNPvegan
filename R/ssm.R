#' Fit a species sampling model
#'
#'
#' @param frequencies A \code{K}-dimensional vector of frequencies
#' @param model Model to fit. Available models are "DP" (Dirichlet Process) and "PY" (Pitman-Yor process)
#' #'
#' @return An object of class "ssm"
#'
#' @export
ssm <- function(frequencies, model) {
  frequencies <- frequencies[frequencies > 0]
  if (model == "DP") {
    fit <- max_EPPF_DP(frequencies)
    out <- list(frequencies = frequencies, param = fit$par, loglik = -fit$objective)
    class(out) <- c("ssm", "DP")
    return(out)
  }

  if (model == "PY") {
    fit <- max_EPPF_PY(frequencies)
    out <- list(frequencies = frequencies, param = fit$par, loglik = -fit$objective)
    class(out) <- c("ssm", "PY")
    return(out)
  }
}


#' Predict method for the DP model
#'
#'
#' Predict method for a DP model, producing either the rarefaction curve or the extrapolation curve.
#'
#' @param object An object of class \code{\link[ssm]{ssm}}.
#' @param newdata A new data frame containing
#' @param ... Further arguments passed to or from other methods.
#'
#' @details The method...
#'
#' @export
#'
predict.DP <- function(object, newdata = NULL, ...) {
  n <- sum(object$frequencies)
  expected_cl_py(1:n, sigma = 0, alpha = object$param)
}

#' Predict method for the PY model
#'
#'
#' Predict method for a PY model, producing either the rarefaction curve or the extrapolation curve.
#'
#' @param object An object of class \code{\link[ssm]{ssm}}.
#' @param newdata A new data frame containing
#' @param ... Further arguments passed to or from other methods.
#'
#' @details The method...
#'
#' @export
#'
predict.PY <- function(object, newdata = NULL, ...) {
  n <- sum(object$frequencies)
  expected_cl_py(1:n, sigma = object$param[2], alpha = object$param[1])
}

#' @export
#'
summary.DP <- function(object, ...) {
  Poch2 <- function(x) x * (x + 1)

  alpha <- object$param[1]
  freq <- object$frequencies
  n <- sum(freq)
  K <- length(freq)

  out <- cbind(
    Abundance = n,
    Richness = K,
    alpha = alpha,
    Coverage = n  / (alpha + n),
    Additional_species = round(extrapolate_cl_py(m = n, n = n, K = K, sigma = 0, alpha = alpha)) - K,
    Gini = 1 - 1 / Poch2(alpha + n) * (alpha + sum(Poch2(freq)))
  )
    cat("Model: Dirichlet Process",
        paste0("\t Abundance: ", out[1]),
        paste0("\t Richness: ", out[2]),
        paste0("\t alpha: ", round(out[3],4)),
        paste0("\t Sample coverage: ", round(out[4],4)),
        paste0("\t Expected species after additional ",  n, " samples: ", out[5]),
        paste0("\t Posterior Gini diversity: ",  round(out[6],4)),
        sep = "\n"
    )
  invisible(out)
}


#' @export
#'
summary.PY <- function(object, ...) {
  Poch2 <- function(x) x * (x + 1)

  alpha <- object$param[1]
  sigma <- object$param[2]
  freq <- object$frequencies
  n <- sum(freq)
  K <- length(freq)

  out <- cbind(
    Abundance = n,
    Richness = K,
    alpha = alpha,
    sigma = sigma,
    Coverage = (n - sigma * K) / (alpha + n),
    Future_species = round(extrapolate_cl_py(m = n, n = n, K = K, sigma = sigma, alpha = alpha)) - K,
    Gini = 1 - 1 / Poch2(alpha + n) * ((1 - sigma) * (alpha + K * sigma) + sum(Poch2(freq - sigma)))
  )

  cat("Model: Dirichlet Process",
      paste0("\t Abundance: ", out[1]),
      paste0("\t Richness: ", out[2]),
      paste0("\t alpha: ", round(c(out[3]),4)),
      paste0("\t sigma: ", round(c(out[4]),4)),
      paste0("\t Sample coverage: ", round(out[5],4)),
      paste0("\t Expected species after additional ",  n, " samples: ", out[6]),
      paste0("\t Posterior Gini diversity: ",  round(out[7],4)),
      sep = "\n"
  )
  invisible(out)
}

#' Plot for the Pitman-Yor model
#' @param object An object of class \code{\link[ssm]{ssm}}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @details The method...
#'
#' @export
#'
plot.PY <- function(object, ...) {
  n <- sum(object$frequencies)
  alpha <- object$param[1]
  sigma <- object$param[2]


  M_l <- as.numeric(table(factor(object$frequencies, levels = 1:n)))
  # P_l <- M_l / sum(M_l)

  idx <- 1:(which.min(M_l) - 1) # which(P_l > 0)

  data_plot <- data.frame(Size = idx, M_l = M_l[idx], Theoretical = expected_m_py(idx, n = n, sigma = sigma, alpha = alpha))
  p <- ggplot(data = data_plot, aes(x = Size, y = M_l)) +
    geom_point() +
    geom_line(aes(y = Theoretical), color = "blue", linetype = "dashed") +
    scale_y_log10() +
    scale_x_log10() +
    theme_bw() +
    xlab("l") +
    ylab(expression(M[l]))
  p
}
