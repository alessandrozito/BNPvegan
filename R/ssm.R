#' Species sampling model
#'
#' @param frequencies A \code{K}-dimensional vector of frequencies
#' @param model Model to fit. Available models are \code{DP} (Dirichlet Process) and \code{PY} (Pitman-Yor process)
#'
#' @return An object of class \code{ssm}
#'
#' @export
ssm <- function(frequencies, model) {
  frequencies <- frequencies[frequencies > 0]
  if (model == "DP") {
    fit <- max_EPPF_DP(frequencies)
    out <- list(frequencies = frequencies, param = fit$par, logLik = -fit$objective)
    class(out) <- c("ssm", "DP")
    return(out)
  }

  if (model == "PY") {
    fit <- max_EPPF_PY(frequencies)
    out <- list(frequencies = frequencies, param = fit$par, logLik = -fit$objective)
    class(out) <- c("ssm", "PY")
    return(out)
  }
}

#' @export
coef.ssm <- function(object, ...) {
  object$param
}

#' @export
logLik.ssm <- function(object, ...) {
  object$logLik
}


#' @export
#'
predict.DP <- function(object, m, ...) {
  alpha <- object$param[1]
  freq <- object$frequencies
  n <- sum(freq)
  K <- length(freq)
  extrapolate_cl_py(m = m, n = n, K = K, sigma = 0, alpha = alpha)
}

#' @export
#'
predict.PY <- function(object, newdata = NULL, ...) {
  alpha <- object$param[1]
  sigma <- object$param[2]
  freq <- object$frequencies
  n <- sum(freq)
  K <- length(freq)
  extrapolate_cl_py(m = m, n = n, K = K, sigma = sigma, alpha = alpha)
}

#' @export
#'
summary.DP <- function(object, ...) {
  alpha <- object$param[1]
  freq <- object$frequencies
  n <- sum(freq)
  K <- length(freq)
  Expected <- round(extrapolate_cl_py(m = n, n = n, K = K, sigma = 0, alpha = alpha)) - K
  Gini <- Gini(object)

  out <- t(c(alpha, object$logLik))
  colnames(out) <- c("alpha", "loglik")

  cat("Model:",
    "\t Dirichlet process",
    "\nQuantities:",
    paste0("\t Abundance: ", n),
    paste0("\t Richness: ", K),
    paste0("\t Estimated sample coverage: ", round(coverage(object), 4)),
    paste0("\t Expected species after additional ", n, " samples: ", Expected + K),
    paste0("\t New expected species after additional ", n, " samples: ", Expected),
    paste0("\t Posterior Gini diversity: ", round(Gini, 4)),
    "\nParameters:",
    paste0("\t ", knitr::kable(out, "simple")),
    sep = "\n"
  )
}


#' @export
#'
summary.PY <- function(object, ...) {
  alpha <- object$param[1]
  sigma <- object$param[2]
  freq <- object$frequencies
  n <- sum(freq)
  K <- length(freq)
  Expected <- round(extrapolate_cl_py(m = n, n = n, K = K, sigma = sigma, alpha = alpha)) - K
  Gini <- Gini(object)

  out <- t(c(alpha, sigma, object$logLik))
  colnames(out) <- c("alpha", "sigma", "loglik")

  cat("Model:",
    "\t Pitman-Yor process",
    "\nQuantities:",
    paste0("\t Abundance: ", n),
    paste0("\t Richness: ", K),
    paste0("\t Estimated sample coverage: ", round(coverage(object), 4)),
    paste0("\t Expected species after additional ", n, " samples: ", Expected + K),
    paste0("\t New expected species after additional ", n, " samples: ", Expected),
    paste0("\t Posterior Gini diversity: ", round(Gini, 4)),
    "\nParameters:",
    paste0("\t ", knitr::kable(out, "simple")),
    sep = "\n"
  )
}


#' Plot for the DP model
#' @param object An object of class \code{\link[ssm]{ssm}}.
#' @param m additional sample to predict. Valid only for \code{type = "extrapolation"}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @details The method...
#'
#' @export
#'
plot.DP <- function(object, type = "rarefaction", m = NULL, ...) {
  n <- sum(object$frequencies)
  alpha <- object$param[1]
  K <- length(object$frequencies)

  if (type == "freq") {
    M_l <- as.numeric(table(factor(object$frequencies, levels = 1:n)))
    # P_l <- M_l / sum(M_l)
    idx <- 1:(which.min(M_l) - 1) # which(P_l > 0)
    data_plot <- data.frame(Size = idx, M_l = M_l[idx], Theoretical = expected_m_py(idx, n = n, sigma = 0, alpha = alpha))
    p <- ggplot(data = data_plot, aes(x = Size, y = M_l)) +
      geom_point() +
      geom_line(aes(y = Theoretical), color = "blue", linetype = "dashed") +
      scale_y_log10() +
      scale_x_log10() +
      theme_bw() +
      xlab("l") +
      ylab(expression(M[l]))
    return(p)
  } else if (type == "coverage") {
    p <- ggplot() +
      xlim(qbeta(0.001, n, alpha), qbeta(0.999, n, alpha)) +
      geom_function(fun = function(x) dbeta(x, n, alpha)) +
      theme_bw() +
      xlab("Sample coverage") +
      ylab("Density")
    return(p)
  } else if (type == "rarefaction") {
    data_plot <- data.frame(n = 1:n, rar = rarefaction(object))
    p <- ggplot(data = data_plot, aes(x = n, y = rar)) +
      geom_line() +
      theme_bw() +
      xlab("n") +
      ylab("Rarefaction")+
      ggplot2::facet_wrap(~"Rarefaction curve")
    return(p)
  } else if (type == "extrapolation") {
    if(is.null(m)){
      m = n
    }
    data_plot <- data.frame(n = c(1:(m+n)), rar = c(rarefaction(object), extrapolation(object, 1:m)))
    p <- ggplot(data = data_plot, aes(x = n, y = rar)) +
      geom_line() +
      geom_vline(xintercept = n, linetype = "dashed") +
      theme_bw() +
      xlab("n") +
      ylab("Rarefaction and extrapolation")+
      ggplot2::facet_wrap(~"Rarefaction and extrapolation curve")
    return(p)
  }
}


#' Plot for the Pitman-Yor model
#' @param object An object of class \code{\link[ssm]{ssm}}.
#' @param m additional sample to predict. Valid only for \code{type = "extrapolation"}. By default, it is equal to the sample size.
#' @param ... Further arguments passed to or from other methods.
#'
#' @details The method...
#'
#' @export
#'
plot.PY <- function(object, type = "rarefaction", m = NULL, ...) {
  n <- sum(object$frequencies)
  K <- length(object$frequencies)
  alpha <- object$param[1]
  sigma <- object$param[2]

  if (type == "freq") {
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
    return(p)
  } else if (type == "coverage") {
    p <- ggplot() +
      xlim(qbeta(0.001, n - sigma * K, alpha + sigma * K), qbeta(0.999, n - sigma * K, alpha + sigma * K)) +
      geom_function(fun = function(x) dbeta(x, n - sigma * K, alpha + sigma * K)) +
      theme_bw() +
      xlab("Sample coverage") +
      ylab("Density")
    return(p)
  } else if (type == "rarefaction") {
    data_plot <- data.frame(n = 1:n, rar = rarefaction(object))
    p <- ggplot(data = data_plot, aes(x = n, y = rar)) +
      geom_line() +
      theme_bw() +
      xlab("n") +
      ylab("Rarefaction")+
      facet_wrap(~"Rarefaction curve")
    return(p)
  } else if (type == "extrapolation") {
    if(is.null(m)){
      m = n
    }
    data_plot <- data.frame(n = c(1:(m+n)), rar = c(rarefaction(object), extrapolation(object, 1:m)))
    p <- ggplot(data = data_plot, aes(x = n, y = rar)) +
      geom_line() +
      geom_vline(xintercept = n, linetype = "dashed") +
      theme_bw() +
      xlab("n") +
      ylab("Rarefaction and extrapolation")+
      facet_wrap(~"Rarefaction and extrapolation curve")
    return(p)
  }
}
