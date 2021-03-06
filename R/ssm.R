#' Species sampling model
#'
#' @param abundances A \code{K}-dimensional vector of abundances
#' @param model Model to fit. Available models are \code{DP} (Dirichlet Process) and \code{PY} (Pitman-Yor process)
#'
#' @return An object of class \code{ssm}
#'
#' @export
ssm <- function(abundances, model) {
  abundances <- abundances[abundances > 0]
  if (model == "DP") {
    fit <- max_EPPF_DP(abundances)
    out <- list(abundances = abundances, param = fit$par, logLik = -fit$objective)
    class(out) <- c("ssm", "DP")
    return(out)
  }

  if (model == "PY") {
    fit <- max_EPPF_PY(abundances)
    out <- list(abundances = abundances, param = fit$par, logLik = -fit$objective)
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
  freq <- object$abundances
  n <- sum(freq)
  K <- length(freq)
  extrapolate_cl_py(m = m, n = n, K = K, sigma = 0, alpha = alpha)
}

#' @export
#'
predict.PY <- function(object, m, ...) {
  alpha <- object$param[1]
  sigma <- object$param[2]
  freq <- object$abundances
  n <- sum(freq)
  K <- length(freq)
  extrapolate_cl_py(m = m, n = n, K = K, sigma = sigma, alpha = alpha)
}

#' @export
#'
summary.DP <- function(object, ...) {
  alpha <- object$param[1]
  freq <- object$abundances
  n <- sum(freq)
  K <- length(freq)
  Expected <- round(extrapolate_cl_py(m = n, n = n, K = K, sigma = 0, alpha = alpha)) - K
  Gini <- Gini(object)

  out <- t(c(alpha, object$logLik))
  colnames(out) <- c("alpha", "loglik")

  cat("Model:",
      "\t Dirichlet process (DP)",
      "\nQuantities:",
      paste0("\t Abundance: ", n),
      paste0("\t Richness: ", K),
      paste0("\t Estimated sample coverage: ", round(coverage(object), 4)),
      paste0("\t Posterior Gini diversity: ", round(Gini, 4)),
      "\nExtrapolations:",
      paste0("\t Expected species after additional ", n, " samples: ", Expected + K),
      paste0("\t New expected species after additional ", n, " samples: ", Expected),
      "\nAsymptotics:",
      paste0("\t Expected species at infinity: Infinite"),
      paste0("\t Standard deviation at infinity: NA"),
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
  freq <- object$abundances
  n <- sum(freq)
  K <- length(freq)
  Expected <- round(extrapolate_cl_py(m = n, n = n, K = K, sigma = sigma, alpha = alpha)) - K
  Gini <- Gini(object)

  out <- t(c(alpha, sigma, object$logLik))
  colnames(out) <- c("alpha", "sigma", "loglik")

  cat("Model:",
    "\t Pitman-Yor process (PY)",
    "\nQuantities:",
    paste0("\t Abundance: ", n),
    paste0("\t Richness: ", K),
    paste0("\t Estimated sample coverage: ", round(coverage(object), 4)),
    paste0("\t Posterior Gini diversity: ", round(Gini, 4)),
    "\nExtrapolations:",
    paste0("\t Expected species after additional ", n, " samples: ", Expected + K),
    paste0("\t New expected species after additional ", n, " samples: ", Expected),
    "\nAsymptotics:",
    paste0("\t Expected species at infinity: Infinite"),
    paste0("\t Standard deviation at infinity: NA"),
    "\nParameters:",
    paste0("\t ", knitr::kable(out, "simple")),
    sep = "\n"
  )
}


#' Plot for the DP model
#'
#' @param x An object of class \code{ssm}.
#' @param m additional sample to predict. Valid only for \code{type = "extrapolation"}.
#' @param type Type of plot. Available options are:  \code{"rarefaction"},\code{"extrapolation"}, \code{"freq"} and \code{"coverage"}
#' @param plot_sample Whether to add the sample-based rarefaction curve. Available for \code{type = "rarefaction"} only
#' @param verbose Monitor the output of the sample-based rarefaction curve. Available for \code{type = "rarefaction"} only
#' @param n_points Number of points to plot
#' @param ... Further arguments passed to or from other methods.
#'
#' @details The method produces summary plots for the Dirichlet output
#'
#' @export
#'
plot.DP <- function(x, type = "rarefaction", plot_sample = TRUE, m = NULL, verbose = TRUE, n_points = 100, ...) {
  n <- sum(x$abundances)
  alpha <- x$param[1]
  K <- length(x$abundances)

  if (type == "freq") {
    M_l <- as.numeric(table(factor(x$abundances, levels = 1:n)))
    # P_l <- M_l / sum(M_l)
    idx <- 1:(which.min(M_l) - 1) # which(P_l > 0)
    data_plot <- data.frame("Size" = idx, "M_l" = M_l[idx], "Theoretical" = expected_m_py(idx, n = n, sigma = 0, alpha = alpha))
    p <- ggplot(data = data_plot, aes_string(x = "Size", y = "M_l")) +
      geom_point() +
      geom_line(aes_string(y = "Theoretical"), color = "blue", linetype = "dashed") +
      scale_y_log10() +
      scale_x_log10() +
      theme_bw() +
      xlab("r") +
      ylab(expression(M[r]))+
      facet_wrap(~"Frequency-of-abundances plot")
    return(p)
  } else if (type == "coverage") {
    p <- ggplot() +
      xlim(qbeta(0.001, n, alpha), qbeta(0.999, n, alpha)) +
      geom_function(fun = function(x) dbeta(x, n, alpha)) +
      theme_bw() +
      xlab("Sample coverage") +
      ylab("Density")+
      facet_grid(~"Posterior sample coverage")
    return(p)
  } else if (type == "rarefaction") {

    if(plot_sample == TRUE){
      # Model-based rarefaction
      rar <- rarefaction(x)
      accum <- rarefaction(as.integer(x$abundances), verbose = verbose)
      df <- data.frame("n" = 1:n, "rar" = rar, "accum" = accum)

      if (nrow(df) > n_points) {
        seqX <- 1:nrow(df)
        seqY <- split(seqX, sort(seqX %% n_points))
        df <- df[unlist(lapply(seqY, function(a) utils::tail(a, 1))), ]
      }

      p <- ggplot(df) +
        geom_point(aes_string(x = "n", y = "accum"), shape = 1) +
        geom_line(aes_string(x = "n", y = "rar"), color = "red", size = 0.9) +
        theme_bw() +
        facet_wrap(~"Rarefaction curve") +
        ylab(expression(K[n]))

    } else {
      # Model-based rarefaction
      rar <- rarefaction(x)
      data_plot <- data.frame("n" = 1:n, "rar" = rarefaction(x))

      p <- ggplot(data = data_plot, aes_string(x = "n", y = "rar")) +
        geom_line(color = "red", size = 0.9) +
        theme_bw() +
        xlab("n") +
        ylab("Rarefaction")+
        facet_wrap(~"Rarefaction curve")+
        ylab(expression(K[n]))
    }
    return(p)

  } else if (type == "extrapolation") {
    if(is.null(m)){
      m = n
    }
    if(plot_sample == TRUE){
      # Model-based rarefaction
      rar <- rarefaction(x)
      accum <- rarefaction(as.integer(x$abundances), verbose = verbose)
      ext <- extrapolation(x, 1:m)
      df <- data.frame("n" = 1:(n+m), "curve" = c(rar,ext) , "accum" = c(accum, rep(NA, length(ext))))

      if (nrow(df) > n_points) {
        seqX <- 1:nrow(df)
        seqY <- split(seqX, sort(seqX %% n_points))
        df <- df[unlist(lapply(seqY, function(a) utils::tail(a, 1))), ]
      }

      p <- ggplot(df) +
        geom_point(aes_string(x = "n", y = "accum"), shape = 1, na.rm=TRUE) +
        geom_line(aes_string(x = "n", y = "curve"), color = "red", size = 0.9) +
        theme_bw() +
        facet_wrap(~"Rarefaction and extrapolation") +
        geom_segment(x = n, xend = n, y = 0, yend = Inf, linetype = "dashed") +
        ylab(expression(K[n]))

    } else {

      data_plot <- data.frame("n" = c(1:(m+n)), "rar" = c(rarefaction(x), extrapolation(x, 1:m)))
      p <- ggplot(data = data_plot, aes_string(x = "n", y = "rar")) +
        geom_line() +
        geom_segment(x = n, xend = n, y = 0, yend = Inf, linetype = "dashed") +
        theme_bw() +
        xlab("n") +
        ylab("Rarefaction and extrapolation")+
        facet_wrap(~"Rarefaction and extrapolation curve")
    }

    return(p)
  }
}


#' Plot for the Pitman-Yor model
#'
#' @param x An object of class \code{ssm}.
#' @param m additional sample to predict. Valid only for \code{type = "extrapolation"}. By default, it is equal to the sample size.
#' @param type Type of plot. Available options are:  \code{"rarefaction"},\code{"extrapolation"}, \code{"freq"} and \code{"coverage"}
#' @param plot_sample Whether to add the sample-based rarefaction curve. Available for \code{type = "rarefaction"} only
#' @param verbose Monitor the output of the sample-based rarefaction curve. Available for \code{type = "rarefaction"} only
#' @param n_points Number of points to plot
#' @param ... Further arguments passed to or from other methods.
#'
#' @details The method produces summary plots for the Pitman-Yor output
#'
#' @export
#'
plot.PY <- function(x, type = "rarefaction", plot_sample = TRUE, m = NULL, verbose = TRUE, n_points = 100, ...) {
  n <- sum(x$abundances)
  K <- length(x$abundances)
  alpha <- x$param[1]
  sigma <- x$param[2]

  if (type == "freq") {
    M_l <- as.numeric(table(factor(x$abundances, levels = 1:n)))
    # P_l <- M_l / sum(M_l)
    idx <- 1:(which.min(M_l) - 1) # which(P_l > 0)
    data_plot <- data.frame("Size" = idx, "M_l" = M_l[idx],
                            "Theoretical" = expected_m_py(idx, n = n, sigma = sigma, alpha = alpha))
    p <- ggplot(data = data_plot, aes_string(x = "Size", y = "M_l")) +
      geom_point() +
      geom_line(aes_string(y = "Theoretical"), color = "blue", linetype = "dashed") +
      scale_y_log10() +
      scale_x_log10() +
      theme_bw() +
      xlab("r") +
      ylab(expression(M[r]))+
      facet_wrap(~"Frequency-of-abundances plot")
    return(p)
  } else if (type == "coverage") {
    p <- ggplot() +
      xlim(qbeta(0.001, n - sigma * K, alpha + sigma * K), qbeta(0.999, n - sigma * K, alpha + sigma * K)) +
      geom_function(fun = function(x) dbeta(x, n - sigma * K, alpha + sigma * K)) +
      theme_bw() +
      xlab("Sample coverage") +
      ylab("Density")+
      facet_grid(~"Posterior sample coverage")
    return(p)
  } else if (type == "rarefaction") {
    if(plot_sample == TRUE){
      # Model-based rarefaction
      rar <- rarefaction(x)
      accum <- rarefaction(as.integer(x$abundances), verbose = verbose)
      df <- data.frame("n" = 1:n, "rar" = rar, "accum" = accum)

      if (nrow(df) > n_points) {
        seqX <- 1:nrow(df)
        seqY <- split(seqX, sort(seqX %% n_points))
        df <- df[unlist(lapply(seqY, function(a) utils::tail(a, 1))), ]
      }

      p <- ggplot(df) +
        geom_point(aes_string(x = "n", y = "accum"), shape = 1) +
        geom_line(aes_string(x = "n", y = "rar"), color = "red", size = 0.9) +
        theme_bw() +
        facet_wrap(~"Rarefaction curve") +
        ylab(expression(K[n]))

    } else {
      # Model-based rarefaction
      rar <- rarefaction(x)
      data_plot <- data.frame("n" = 1:n, "rar" = rarefaction(x))

      p <- ggplot(data = data_plot, aes_string(x = "n", y = "rar")) +
        geom_line(color = "red", size = 0.9) +
        theme_bw() +
        xlab("n") +
        ylab("Rarefaction")+
        facet_wrap(~"Rarefaction curve")+
        ylab(expression(K[n]))
    }
    return(p)
  } else if (type == "extrapolation") {
    if(is.null(m)){
      m = n
    }
    if(plot_sample == TRUE){
      # Model-based rarefaction
      rar <- rarefaction(x)
      accum <- rarefaction(as.integer(x$abundances), verbose = verbose)
      ext <- extrapolation(x, 1:m)
      df <- data.frame("n" = 1:(n+m), "curve" = c(rar,ext) , "accum" = c(accum, rep(NA, length(ext))))

      if (nrow(df) > n_points) {
        seqX <- 1:nrow(df)
        seqY <- split(seqX, sort(seqX %% n_points))
        df <- df[unlist(lapply(seqY, function(a) utils::tail(a, 1))), ]
      }

      p <- ggplot(df) +
        geom_point(aes_string(x = "n", y = "accum"), shape = 1, na.rm=TRUE) +
        geom_line(aes_string(x = "n", y = "curve"), color = "red", size = 0.9) +
        theme_bw() +
        facet_wrap(~"Rarefaction and extrapolation") +
        geom_segment(x = n, xend = n, y = 0, yend = Inf, linetype = "dashed") +
        ylab(expression(K[n]))

    } else {

      data_plot <- data.frame("n" = c(1:(m+n)), "rar" = c(rarefaction(x), extrapolation(x, 1:m)))
      p <- ggplot(data = data_plot, aes_string(x = "n", y = "rar")) +
        geom_line() +
        geom_segment(x = n, xend = n, y = 0, yend = Inf, linetype = "dashed") +
        theme_bw() +
        xlab("n") +
        ylab("Rarefaction and extrapolation")+
        facet_wrap(~"Rarefaction and extrapolation curve")
    }
    return(p)
  }
}
