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

#' Summary for the DP model
#'
#'
#' Gini heterogeneity
#'
#' @param object An object of class \code{\link[ssm]{ssm}}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @details The method...
#'
#' @export
#'
summary.DP <- function(object, ...) {
  tab <- cbind(Simpson_Emp = Simpson(object$frequencies),
               Shannon_Emp = Simpson(object$frequencies),
               Simpson_MB  = 1 / (object$param[1] + 1))
  knitr::kable(tab)
}


#' Heterogeneity for the PY model
#'
#'
#' Gini heterogeneity
#'
#' @param object An object of class \code{\link[ssm]{ssm}}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @details The method...
#'
#' @export
#'
summary.PY <- function(object, ...) {
  tab <- cbind(Abundance   = sum(object$frequencies),
               Richness   = length(object$frequencies),
               Simpson_Emp = Simpson(object$frequencies),
               Shannon_Emp = Shannon(object$frequencies),
               Rare_Emp = mean(object$frequencies == 1),
               Simpson_MB  = (1 - object$param[2]) / (object$param[1] + 1))
  tab
}

#' @param object An object of class \code{\link[ssm]{ssm}}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @details The method...
#'
#' @export
#'
plot.PY <- function(object,...) {
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
