#' Fit the three parameter log-logistic as in Zito et al. (2020)
#'
#' @param discoveries Sequence of discovery indicators
#'
#' @return A list, including the parameters
#' @export
#'
#' @examples
#' fit_LL3(c(1, 0, 1, 1, 0, 1, 0, 0, 0))
fit_LL3 <- function(discoveries) {

  # Step 1 - introduce useful quantities
  N <- length(discoveries)
  X <- as.matrix(cbind(1, log(1:(N - 1)), c(1:(N - 1)))) # Predictors are 1,log(n),n

  # Step 2 - Fit the logistic regression via fastglm
  glmfit <- fastglm::fastglmPure(y = discoveries[-1], x = X, family = binomial("logit"), method = 3)
  pars <- c("alpha" = exp(glmfit$coefficients[1]), "sigma" = 1 + glmfit$coefficients[2], "phi" = exp(glmfit$coefficients[3]))

  # Step 4 - Check if the parameters are in bound
  if (pars[3] > 1 | pars[2] >= 1) {
    # Need to re-estimate the model using the two parameter log-logistic.
    out <- fit_LL2(discoveries)
    out$parameters <- c(out$parameters, "phi" = 1)
    out$Out_of_bound <- TRUE

    return(out)
  } else {
    return(list("discoveries" = discoveries, "parameters" = pars, "Out_of_bound" = FALSE))
  }
}
