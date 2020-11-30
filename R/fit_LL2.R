#' Fit the two parameter log-logistic as in Zito et al. (2020)
#'
#' @param discoveries Sequence of discoveries indicators
#'
#' @return Parameters of the log logistic distribution
#' @export
#'
#' @examples
#' fit_LL2(c(1, 0, 1, 1, 0, 1, 0, 0, 0))
fit_LL2 <- function(discoveries) {

  # Step 1 - introduce useful quantities
  N <- length(discoveries)
  X <- as.matrix(cbind(1, log(1:(N - 1))))

  # Step 2 - Fit the logistic regression via fastglm
  glmfit <- fastglm::fastglmPure(y = discoveries[-1], x = X, family = binomial("logit"), method = 3)
  pars <- c("alpha" = exp(glmfit$coefficients[1]), "sigma" = 1 + glmfit$coefficients[2])

  # Step 3 - Check if the parameters are in bound
  if (pars[2] >= 1) {
    # Need to re-estimate the model using the one parameter log-logistic. This is equivalent to the Dirichlet process
    glmfit <- fastglm::fastglmPure(y = discoveries[-1], x = as.matrix(X[, 1]), family = binomial("logit"), method = 3, offset = -X[, 2])
    pars <- c("alpha" = exp(glmfit$coefficients[1]), "sigma" = 0)
    Out_of_bound <- TRUE
  } else {
    Out_of_bound <- FALSE
  }
  # Return the values needed
  return(list("discoveries" = discoveries, "parameters" = pars, "Out_of_bound" = Out_of_bound))
}
