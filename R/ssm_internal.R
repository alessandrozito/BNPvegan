#' @import stats ggplot2
#' @importFrom Rcpp evalCpp sourceCpp
#' @importFrom vegan diversity
#' @useDynLib BNPvegan
#'

logEPPF_PY <- function(alpha, sigma, abundances) {
  if (any(sigma < 0, alpha <= -sigma + 1e-04)) {
    return(-Inf)
  }

  # Sample size
  n <- sum(abundances)
  # Number of abundances
  K <- length(abundances)
  # Loglikelihood
  loglik <- sum(log(alpha + 1:(K - 1) * sigma)) - lgamma(alpha + n) + lgamma(alpha + 1) + sum(lgamma(abundances - sigma)) - K * lgamma(1 - sigma)

  loglik
}

logEPPF_DP <- function(alpha, abundances) {

  # Sample size
  n <- sum(abundances)
  # Number of abundances
  K <- length(abundances)

  # Loglikelihood
  loglik <- K * log(alpha) - lgamma(alpha + n) + lgamma(alpha)

  loglik
}

max_EPPF_PY <- function(abundances) {
  start <- c(1, 0.5) # Initialization of the maximization algorithm
  out <- nlminb(
    start = start,
    function(param) -logEPPF_PY(alpha = param[1], sigma = param[2], abundances = abundances),
    lower = c(-Inf, 1e-16), upper = c(Inf, 1 - 1e-10)
  )
  return(out)
}

max_EPPF_DP <- function(abundances) {
  start <- 1 # Initialization of the maximization algorithm
  out <- nlminb(
    start = start,
    function(param) -logEPPF_DP(alpha = param, abundances = abundances),
    lower = 1e-10, upper = Inf
  )
  return(out)
}
