logEPPF_PY <- function(theta, alpha, frequencies) {
  if (any(alpha < 0, theta <= - alpha + 1e-04)) {
    return(-Inf)
  }

  # Sample size
  n <- sum(frequencies)
  # Number of frequencies
  K <- length(frequencies)
  # Loglikelihood
  loglik <- sum(log(theta + 1:(K - 1) * alpha)) - lgamma(theta + n) + lgamma(theta + 1) + sum(lgamma(frequencies - alpha)) - K * lgamma(1 - alpha)

  loglik
}

logEPPF_DP <- function(theta, frequencies) {

  # Sample size
  n <- sum(frequencies)
  # Number of frequencies
  K <- length(frequencies)

  # Loglikelihood
  loglik <- K * log(theta) - lgamma(theta + n) + lgamma(theta)

  loglik
}

max_EPPF_PY <- function(frequencies) {
  start <- c(1, 0.5) # Initialization of the maximization algorithm
  out <- nlminb(
    start = start,
    function(param) -logEPPF_PY(theta = param[1], alpha = param[2], frequencies = frequencies),
    lower = c(-Inf, 1e-16), upper = c(Inf, 1 - 1e-10)
  )
  return(out)
}

max_EPPF_DP <- function(frequencies) {
  start <- 1 # Initialization of the maximization algorithm
  out <- nlminb(
    start = start,
    function(param) -logEPPF_DP(theta = param, frequencies = frequencies),
    lower = 1e-10, upper = Inf
  )
  return(out)
}
