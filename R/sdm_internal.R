logLik_LL3 <- function(d, beta_0, beta_1, beta_2) {
  # Compute the log-likelihood for a sample from the posterior
  beta <- c(beta_0, beta_1, beta_2)
  d <- d[-1] # Remove the first discovery
  n <- length(d)
  X <- cbind(1, log(1:n), c(1:n))
  eta <- c(X %*% beta) # odds
  loglik <- sum(d * eta - log(1 + exp(eta)))
  return(-loglik)
}

max_logLik_LL3 <- function(d){
  start <- c(1, 0, -1)
  out <- nlminb(start = start,
                objective = function(par) logLik_LL3(d, beta_0 = par[1], beta_1 = par[2], par[3]),
                upper = c(Inf, 1e-16, 0)
  )
  return(out)
}

prob_LL3 <- function(n, alpha, sigma, phi){
  return(alpha*phi^n/(alpha*phi^n + n^(1-sigma)))
}

expected_Kinf <- function(alpha, sigma, phi){
  if(phi==1){
    if(sigma >= 0){
      E_KInf <- Inf
    } else {
      E_KInf <- ceiling(alpha^(1/(1-sigma))*pi /((1-sigma)*sin(pi/(1-sigma))))
    }
  } else {
    E_KInf <- ceiling(integrate(f = prob_LL3, alpha = alpha, sigma = sigma, phi = phi, lower = 0, upper = Inf)$value)
  }
  return(E_KInf)
}


expected_rarefaction <- function(N, alpha, sigma, phi){
  n <- c(1:N)-1
  cumsum(prob_LL3(n, alpha, sigma, phi))
}





