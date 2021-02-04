# Loglikelihood for LL3
logLik_LL3 <- function(d, beta_0, beta_1, beta_2, X) {
  # Compute the log-likelihood for a sample from the posterior
  beta <- c(beta_0, beta_1, beta_2)
  eta <- c(X %*% beta) # odds
  loglik <- sum(d * eta - log(1 + exp(eta)))
  return(-loglik)
}

# Loglikilhood maximizer for LL3
max_logLik_LL3 <- function(d, X) {
  start <- c(1, 0, -1)
  out <- nlminb(
    start = start,
    objective = function(par) logLik_LL3(d, beta_0 = par[1], beta_1 = par[2], par[3], X = X),
    upper = c(Inf, -1e-7, -1e-7)
  )
  return(out)
}

# Loglikelihood of the Weibull
Weibull_logLik <- function(d, phi, lambda) {
  # Prob_new = phi^(t-1)^lambda
  t_one <- which(d == 1)
  t_one <- t_one[-1] # <=== the first observation is always a 1, needs to be removed
  t_zero <- which(d == 0)
  logLik <- log(phi) * sum((t_one - 1)^lambda) + sum(log(1 - exp(log(phi) * ((t_zero - 1)^lambda))))
  return(-logLik)
}

# Loglikilhood maximizer for the Weibull
max_logLik_Weibull <- function(d) {
  start <- c(0.5, 1)
  out <- nlminb(
    start = start,
    objective = function(par) Weibull_logLik(d, phi = par[1], lambda = par[2]),
    lower = c(0, 0), upper = c(1 - 9e-7, Inf)
  )
  return(out)
}



# Discovery probability of LL3
prob_LL3 <- function(n, alpha, sigma, phi) {
  return(alpha * phi^n / (alpha * phi^n + n^(1 - sigma)))
}


prob_LL3_squared <- function(n, alpha, sigma, phi) {
  p <- prob_LL3(n, alpha, sigma, phi)
  return(p^2)
}

# Discovery probability for Weibull
prob_Weibull <- function(n, phi, lambda) {
  nl <- n^lambda
  return(phi^nl)
}

prob_Weibull_squared <- function(n, phi, lambda) {
  p <- prob_Weibull(n, phi, lambda)
  return(p^2)
}


# Compute the asymptotic moments for k_inf
moments_Kinf <- function(par, n, k, model) {
  # k = number of species already detected
  # n = number of samples observed

  if (model == "LL3") {
    alpha <- par[1]
    sigma <- par[2]
    phi <- par[3]
    if (phi == 1) {
      if (sigma >= 0) {
        E_KInf <- Inf
        Var_KInf <- Inf
      } else {
        E_KInf <- alpha^(1 / (1 - sigma)) * pi / ((1 - sigma) * sin(pi / (1 - sigma)))
        Var_KInf <- E_KInf - cubature::cubintegrate(f = prob_LL3_squared, alpha = alpha, sigma = sigma, phi = phi, lower = 0, upper = Inf)$integral + 0.5
      }
    } else {
      E_KInf <- cubature::cubintegrate(f = prob_LL3, alpha = alpha, sigma = sigma, phi = phi, lower = n, upper = Inf)$integral + k
      Var_KInf <- E_KInf - cubature::cubintegrate(f = prob_LL3_squared, alpha = alpha, sigma = sigma, phi = phi, lower = n, upper = Inf)$integral + 0.5
    }
  } else if (model == "Weibull") {
    phi <- par[1]
    lambda <- par[2]
    if (n == 0 & k == 0) {
      E_KInf <- gamma(1 + 1 / lambda) * (-1 / log(phi))^(1 / lambda)
      Var_KInf <- E_KInf - cubature::cubintegrate(f = prob_Weibull_squared, phi = phi, lambda = lambda, lower = 0, upper = Inf)$integral + 0.5
    } else {
      E_KInf <- cubature::cubintegrate(f = prob_Weibull, phi = phi, lambda = lambda, lower = n, upper = Inf)$integral + k
      Var_KInf <- E_KInf - cubature::cubintegrate(f = prob_Weibull_squared, phi = phi, lambda = lambda, lower = n, upper = Inf)$integral + 0.5
    }
  }

  return(c("E_KInf" = ceiling(E_KInf), "sd_KInf" = sqrt(Var_KInf)))
}


# Compute the expected rarefaction (cumulative sum of the discovery probabilities)
expected_rarefaction <- function(N, alpha, sigma, phi) {
  n <- c(1:N) - 1
  cumsum(prob_LL3(n, alpha, sigma, phi))
}
