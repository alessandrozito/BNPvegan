logLik_LL3 <- function(d, beta_0, beta_1, beta_2, X) {
  # Compute the log-likelihood for a sample from the posterior
  beta <- c(beta_0, beta_1, beta_2)
  eta <- c(X %*% beta) # odds
  loglik <- sum(d * eta - log(1 + exp(eta)))
  return(-loglik)
}

max_logLik_LL3 <- function(d, X){
  start <- c(1, 0, -1)
  out <- nlminb(start = start,
                objective = function(par) logLik_LL3(d, beta_0 = par[1], beta_1 = par[2], par[3], X=X),
                upper = c(Inf, -1e-16, -1e-7)
  )
  return(out)
}

prob_LL3 <- function(n, alpha, sigma, phi){
  return(alpha*phi^n/(alpha*phi^n + n^(1-sigma)))
}

prob_LL3_squared <- function(n, alpha, sigma, phi){
  p <- prob_LL3(n, alpha, sigma, phi)
  return(p^2)
}

moments_Kinf <- function(alpha, sigma, phi){
  if(phi==1){
    if(sigma >= 0){
      E_KInf <- Inf
      Var_KInf <- Inf
    } else {
      E_KInf <- alpha^(1/(1-sigma))*pi /((1-sigma)*sin(pi/(1-sigma)))
      Var_KInf <- E_KInf - integrate(f = prob_LL3_squared, alpha = alpha, sigma = sigma, phi = phi, lower = 0, upper = Inf)$value + 0.5
    }
  } else {
    E_KInf <- integrate(f = prob_LL3, alpha = alpha, sigma = sigma, phi = phi, lower = 0, upper = Inf)$value
    Var_KInf <- E_KInf - integrate(f = prob_LL3_squared, alpha = alpha, sigma = sigma, phi = phi, lower = 0, upper = Inf)$value + 0.5
  }
  return(c("E_KInf" = ceiling(E_KInf), "sd_KInf" = sqrt(Var_KInf)))
}


expected_rarefaction <- function(N, alpha, sigma, phi){
  n <- c(1:N)-1
  cumsum(prob_LL3(n, alpha, sigma, phi))
}


#fit_true <- logit_regression(y = LL3out$d[-1], X = X)
#c(exp(fit_true$par[1]), 1 + fit_true$par[2], exp(fit_true$par[3]))
#loglik <- max(fit_true$Convergence[,2])

#plot(cumsum(fit$discoveries), type = "l", col = "red", lwd=2, lty=1)
#lines(cumsum(LL3out$d),lwd=2, lty=1, col = "blue")

#hist(fit$resampling_output$param[,4])
#abline(v = loglik,col = "blue", lwd=3, lty=2)
#abline(v = fit$loglik, col="red", lwd=3, lty=2)




