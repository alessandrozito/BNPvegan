expected_cl_py <- function(n, sigma, alpha) {
  n <- as.integer(n)
  if (sigma == 0) {
    out <- alpha * sum(1 / (alpha - 1 + 1:n))
  } else {
    out <- 1 / sigma * exp(lgamma(alpha + sigma + n) - lgamma(alpha + sigma) - lgamma(alpha + n) + lgamma(alpha + 1)) - alpha / sigma
  }

  return(out)
}
#' @export
expected_cl_py <- Vectorize(expected_cl_py, vectorize.args = "n")


expected_m_dp <- function(m, n, alpha) {
  out <- log(alpha) + lchoose(n, m) + lgamma(m) + lgamma(alpha + n - m) - lgamma(alpha + n)
  exp(out)
}
expected_m_dp <- Vectorize(expected_m_dp, vectorize.args = "m")

expected_m_py <- function(m, n, sigma, alpha) {
  out <- log(alpha) + lchoose(n, m) + lgamma(m - sigma) - lgamma(1 - sigma) - lgamma(alpha + n) + lgamma(alpha) + lgamma(alpha + sigma + n - m) - lgamma(alpha + sigma)
  exp(out)
}

#' @export
expected_m_py <- Vectorize(expected_m_py, vectorize.args = "m")


extrapolate_cl_py <- function(m, K, n, sigma, alpha) {
  n <- as.integer(n)
  #if (sigma == 0) {
  #  out <- alpha * sum(1 / (alpha - 1 + 1:n))
  #} else {
    out <- (K + alpha / sigma) * (exp(lgamma(alpha + n +  sigma + m) - lgamma(alpha + n + sigma) - lgamma(alpha + n + m) + lgamma(alpha + n)) - 1)
  #}

  return(out)
}

#' @export
extrapolate_cl_py <- Vectorize(extrapolate_cl_py, vectorize.args = "m")
