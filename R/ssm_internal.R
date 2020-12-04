expected_cl_py <- function(n, sigma, alpha) {
  n <- as.integer(n)
  if (sigma == 0) {
    out <- alpha * sum(1 / (alpha - 1 + 1:n))
  } else {
    out <- 1 / sigma * exp(lgamma(alpha + sigma + n) - lgamma(alpha + sigma) - lgamma(alpha + n) + lgamma(alpha + 1)) - alpha / sigma
  }

  return(out)
}
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
expected_m_py <- Vectorize(expected_m_py, vectorize.args = "m")

frequency_check_PY <- function(frequencies) {
  n <- sum(frequencies)
  fit_PY <- max_EPPF_PY(frequencies)
  alpha <- fit_PY$par[1]
  sigma <- fit_PY$par[2]


  M_l <- as.numeric(table(factor(frequencies, levels = 1:n)))
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
