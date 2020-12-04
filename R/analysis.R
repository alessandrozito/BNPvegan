expected_cl_py <- function(n, alpha, theta) {
  n <- as.integer(n)
  if (alpha == 0) {
    out <- theta * sum(1 / (theta - 1 + 1:n))
  } else {
    out <- 1 / alpha * exp(lgamma(theta + alpha + n) - lgamma(theta + alpha) - lgamma(theta + n) + lgamma(theta + 1)) - theta / alpha
  }

  return(out)
}

expected_m_dp <- function(m, n, theta) {
  out <- log(theta) + lchoose(n, m) + lgamma(m) + lgamma(theta + n - m) - lgamma(theta + n)
  exp(out)
}
expected_m_dp <- Vectorize(expected_m_dp, vectorize.args = "m")

expected_m_py <- function(m, n, alpha, theta) {
  out <- log(theta) + lchoose(n, m) + lgamma(m - alpha) - lgamma(1 - alpha) - lgamma(theta + n) + lgamma(theta) + lgamma(theta + alpha + n - m) - lgamma(theta + alpha)
  exp(out)
}
expected_m_py <- Vectorize(expected_m_py, vectorize.args = "m")

frequency_check_PY <- function(frequencies) {
  n <- sum(frequencies)
  fit_PY <- max_EPPF_PY(frequencies)
  theta <- fit_PY$par[1]
  alpha <- fit_PY$par[2]


  M_l <- as.numeric(table(factor(frequencies, levels = 1:n)))
  # P_l <- M_l / sum(M_l)

  idx <- 1:(which.min(M_l) - 1) # which(P_l > 0)

  data_plot <- data.frame(Size = idx, M_l = M_l[idx], Theoretical = expected_m_py(idx, n = n, alpha = alpha, theta = theta))
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
