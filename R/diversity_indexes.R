Simpson <- function(frequencies) {
  freq_rel <- frequencies / sum(frequencies)
  out <- sum(freq_rel^2)
}

Shannon <- function(frequencies) {
  freq_rel <- frequencies / sum(frequencies)
  -sum(freq_rel * log(freq_rel))
}

Gini <- function(frequencies) {
  1 - Simpson(frequencies)
}

Gini_norm <- function(frequencies) {
  K <- length(frequencies)
  G <- Gini(frequencies)
  G * K / (K - 1)
}

Shannon_norm <- function(frequencies) {
  K <- length(frequencies)
  H <- Shannon(frequencies)
  H / log(K)
}

unique_ratio <- function(frequencies) {
  mean(frequencies == 1)
}
