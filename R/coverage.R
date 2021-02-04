#' @export
coverage <- function(x, ...) {
  UseMethod("coverage", x)
}

#' @export
coverage.numeric <- function(frequencies) {
  n <- sum(frequencies)
  m1 <- sum(frequencies == 1)
  1 - m1 / n
}

#' @export
coverage.DP <- function(object, ...) {
  alpha <- object$param[1]
  freq <- object$frequencies
  n <- sum(freq)
  K <- length(freq)

  n / (alpha + n)
}

#' @export
coverage.PY <- function(object, ...) {
  alpha <- object$param[1]
  sigma <- object$param[2]
  freq <- object$frequencies
  n <- sum(freq)
  K <- length(freq)

  (n - sigma * K) / (alpha + n)
}
