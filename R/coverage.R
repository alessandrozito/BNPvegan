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


#' @export
rcoverage <- function(x, ...) {
  UseMethod("rcoverage", x)
}

#' @export
rcoverage.DP <- function(object, R = 1000, ...) {
  alpha <- object$param[1]
  freq <- object$frequencies
  n <- sum(freq)

  rbeta(R, n, alpha)
}

#' @export
rcoverage.PY <- function(object, R = 1000, ...) {
  alpha <- object$param[1]
  sigma <- object$param[2]
  freq <- object$frequencies
  n <- sum(freq)
  K <- length(freq)

  rbeta(R, n - sigma*K, alpha + sigma*K)
}
