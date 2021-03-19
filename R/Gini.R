#' @export
Gini <- function(x, ...) {
  UseMethod("Gini", x)
}

#' @export
Gini.numeric <- function(object, ...) {
  freq_rel <- frequencies / sum(frequencies)
  out <- 1 - sum(freq_rel^2)
  out
}

#' @export
#'
Gini.DP <- function(object, ...) {
  Poch2 <- function(x) x * (x + 1)

  alpha <- object$param[1]
  freq <- object$frequencies
  n <- sum(freq)
  out <- 1 - 1 / Poch2(alpha + n) * (alpha + sum(Poch2(freq)))
  out
}

#' @export
#'
Gini.PY <- function(object, ...) {
  Poch2 <- function(x) x * (x + 1)

  alpha <- object$param[1]
  sigma <- object$param[2]
  freq <- object$frequencies
  n <- sum(freq)
  K <- length(freq)

  out <- 1 - 1 / Poch2(alpha + n) * ((1 - sigma) * (alpha + K * sigma) + sum(Poch2(freq - sigma)))
  out
}


