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

#' @export
diversity <- function(x, ...) {
  UseMethod("diversity", x)
}



#' @export
#'
diversity.data.frame <- function(dataset) {
  if (is.null(dim(dataset))) dataset <- matrix(dataset, nrow = 1)
  tab <- matrix(0, nrow(dataset), 4)
  rownames(tab) <- rownames(dataset)
  colnames(tab) <- c( "Gini", "Normalized Gini", "Entropy", "Normalized entropy")
  for (i in 1:nrow(dataset)) {
    frequencies <- dataset[i, ]
    frequencies <- frequencies[frequencies > 0]
    tab[i, ] <- c(
      Gini(frequencies),
      Gini_norm(frequencies),
      Shannon(frequencies),
      Shannon_norm(frequencies),
    )
  }
  tab
}

#' @export
#'
diversity.numeric <- function(frequencies) {
  tab <- matrix(0, 1, 4)
  colnames(tab) <- c("Gini", "Normalized Gini", "Entropy", "Normalized entropy")
  tab[1, ] <- c(
      Gini(frequencies),
      Gini_norm(frequencies),
      Shannon(frequencies),
      Shannon_norm(frequencies),
    )
  tab
}

