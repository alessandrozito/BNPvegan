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


#' Diversity measures
#'
#'
#' @param frequencies A \code{K}-dimensional vector of frequencies
#' #'
#' @return Collection of diversity indices
#'
#' @export
#'
diversity <- function(dataset) {
  if (is.null(dim(dataset))) dataset <- matrix(dataset, nrow = 1)
  tab <- matrix(0, nrow(dataset), 6)
  rownames(tab) <- rownames(dataset)
  colnames(tab) <- c("Simpson", "Gini", "Normalized Gini", "Entropy", "Normalized entropy", "Rare species ratio")
  for (i in 1:nrow(dataset)) {
    frequencies <- dataset[i, ]
    frequencies <- frequencies[frequencies > 0]
    tab[i, ] <- c(
      Simpson(frequencies),
      Gini(frequencies),
      Gini_norm(frequencies),
      Shannon(frequencies),
      Shannon_norm(frequencies),
      unique_ratio(frequencies)
    )
  }
  print(knitr::kable(tab))
  invisible(tab)
}
