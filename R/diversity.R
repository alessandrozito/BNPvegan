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

diversities <- function(dataset) {
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


#' Diversity measures
#'
#'
#' @param frequencies A \code{K}-dimensional vector of frequencies
#' @param rel a
#' @param plot a
#'
#' @return Collection of diversity indices
#'
#' @export
#'

freq_of_freq <- function(dataset, rel = FALSE, plot = FALSE) {
  if (is.null(dim(dataset))) dataset <- matrix(dataset, nrow = 1)

  K <- max(dataset)
  tab <- matrix(0, nrow(dataset), K)
  rownames(tab) <- rownames(dataset)
  colnames(tab) <- 1:K
  for (i in 1:nrow(dataset)) {
    frequencies <- dataset[i, ]
    frequencies <- factor(frequencies[frequencies > 0], levels = 1:K)
    tab[i, ] <- as.numeric(table(frequencies))
    if (rel) tab[i, ] <- tab[i, ] / length(frequencies)
  }
  tab <- t(tab)

  if (plot) {
    data_plot <- reshape2::melt(tab)
    data_plot$value[data_plot$values == 0] <- NA
    colnames(data_plot)[2] <- "Sample"

    p <- ggplot(data = data_plot, aes(x = Var1, y = value, col = Sample)) +
      geom_point(size = 0.8, alpha = 0.5) +
      geom_smooth(size = 0.8, se = F) +
      scale_y_log10() +
      scale_x_log10() +
      theme_bw() +
      xlab("Frequency") +
      ylab("Frequency of frequencies")
    return(p)
  }
  tab
}
