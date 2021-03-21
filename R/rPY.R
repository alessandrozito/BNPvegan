#' Random samples from a Pitman-Yor process
#'
#' @param size Number of observations to draw
#' @param alpha Centrality parameter
#' @param sigma Strenght parameter
#'
#' @return  A vector of frequencies
#' @export
#'
#' @examples
rPY <- function(size, alpha, sigma){
  if(size == 1){
    return(1)
  }
  # Sample a Pitman-Yor of a given size
  counts <- c(1)
  for (n in 2:size) {
    K <- length(counts)
    prob_new <- (alpha + K * sigma) / (n - 1 + alpha)
    prob_old <- (counts - sigma) / (n - 1 + alpha)

    # Sample the observation
    x <- sample(x = c(0, 1:K), size = 1, prob = c(prob_new, prob_old), replace = FALSE)

    if (x == 0) {
      counts <- c(counts, 1)
    } else {
      counts[x] <- counts[x] + 1
    }
  }
  return(counts)
}


#' Random samples from a Pitman-Yor process
#'
#' @param size Number of observations to draw
#' @param alpha Centrality parameter
#'
#' @return A vector of frequencies
#' @export
rDP <- function(size, alpha){
  return(rPY(size = size, alpha = alpha, sigma = 0))
}
