#' Random sampler from ad Dirichlet distribution
#'
#' @param size Number of samples to draw
#' @param alpha Centrality parameter. Has to be a vector
#'
#' @return Samples from a dirichlet distribution
rdirichlet <- function(size, alpha){
  # Draw the beta distributions
  n <- length(alpha)
  if(n==1){
    cat("Specify a longer vector parameters alpha")
  } else {
    if(size == 1){
      rdir <- rgamma(n = n, shape = alpha)
      rdir <- rdir/sum(rdir)
    } else {
      # sample a whole matrix
      rdir <- matrix(rgamma(n = n*size, shape = alpha), ncol = n, nrow = size, byrow = TRUE)
      rdir <- rdir/rowSums(rdir)
    }
    return(rdir)
  }
}
