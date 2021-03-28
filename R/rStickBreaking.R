#' Random Sampler for a Stick breaking construction of a Pitman-Yor
#'
#' @param alpha Centrality parameter of the Pitman-Yor
#' @param sigma Tail parameter of the Pitman-Yor
#' @param truncation_method Howe to truncate the stick breaking.
#'        Available options are \code{"dynamic"} and \code{"static"}.
#'        The first one stops sampling the from the distribution when the sampled weights
#'        sum is higher than the \code{prob_cutoff}. The second instead performs a
#'        hard truncation at the specified \code{static_cutoff}.
#' @param prob_cutoff Probability cutoff to stop the sampler. Valid for \code{truncation_method = "dynamic"}
#' @param static_cutoff Integer cutoff to stop the sampler. Valid for \code{truncation_method = "static"}
#'
#' @return A vector of weights representig the Pitman-Yor weights. Notice that the last element in the vector represents the aggregate weight after the truncation.
rStickBreaking <- function(alpha, sigma, truncation_method = "dynamic", prob_cutoff = 0.99, static_cutoff = 1000){
  if(truncation_method == "dynamic"){
    h <- 1
    betas <- rbeta(n = 1, 1 - sigma, alpha + h*sigma)
    weights <- betas
    while(sum(weights) < prob_cutoff){
      h <- h + 1
      betas <- c(betas, rbeta(n = 1, 1 - sigma, alpha + h*sigma))
      weights <- c(weights, betas[h] * prod(1 - betas[1:(h-1)]))
    }
    # Add the remaining weights
    weights <- c(weights, 1 - sum(weights))
  } else if (truncation_method == "static"){
    h <- c(1:static_cutoff)
    betas <- rbeta(n = static_cutoff, 1 - sigma, alpha + h*sigma)
    weights <- c(betas[1], rep(0, static_cutoff))
    for(i in 2:static_cutoff){
      weights[i] <- betas[i] * prod(1-betas[1:(i-1)])
    }
    # Add the remaining weights
    weights[static_cutoff + 1] <- 1 - sum(weights[1:static_cutoff])
  }
  return(weights)
}
