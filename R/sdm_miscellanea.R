#' Obtain the sequence of new discoveries form
#'
#' @param sequence A sequence of species in order of appearence
#'
#' @return A sequence of binary variables
extract_discoveries <- function(sequence) {
  as.numeric(!duplicated(sequence))
}

#' Build a random sequence of observed species form a vector of counts
#'
#' @param abundances Vector of counts of the species observed
#'
#' @return vector
sample_sequence <- function(abundances) {

  # Extract the crude sequence
  sequence <- rep(1:length(abundances), times = abundances)

  # Randomize it
  sequence <- sample(sequence, size = sum(abundances), replace = FALSE)

  return(sequence)
}

# THIS FUNCTION IS HIGHLY INEFFICIENT AND SHOULD BE REWRITTEN

rLL3 <- function(L, alpha, sigma, phi) {
  species <- c()
  D <- c()
  for (n in 0:(L - 1)) {
    d <- rbinom(1, 1, prob_LL3(n, alpha, sigma, phi))
    D <- c(D, d)
    if (d == 1) {
      species <- c(species, as.character(n))
    } else {
      species <- c(species, sample(species, 1, replace = TRUE))
    }
  }
  return(list("d" = D, "abundances" = unname(table(species))))
}
