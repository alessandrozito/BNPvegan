#' Obtain the sequence of new discoveries form
#'
#' @param sequence A sequence of species in order of appearence
#'
#' @return A sequence of binary variables
#' @export
#'
#' @examples
#' Species <-  c("SpeciesA", "SpeciesB", "SpeciesB", "SpeciesC")
#' extract_discoveries(Species)
extract_discoveries <- function(sequence) {
  species_appeared <- NULL
  m <- length(sequence)
  discoveries <- rep(0, m)
  for (i in 1:m) {
    if (!(sequence[i] %in% species_appeared)) {
      discoveries[i] <- 1
      species_appeared <- c(species_appeared, sequence[i])
    }
  }
  return(discoveries)
}


#' Build a random sequence of observed species form a vector of counts
#'
#' @param frequencies Vector of counts of the species observed
#'
#' @return vector
#' @export
#'
#' @examples
#' sample_sequence(c(1, 2, 45, 7))
sample_sequence <- function(frequencies) {

  n <- sum(frequencies)
  K <- length(frequencies)

  # Extract the crude sequence
  sequence <- rep(1:K, each=frequencies)

  # Randomize it
  sequence <- sample(sequence, size = n, replace = FALSE)

  return(sequence)
}
