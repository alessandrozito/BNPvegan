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
#' @param frequencies Vector of counts of the species observed
#'
#' @return vector
sample_sequence <- function(frequencies) {

  # Extract the crude sequence
  sequence <- rep(1:length(frequencies), times=frequencies)

  # Randomize it
  sequence <- sample(sequence, size = sum(frequencies), replace = FALSE)

  return(sequence)
}
