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
#' @param species_counts Vector of counts of the species observed
#' @param seed Optional seed. Default is 42
#'
#' @return vector
#' @export
#'
#' @examples
#' build_random_sequences(c(1, 2, 45, 7))
extract_sequence <- function(species_counts, seed = 42) {

  # Turn the count function into an integer list
  species_counts <- round(species_counts)

  # Extract the names of the species in the sequence
  if (is.null(names(names(species_counts)))) {
    names(species_counts) <- paste0("sp", as.character(c(1:length(species_counts))))
  }
  species_names <- names(species_counts)

  # Filter out the observations with counts equal to 0
  species_names <- species_names[species_counts > 0]
  species_counts <- species_counts[species_counts > 0]

  # Extract the crude sequence
  sequence <- c()
  for (j in 1:length(species_counts)) {
    subseq <- rep(x = species_names[j], times = as.numeric(species_counts[j]))
    sequence <- c(sequence, subseq)
  }

  # Randomize it
  set.seed(seed)
  sequence <- sample(sequence, size = length(sequence), replace = FALSE)

  return(sequence)
}
