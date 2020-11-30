#' Obtain the sequence of new discoveries form
#'
#' @param sequence A sequence of species in order of appearence
#'
#' @return A sequence of binary variables
#' @export
#'
#' @examples
#' extract_discoveries(c("sp1", "sp1", "sp3", "sp4"))
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
