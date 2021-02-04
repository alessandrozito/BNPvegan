#' @export
coverage <- function(frequencies){
  n <- sum(n)
  m1 <- sum(frequencies == 1)
  1 - m1 / n
}
