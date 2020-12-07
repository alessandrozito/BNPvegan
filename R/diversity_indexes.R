
Simpson <- function(frequencies){
  freq_rel <- frequencies / sum(frequencies)
  out <- sum(freq_rel^2)
}


Shannon <- function(frequencies){
  freq_rel <- frequencies / sum(frequencies)
  - sum(freq_rel*log(freq_rel))
}
