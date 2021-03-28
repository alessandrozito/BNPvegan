#' Frequency of abundances
#' @param object An object of class \code{numeric} or \code{ssm}
#' @param ... Additional parameters
#' @export
freq_abundance <- function(object, ...) {
  UseMethod("freq_abundance", object)
}

#' Frequency of abundances for a vector of species counts
#'
#' @param object An object of class \code{numeric} containing species abundances
#' @param r frequencies of abundances
#' @param ... Additional parameters
#'
#' @export
#' @details Return the table of the species counts
freq_abundance.numeric <- function(object, r, ...) {
  res <- rep(0, length(r))
  names(res) <- r
  for(j in 1:length(r)){
    res[j] <- sum(object == r[j])
  }
  return(res)
}

#' Model-based expected frequency of abundances for the Dirichlet process
#'
#' @param object An object of class \code{ssm, DP}
#' @param r frequencies of abundances
#' @param ... Additional parameters
#'
#' @export
#' @details Specify the model based expected frequency of abundances for a
#' Dirichlet process
#' @examples r<- c(1:10)
#' fit <- ssm(fungalOTU, model = "DP")
#' freq_abundance(fit, r)
freq_abundance.DP <- function(object, r, ...) {
  tab <- expected_m_dp(r = r, n = sum(object$abundances), alpha = coef(object)[1])
  names(tab) <- r
  return(tab)
}

#' Model-based expected frequency of abundances for the Pirman-Yor process
#'
#' @param object An object of class \code{ssm, PY}
#' @param r frequencies of abundances
#' @param ... Additional parameters
#'
#' @export
#' @details Specify the model based expected frequency of abundances for a
#' Pitman-Yor process
#' @examples r<- c(1:10)
#' fit <- ssm(fungalOTU, model = "PY")
#' freq_abundance(fit, r)
freq_abundance.PY <- function(object, r, ...) {
  tab <- expected_m_py(r = r, n = sum(object$abundances), alpha = coef(object)[1], sigma = coef(object)[2])
  names(tab) <- r
  return(tab)
}

