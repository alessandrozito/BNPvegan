#' Fit the sequential discovery model of Zito et al. (2020+)
#'
#' @param frequencies Counts of the species
#' @param n_resamples number of accumulation curves to generate
#' @param verbose Whether to monitor the function or not
#'
#' @details This function.
#' @export
#'
sdm <- function(frequencies, n_resamples = 500L, verbose = TRUE){

  # Initialize an empty matrix for the parameters
  param <- matrix(NA, nrow = n_resamples, ncol = 3)
  colnames(param) <- c("alpha","sigma", "phi")
  # List to store the position of the discoveries in the list
  discoveries_indexes <- matrix(NA, nrow = n_resamples, ncol = length(frequencies))
  verbose_step <- round(n_resamples/10)

  for(i in 1:n_resamples){
    # Obtain the discovery sequence
    seqD <- sample_sequence(frequencies)
    d <- extract_discoveries(seqD)
    discoveries_indexes[i, ] <- which(d==1)

    # Fit the three-parameter log-logistic
    fit <- max_logLik_LL3(d)
    param[i,] <- c(exp(fit$par[1]), 1 + fit$par[2], exp(fit$par[3]))

    # Monitor output
    if(verbose){
      if(i%%verbose_step==0){
        cat(paste0("Number of accumulation curves fitted: ", i, " [",round(100*i/n_resamples) ,'%]'), sep ="\n")
      }
    }
  }
  out <- list(frequencies = frequencies,
              n_resamples = n_resamples,
              discoveries_indexes = discoveries_indexes,
              param = param,
              loglik = -fit$objective)
  class(out) <- "sdm"
  return(out)
}


#' Summary for a species discovery model
#'
#' @param object An object of class \code{\link[sdm]{sdm}}.
#' @param ... additional parameters
#'
#' @details A function to summarize the three parameter log-logistic output
#' @export
summary.sdm <- function(object, plot = TRUE, ...){

  # Sample abundance and richness
  abundance <- sum(object$frequencies)
  richness <- length(object$frequencies)

  # Count the number of divergent accumulation curves
  n_div <- sum(object$param[,2]>0 & object$param[,3]==1)

  # Summary of the parameters
  pars_tab <- matrix(NA, nrow = 3, ncol = 6)
  colnames(pars_tab) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
  rownames(pars_tab) <- c("alpha", "sigma", "phi")
  pars_tab[1,-4] <- unname(quantile(object$param[,1]))
  pars_tab[2,-4] <- unname(quantile(object$param[,2]))
  pars_tab[3,-4] <- unname(quantile(object$param[,3]))
  pars_tab[1,4] <- mean(object$param[,1])
  pars_tab[2,4] <- mean(object$param[,2])
  pars_tab[3,4] <- mean(object$param[,3])



  # Compute asymptotic richness
  EK <- apply(object$param, 1, FUN = function(p) expected_Kinf(alpha = p[1], sigma = p[2], phi = p[3]))
  asymp_tab <-rbind(summary(EK[EK<Inf]), (summary(EK[EK<Inf])-richness)/richness)
  rownames(asymp_tab) <- c("Asymptotic Richness", '% increase')

  if(plot == TRUE){
    param_plot <- ggplot2::ggplot(data = tidyr::gather(data.frame(object$param)))+
      ggplot2::geom_histogram(ggplot2::aes(x=value),alpha=0.45, bins = 30,color = "black")+
      ggplot2::theme_bw()+
      ggplot2::facet_wrap(~key, scales = "free")+
      ggplot2::xlab("")
    richness_plot <- ggplot2::ggplot(data.frame("Asy_richness" = EK[EK<Inf]))+
      ggplot2::geom_histogram(ggplot2::aes(x=Asy_richness), alpha = 0.45, bins = 30,color = "black")+
      ggplot2::theme_bw()+
      ggplot2::facet_wrap(~"Asymptotic Richness")+
      ggplot2::xlab("N. species")+
      ggplot2::geom_segment(x=richness, xend = richness, y = 0, yend = Inf, color = "red", linetype = "dashed")
  } else {
    param_plot<- NULL
    richness_plot<-NULL
  }
  # Print the summary
  cat("Model:",
      "\t Three-parameter log-logistic (LL3)",
      "\nQuantities:",
      paste0("\t Abundance: ", abundance),
       paste0("\t Richness: ", richness),
     paste0("\t Number of divergent accumulation curves: ", n_div, " out of ", object$n_resamples, " sampled [", round(100*n_div/object$n_resamples,2),'%]'),
     "\nEst. Richness for non-divergent accumulation curves:",
     knitr::kable(round(asymp_tab,2), "simple"),
     "\nParameters:",
     knitr::kable(pars_tab, "simple"),
     sep= "\n" )
  print(param_plot)
  print(richness_plot)

  # Output
  out <-list(Abundance = abundance,
             Richness = richness,
             n_resamples = object$n_resamples,
             n_div = n_div,
             Asymp_richness = EK,
             Asymp_richness_summary = asymp_tab,
             Asymp_richness_plot = richness_plot,
             param_plot = param_plot,
             param_summary = pars_tab)
  class(out) <- union("summary", class(object))
  return(invisible(out))
}

#' Print method for the summary
#'
#' @param x object of class \code{\link[summary.sdm]{summary.sdm}}
#' @param ... other parameters
#' @export
print.summary.sdm <- function(x, ...){
  cat("Model:",
      "\t Three-parameter log-logistic (LL3)",
      "\nQuantities:",
      paste0("\t Abundance: ", x$Abundance),
      paste0("\t Richness: ", x$Richness),
      paste0("\t Number of divergent accumulation curves: ", x$n_div, " out of ", x$n_resamples, " sampled [", round(100*x$n_div/x$n_resamples,2),'%]'),
      "\nEst. Richness for non-divergent accumulation curves:",
      knitr::kable(round(x$Asymp_richness_summary,2), "simple"),
      "\nParameters:",
      knitr::kable(x$param_summary, "simple"),
      sep= "\n" )
  #print(x$param_plot)
  #print(x$Asymp_richness_plot)
  invisible(x)
}


#' Compute rarefaction curves
#'
#' @param object An object of class \code{\link[sdm]{sdm}}.
#' @param newdata a new dataset
#' @param ... other parameters

#' @export
plot.sdm <- function(object, newdata = NULL, ...){
  # Step 1 - compute the average accumulation curve by averaging across re-samples
  N <- sum(object$frequencies)
  K_mat <- matrix(0, nrow = object$n_resamples, ncol = N)
  for(i in 1:object$n_resamples){
    K_mat[i,object$discoveries_indexes[i,]] <- 1
    K_mat[i,] <- cumsum(K_mat[i,])
  }
  avg_accumulation <- colMeans(K_mat)

  # Step 2 - compute the average rarefaction curve
  avg_rarefaction <- rowMeans(apply(object$param, 1, function(p) expected_rarefaction(N = N, alpha = p[1], sigma = p[2], phi = p[3])))
  df <- data.frame("n" = c(1:N), "accum"= avg_accumulation, "rar" = avg_rarefaction)

  if(nrow(df)>100){
    # Plot one every 100 points
    df <- df[seq(1, nrow(df), 100),]
  }

  ggplot2::ggplot(df)+
    ggplot2::geom_point(ggplot2::aes(x = n, y = accum), shape = 1)+
    ggplot2::geom_line(ggplot2::aes(x = n, y = rar), color = "red")+
    ggplot2::theme_bw()+
    ggplot2::facet_wrap(~"Average rarefaction")+
    ggplot2::ylab("Number of species")
}

