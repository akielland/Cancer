###############################################################################
## create outputs for report 2 For Ridge where there is no feature selection ##
###############################################################################

library(ggplot2)

report02_out_ridge <- function(object){
  cor_vec <- object$cor_vec
  n_models <- length(cor_vec)
  cat("number of models fitted:", n_models, "\n")
  fr <- sum(is.na(cor_vec))/n_models
  cat("Fraction of model fits with no selected genes:", fr, "\n")
  cat("\n")
  
  ## Correlations
  cat("CORRELATIONS RESULTS", "\n")
  cor_mean <- mean(na.omit(cor_vec))
  cat("Mean:", cor_mean, "\n") 
  cor_median <- median(cor_vec, na.rm=TRUE)
  cat("Median:", cor_median, "\n")
  cor_var <- var(cor_vec, na.rm=TRUE)
  cat("Variance:", cor_var, "\n")
  cor_sd <- sd(na.omit(cor_vec)) 
  cat("st.dev.:", cor_sd, "\n")
  
  # Histogram correlations
  correlations_finite <- cor_vec[is.finite(cor_vec)]
  cor_df <- data.frame(correlation = correlations_finite)
  
  h <- ggplot(cor_df, aes(x=correlation)) +
    geom_histogram(bins = 30, color = "black", fill = "white") +
    xlab("Correlation") +
    ylab("Frequency") +
    ggtitle("Histogram of Correlation Values")
  show(h)
  
  ## MSE
  cat("MSE RESULTS", "\n")
  MSE_vec <- object$MSE_vec
  MSE_mean <- mean(na.omit(MSE_vec))
  cat("Mean:", MSE_mean, "\n") 
  MSE_median <- median(MSE_vec, na.rm=TRUE)
  cat("Median:", MSE_median, "\n")
  MSE_var <- var(MSE_vec, na.rm=TRUE)
  cat("Variance:", MSE_var, "\n")
  MSE_sd <- sd(na.omit(MSE_vec)) 
  cat("st.dev.:", MSE_sd)
  
  
  # Histogram MSE
  MSE_df <- data.frame(MSE = MSE_vec)
  
  h <- ggplot(MSE_df, aes(x=MSE)) +
    geom_histogram(bins = 30, color = "black", fill = "white") +
    xlab("MSE") +
    ylab("Frequency") +
    ggtitle("Histogram of MSE Values")
  show(h)
  
  out <- list(cor_mean=cor_mean,
              cor_median=cor_median,
              cor_var=cor_var,
              cor_sd=cor_sd,
              MSE_mean=MSE_mean,
              MSE_median=MSE_median,
              MSE_var=MSE_var,
              MSE_sd=MSE_sd)
  return(out)
}
