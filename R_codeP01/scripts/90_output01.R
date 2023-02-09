#################################
## create outputs for report 2 ##
#################################

library(ggplot2)

report02_out <- function(object){
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
  
  ## Most prevalent Features
  # Order the features based on their selection frequency
  coef_matrix <- object$coef_matrix
  frequency <- data.frame(Feature = colnames(coef_matrix), Frequency = colSums(coef_matrix != 0) / (n_models))
  frequency <- frequency[order(frequency$Frequency, decreasing = TRUE),]
  
  # Bar plot of the selection frequency of the features
  n_best <- 50
  h <- ggplot(frequency[1:n_best ,], aes(x = Frequency, y = reorder(Feature, Frequency))) +
    geom_bar(stat = "identity") +
    xlab("Selection Frequency") +
    ylab("Features") +
    ggtitle("Selection Frequency of Features") +
    theme(axis.text.y = element_text(angle = 0, hjust = 0))
  show(h)
  
  # extract the top features
  perc_best <- 0.5 # how often they where selected in percentage
  top_features_with_index = which(colSums(coef_matrix != 0) >= perc_best * n_models)
  top_feature_names = colnames(coef_matrix)[top_features_with_index]
  cat("\n")
  cat("Features selected 50% or more times:", "\n")
  cat(top_feature_names, "\n")
  
  num_features_to_keep <- 20 
  # count the frequency of each feature in coef_matrix
  counts <- colSums(coef_matrix != 0)
  # sort the features based on their frequency
  sorted_features <- names(sort(counts, decreasing = TRUE))
  cat("Top 20 featrues:", "\n")
  print(sorted_features[1:num_features_to_keep])
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
