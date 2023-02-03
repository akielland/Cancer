###########################################################
## Repeated cross-validation for testing the Ensembel models ##
###############################################################
##
## Test against the fold
## Cross-validation (5-fold) used for training/tuning lambda and selecting features
## output: selected genes; proliferation.score correlation; SEM

## Post Lasso model is made at the bottom
library(stats)
library(caret)

library(glmnet)


library(ggplot2)

###########################################################
##  Repeated cross validation ##
###########################################################

pca_results <- function(fm, data){
  X_ <- model.matrix(fm, data = data)[,-1]  
  pca_results <- prcomp(X_, scale = TRUE) # create PCA object
  # Get the number of components that explain at least X % of the variance
  n_components <- sum(pca_results$sdev^2 / sum(pca_results$sdev^2) >= 0.1)
  # Extract the first n_components principal components
  df_first_pca <- data.frame(pca_results$x[, 1:n_components])
  return(df_first_pca)
}



# Function: repeated k-fold cross validation
lasso_rep_cv <- function(df_data, func=lasso_sample, folds=5, repeats=200, method="pearson") {
  n_models <- repeats * folds
  print(n_models)
  cor_vec <- rep(NA, n_models)
  MSE_vec <- rep(NA, n_models)
  #coef_matrix <- matrix(NA, nrow = n_models, ncol = ncol(df_data[,-1]))
  #colnames(coef_matrix) <- colnames(df_data[, -1])
  
  count <- 1
  # Repeat the cross-validation process
  for (i in 1:repeats) {
    # Create the folds for evaluating the performance
    kf <- caret::createFolds(df_data[,1], k = folds, list = TRUE, returnTrain = TRUE)
    # Loop through the folds
    for (j in 1:folds) {
      # Get the training and testing data
      train_data <- df_data[kf[[j]],]
      test_data <- df_data[-kf[[j]],]
      
      immune_inf <- pca_results(fm_immune_inf, data=prolif_771genes)
      prolif <- pca_results(fm_prolif, data=prolif_771genes)
      ...
      
      fm <- Y ~ immune_inf + prolif + ...

      # Fit the function on the training data and get results
      fit <- func(train_data)      
      fit <- lm(fm, data = train_data)
      fit <- XGBoost_sample(fm, train_data)
      fit <- lasso_sample(fm, train_data)
      

     # coef_matrix[coef_matrix_row_index, ] <- coef(fit, s = "lambda.min")[-1]
      
      pred = predict(fit, newx = as.matrix(test_data)[,-1], type = "response", s = "lambda.min")  
      cor_vec[coef_matrix_row_index]  <- suppressWarnings(cor(pred, test_data[,1], method=method))
      MSE_vec[coef_matrix_row_index] <- mean((pred - test_data[,1])^2)
      
      cat(coef_matrix_row_index, "")
      coef_matrix_row_index <- coef_matrix_row_index + 1
    }
  }
  return(list(cor_vec=cor_vec, coef_matrix=coef_matrix, MSE_vec=MSE_vec))
}