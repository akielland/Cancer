###################################################################
## Repeated cross-validation for testing the Stacking - Ensemble models ##
###################################################################
##
## Test against the fold
## Cross-validation (5-fold) used for training/tuning lambda and selecting features
## output: selected genes; proliferation.score correlation; SEM

library(caret)
library(glmnet)
library(ggplot2)

###########################################################
##  Repeated cross validation ##
###########################################################




t <- lasso_rep_cv(prolif_771genes, folds=2, repeats=1, method="pearson")

# Function: repeated k-fold cross validation
lasso_rep_cv <- function(df_data, folds=5, repeats=1, method="pearson") {
  n_models <- repeats * folds
  print(n_models)
  cor_vec <- rep(NA, n_models)
  MSE_vec <- rep(NA, n_models)
  coef_matrix <- matrix(NA, nrow = n_models, ncol = ncol(df_data[,-1]))
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

      # PCA's:
      immune_inf <- pca_results(char_immune_inf, data=train_data, percentage=0.9)
      prolif <- pca_results(char_prolif, data=train_data, percentage=0.9)
      ER_singaling <- pca_results(char_ER_signaling, data=train_data, percentage=0.9)
      
      train_PCA <- data.frame(Y=train_data[1], immune_inf, prolif, ER_singaling)
      
      # Fit the function on the training data and get results
      # model_fit <- lm(Y ~ ., data = train_PCA)
      # 
      # fit <- func(train_data)      
      # fit <- lm(fm, data = train_data)
      # fit <- XGBoost_sample(fm, train_data)
      # 
      pred <- lasso_block(train_PCA, test_data)
      # pred = predict(fit, newx = as.matrix(test_data)[,-1], type = "response", s = "lambda.min")  
    
     # coef_matrix[coef_matrix_row_index, ] <- coef(fit, s = "lambda.min")[-1]

      cor_vec[count]  <- suppressWarnings(cor(pred, test_data[,1], method=method))
      # MSE_vec[coef_matrix_row_index] <- mean((pred - test_data[,1])^2)
      # 
      # cat(coef_matrix_row_index, "")
      # coef_matrix_row_index <- coef_matrix_row_index + 1
    }
  }
  return(cor_vec)
  #return(list(cor_vec=cor_vec, coef_matrix=coef_matrix, MSE_vec=MSE_vec))
}
