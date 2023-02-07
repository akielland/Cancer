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

library(caret)
library(glmnet)

# Load the data
data(iris)

# Number of repetitions of k-fold cross-validation
reps <- 10

# Number of folds in k-fold cross-validation
folds <- 5

# Create a vector to store the MSE values
mse_values <- rep(NA, reps)

# Repeat k-fold cross-validation
for (i in 1:reps) {
  # Split the data into training and testing sets using k-fold cross-validation
  set.seed(i)
  cv_folds <- createFolds(iris$Species, k = folds, returnTrain = TRUE)
  cv_indices <- lapply(cv_folds, function(x) which(x))
  
  # Initialize a list to store the fitted models for each fold
  cv_models <- list()
  
  # Loop over each fold
  for (j in 1:folds) {
    # Extract the training and testing data for this fold
    train <- iris[cv_indices[[j]], ]
    test <- iris[-cv_indices[[j]], ]
    
    # Fit models using different algorithms
    models <- list(lm = train(Species ~ ., data=train, method="lm"),
                   rf = train(Species ~ ., data=train, method="rf"),
                   svm = train(Species ~ ., data=train, method="svmRadial"))
    
    # Make predictions using the fitted models
    preds <- lapply(models, function(x) predict(x, newdata=test))
    
    # Stack the predictions
    preds <- do.call(cbind, preds)
    
    # Fit a meta-model to the stacked predictions
    meta_model <- train(Species ~ ., data=data.frame(preds), method="glmnet")
    
    # Make final predictions using the meta-model
    final_preds <- predict(meta_model, newx=preds)
    
    # Store the fitted meta-model for this fold
    cv_models[[j]] <- meta_model
    
    # Calculate the MSE on the hold-out set
    mse_values[i] <- mse_values[i] + mean((final_preds - test$Species)^2) / folds
  }
}

# Calculate the average MSE over the repetitions
mean(mse_values)



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
