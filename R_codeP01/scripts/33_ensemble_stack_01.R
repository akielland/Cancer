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




t <- lasso_rep_cv(prolif_771genes, folds=2, repeats=2, method="pearson")

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
    kf_1 <- caret::createFolds(df_data[,1], k = folds, list = TRUE, returnTrain = TRUE)
    # Loop through the folds
    for (j in 1:folds) {
      # Get the training and testing data
      train_data <- df_data[kf_1[[j]],]
      test_data <- df_data[-kf_1[[j]],]

      # Split into training and prediction set
      ind <- sample.int(2, nrow(train_data), replace=TRUE, prob=c(0.8, 0.2))
      train0 <- train_data[ind==1,]
      pred0 <- train_data[ind==2,]
      
      # Partition feature space for lasso
      X_prolif <- as.matrix(train0[c(char_prolif)])
      X_immune_inf <- as.matrix(train0[c(char_immune_inf)])
      X_ER_singaling <- as.matrix(train0[c(char_ER_signaling)])
    
      fit_prolif <- cv.glmnet(X_prolif, train0$Y, nfolds=5)
      fit_immune_inf <- cv.glmnet(X_immune_inf, train0$Y, nfolds=5)
      fit_ER_singaling <- cv.glmnet(X_ER_singaling, train0$Y, nfolds=5)
      
      pred_prolif = predict(fit_prolif, newx = as.matrix(pred0[c(char_prolif)]), type = "response", s = "lambda.min")
      pred_immune = predict(fit_immune_inf, newx = as.matrix(pred0[c(char_immune_inf)]), type = "response", s = "lambda.min")
      pred_ER_signal = predict(fit_ER_singaling, newx = as.matrix(pred0[c(char_ER_signaling)]), type = "response", s = "lambda.min")
     
       # stack predictions
      preds <- cbind(pred_prolif, pred_immune, pred_ER_signal)
      colnames(preds) <- c("pred_prolif", "pred_immune", "pred_ER_signal")
      Y_preds <- cbind("Y"=pred0$Y, preds)
      meta_model <- lm(Y ~ pred_prolif+pred_immune+pred_ER_signal, data = as.data.frame(Y_preds))
      
      pred <- predict(meta_model, newx=test_data[-1])

   # coef_matrix[coef_matrix_row_index, ] <- coef(fit, s = "lambda.min")[-1]

      cor_vec[count]  <- suppressWarnings(cor(pred, test_data[,1], method=method))
    # MSE_vec[coef_matrix_row_index] <- mean((pred - test_data[,1])^2)
    # 
      cat(count, "")
      count <- count + 1
    }
  }
  return(cor_vec)
  #return(list(cor_vec=cor_vec, coef_matrix=coef_matrix, MSE_vec=MSE_vec))
}


library(caret)
library(glmnet)

# Load the data
data(iris)

reps <- 10
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
