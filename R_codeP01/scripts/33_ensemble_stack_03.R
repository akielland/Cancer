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


# Name list og gene signature sets 
char_list <- list(imm_inf_ = char_immune_inf,
                  prolif_ = char_prolif,
                  ER_sing_ = char_ER_signaling,
                  anti_pres_ = char_antigen_present,
                  angiogen_ = char_angiogenesis
)

# function to calculate new feature set from level 0 learners
level_0_model <- function(named_list, train_data){
  pred_y_L0 <- data.frame(Y=train_data$Y)
  fits <- list()
  
  for (key in names(named_list)) {
    char_matrix <- named_list[[key]]
    X_ <- as.matrix(train_data[c(char_matrix)])
    fit_ <- cv.glmnet(X_, train_data$Y, nfolds=5)
    pred_ = predict(fit_, newx = X_, type = "response", s = "lambda.min")
    colnames(pred_) <- key
    pred_y_L0 <- cbind(pred_y_L0, pred_)
    
    fits[[key]] <- fit_

  }
  return(list(pred_y_L0=pred_y_L0, fits=fits))
}
t <-level_0(char_list, prolif_771genes)
class(t)

# function to calculate predictions of the level 0 models
level_0_test <- function(named_list, test_data, fits){
  pred_test_L0 <- data.frame(Y=test_data$Y)
  
  for (key in names(named_list)) {
    char_matrix <- named_list[[key]]
    X_ <- as.matrix(test_data[c(char_matrix)])
    fit_ <- fits[[key]]
    test_pred = predict(fit_, newx = X_, type = "response", s = "lambda.min")
    colnames(test_pred) <- key
    pred_test_L0[key] <- test_pred
  }
  return(pred_test_L0)
}

t1 <- level_0_test(char_list, prolif_771genes, t$fit)


# Function: repeated k-fold cross validation
lasso_rep_cv <- function(df_data, named_list, folds=5, repeats=1, method="pearson") {
  n_models <- repeats * folds
  print(n_models)
  cor_vec <- rep(NA, n_models)
  MSE_vec <- rep(NA, n_models)
  coef_matrix <- matrix(NA, nrow = n_models, ncol = length(char_list))
  colnames(coef_matrix) <- names(named_list)
  print(coef_matrix)
  
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

      # # Split into training and prediction set
      # ind <- sample.int(2, nrow(train_data), replace=TRUE, prob=c(0.8, 0.2))
      # train0 <- train_data[ind==1,]
      # pred0 <- train_data[ind==2,]
      
      # Partition feature space; run lasso and stack in-sample predictions
      predictions <- level_0_model(char_list, train_data)
      pred_y_L0 <- predictions$pred_y_L0
      fits <- predictions$fits
      # fit lasso meta-model
     # meta_model <- cv.glmnet(as.matrix(pred_y_L0[,-1]), train_data$Y, nfolds=5)
      meta_model <- cv.glmnet(as.matrix(pred_y_L0[,-1]), train_data$Y, nfolds=5, alpha=0.5)
   
      ####################################################################
      ## test part
      ####################################################################
      # test data prediction from level 0 
      pred_y_test_L0 <- level_0_test(named_list, test_data, fits)
      # predict meta_model with test-data
      test_pred <- predict(meta_model, newx = as.matrix(pred_y_test_L0[,-1]), type = "response", s = "lambda.min")
     
      coef_matrix[count, ] <- coef(meta_model, s = "lambda.min")[-1]
      cor_vec[count]  <- suppressWarnings(cor(test_pred, test_data[,1], method=method))
      MSE_vec[count] <- mean((test_pred - test_data[,1])^2)
    
      cat(count, "")
      count <- count + 1
    }
  }
  #return(cor_vec)
  return(list(cor_vec=cor_vec, coef_matrix=coef_matrix, MSE_vec=MSE_vec))
}

t3 <- lasso_rep_cv(prolif_771genes, char_list, folds=5, repeats=1, method="pearson")

reps <- 2
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
