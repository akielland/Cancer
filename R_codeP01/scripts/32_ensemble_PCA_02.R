###################################################################
## Repeated cross-validation for testing the PCA-Ensemble models ##
###################################################################
##THERE IS NEWER VERSION OF THIS CODE (but this workes)
## Test against the fold
## Cross-validation (5-fold) used for training/tuning lambda and selecting features
## output: selected genes; proliferation.score correlation; SEM

## Post Lasso model is made at the bottom
library(stats)
library(caret)
library(glmnet)
library(ggplot2)
library(dplyr)
library(stringr)


###########################################################
##  Repeated cross validation ##
###########################################################

pca_components<- function(name, characters, data, percentage=0.9){
  sign_data <-  data[, characters]
  pca_results <- prcomp(sign_data, scale = TRUE, center=TRUE) # create PCA object
  # Get the number of components that explain at least 90 % of the variance
  n_components <- sum(pca_results$sdev^2 / sum(pca_results$sdev^2) >= (1-percentage))
  # Extract the first n_components principal components
  df_important_pca <- data.frame(pca_results$x[, 1:n_components])
  colnames(df_important_pca) <- c(paste0(name, 1:n_components))
  return(df_important_pca)
}

# Lasso base function returning a trained object
lasso_sample <- function(train_data){
  X_ <- as.matrix(train_data |> select(-1))
  Y_ <- as.matrix(train_data |>  select(1))
  # NB: NEED TO LOOK INTO the SETTINGS of Lasso
  # fit.cv <- cv.glmnet(X_, Y_, family = "gaussian", alpha = 1, standardize = TRUE, nlambda = 100, nfolds = 5)
  fit.cv <- cv.glmnet(X_, Y_, nfolds = 5)
  return(fit.cv)
}

# PCA's:

char_list <- list(imm_inf_ = char_immune_inf,
                  prolif_ = char_prolif,
                  ER_sing_ = char_ER_signaling,
                  anti_pres_ = char_antigen_present,
                  angiogen_ = char_angiogenesis
                  )


pca_as_features <- function(named_list, df_data, percentage=0.9){
  df_Y_PCA <- data.frame(Y=df_data[,1])
  for (key in names(named_list)) {
    char_matrix <- named_list[[key]]
    df_ <- pca_components(key, char_matrix, data=df_data, percentage=0.9)
    df_Y_PCA <- cbind(df_Y_PCA, df_)

  }
  return(df_Y_PCA)
}
t <-pca_as_features(char_list, prolif_771genes)



# Function: repeated k-fold cross validation
lasso_PCA_rep_cv <- function(df_data, folds=5, repeats=1, method="pearson") {
  n_models <- repeats * folds
  print(n_models)

  # PCA's:
  immune_inf <- pca_results("imm_inf_", char_immune_inf, data=df_data, percentage=0.9)
  prolif <- pca_results("prolif_", char_prolif, data=df_data, percentage=0.9)
  ER_singaling <- pca_results("ER_sing_", char_ER_signaling, data=df_data, percentage=0.9)
  antigen_present<- pca_results("anti_pres_", char_antigen_present, data=df_data, percentage=0.9)
  angiogenesis<- pca_results("angiogen_", char_angiogenesis, data=df_data, percentage=0.9)

  
  # make df with response Y and the PCA's from the different signatures
  df_Y_PCA <- data.frame(Y=df_data[,1], prolif, immune_inf, ER_singaling, antigen_present, angiogenesis)
  
  cor_vec <- rep(NA, n_models)
  MSE_vec <- rep(NA, n_models)
  coef_matrix <- matrix(NA, nrow = n_models, ncol = ncol(df_Y_PCA[,-1]))
  colnames(coef_matrix) <- colnames(df_Y_PCA[,-1])
  
  coef_matrix_row_index <- 1
  
  # Repeat the cross-validation process
  for (i in 1:repeats) {
    # Create the folds for evaluating the performance
    kf <- caret::createFolds(df_data[,1], k = folds, list = TRUE, returnTrain = TRUE)
    # Loop through the folds
    for (j in 1:folds) {
      # Get the training and testing data
      train_data <- df_Y_PCA[kf[[j]],]
      test_data <- df_Y_PCA[-kf[[j]],]
  
      # model_fit <- lm(Y ~ ., data = train_PCA)
      # Fit the function on the training data and get results
      fit <- lasso_sample(train_data)
      
      coef_matrix[coef_matrix_row_index, ] <- coef(fit, s = "lambda.min")[-1]
      
      pred = predict(fit, newx = as.matrix(test_data)[,-1], type = "response", s = "lambda.min")  
      cor_vec[coef_matrix_row_index]  <- suppressWarnings(cor(pred, test_data[,1], method=method))
      MSE_vec[coef_matrix_row_index] <- mean((pred - test_data[,1])^2)
      
      cat(coef_matrix_row_index, "")
      coef_matrix_row_index <- coef_matrix_row_index + 1
    }
  }
  return(list(cor_vec=cor_vec, coef_matrix=coef_matrix, MSE_vec=MSE_vec))
}


t1 <- lasso_PCA_rep_cv(prolif_771genes, folds=5, repeats=100, method="pearson")


