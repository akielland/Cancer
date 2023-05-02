







# This function return the PCs based on percentage of variance explained from a given signature gene set
PCs <- function(name, characters, data, percentage = 0.9) {
  sign_data <- data[, characters]
  pca_results <- prcomp(sign_data, scale = TRUE, center = TRUE)  # create PCA object
  # Get the number of components that explain at least 90 % of the variance
  n_components <- sum(pca_results$sdev^2 / sum(pca_results$sdev^2) >= (1 - percentage))
  # Extract the first n_components principal components
  df_important_pca <- data.frame(pca_results$x[, 1:n_components])
  colnames(df_important_pca) <- c(paste0(name, 1:n_components))
  return(df_important_pca)
}

# Create df from the PCs and interaction terms
pca_as_features_interaction_combined <- function(named_list, df_data, percentage = 0.9) {
  df <- data.frame(matrix(nrow = nrow(df_data), ncol = 0))
  
  for (key in names(named_list)) {
    char_matrix <- named_list[[key]]
    df_ <- PCs(key, char_matrix, data = df_data, percentage)
    df <- cbind(df, df_)
  }
  
  # add interaction terms to dataframe with features
  interaction_matrix <- model.matrix(~ .^2 - 1, data = df)
  colnames(interaction_matrix) <- sub(":", "*", colnames(interaction_matrix))  # replace ':' with '*' in column names
  
  df_Y <- data.frame(Y = df_data[, 1])
  df_Y_PCA_interaction_combined <- cbind(df_Y, interaction_matrix)
  return(df_Y_PCA_interaction_combined)
}



t2 <- pca_as_features_interaction_combined(char_list_5, ROR_prolif_771genes)
colnames(t2)



# Function: repeated k-fold cross validation: NEED TO select if interaction terms or not
PCA_rep_cv_02 <- function(df_data, char_list, alpha, folds=5, repeats=1, interactions=FALSE, method="pearson") {
  n_models <- repeats * folds
  print(n_models)
  print(alpha)
  
  # make df with response Y and the PCA's from the different signatures
  if (interactions){
    df_Y_PCA <-pca_as_features_interaction_combined(char_list, df_data)
    print(colnames(df_Y_PCA))
  } else {
    df_Y_PCA <- pca_as_features(char_list, df_data, percentage)
    print(colnames(df_Y_PCA))
  }
  
  cor_vec <- rep(NA, n_models)
  MSE_vec <- rep(NA, n_models)
  coef_matrix <- matrix(NA, nrow = n_models, ncol = ncol(df_Y_PCA[,-1]))
  colnames(coef_matrix) <- colnames(df_Y_PCA[,-1])
  row_index <- 1
  
  # Repeat the cross-validation process
  for (i in 1:repeats) {
    # Create the folds for evaluating the performance
    kf <- caret::createFolds(df_data[,1], k = folds, list = TRUE, returnTrain = TRUE)
    # Loop through the folds
    for (j in 1:folds) {
      # Get the training and testing data
      train_data <- df_Y_PCA[kf[[j]],]
      test_data <- df_Y_PCA[-kf[[j]],]
      
      # # Fit the function on the training data and get results
      fit <- elastic_net_sample(train_data, alpha)
      # if (model == "lasso") {
      #   fit <- lasso_sample(train_data)
      # } else if (model == "elastic_net") {
      #   fit <- elastic_net_sample(train_data)
      # }
      
      
      coef_matrix[row_index, ] <- coef(fit, s = "lambda.min")[-1]
      
      pred = predict(fit, newx = as.matrix(test_data)[,-1], type = "response", s = "lambda.min")  
      cor_vec[row_index]  <- suppressWarnings(cor(pred, test_data[,1], method=method))
      MSE_vec[row_index] <- mean((pred - test_data[,1])^2)
      
      if(!row_index %% 10) cat(row_index, "")
      row_index <- row_index + 1
    }
  }
  return(list(cor_vec=cor_vec, coef_matrix=coef_matrix, MSE_vec=MSE_vec))
}

percentage = 0.9
t1 <- PCA_rep_cv_02(prolif_771genes, char_list_5, alpha=0.5, folds=5, repeats=5, interactions=TRUE, method="pearson")
mean(t1$cor, na.rm=TRUE)
t1$coef_matrix
t1 <- PCA_rep_cv_02(ROR_prolif_771genes, char_list_10, alpha=0, folds=5, repeats=5, interactions=FALSE, method="pearson")
mean(t1$cor, na.rm=TRUE)
t1$coef_matrix
t1 <- PCA_rep_cv_02(ROR_prolif_771genes, char_list_10, alpha=0, folds=5, repeats=5, interactions=TRUE, method="pearson")
mean(t1$cor, na.rm=TRUE)
t1$coef_matrix


