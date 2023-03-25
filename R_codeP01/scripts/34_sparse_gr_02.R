############################################################################
## Repeated cross-validation for testing Sparse-group Lasso with SPARSGEL ##
############################################################################
##
## Test against the fold
## Cross-validation (5-fold) used for training/tuning lambda and selecting features
## output: selected genes; proliferation.score correlation; SEM

char_list <- list(imm_inf_ = char_immune_inf,
                  prolif_ = char_prolif,
                  ER_sing_ = char_ER_signaling,
                  anti_pres_ = char_antigen_present,
                  angiogen_ = char_angiogenesis
                  )


prepare_group_index_vector <- function(df, dic_groups) {
  # Input
  # df = dataframe containing with predictors and response variable in first column
  # col_groups = a list of character vectors, where each vector contains the column names belonging to the same group
  # Output: group assignment vector with 0 where variables where not assign to a specific group
  
  # Initialize the group assignment vector
  groups <- rep(0, ncol(df)-1)
  
  # Remove the response variable from the dataframe
  predictor_df <- df[, -1]
  
  # Assign group numbers to predictors within choosen groups as defined by col_groups
  group_number <- 1
  for (group_cols in dic_groups) {
    for (col in group_cols) {
      col_index <- which(colnames(predictor_df) == col)
      groups[col_index] <- group_number
    }
    group_number <- group_number + 1
  }
  
  return(groups)
}


# Sparse Group Lasso function with grouped predictors a arguments
sparse_group_lasso_sparsegl <- function(df, groups) {
  # Input:
  # df = dataframe with predictors and the response variable in first column
  # groups = group assignment vector with zeros for ungrouped predictors
  
  # Assign individual groups to ungrouped predictors
  ungrouped <- which(groups == 0)
  if (length(ungrouped) > 0) {
    max_group <- max(groups)
    for (i in ungrouped) {
      max_group <- max_group + 1
      groups[i] <- max_group
    }
  }
  
  # Extract predictor matrix and response vector; assuming response variable in first column
  X <- as.matrix(df[, -1])
  X <- scale(X)
  y <- df[, 1]
  # Run cross-validated sparse group lasso
  cv_model <- cv.sparsegl(X, y, groups, nfolds=5)
  
  optimal_lambda <- cv_model$lambda.min
  coef <- coef(cv_model, s = optimal_lambda)
  
  return(list(model=cv_model, coef=coef, optimal_lambda=optimal_lambda))
}
set.seed(123)
gr_indices <- prepare_group_index_vector(ROR_prolif_771genes, char_list)
result <- sparse_group_lasso_sparsegl(ROR_prolif_771genes, gr_indices)

dense_coef <- as.matrix(result$coef) # Convert the sparse matrix to a dense matrix
non_zero_coef <- dense_coef[dense_coef != 0]  # Get non-zero coefficients and their corresponding names
names(non_zero_coef) <- rownames(result$coef)[dense_coef != 0]  # Print non-zero coefficients with their corresponding names
print(non_zero_coef)



# Function: repeated k-fold cross validation: 
sparsgel_rep_cv <- function(df_data, groups, folds=5, repeats=1, method="pearson") {
  n_models <- repeats * folds
  print(n_models)
  
  cor_vec <- rep(NA, n_models)
  MSE_vec <- rep(NA, n_models)
  coef_matrix <- matrix(NA, nrow = n_models, ncol = ncol(df_data[,-1]))
  colnames(coef_matrix) <- colnames(df_data[,-1])
  row_index <- 1
  
  # Repeat the cross-validation process
  for (i in 1:repeats) {
    # Create the folds for evaluating the performance
    kf <- caret::createFolds(df_data[,1], k = folds, list = TRUE, returnTrain = TRUE)
    # Loop through the folds
    for (j in 1:folds) {
      # Get the training and testing data
      train_data <- df_data[kf[[j]],]
      test_data <- df_data[-kf[[j]],]
  
      # # Fit the function on the training data and get results
      out <- sparse_group_lasso_sparsegl(train_data, groups)
      
      model <- out$model
      coef_matrix[row_index, ] <- out$coef[-1]   # remove intercept
      
      pred = predict(model, newx = as.matrix(test_data)[,-1], type = "response", s = "lambda.min")  
      cor_vec[row_index]  <- suppressWarnings(cor(pred, test_data[,1], method=method))
      MSE_vec[row_index] <- mean((pred - test_data[,1])^2)
      
      if(!row_index %% 10) cat(row_index, "")
      row_index <- row_index + 1
    }
  }
  return(list(cor_vec=cor_vec, MSE_vec=MSE_vec, coef_matrix=coef_matrix))
}


gr_indices <- prepare_group_index_vector(ROR_prolif_771genes, char_list)
result <- sparse_group_lasso_sparsegl(ROR_prolif_771genes, gr_indices)

dense_coef <- as.matrix(result$coef) # Convert the sparse matrix to a dense matrix
non_zero_coef <- dense_coef[dense_coef != 0]  # Get non-zero coefficients and their corresponding names
names(non_zero_coef) <- rownames(result$coef)[dense_coef != 0]  # Print non-zero coefficients with their corresponding names
print(non_zero_coef)


t1 <- sparsgel_rep_cv(ROR_prolif_771genes, groups=gr_indices, folds=5, repeats=20, method="pearson")
mean(t1$cor_vec, na.rm=TRUE)
mean(t1$MSE_vec, na.rm=TRUE)
t1$coef_matrix
dim(t1$coef_matrix)
sd(t1$cor_vec, na.rm=TRUE)
sd(t1$MSE_vec, na.rm=TRUE)


# RUN: sgl_c_771_RORprolif
set.seed(123)
sgl_c_771_RORprolif <- sparsgel_rep_cv(ROR_prolif_771genes, groups=gr_indices, folds=5, repeats=200, method="pearson")
head(sgl_c_771_RORprolif$coef_matrix)[,1:8]
save(sgl_c_771_RORprolif, file="/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/sgl_c_771_RORprolif.RData")
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/sgl_c_771_RORprolif.RData")
mean(sgl_c_771_RORprolif$cor_vec, na.rm=T)
sd(sgl_c_771_RORprolif$cor_vec)
