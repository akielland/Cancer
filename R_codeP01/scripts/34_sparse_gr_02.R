###################################################################
## Repeated cross-validation for testing Sparse-group Lasso with SPARSGEL ##
###################################################################
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



prepare_group_index_vector <- function(df, col_groups) {
  # Input
  # df = dataframe containing with predictors and response variable in first column
  # col_groups = a list of character vectors, where each vector contains the column names belonging to the same group
  # Output: group assignment vector with 0 where variables where not assign to a specific group
  
  # Initialize the group assignment vector
  group <- rep(0, ncol(df)-1)
  
  # Remove the response variable from the dataframe
  predictor_df <- df[, -1]
  
  # Assign group numbers to predictors within choosen groups as defined by col_groups
  group_number <- 1
  for (group_cols in col_groups) {
    for (col in group_cols) {
      col_index <- which(colnames(predictor_df) == col)
      group[col_index] <- group_number
    }
    group_number <- group_number + 1
  }
  
  return(group)
}


# Sparse Group Lasso function with grouped and ungrouped predictors
sparse_group_lasso_sparsegl <- function(df, group) {
  # Input:
  # df = dataframe containing with predictors and response variable in first column
  # groups = group assignment vector with zeros for ungrouped predictors
  
  # Assign individual groups to ungrouped predictors
  ungrouped <- which(group == 0)
  if (length(ungrouped) > 0) {
    max_group <- max(group)
    for (i in ungrouped) {
      max_group <- max_group + 1
      group[i] <- max_group
    }
  }
  print(length(group))
  # print(group)
  
  # Extract predictor matrix and response vector; assuming response variable in first column
  X <- as.matrix(df[, -1])
  X <- scale(X)
  y <- df[, 1]
  # Run cross-validated sparse group lasso
  cv_model <- cv.sparsegl(X, y, group, nfolds=5)
  
  optimal_lambda <- cv_model$lambda.min
  coef <- coef(cv_model, s = optimal_lambda)
  
  return(list(coef=coef, optimal_lambda=optimal_lambda))
}

gr_indices <- prepare_group_index_vector(ROR_prolif_771genes, char_list)

# Run sparse group lasso
result <- sparse_group_lasso_sparsegl(ROR_prolif_771genes, gr_indices)


dense_coef <- as.matrix(result$coef) # Convert the sparse matrix to a dense matrix
non_zero_coef <- dense_coef[dense_coef != 0]  # Get non-zero coefficients and their corresponding names
names(non_zero_coef) <- rownames(result$coef)[dense_coef != 0]  # Print non-zero coefficients with their corresponding names
print(non_zero_coef)



# Function: repeated k-fold cross validation: 
sparsgel_rep_cv <- function(df_data, alpha, folds=5, repeats=1, interactions=FALSE, method="pearson") {
  n_models <- repeats * folds
  print(n_models)
  print(alpha)

  # make df with response Y and the PCA's from the different signatures
  if (interactions){
    df_Y_PCA <-pca_as_features_interaction(char_list, df_data)
  } else {
    df_Y_PCA <- pca_as_features(char_list, df_data, percentage)
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
t1 <- PCA_rep_cv(prolif_771genes, alpha=0.5, folds=5, repeats=5, interactions=TRUE, method="pearson")
mean(t1$cor, na.rm=TRUE)
t1$coef_matrix

# Ridge
# RUN: r_pca_c_771_RORprolif
set.seed(123)
percentage = 0.9
r_pca_c_771_RORprolif <- PCA_rep_cv(ROR_prolif_771genes, alpha=0, folds=5, repeats=200, interactions=FALSE, method="pearson")
head(r_pca_c_771_RORprolif$coef_matrix)[,1:8]
save(r_pca_c_771_RORprolif, file="/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/r_pca_c_771_RORprolif.RData")
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/r_pca_c_771_RORprolif.RData")
mean(r_pca_c_771_RORprolif$cor_vec, na.rm=T)
sd(r_pca_c_771_RORprolif$cor_vec)

# RUN: r_pca_c_interact_771_RORprolif
set.seed(123)
r_pca_c_interact_771_RORprolif <- PCA_rep_cv(ROR_prolif_771genes, alpha=0, folds=5, repeats=200, interactions=TRUE, method="pearson")
head(r_pca_c_interact_771_RORprolif$coef_matrix)[,1:8]
save(r_pca_c_interact_771_RORprolif, file="/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/r_pca_c_interact_771_RORprolif.RData")
mean(r_pca_c_interact_771_RORprolif$cor_vec, na.rm=T)

# Lasso
# RUN: l_pca_c_771_RORprolif
set.seed(123)
percentage = 0.9
l_pca_c_771_RORprolif <- PCA_rep_cv(ROR_prolif_771genes, alpha=1, folds=5, repeats=200, interactions=FALSE, method="pearson")
head(l_pca_c_771_RORprolif$coef_matrix)[,1:8]
save(l_pca_c_771_RORprolif, file="/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/l_pca_c_771_RORprolif.RData")
mean(l_pca_c_771_RORprolif$cor_vec, na.rm=T)

# RUN: l_pca_c_interact_771_RORprolif
set.seed(123)
l_pca_c_interact_771_RORprolif <- PCA_rep_cv(ROR_prolif_771genes, alpha=1, folds=5, repeats=200, interactions=TRUE, method="pearson")
head(l_pca_c_interact_771_RORprolif$coef_matrix)[,1:8]
save(l_pca_c_interact_771_RORprolif, file="/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/l_pca_c_interact_771_RORprolif.RData")
mean(l_pca_c_interact_771_RORprolif$cor_vec, na.rm=T)

# Elastic
# RUN: en_pca_c_771_RORprolif
set.seed(123)
percentage = 0.9
en_pca_c_771_RORprolif <- PCA_rep_cv(ROR_prolif_771genes, alpha=0.5, folds=5, repeats=200, interactions=FALSE, method="pearson")
head(en_pca_c_771_RORprolif$coef_matrix)[,1:8]
save(en_pca_c_771_RORprolif, file="/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/en_pca_c_771_RORprolif.RData")
mean(en_pca_c_771_RORprolif$cor_vec, na.rm=T)

# RUN: en_pca_c_interact_771_RORprolif
set.seed(123)
en_pca_c_interact_771_RORprolif <- PCA_rep_cv(ROR_prolif_771genes, alpha=0.5, folds=5, repeats=200, interactions=TRUE, method="pearson")
head(en_pca_c_interact_771_RORprolif$coef_matrix)[,1:8]
save(en_pca_c_interact_771_RORprolif, file="/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/en_pca_c_interact_771_RORprolif.RData")
mean(en_pca_c_interact_771_RORprolif$cor_vec, na.rm=T)

