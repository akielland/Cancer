###################################################################
## Repeated cross-validation for testing Sparse-group Lasso with SPARSGEL ##
###################################################################
##
## Test against the fold
## Cross-validation (5-fold) used for training/tuning lambda and selecting features
## output: selected genes; proliferation.score correlation; SEM


# PCA's:
char_list <- list(imm_inf_ = char_immune_inf,
                  prolif_ = char_prolif,
                  ER_sing_ = char_ER_signaling,
                  anti_pres_ = char_antigen_present,
                  angiogen_ = char_angiogenesis
)

# This function return only the best PC from a given signature gene set
PC <- function(name, characters, data){
  sign_data <-  data[, characters]
  pca_results <- prcomp(sign_data, scale = TRUE, center=TRUE) # create PCA object
  # Extract the first principal component
  df_important_pca <- data.frame(pca_results$x[, 1])
  colnames(df_important_pca) <- paste0(name)
  return(df_important_pca)
}

# This function return the PCs based on percentage of variance explained from a given signature gene set
PCs <- function(name, characters, data, percentage=0.9){
  sign_data <-  data[, characters]
  pca_results <- prcomp(sign_data, scale = TRUE, center=TRUE) # create PCA object
  # Get the number of components that explain at least 90 % of the variance
  n_components <- sum(pca_results$sdev^2 / sum(pca_results$sdev^2) >= (1-percentage))
  # Extract the first n_components principal components
  df_important_pca <- data.frame(pca_results$x[, 1:n_components])
  colnames(df_important_pca) <- c(paste0(name, 1:n_components))
  return(df_important_pca)
}

# Create df from the PCs
pca_as_features <- function(named_list, df_data, percentage=0.9){
  df_Y_PCA <- data.frame(Y=df_data[,1])
  df <- data.frame(matrix(nrow = nrow(df_data), ncol = 0))
  for (key in names(named_list)) {
    char_matrix <- named_list[[key]]
    
    df_ <- PCs(key, char_matrix, data=df_data, percentage)
    df <- cbind(df, df_)
  }
  df_Y_PCA <- cbind(df_Y_PCA, df)
  return(df_Y_PCA)
}

  




# Elastic Net function returning a trained object
elastic_net_sample <- function(train_data, alpha){
  X_ <- as.matrix(train_data |> select(-1))
  Y_ <- as.matrix(train_data |>  select(1))
  # Fit the model with elastic net regularization
  fit.en <- glmnet(X_, Y_, family = "gaussian", alpha=alpha, standardize = TRUE, nlambda = 100, nfolds = 5)
  # Perform cross-validation to select the best model
  fit.cv <- cv.glmnet(X_, Y_, lambda = fit.en$lambda, alpha=alpha, nfolds = 5)
  return(fit.cv)
}



# Install and load the sparsgel package
if (!requireNamespace("sparsgel", quietly = TRUE)) {
  install.packages("sparsgel")
}
library(sparsgel)





# Function: repeated k-fold cross validation: NEED TO select if interaction terms or not
PCA_rep_cv <- function(df_data, alpha, folds=5, repeats=1, interactions=FALSE, method="pearson") {
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

