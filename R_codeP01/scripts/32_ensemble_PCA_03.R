###################################################################
## Repeated cross-validation for testing the PCA-Ensemble models ##
###################################################################
##
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

# create df with interaction terms between PCs
pca_as_features_interaction <- function(named_list, df_data){
  df <- data.frame(matrix(nrow = nrow(df_data), ncol = 0))
  
  for (key in names(named_list)) {
    char_matrix <- named_list[[key]]
    df_ <- PC(key, char_matrix, data=df_data)
    df <- cbind(df, df_)
  }
  
  # add interaction terms to dataframe with features
  interaction_matrix <- model.matrix(~ .^2 - 1, data = df)
  colnames(interaction_matrix) <- sub(":", "*", colnames(interaction_matrix)) # replace ':' with '*' in column names
  
  df_Y <- data.frame(Y=df_data[,1])
  df_Y_PCA_interaction <- cbind(df_Y, interaction_matrix)
  return(df_Y_PCA_interaction)
}

t <-pca_as_features_interaction(char_list, prolif_771genes)
t <-pca_as_features(char_list, prolif_771genes)


# Lasso function returning a trained object
lasso_sample <- function(train_data){
  X_ <- as.matrix(train_data |> select(-1))
  Y_ <- as.matrix(train_data |>  select(1))
  # Fit the model with lasso regularization
  fit.lasso <- glmnet(X_, Y_, family = "gaussian", alpha = 1, standardize = TRUE, nlambda = 100)
  # Perform cross-validation to select the best model
  fit.cv <- cv.glmnet(X_, Y_, lambda = fit.lasso$lambda, alpha = 1, nfolds = 5)
  return(fit.cv)
}


# Elastic Net function returning a trained object
elastic_net_sample <- function(train_data){
  X_ <- as.matrix(train_data |> select(-1))
  Y_ <- as.matrix(train_data |>  select(1))
  # Fit the model with elastic net regularization
  fit.en <- glmnet(X_, Y_, family = "gaussian", alpha = 0.5, standardize = TRUE, nlambda = 100, nfolds = 5)
  # Perform cross-validation to select the best model
  fit.cv <- cv.glmnet(X_, Y_, lambda = fit.en$lambda, alpha = 0.5, nfolds = 5)
  return(fit.cv)
}


# Function: repeated k-fold cross validation: NEED TO select if interaction terms or not
lasso_PCA_rep_cv <- function(df_data, folds=5, repeats=1, interactions=FALSE, method="pearson", model) {
  n_models <- repeats * folds
  print(n_models)

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
  
      # Fit the function on the training data and get results
      if (model == "lasso") {
        fit <- lasso_sample(train_data)
      } else if (model == "elastic_net") {
        fit <- elastic_net_sample(train_data)
      }
    
      coef_matrix[row_index, ] <- coef(fit, s = "lambda.min")[-1]
      
      pred = predict(fit, newx = as.matrix(test_data)[,-1], type = "response", s = "lambda.min")  
      cor_vec[row_index]  <- suppressWarnings(cor(pred, test_data[,1], method=method))
      MSE_vec[row_index] <- mean((pred - test_data[,1])^2)
      
      cat(row_index, "")
      row_index <- row_index + 1
    }
  }
  return(list(cor_vec=cor_vec, coef_matrix=coef_matrix, MSE_vec=MSE_vec))
}


t1 <- lasso_PCA_rep_cv(prolif_771genes, folds=5, repeats=50, interactions=TRUE, method="spearman", model="lasso")
mean(t1$cor, na.rm=TRUE)
t1$coef_matrix

percentage = 0.9
# RUN: pc_c_obj_771_prolif
set.seed(123)
pc_c_obj_771_prolif <- lasso_PCA_rep_cv(prolif_771genes, folds=5, repeats=200, method="pearson")
head(pc_c_obj_771_prolif$coef_matrix)[,1:8]
save(pc_c_obj_771_prolif, file="/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/pc_c_obj_771_prolif.RData")
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/pc_c_obj_771_prolif.RData")
mean(pc_c_obj_771_prolif$cor_vec, na.rm=T)
sd(pc_c_obj_771_prolif$cor_vec)

# RUN: pc_c_obj_771_RORprolif
set.seed(123)
pc_c_obj_771_RORprolif <- lasso_PCA_rep_cv(ROR_prolif_771genes, folds, repeats, method="pearson")
head(pc_c_obj_771_RORprolif$coef_matrix)[,1:8]
save(pc_c_obj_771_RORprolif, file="/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/pc_c_obj_771_RORprolif.RData")
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/pc_c_obj_771_RORprolif.RData")
mean(pc_c_obj_771_RORprolif$cor_vec, na.rm=T)
sd(pc_c_obj_771_RORprolif$cor_vec)

# RUN: l_pc_c_obj_interaction_771_RORprolif
set.seed(123)
folds=5
repeats=200
model="lasso"
l_pc_c_obj_interaction_771_RORprolif <- lasso_PCA_rep_cv(prolif_771genes, folds, repeats, interactions=TRUE, method="pearson", model)
head(l_pc_c_obj_interaction_771_RORprolif$coef_matrix)[,1:8]
save(l_pc_c_obj_interaction_771_RORprolif, file="/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/l_pc_c_obj_interaction_771_RORprolif.RData")
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/l_pc_c_obj_interaction_771_RORprolif.RData")
mean(l_pc_c_obj_interaction_771_RORprolif$cor_vec, na.rm=T)
sd(l_pc_c_obj_interaction_771_RORprolif$cor_vec)

# RUN: en_pc_c_obj_interaction_771_RORprolif
set.seed(123)
folds=5
repeats=200
model="elstic_net"
en_pc_c_obj_interaction_771_RORprolif <- lasso_PCA_rep_cv(prolif_771genes, folds, repeats, interactions=TRUE, method="pearson", model)
head(en_pc_c_obj_interaction_771_RORprolif$coef_matrix)[,1:8]
save(en_pc_c_obj_interaction_771_RORprolif, file="/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/en_pc_c_obj_interaction_771_RORprolif.RData")
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/en_pc_c_obj_interaction_771_RORprolif.RData")
mean(en_pc_c_obj_interaction_771_RORprolif$cor_vec, na.rm=T)
sd(en_pc_c_obj_interaction_771_RORprolif$cor_vec)
