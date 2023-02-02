#################################################################
##  XGBoost on Repeated cross-validation ##
#################################################################

## Boosting using 5-fold cross validation to evaluate model
## Base learners:
## - btree
## - (linear)
## - (spline)
## 
## Test against folds
## output: - vector with correlations 
##         - matrix with features
##         - vector with SEM
##         - (vector with number of rounds used in each model)

library(xgboost)
library(Matrix)
library(data.table)
library(caret)
library(ggplot2)

# structure of base learners
fm01 <- Y ~ CCND1 + CCNE1 + CDKN1A + ESR1 + MYC + RB1
fm05 <- as.formula(paste("Y", paste(genes, collapse="+"), sep=" ~ "))

XGBoost_sample <- function(fm, df_train) {
  # Convert dataframe to dataMatrix to DMatrix object
  train_sparse = sparse.model.matrix(object = fm, data = df_train)
  dtrain = xgb.DMatrix(data = train_sparse, label = df_train$Y)

  # param <- list(eta = 0.1, max_depth = 5)
   param <- list(booster = "gbtree", 
                objective = "reg:squarederror", 
                eval_metric = "rmse", 
                # nthread = 4, 
                max_depth = 1, 
                min_child_weight = 1,
                tree_method = "hist",   # tree_methods=  exact, approx and hist, depending on the dataset size, memory and time constraints
                #eta_decay=0.1,
                eta=0.3)
 
  cv <- xgb.cv(data = dtrain, 
               params = param, 
               nrounds = 100, 
               nfold = 5, 
               #showsd = T, stratified = T,
               print_every_n = 10,
               early_stopping_rounds = 5,
               verbose = FALSE)
  
  n_rounds <- cv$best_iteration
  
  xgb_model = xgb.train(data = dtrain,
                        # The following is needed if to monitor training and validation error
                        # watchlist = list(train = dtrain, val = dvalidation),
                        params = param,
                        nrounds = n_rounds,
                        verbose = FALSE)
  
  return(xgb_model)
}

XGBoost_sample(fm01, prolif_6genes)

# Function: repeated k-fold cross validation
residuals_xgb_rep_cv = function(fm, df_data, pred_mech, additive=TRUE, folds=5, repeats=200, method="pearson"){
  n_models <- repeats * folds
  print(n_models)
  
  cor_vec <- rep(NA, n_models)
  MSE_vec <- rep(NA, n_models)
  n_rounds_vec <- integer(length = n_models)
  
  coef_matrix <- matrix(NA, nrow = n_models, ncol = ncol(df_data)-1)
  colnames(coef_matrix) <- colnames(df_data[, -1])
  
  coef_matrix_row_index <- 1
  
  # calculating the residuals
  if (additive==TRUE){res <- df_data$Y - pred_mech}
  else{res <- df_data$Y / pred_mech}
  
  # create df with training data and update the training data with residuals instead of Y
  df_train_res <- df_data
  df_train_res$Y <- res
  
  # Repeat the cross-validation process
  for (i in 1:repeats) {
    # Create the folds for evaluating the performance
    kf <- caret::createFolds(df_data[,1], k = folds, list = TRUE, returnTrain = TRUE)
    # Loop through the folds
    for (j in 1:folds) {
      # Get the training and testing data: here picked from two different df
      train_data_res <- df_train_res[kf[[j]],]
      test_data <- df_data[-kf[[j]],]
      test_pred_mech <- pred_mech[-kf[[j]]]
      
      # Fit the function on the training data and get results
      xgb_model <- XGBoost_sample(fm, train_data_res)
    
      # Extract feature importance
      feature_importance <- data.table(xgb.importance(colnames(df_data), model = xgb_model))
      # print(feature_importance[,1:2])
      # Sum up the feature importance
      coef_matrix[coef_matrix_row_index, feature_importance$Feature] <- 1
      # print(coef_matrix)

      # Test set converted to DMatrix object
      test_sparse = sparse.model.matrix(object = fm, data = test_data)
      d_test = xgb.DMatrix(data = test_sparse, label = test_data$Y)
      
      pred_res = predict(object = xgb_model, newdata = d_test)
      
      if (additive==TRUE){pred <- test_pred_mech + pred_res}
      else{pred <- test_pred_mech * pred_res}
      
      cor_vec[coef_matrix_row_index] <- suppressWarnings(cor(pred, test_data$Y, method = method))
      MSE_vec[coef_matrix_row_index] <- mean((test_data$Y - pred)^2)
      
      cat(coef_matrix_row_index, "")
      coef_matrix_row_index <- coef_matrix_row_index + 1
    }
  }  
  return(list(cor_vec=cor_vec, MSE_vec=MSE_vec, coef_matrix=coef_matrix, n_rounds_vec=n_rounds_vec))
}

bb_object_t <- residuals_xgb_rep_cv(fm01, prolif_6genes, pred_mech, additive=TRUE, folds=5, repeats=2, method="pearson")
bb_object_t


# Set repeats and folds of the cross-validations
repeats = 200
folds = 5

# RUN: xc_obj_6_prolif
set.seed(123)
xc_obj_6res_prolif <- residuals_xgb_rep_cv(fm01, prolif_6genes, pred_mech, additive=TRUE, folds, repeats, method="pearson")
head(xc_obj_6res_prolif$coef_matrix)[,1:6]
save(xc_obj_6res_prolif, file="/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/xc_obj_6res_prolif.RData")
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instancesxc_obj_6res_prolif.RData")
mean(xc_obj_6res_prolif$cor_vec)
sd(xc_obj_6res_prolif$cor_vec)

# RUN: xc_obj_6_RORprolif
set.seed(123)
xc_obj_6_RORprolif <- lasso_rep_cv(RORprolif_6genes, func=lasso_sample, folds, repeats, method="pearson")
head(xc_obj_6_RORprolif$coef_matrix)[,1:6]
save(xc_obj_6_RORprolif, file="/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/xc_obj_6_RORprolif.RData")
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/xc_obj_6_RORprolif.RData")
mean(xc_obj_6_RORprolif$cor_vec)
sd(xc_obj_6_RORprolif$cor_vec)


# RUN: xc_obj_771_prolif
set.seed(123)
xc_obj_771res_prolif <- residuals_xgb_rep_cv(fm05, prolif_771genes, pred_mech, additive=TRUE, folds, repeats, method="pearson")
head(xc_obj_771res_prolif$coef_matrix)[,1:8]
save(xc_obj_771res_prolif, file="/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/xc_obj_771res_prolif.RData")
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/xc_obj_771res_prolif.RData")
mean(xc_obj_771res_prolif$cor_vec)
sd(xc_obj_771res_prolif$cor_vec)

# RUN: xc_obj_771_RORprolif
set.seed(123)
xc_obj_771_RORprolif <- lasso_rep_cv(RORprolif_771genes, func=lasso_sample, folds, repeats, method="pearson")
head(xc_obj_771_RORprolif$coef_matrix)[,1:8]
save(xc_obj_771_RORprolif, file="/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/xc_obj_771_RORprolif.RData")
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/xc_obj_771_RORprolif.RData")
mean(xc_obj_771_RORprolif$cor_vec)
sd(xc_obj_771_RORprolif$cor_vec)


# RUN: xc_obj_nodes_prolif
set.seed(123)
xc_obj_nodes_res_prolif <- residuals_xgb_rep_cv(fm05, prolif_nodes, pred_mech, additive=TRUE, folds, repeats, method="pearson")
head(xc_obj_nodes_prolif$coef_matrix)[,1:8]
save(xc_obj_nodes_prolif, file="/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/xc_obj_nodes_prolif.RData")
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/xc_obj_nodes_prolif.RData")
mean(xc_obj_nodes_prolif$cor_vec)
sd(xc_obj_nodes_prolif$cor_vec)

# RUN: xc_obj_nodes_RORprolif
set.seed(123)
xc_obj_nodes_RORprolif <- lasso_rep_cv(RORprolif_nodes, func=lasso_sample, folds, repeats, method="pearson")
head(xc_obj_nodes_RORprolif$coef_matrix)[,1:8]
save(xc_obj_nodes_RORprolif, file="/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/xc_obj_nodes_RORprolif.RData")
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/xc_obj_nodes_RORprolif.RData")
mean(xc_obj_nodes_RORprolif$cor_vec)
sd(xc_obj_nodes_RORprolif$cor_vec)

###################




histogram(bb_object$cor_vec, breaks = 99,
          xlab = "Proliferation score", 
          main = "Boosting (trail 1, arm Letro+Ribo)")

# Various objects
save(bb_object, file="bb_object_6genes.RData")
save(bb_object, file="bb_object_nodes01.RData")
save(bb_object, file="bb_object_AllGenes01.RData")
save(bb_object, file="bb_object_ALLGenes_ROR.RData")

# stored object can be loaded to save time
load("bb_object_6Genes.RData")
load("bb_object_nodes01.RData")
load("bb_object_AllGenes01.RData")
load("bb_object_ALLGenes_ROR.RData")




# Convert the dataframe to a DMatrix object - This dosen't seems to work
dtrain <- xgb.DMatrix(data = as.matrix(Proliferation_6genes[, -1]), label = Proliferation_6genes$Y)
dtrain <- xgb.DMatrix(data = data.matrix(Proliferation_6genes[,-1]), label = Proliferation_6genes$Y)

model_formula = fm01
# Train set: 
train_sparse = sparse.model.matrix(object = model_formula, data = Proliferation_6genes)
dtrain = xgb.DMatrix(data = train_sparse, label = Proliferation_6genes$Y)
class(train_sparse)
# Validation set
validation_sparse = sparse.model.matrix(object = model_formula, data = Proliferation_6genes)
dvalidation = xgb.DMatrix(data = validation_sparse, label = Proliferation_6genes$Y)
# Test set
test_sparse = sparse.model.matrix(object = model_formula, data = test_df)
dtest = xgb.DMatrix(data = test_sparse, label = test_df$y)

# Hyper parameters:
# eta: the learning rate. Comparable to step size in gradient descent algorithms
# max_depth: the max depth of each tree. To get stumps: set tree_depth = 1
xgb_params = list(eta = 0.1, max_depth = 1)
params <- list(objective = "binary:logistic",
               eta = 0.3,
               max_depth = 2)
# THIS IS WHAT YOU NEED TO UNDERSTAND:
# If we take small steps (low eta), we require many trees (nrounds) to get a good result.
# If we make deep trees then each tree is better => probably need fewer trees to get a good result

# How many sequential trees should we train? Number of iterations
N_ROUNDS = 1000

start_time = Sys.time()
xgb_model = xgb.train(data = dtrain,
                      # The following is needed if to monitor training and validation error
                      watchlist = list(train = dtrain, 
                                       val = dvalidation),
                      params = xgb_params,
                      nrounds = N_ROUNDS,
                      verbose = FALSE)
end_time = Sys.time()
cat("It took ", difftime(end_time, start_time, units = "secs"), " to train the model. \n")

# Extract and plot train and validation error: 
training_error = data.frame(iteration = xgb_model$evaluation_log$iter,
                            rmse = xgb_model$evaluation_log$train_rmse, 
                            type = "train_rmse")
validation_error = data.frame(iteration = xgb_model$evaluation_log$iter,
                              rmse = xgb_model$evaluation_log$val_rmse, 
                              type = "val_rmse")

both_errors = rbind(training_error, validation_error)

ggplot(both_errors, aes(x = iteration, y = rmse, col = type)) + 
  geom_line()


boost_m01 = glmboost(fm01, data = Proliferation_6genes)
coef(boost_m01, which = "")
par(mfrow=c(1,2))
plot(boost_m01, off2int=TRUE)
plot(boost_m01)

pred_boost01 = predict(boost_m01, Proliferation_6genes)

cor(pred_boost01, Proliferation_6genes$Y)
plot(pred_boost01, Proliferation_6genes$Y)
abline(lm(Proliferation_6genes$Y ~ pred_boost01))


# smooth - P-spline as base learner
spline01 <- gamboost(fm01, data = Proliferation_6genes,
                     baselearner = "bbs", # dfbase=4
                     control = boost_control(mstop = 40))

coef(spline01, which = "")
par(mfrow=c(2,3))
plot(spline01, off2int=TRUE)

cols_chosen = unique(spline01$xselect())
colnames(Proliferation_6genes)[cols_chosen]


pred_spline01 = predict(spline01, Proliferation_6genes)

cor(pred_spline01, Proliferation_6genes$Y)
par(mfrow=c(1,1))
plot(pred_spline01, Proliferation_6genes$Y)
abline(lm(Proliferation_6genes$Y ~ pred_spline01))

# boosting with stumps
stump01 = gamboost(fm01, data = Proliferation_6genes, 
                   family=Gaussian(), 
                   baselearner='btree',
                   boost_control(mstop = 100))

cols_chosen = unique(stump01$xselect())
colnames(Proliferation_6genes)[cols_chosen+2]

par(mfrow=c(2,3))
plot(stump01)


pred_boost01 = predict(boost_m01, Proliferation_6genes)

cor(pred_boost01, Proliferation_6genes$Y)
plot(pred_boost01, Proliferation_6genes$Y)
abline(lm(Proliferation_6genes$Y ~ pred_boost01))
