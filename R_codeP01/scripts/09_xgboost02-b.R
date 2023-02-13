#################################################################
##  XGBoost on Bootstrap samples - using original sample as test set ##
#################################################################

## Boosting using bootstrap sample to train 1000 models
## Baselearner:
## - linear
## - spline
## - btree
## 
## Test against original data sample
## output: proliferation.score correlation; SEM; Coefficient of genes

library(xgboost)
library(Matrix)
library(data.table)


# structure of base learners
fm01 <- Y ~ CCND1 + CCNE1 + CDKN1A + ESR1 + MYC + RB1
fm05 <- as.formula(paste("Y", paste(all_genes, collapse="+"), sep=" ~ "))

XGBoost_sample <- function(fm, df_train) {
  # Convert the dataMatrix to a DMatrix object
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
  #return(list(cor=cor, features_names=features_names, n_rounds=n_rounds, importance=importance, MSE=MSE))
}

t_ <- XGBoost_sample(fm01, Proliferation_6genes)

# 1000 bootstrap fits
xgboost_boot = function(fm, df_train, df_test, method="pearson", n_bootstraps=1000){
  # run many bootstraps
  # output: - vector with correlations 
  #         - matrix with features
  #         - vector with SEM
  n <- nrow(df_train)
  cor_vec <- rep(NA, n_bootstraps)
  MSE_vec <- rep(NA, n_bootstraps)
  n_rounds_vec <- integer(length = n_bootstraps)
  
  coef_matrix <- matrix(NA, nrow = n_bootstraps, ncol = ncol(df_train)-1)
  colnames(coef_matrix) <- colnames(df_train[, -1])
  
  # Test set Converted DMatrix object
  validation_sparse = sparse.model.matrix(object = fm, data = df_test)
  d_validation = xgb.DMatrix(data = validation_sparse, label = df_test$Y)
  
  for (i in c(1:n_bootstraps)) {
    int <- sample.int(n, size = n, replace = TRUE)
    train_data = df_train[int,]
    
    xgb_model <- XGBoost_sample(fm, train_data)
    
    # Extract feature importance
    feature_importance <- data.table(xgb.importance(colnames(df_train), model = xgb_model))

    # Sum up the feature importance
    coef_matrix[i, feature_importance$Feature] <- 1

    pred = predict(object = xgb_model, newdata = d_validation)
    cor_vec[i] <- suppressWarnings(cor(pred, df_test$Y, method = method))
    MSE_vec[i] <- mean((df_test$Y - pred)^2)
    # cor_vec[i]  <- suppressWarnings(cor(pred, df_test[,1], method = method))
    # MSE_vec[i] <- mean((pred - df_test[,1])^2)        
    cat(i, "")
  }  
  #return(list(cor_vec=cor_vec, coef_matrix=coef_matrix, MSE_vec=MSE_vec, n_rounds_vec=n_rounds_vec))
  return(list(cor_vec=cor_vec, MSE_vec=MSE_vec, n_rounds_vec=n_rounds_vec, coef_matrix=coef_matrix))
}


# 1000 bootstrap fits
# RUN: xg_b_obj_771_p
set.seed(123)
xg_b_obj_771_p <- xgboost_boot(fm=fm05, prolif_771genes, prolif_771genes, method="pearson", n_bootstraps=1000)
head(xg_b_obj_771_p$coef_matrix)[,1:6]
save(xg_b_obj_771_p, file="/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/xg_b_obj_771_p.RData")
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/xg_b_obj_771_p.RData")

# RUN: e_b_obj_771_ROR_p
set.seed(123)
e_b_obj_771_ROR_p <- xgboost_boot(fm=fm05, ROR_prolif_771genes, ROR_prolif_771genes, method="pearson", n_bootstraps=1000)
head(e_b_obj_771_ROR_p$coef_matrix)[,1:6]
save(e_b_obj_771_ROR_p, file="/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/e_b_obj_771_ROR_p.RData")
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/e_b_obj_771_ROR_p.RData")



set.seed(123)
bb_object_t <- XGboost_boot(fm01, Proliferation_6genes, Proliferation_6genes, n_bootstraps=10)
bb_object_t
bb_object <- boost_boot(fm01, Proliferation_6genes, Proliferation_6genes, n_bootstraps=1000)
bb_object <- boost_boot(fm05, dfA03, dfA03, n_bootstraps=1000)

bb_object_t

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
