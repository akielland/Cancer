#################################################################
##  Bloks for ensamble models ##
#################################################################
## input: - train data
##        - formula defining the signatures
## output: perdition Y

library(xgboost)
library(Matrix)
library(data.table)
library(caret)
library(ggplot2)



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



