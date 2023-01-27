# Here I added the following parameters:
#   max_depth = 1: This tells XGBoost to use a tree with a depth of 1 (i.e., a stump).
#   min_child_weight = 1: This tells XGBoost to use a single decision tree as the base learner
#   tree_method = "hist": This tells XGBoost to use histogram-based algorithm for decision tree learning.
#   You can use different tree_methods like exact, approx and hist, depending on the dataset size, memory and time constraints.
# 



# Load the training and test datasets
train <- read.csv("train.csv")
test <- read.csv("test.csv")

# Define the label and feature columns
label_column <- "target"
feature_columns <- setdiff(colnames(train), label_column)

# Define the xgboost parameters
param <- list(booster = "gbtree", 
              objective = "reg:squarederror", 
              eval_metric = "rmse", 
              nthread = 4, 
              max_depth = 1, 
              min_child_weight = 1,
              tree_method = "hist")

# Define the xgboost train control with large number of rounds and early stopping
train_control <- xgb.cv(data = as.matrix(train[, feature_columns]), 
                        label = train[, label_column], 
                        nrounds = 10000,
                        nfold = 5, 
                        early_stopping_rounds = 10, 
                        verbose = FALSE, 
                        print_every_n = 10, 
                        seed = 123, 
                        params = param)

# Select the optimal nrounds
best_nrounds <- train_control$best_ntreelimit

# Train the final model using the optimal nrounds
model <- xgboost(data = as.matrix(train[, feature_columns]), 
                 label = train[, label_column], 
                 nrounds = best_nrounds, 
                 nthread = 4, 
                 verbose = FALSE, 
                 params = param)

# Make predictions on the test dataset
predictions <- predict(model, as.matrix(test[, feature_columns]))






# Load the dataframe
data(mtcars)

# Convert the dataframe to a DMatrix object
dtrain <- xgb.DMatrix(data = as.matrix(mtcars[, -1]), label = mtcars$mpg)

dtrain <- xgb.DMatrix(data = data.matrix(mtcars[,-1]), label = mtcars$mpg)





# Create the dataset
data(agaricus.train, package='xgboost')
dtrain <- xgb.DMatrix(agaricus.train$data, label = agaricus.train$label)

# Set the xgboost parameters
params <- list(objective = "binary:logistic",
               eta = 0.3,
               max_depth = 2)

# Find the best iteration using the cv function
cv <- xgb.cv(params = params, data = dtrain, nrounds = 10, nfold = 5, showsd = T, stratified = T)
best_iteration <- which.min(cv$evaluation_log$test_error_mean)

# Train the model with the best iteration
xgb.train(params = params, data = dtrain, nrounds = best_iteration)

#######################################################################
# Load necessary libraries
library(xgboost)

# Load the training and test datasets
train <- read.csv("train.csv")
test <- read.csv("test.csv")

# Define the label and feature columns
label_column <- "target"
feature_columns <- setdiff(colnames(train), label_column)

# Define the xgboost parameters
param <- list(booster = "gbtree", eval_metric = "auc", nthread = 4)

# Define the xgboost train control with large number of rounds and early stopping
train_control <- xgb.cv(data = as.matrix(train[, feature_columns]), 
                        label = train[, label_column], 
                        nrounds = 10000,
                        nfold = 5, 
                        early_stopping_rounds = 10, 
                        verbose = FALSE, 
                        print_every_n = 10, 
                        seed = 123, 
                        params = param)

# Select the optimal nrounds
best_nrounds <- train_control$best_ntreelimit

# Train the final model using the optimal nrounds
model <- xgboost(data = as.matrix(train[, feature_columns]), 
                 label = train[, label_column], 
                 nrounds = best_nrounds, 
                 nthread = 4, 
                 verbose = FALSE, 
                 params = param)

# Make predictions on the test dataset
predictions <- predict(model, as.matrix(test[, feature_columns]))

# Perform any post-processing or evaluation on the predictions




library(glmnet)

# load data
data(mtcars)
x = model.matrix(mpg ~ ., mtcars)[,-1]
y = mtcars$mpg

# specify a range of alpha values to search
alphas = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)

# specify lambda values
lambdas = 10^seq(3, -3, length.out = 100)

# run elastic net with cross-validation and alpha search
fit = cv.glmnet(x, y, alpha = alphas, lambda = lambdas, nfolds = 5)

# find the row with the minimum cross-validated error
idx = which.min(fit$cvm)

# extract the corresponding alpha and lambda values
best_alpha = fit$alpha[idx]
best_lambda = fit$lambda[idx]

# print the best alpha and lambda values
print(paste("Best alpha:", best_alpha))
print(paste("Best lambda:", best_lambda))

____________
min_error = min(fit$cvm)
idx = which(fit$cvm == min_error)
best_alpha = fit$alpha[idx]
best_lambda = fit$lambda[idx]



# Create a dataset for demonstration purposes
set.seed(123)
x <- matrix(rnorm(200*5), ncol = 5)
y <- rnorm(200)

# Fit the model using mboost
model <- mboost(x, y, control = boost_control(mstop = 5, base_learners = "stump"))
model <- gamboost(x, y, control = boost_control(mstop = 5, base_learners = "stump"))
model <- gamboost(y~.,data = df,control = boost_control(mstop = 5, base_learners = "stump"))
# Print the model
print(model)