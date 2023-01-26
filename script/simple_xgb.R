library(xgboost)
library(Matrix)


# To make it fully self-contained, let me simulate some data: 

df = data.frame(X1 = rnorm(n = 10000, mean = 50, sd = 20), 
                X2 = floor(runif(10000, min = 1, max = 10)))

df$y = 0.07*df$X1 + 0.1*df$X2 + 0.1*sin(df$X1*df$X2) + rnorm(10000, mean = 0, sd = 1)



# Split into training and test set
train_df = df[1:5000,]
validation_df = df[5001:8000,]
test_df = df[8001:10000,]

### XGBOOST
# XGBoost must get the data in a so-called xgb.DMatrix object;
# We must do a two-step transformation: from normal data frame to sparse model matrix to xgb.DMatrix
# I don't really understand all the details of this, and it is not super important

# First, specify the formula for the prediction model: 
model_formula = as.formula(y ~ X1 + X2)


# Train set: 
train_sparse = sparse.model.matrix(object = model_formula, data = train_df)
dtrain = xgb.DMatrix(data = train_sparse,label = train_df$y)

# Validation set
validation_sparse = sparse.model.matrix(object = model_formula, data = validation_df)
dvalidation = xgb.DMatrix(data = validation_sparse,label = validation_df$y)

# Test set
test_sparse = sparse.model.matrix(object = model_formula, data = test_df)
dtest = xgb.DMatrix(data = test_sparse, label = test_df$y)


# Hyper parameters:
# eta: the learning rate. Comparable to step size in gradient descent algorithms
# max_depth: the max depth of each tree. To get stumps: set tree_depth = 1

xgb_params = list(eta = 0.1, max_depth = 5)

# How many sequential trees should we train? Numbre of iterations
NROUNDS = 1000

# THIS IS WHAT YOU NEED TO UNDERSTAND:
# If we take small steps (low eta), we require many trees (nrounds) to get a good result.
# If we make deep trees then each tree is better => probably need fewer trees to get a good result

# Now we are ready to train the actual model:
start_time = Sys.time()
xgb_model = xgb.train(data = dtrain,
                      # The following is needed if to monitor training and validation error
                      watchlist = list(train = dtrain, 
                                       val = dvalidation),
                      params = xgb_params,
                      nrounds = NROUNDS,
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

# Feature importance:
xgb.importance(model = xgb_model)

# verbose = TRUE will print out stuff from the training process

# We can make new predictions with the model like this: 
xgb_predictions = predict(object = xgb_model, 
                          newdata = dtest)

# Evaluate
xgb_error = xgb_predictions - test_df$SalePric
xgb_rmse = sqrt(mean(xgb_error^2))
