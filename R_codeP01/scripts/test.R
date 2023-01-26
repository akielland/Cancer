


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