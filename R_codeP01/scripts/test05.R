library(caret)
library(glmnet)
library(gbm)

# Split the data into training and testing sets
set.seed(123)
index <- createDataPartition(mtcars$mpg, p = 0.7, list = FALSE)
train_data <- mtcars[index, ]
test_data <- mtcars[-index, ]

# Fit the linear regression model
linear_model <- train(mpg ~ ., data = train_data, method = "lm", preProcess = c("center", "scale"))

# Fit the boosting model
boosting_model <- train(mpg ~ ., data = train_data, method = "gbm", preProcess = c("center", "scale"))

# Fit the lasso model
lasso_model <- train(mpg ~ ., data = train_data, method = "glmnet", preProcess = c("center", "scale"))

# Predict on test data
linear_prediction <- predict(linear_model, newdata = test_data)
boosting_prediction <- predict(boosting_model, newdata = test_data)
lasso_prediction <- predict(lasso_model, newdata = test_data)

# Calculate the MSE
linear_mse <- mean((linear_prediction - test_data$mpg)^2)
boosting_mse <- mean((boosting_prediction - test_data$mpg)^2)
lasso_mse <- mean((lasso_prediction - test_data$mpg)^2)

# Print the MSE values
cat("Linear MSE:", linear_mse, "\n")
cat("Boosting MSE:", boosting_mse, "\n")
cat("Lasso MSE:", lasso_mse, "\n")


library(caret)
library(glmnet)
library(gbm)

# Split the data into training and testing sets
set.seed(123)
index <- createDataPartition(mtcars$mpg, p = 0.7, list = FALSE)
train_data <- mtcars[index, ]
test_data <- mtcars[-index, ]

# Fit the linear regression model
linear_model <- train(mpg ~ ., data = train_data, method = "lm", preProcess = c("center", "scale"))

# Fit the boosting model
boosting_model <- train(mpg ~ ., data = train_data, method = "gbm", preProcess = c("center", "scale"))

# Fit the lasso model
lasso_model <- train(mpg ~ ., data = train_data, method = "glmnet", preProcess = c("center", "scale"))

# Predict on test data
linear_prediction <- predict(linear_model, newdata = test_data)
boosting_prediction <- predict(boosting_model, newdata = test_data)
lasso_prediction <- predict(lasso_model, newdata = test_data)

# Calculate the MSE
linear_mse <- mean((linear_prediction - test_data$mpg)^2)
boosting_mse <- mean((boosting_prediction - test_data$mpg)^2)
lasso_mse <- mean((lasso_prediction - test_data$mpg)^2)

# Print the MSE values
cat("Linear MSE:", linear_mse, "\n")
cat("Boosting MSE:", boosting_mse, "\n")
cat("Lasso MSE:", lasso_mse, "\n")
