library(caret)
set.seed(123)
mtcars_pca_combined <- cbind(mtcars_pca_numeric, mtcars_pca_categorical)
mtcars_pca_combined$mpg <- as.numeric(mtcars_pca_combined$mpg)

# Splitting the data into training and test sets
splitIndex <- createDataPartition(mtcars_pca_combined$mpg, p = 0.7, list = FALSE)
train_data_pca_combined <- mtcars_pca_combined[ splitIndex,]
test_data_pca_combined <- mtcars_pca_combined[-splitIndex,]

# Define the 3 models
model1 <- lm(mpg ~ ., data = train_data_pca_combined)
model2 <- train(mpg ~ ., data = train_data_pca_combined, method = "gbm", verbose = FALSE)
model3 <- train(mpg ~ ., data = train_data_pca_combined, method = "glmnet", tuneLength = 5, verbose = FALSE)

# 10-fold cross-validation
control <- trainControl(method = "cv", number = 10)
model1_cv <- train(mpg ~ ., data = train_data_pca_combined, method = "lm", trControl = control)
model2_cv <- train(mpg ~ ., data = train_data_pca_combined, method = "gbm", trControl = control, verbose = FALSE)
model3_cv <- train(mpg ~ ., data = train_data_pca_combined, method = "glmnet", trControl = control, tuneLength = 5, verbose = FALSE)

# Evaluate the models on the test set
model1_test_pred <- predict(model1, newdata = test_data_pca_combined)
model2_test_pred <- predict(model2, newdata = test_data_pca_combined)
model3_test_pred <- predict(model3, newdata = test_data_pca_combined)

# Calculate the mean squared error (MSE) of each model on the test set
model1_test_mse <- mean((model1_test_pred - test_data_pca_combined$mpg)^2)
model2_test_mse <- mean((model2_test_pred - test_data_pca_combined$mpg)^2)
model3_test_mse <- mean((model3_test_pred - test_data_pca_combined$mpg)^2)

# Compare the MSE of each model
cat("Linear Regression MSE:", model1_test_mse, "\n")
cat("Gradient Boosting MSE:", model2_test_mse, "\n")
cat("Lasso Regression MSE:", model3_test)
    