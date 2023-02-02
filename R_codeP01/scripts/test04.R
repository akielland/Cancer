library(caret)
library(glmnet)
library(gbm)

# Split the data into training and testing sets
set.seed(123)
index <- createDataPartition(mtcars$mpg, p = 0.7, list = FALSE)
train_data <- mtcars[index, ]
test_data <- mtcars[-index, ]

# Numeric features PCA
train_pca_numeric <- prcomp(train_data[, c("disp", "hp", "drat", "wt", "qsec")], center = TRUE, scale. = TRUE)
test_pca_numeric <- predict(train_pca_numeric, newdata = test_data[, c("disp", "hp", "drat", "wt", "qsec")])
train_pca_numeric_components <- data.frame(pc1 = train_pca_numeric$x[, 1], pc2 = train_pca_numeric$x[, 2], pc3 = train_pca_numeric$x[, 3], pc4 = train_pca_numeric$x[, 4], pc5 = train_pca_numeric$x[, 5])
test_pca_numeric_components <- data.frame(pc1 = test_pca_numeric[, 1], pc2 = test_pca_numeric[, 2], pc3 = test_pca_numeric[, 3], pc4 = test_pca_numeric[, 4], pc5 = test_pca_numeric[, 5])

# Categorical features PCA
train_pca_categorical <- model.matrix(mpg ~ . -1, data = train_data[, c("cyl", "vs", "am", "gear", "carb")])
test_pca_categorical <- model.matrix(mpg ~ . -1, data = test_data[, c("cyl", "vs", "am", "gear", "carb")])

# Combine numeric and categorical PCA components
train_pca_combined <- cbind(train_pca_numeric_components, train_pca_categorical)
test_pca_combined <- cbind(test_pca_numeric_components, test_pca_categorical)

# Fit the linear regression model
linear_model <- train(mpg ~ ., data = train_data, method = "lm", preProcess = c("center", "scale"))

# Fit the boosting model
boosting_model <- train(mpg ~ ., data = train_data, method = "gbm", preProcess = c("center", "scale"))

# Fit the lasso model
lasso_model <- train(mpg ~ ., data = train_data, method = "glmnet", preProcess = c("center", "scale"))

# Predict on test data
linear_prediction <- predict(linear_model, newdata = test_data)
boosting_prediction <- predict(boosting_model, newdata = test_data)
lasso_pred





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



## Lasso
# Load required libraries
library(stats)
library(ggplot2)
library(caret)
library(glmnet)

# Load sample data
data(mtcars)

# Split the data into two subsets
mtcars_numeric <- mtcars[, c("mpg", "disp", "hp", "drat", "wt")]
mtcars_categorical <- mtcars[, c("am", "gear", "carb")]

# Perform PCA on the numeric subset
pca_results_numeric <- princomp(mtcars_numeric, cor = TRUE)

# Get the number of components that explain at least 50% of the variance
n_components_numeric <- sum(pca_results_numeric$sdev^2 / sum(pca_results_numeric$sdev^2) >= 0.5)

# Extract the first n_components_numeric principal components
mtcars_pca_numeric <- data.frame(pc1_numeric = pca_results_numeric$scores[, 1:n_components_numeric])

# Perform PCA on the categorical subset
pca_results_categorical <- princomp(mtcars_categorical, cor = TRUE)

# Get the number of components that explain at least 50% of the variance
n_components_categorical <- sum(pca_results_categorical$sdev^2 / sum(pca_results_categorical$sdev^2) >= 0.5)

# Extract the first n_components_categorical principal components
mtcars_pca_categorical <- data.frame(pc1_categorical = pca_results_categorical$scores[, 1:n_components_categorical])

# Combine the numeric and categorical PCA results and use them as predictors in a Lasso model
mtcars_pca_combined <- cbind(mtcars_pca_numeric, mtcars_pca_categorical, mpg = mtcars$mpg)
model_fit <- cv.glmnet(as.matrix(mtcars_pca_combined[, 1:n_components_numeric + n_components_categorical]), mtcars_pca_combined$mpg, alpha = 1, family = "gaussian")

# Plot the cross-validation curve
plot(model_fit)



## boosting
# Load required libraries
library(stats)
library(ggplot2)
library(caret)
library(xgboost)

# Load sample data
data(mtcars)

# Split the data into two subsets
mtcars_numeric <- mtcars[, c("mpg", "disp", "hp", "drat", "wt")]
mtcars_categorical <- mtcars[, c("am", "gear", "carb")]

# Perform PCA on the numeric subset
pca_results_numeric <- princomp(mtcars_numeric, cor = TRUE)

# Get the number of components that explain at least 50% of the variance
n_components_numeric <- sum(pca_results_numeric$sdev^2 / sum(pca_results_numeric$sdev^2) >= 0.5)

# Extract the first n_components_numeric principal components
mtcars_pca_numeric <- data.frame(pc1_numeric = pca_results_numeric$scores[, 1:n_components_numeric])

# Perform PCA on the categorical subset
pca_results_categorical <- princomp(mtcars_categorical, cor = TRUE)

# Get the number of components that explain at least 50% of the variance
n_components_categorical <- sum(pca_results_categorical$sdev^2 / sum(pca_results_categorical$sdev^2) >= 0.5)

# Extract the first n_components_categorical principal components
mtcars_pca_categorical <- data.frame(pc1_categorical = pca_results_categorical$scores[, 1:n_components_categorical])

# Combine the numeric and categorical PCA results and use them as predictors in a XGBoost model
mtcars_pca_combined <- cbind(mtcars_pca_numeric, mtcars_pca_categorical, mpg = mtcars$mpg)
dtrain <- xgb.DMatrix(data = as.matrix(mtcars_pca_combined[, 1:n_components_numeric + n_components_categorical]), label = mtcars_pca_combined$mpg)
model_fit <- xgboost(data = dtrain, nrounds = 100, objective = "reg:linear", eval_metric = "rmse")

# Print summary of the model
print(model_fit)

# Load required libraries
library(stats)
library(ggplot2)
library(caret)
library(lightgbm)

# Load sample data
data(mtcars)

# Split the data into two subsets
mtcars_n




## linear model
library(stats)
library(ggplot2)
library(caret)

# Load sample data
data(mtcars)

# Split the data into two subsets
mtcars_numeric <- mtcars[, c("mpg", "disp", "hp", "drat", "wt")]
mtcars_categorical <- mtcars[, c("am", "gear", "carb")]

# Perform PCA on the numeric subset
pca_results_numeric <- princomp(mtcars_numeric, cor = TRUE)

# Get the number of components that explain at least 50% of the variance
n_components_numeric <- sum(pca_results_numeric$sdev^2 / sum(pca_results_numeric$sdev^2) >= 0.5)

# Extract the first n_components_numeric principal components
mtcars_pca_numeric <- data.frame(pc1_numeric = pca_results_numeric$scores[, 1:n_components_numeric])

# Perform PCA on the categorical subset
pca_results_categorical <- princomp(mtcars_categorical, cor = TRUE)

# Get the number of components that explain at least 50% of the variance
n_components_categorical <- sum(pca_results_categorical$sdev^2 / sum(pca_results_categorical$sdev^2) >= 0.5)

# Extract the first n_components_categorical principal components
mtcars_pca_categorical <- data.frame(pc1_categorical = pca_results_categorical$scores[, 1:n_components_categorical])

# Combine the numeric and categorical PCA results and use them as predictors in a linear regression model
mtcars_pca_combined <- cbind(mtcars_pca_numeric, mtcars_pca_categorical, mpg = mtcars$mpg)
model_fit <- lm(mpg ~ ., data = mtcars_pca_combined)

# Print summary of the model
summary(model_fit)





# Load required libraries
library(prcomp)
library(ggplot2)
library(caret)

# Load sample data
data(mtcars)

# Split the data into two subsets
mtcars_numeric <- mtcars[, c("mpg", "disp", "hp", "drat", "wt")]
mtcars_categorical <- mtcars[, c("am", "gear", "carb")]

# Perform PCA on the numeric subset
pca_results_numeric <- prcomp(mtcars_numeric, scale = TRUE)

# Get the number of components that explain at least 50% of the variance
n_components_numeric <- sum(pca_results_numeric$sdev^2 / sum(pca_results_numeric$sdev^2) >= 0.5)

# Extract the first n_components_numeric principal components
mtcars_pca_numeric <- data.frame(pc1 = pca_results_numeric$x[, 1:n_components_numeric])



# Perform PCA on the categorical subset
pca_results_categorical <- prcomp(mtcars_categorical, scale = TRUE)

# Get the number of components that explain at least 50% of the variance
n_components_categorical <- sum(pca_results_categorical$sdev^2 / sum(pca_results_categorical$sdev^2) >= 0.5)

# Extract the first n_components_categorical principal components
mtcars_pca_categorical <- data.frame(pc1 = pca_results_categorical$x[, 1:n_components_categorical])

# Combine the numeric and categorical PCA results and use them as predictors in a linear regression model
mtcars_pca_combined <- cbind(mtcars_pca_numeric, mtcars_pca_categorical)
model_fit <- lm(mpg ~ ., data = mtcars_pca_combined)

# Print summary of the model
summary(model_fit)





# Load required libraries
library(prcomp)
library(ggplot2)
library(caret)

# Load sample data
data(mtcars)

# Perform PCA on the data
pca_results <- prcomp(mtcars, scale = TRUE)

# Get the number of components that explain at least 50% of the variance
n_components <- sum(pca_results$sdev^2 / sum(pca_results$sdev^2) >= 0.5)

# Extract the first n_components principal components and use them as predictors in a linear regression model
mtcars_pca <- data.frame(pc1 = pca_results$x[, 1:n_components], mpg = mtcars$mpg)

model_fit <- lm(mpg ~ ., data = mtcars_pca)

# Print summary of the model
summary(model_fit)

# Plot the regression results
ggplot(mtcars_pca, aes(x = pc1, y = mpg)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("PCA-Regression Results")






# Load required libraries
library(prcomp)
library(ggplot2)
library(caret)

# Load sample data
data(mtcars)

# Perform PCA on the data
pca_results <- prcomp(mtcars, scale = TRUE)

# Extract the first two principal components and use them as predictors in a linear regression model
mtcars_pca <- data.frame(pc1 = pca_results$x[, 1], pc2 = pca_results$x[, 2], mpg = mtcars$mpg)

model_fit <- lm(mpg ~ pc1 + pc2, data = mtcars_pca)

# Print summary of the model
summary(model_fit)

# Plot the regression results
ggplot(mtcars_pca, aes(x = pc1, y = pc2, color = mpg)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("PCA-Regression Results")
