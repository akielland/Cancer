library(caret)
library(glmnet)

data("mtcars")

# Split data into training and testing sets
set.seed(123)
train_index <- createDataPartition(mtcars$mpg, p = 0.7, list = FALSE)
train_data <- mtcars[train_index, ]
test_data <- mtcars[-train_index, ]

# Perform PCA on the training data
train_data_numeric <- train_data[, c("disp", "hp", "drat", "wt", "qsec")]
pca_numeric <- prcomp(train_data_numeric, center = TRUE, scale. = TRUE)
train_data_pca_numeric <- data.frame(pca_numeric$x[,1:2])
colnames(train_data_pca_numeric) <- c("PC1", "PC2")

# Perform PCA on the categorical data
train_data_categorical <- train_data[, c("cyl", "gear", "am", "carb")]
train_data_categorical <- sapply(train_data_categorical, as.numeric)
train_data_categorical <- data.frame(train_data_categorical)
pca_categorical <- prcomp(train_data_categorical, center = TRUE, scale. = TRUE)
train_data_pca_categorical <- data.frame(pca_categorical$x[,1:2])
colnames(train_data_pca_categorical) <- c("PC1_cat", "PC2_cat")

# Combine the PCA components from the numeric and categorical data
train_data_pca_combined <- cbind(train_data_pca_numeric, train_data_pca_categorical)

# Train and evaluate a linear regression model using the PCA components
set.seed(456)
model_linear <- train(mpg ~ ., data = train_data_pca_combined, method = "lm", trControl = trainControl(method = "cv", number = 10), tuneLength = 5)

# Train and evaluate a boosting model using the PCA components
set.seed(456)
model_boosting <- train(mpg ~ ., data = train_data_pca_combined, method = "gbm", trControl = trainControl(method = "cv", number = 10), tuneLength = 5)

# Train and evaluate a lasso model using the PCA components
set.seed(456)
model_lasso <- train(mpg ~ ., data = train_data_pca_combined, method = "glmnet", trControl = trainControl(method = "cv", number = 10), tuneLength = 5)

# Compare the models
model_comparison <- resamples(list(linear = model_linear, boosting = model_boosting, lasso = model_lasso))
summary(model_comparison)




library(caret)
library(glmnet)

data("mtcars")

# Split data into training and testing sets
set.seed(123)
train_index <- createDataPartition(mtcars$mpg, p = 0.7, list = FALSE)
train_data <- mtcars[train_index, ]
test_data <- mtcars[-train_index, ]

# Perform PCA on the training data
train_data_numeric <- train_data[, c("disp", "hp", "drat", "wt", "qsec")]
pca_numeric <- prcomp(train_data_numeric, center = TRUE, scale. = TRUE)
train_data_pca_numeric <- data.frame(pca_numeric$x[,1:2])
colnames(train_data_pca_numeric) <- c("PC1", "PC2")

# Perform PCA on the categorical data
train_data_categorical <- train_data[, c("cyl", "gear", "am", "carb")]
train_data_categorical <- sapply(train_data_categorical, as.numeric)
train_data_categorical <- data.frame(train_data_categorical)
pca_categorical <- prcomp(train_data_categorical, center = TRUE, scale. = TRUE)
train_data_pca_categorical <- data.frame(pca_categorical$x[,1:2])
colnames(train_data_pca_categorical) <- c("PC1_cat", "PC2_cat")

# Combine the PCA components from the numeric and categorical data
train_data_pca_combined <- cbind(train_data_pca_numeric, train_data_pca_categorical)

# Train and evaluate a linear regression model using the PCA components
set.seed(456)
model_linear <- train(mpg ~ ., data = train_data_pca_combined, method = "lm", trControl = trainControl(method = "cv", number = 10), tuneLength = 5)

# Train and evaluate a boosting model using the PCA components
set.seed(456)
model_boosting <- train(mpg ~ ., data = train_data_pca_combined, method = "gbm", trControl = trainControl(method = "cv", number = 10), tuneLength = 5)

# Train and evaluate a lasso model using the PCA components
set.seed(456)
model_lasso <- train(mpg ~ ., data = train_data_pca_combined, method = "glmnet", trControl = trainControl(method = "cv", number = 10), tuneLength = 5)

# Compare the models
model_comparison <- resamples(list(linear = model_linear, boosting = model_boosting, lasso = model_lasso))
summary(model_comparison)



library(caret)
library(glmnet)

# Load the mtcars dataset
data(mtcars)

# Split the data into training and testing sets
set.seed(123)
indx <- createDataPartition(mtcars$mpg, p = 0.8, list = FALSE)
train_data <- mtcars[indx, ]
test_data <- mtcars[-indx, ]

# Perform PCA on the numeric variables
numeric_vars <- c("disp", "hp", "drat", "wt", "qsec")
pca_numeric <- prcomp(train_data[, numeric_vars], scale = TRUE)
pca_numeric_train <- as.data.frame(predict(pca_numeric, train_data[, numeric_vars]))
pca_numeric_test <- as.data.frame(predict(pca_numeric, test_data[, numeric_vars]))

# Perform PCA on the categorical variables
categorical_vars <- c("vs", "am", "gear", "carb")
pca_categorical <- prcomp(train_data[, categorical_vars], scale = TRUE)
pca_categorical_train <- as.data.frame(predict(pca_categorical, train_data[, categorical_vars]))
pca_categorical_test <- as.data.frame(predict(pca_categorical, test_data[, categorical_vars]))

# Combine the PCA components
train_data_pca <- cbind(pca_numeric_train, pca_categorical_train)
train_data_pca$mpg <- train_data$mpg
test_data_pca <- cbind(pca_numeric_test, pca_categorical_test)
test_data_pca$mpg <- test_data$mpg

# Fit the models
model_linear <- lm(mpg ~ ., data = train_data_pca)
model_boosting <- train(mpg ~ ., data = train_data_pca, method = "gbm", verbose = FALSE)
model_lasso <- cv.glmnet(as.matrix(train_data_pca[, -ncol(train_data_pca)]), train_data_pca$mpg, alpha = 1)

# Make predictions
preds_linear <- predict(model_linear, newdata = test_data_pca)
preds_boosting <- predict(model_boosting, newdata = test_data_pca)
preds_lasso <- predict(model_lasso, s = "lambda.min", newx = as.matrix(test_data_pca[, -ncol(test_data_pca)]))

# Evaluate the models
rmse_linear <- sqrt(mean((test_data_pca$mpg - preds_linear)^2))
rmse_boosting <- sqrt(mean((test_data_pca$mpg - preds_boosting)^2))
rmse_lasso <- sqrt(
)





library(caret)
library(glmnet)
library(xgboost)

data(mtcars)
mtcars_pca_numeric <- prcomp(mtcars[, c("mpg", "disp", "hp", "drat", "wt", "qsec")], center=TRUE, scale=TRUE)$x[,1:2]
mtcars_pca_categorical <- model.matrix(~factor(mtcars$am) - 1)
mtcars_pca_combined <- cbind(mtcars_pca_numeric, mtcars_pca_categorical)

set.seed(123)
cv_control <- trainControl(method="cv", number=10, returnResamp="all", savePredictions=TRUE)

model_linear <- train(mpg ~ ., data = mtcars_pca_combined, method = "lm", trControl = cv_control)
model_boosting <- train(mpg ~ ., data = mtcars_pca_combined, method = "xgbLinear", trControl = cv_control)
model_lasso <- train(mpg ~ ., data = mtcars_pca_combined, method = "glmnet", trControl = cv_control)

results <- resamples(list(Linear=model_linear, Boosting=model_boosting, Lasso=model_lasso))
summary(results)
