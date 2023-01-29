library(glmnet)
library(caret)

# Define the number of repeats and folds
repeats = 200
folds = 5

# Load your data into a variable called "data"
# Make sure the response variable is in the first column and the feature variables are in the remaining columns
data = Proliferation_6genes
  


# Initialize variables to store the results
# correlations = c()
# selected_features = matrix(0, nrow = repeats * folds, ncol = ncol(data) - 1)

# Initialize lists to store the results
correlations <- list()
selected_features <- list()


# Repeat the cross-validation process
for(i in 1:repeats) {
  # Create the folds
  cv_split = createFolds(data[,1], k = folds, list = TRUE, returnTrain = TRUE)
  
  # Fit the Lasso model and collect the selected features
  fit = cv.glmnet(data[-1], data[,1], family = "gaussian", alpha = 1, standardize = TRUE, type.measure = "mse", foldid = cv_split, nlambda = 100)
  selected_features[[i]] = (coef(fit, s = "lambda.min") != 0)
  
  # Collect the correlations
  predictions <- predict(fit, newx = data[-1])
  correlations[[i]] <- cor(predictions, data[,1])
}





# Fit the Lasso model and collect the results for each repeat and fold
for (i in 1:length(cv_split)) {
  train_data = data[cv_split[[i]],]
  test_data = data[-cv_split[[i]],]
  
  # Fit the Lasso model on the training data
  lasso_fit = cv.glmnet(x = as.matrix(train_data[,-1]), y = train_data[,1], family = "gaussian", alpha = 1)
  
  # Make predictions on the test data
  test_predictions = predict(lasso_fit, newx = as.matrix(test_data[,-1]))
  
  # Calculate the correlation between the predicted and actual test responses
  corr = cor(test_predictions, test_data[,1])
  correlations = c(correlations, corr)
  
  # Collect the selected features from the Lasso model
  selected_features[i, colSums(abs(coef(lasso_fit)) > 0) > 0] = selected_features[i, colSums(abs(coef(lasso_fit)) > 0) > 0] + 1
}
