library(xgboost)
data(mtcars)

# Split the data into training and test sets
train_index <- sample(1:nrow(mtcars), nrow(mtcars)/2)
train <- mtcars[train_index, ]
test <- mtcars[-train_index, ]

# Convert the data into a matrix format required by xgboost
train_matrix <- xgb.DMatrix(data = as.matrix(train[, -1]), label = train[, 1])
test_matrix <- xgb.DMatrix(data = as.matrix(test[, -1]))

# Set the number of iterations
iterations <- 10

# Initialize the matrix to store feature selection results
feature_selection_matrix <- matrix(0, nrow = iterations, ncol = ncol(train) - 1)
colnames(feature_selection_matrix) <- colnames(train[, -1])

# Train the model in each iteration and extract the selected features
for (i in 1:iterations) {
  
  # Train the model
  params <- list(booster = "gbtree", objective = "reg:linear", nrounds = 2, max_depth = 1)
  model <- xgb.train(params = params, data = train_matrix, nrounds = 2)
  
  # Extract the selected features
  feature_importance <- xgb.importance(colnames(train_matrix), model = model)
  selected_features <- feature_importance[, 1][feature_importance[, 3] > 0]
  
  # Store the selected features in the feature selection matrix
  for (feature in selected_features) {
    feature_selection_matrix[i, feature] <- 1
  }
}


library(data.table)
library(xgboost)

# Load the data
train <- iris[,1:4]
train <- as.matrix(train)
labels <- iris[,5]

# Initialize matrix to store the feature importance
num_iterations <- 5
num_features <- ncol(train)
importance_matrix <- matrix(0, nrow=num_iterations, ncol=num_features)
colnames(importance_matrix) <- colnames(train)

# Loop through each iteration
for (i in 1:num_iterations) {
  # Fit XGBoost model
  model <- xgboost(data = train, label = labels, nrounds = 50, max_depth = 1, objective = "reg:squarederror",
                   verbose = FALSE, eta=0.05)
  print(importance_matrix)
  # Extract feature importance
  feature_importance <- data.table(xgb.importance(colnames(train), model = model))
  print(class(feature_importance))
  #feature_importance <- feature_importance[,1:2]
  
  # Sum up the feature importance
  importance_matrix[i, feature_importance$Feature] <- 1
}

# Check the resulting feature importance matrix
print(importance_matrix)









library(xgboost)

# Set the number of iterations and number of features
iterations <- 5
num_features <- 10

# Initialize a matrix to store feature selections in each iteration
selected_features_mat <- matrix(ncol = num_features, nrow = iterations)

for (i in 1:iterations) {
  # Train an xgboost model with stumps as base learners
  xgb_model <- xgboost(data = train_data, label = train_labels,
                       booster = "gbtree", eta = 1, max_depth = 1, 
                       objective = "reg:squarederror", nrounds = 10)
  
  # Select top `num_features` features
  selected_features <- head(colnames(train_data)[order(xgb_model$importance, decreasing = T)], num_features)
  
  # Store the selected features in the `selected_features_mat` matrix
  selected_features_mat[i,] <- selected_features
}



library(xgboost)

# Set the number of iterations and number of features
iterations <- 2
num_features <- 3

train_data <- Proliferation_6genes
train_sparse = sparse.model.matrix(object = fm01, data = train_data)
dtrain = xgb.DMatrix(data = train_sparse, label = train_data$Y)

# Initialize a matrix to store feature selections in each iteration
selected_features_mat <- matrix(ncol = num_features, nrow = iterations)

for (i in 1:iterations) {
  # Train an xgboost model using a formula
  xgb_model <- xgboost(data = dtrain,
                       booster = "gbtree", eta = 1, max_depth = 1, eval_metric = "rmse", tree_method = "hist",
                       objective = "reg:squarederror", nrounds = 8)
  
  print(class(xgb.importance(model = xgb_model)))
  
  # colnames(train_sparse)[order(xgb.importance(model = xgb_model))]
  # 
  # # Select top `num_features` features
  # selected_features <- head(colnames(train_sparse)[order(xgb_model$importance, decreasing = T)], num_features)
  # 
  # # Store the selected features in the `selected_features_mat` matrix
  # selected_features_mat[i,] <- selected_features
}

