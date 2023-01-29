library(caret)
library(ggplot2)

# Custom function for repeated k-fold cross validation
repeated_kfold_cv <- function(data, func=lasso_sample, folds=5, repeats=100) {
  
  correlations <- numeric(repeats * folds)
  print(length(correlations))
  SEM_vec <- rep(NA, repeats * folds)
  coef_matrix <- matrix(NA, nrow = repeats*folds, ncol = ncol(data[,-1]))
  colnames(coef_matrix) <- colnames(data[, -1])
  coef_matrix_row_index <- 1
  
  # Repeat the cross-validation process
  for (i in 1:repeats) {
    # Create the folds for evaluating the performance
    kf <- caret::createFolds(data[,1], k = folds, list = TRUE, returnTrain = TRUE)
    # Loop through the folds
    for (j in 1:folds) {
      # Get the training and testing data
      train_data <- data[kf[[j]],]
      test_data <- data[-kf[[j]],]

      # Fit the function on the training data and get results
      fit <- func(train_data)
      
      coef_matrix[coef_matrix_row_index, ] <- coef(fit, s = "lambda.min")[-1]
      
      pred = predict(fit, newx = as.matrix(test_data)[,-1], type = "response", s = "lambda.min")  
      correlations[coef_matrix_row_index]  <- suppressWarnings(cor(pred, test_data[,1]))
      MSE_vec <- mean((pred - test_data[,1])^2)
      
      cat(coef_matrix_row_index, "")
      coef_matrix_row_index <- coef_matrix_row_index + 1
    }
  }
  return(list(correlations=correlations, coef_matrix=coef_matrix, MSE_vec=MSE_vec))
}
# Set repeats and folds of the cross-validations
repeats = 100
folds = 5

lasso_k_ob  <- repeated_kfold_cv(Proliferation_ALLgenes, lasso_sample, folds, repeats)
head(lasso_k_ob$coef_matrix)[,1:10]
mean(na.omit(lasso_k_ob$correlations))

# Order the features based on their selection frequency
frequency <- data.frame(Feature = colnames(lasso_k_ob$coef_matrix), Frequency = colSums(lasso_k_ob$coef_matrix != 0) / (repeats * folds))
frequency <- frequency[order(frequency$Frequency, decreasing = TRUE),]
frequency[1:3, 1:2]


# Create a bar plot of the selection frequency of the features
ggplot(frequency[1:50,], aes(x = Frequency, y = reorder(Feature, Frequency))) +
  geom_bar(stat = "identity") +
  xlab("Selection Frequency") +
  ylab("Features") +
  ggtitle("Selection Frequency of Features") +
  theme(axis.text.y = element_text(angle = 0, hjust = 0))

# Create a histogram ver correlations
correlations_finite <- lasso_k_ob$correlations[is.finite(lasso_k_ob$correlations)]
corr_df <- data.frame(correlation = correlations_finite)

ggplot(corr_df, aes(x=correlation)) +
  geom_histogram(bins = 30, color = "black", fill = "white") +
  xlab("Correlation") +
  ylab("Frequency") +
  ggtitle("Histogram of Correlation Values")

# Create a histogram ver correlations
correlations_finite <- lasso_k_ob$correlations[is.finite(lasso_k_ob)]
corr_df <- data.frame(correlation = correlations_finite)

ggplot(corr_df, aes(x=correlation)) +
  geom_histogram(bins = 30, color = "black", fill = "white") +
  xlab("Correlation") +
  ylab("Frequency") +
  ggtitle("Histogram of Correlation Values")
