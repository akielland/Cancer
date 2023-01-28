library(caret)


# Custom function for repeated k-fold cross validation
repeated_kfold_cv <- function(data, response=Y, func=lasso_sample, k=5, repeats=100) {
  # Initialize empty lists to store results
  results_list <- list()
  
  # Loop for the number of repeats
  for (i in 1:repeats) {
    # Perform k-fold cross validation
    print(data$response)
    kf <- caret::createFolds(data$Y, k = k, list = TRUE, returnTrain = TRUE)
    # Loop through the folds
    for (j in 1:k) {
      # Get the training and testing data
      train_data <- data[kf[[j]],]
      test_data <- data[-kf[[j]],]
      test_data <- as.matrix(test_data |> select(-Y))
      
      # Fit the function on the training data and get the results
      fit.cv <- func(train_data)
      co <- coef(fit.cv, s = "lambda.min")
      #co <- (co != 0)
      co <- co[-which(rownames(co) == "(Intercept)"), ]
      inds <- which(co != 0)
      #co <- colnames(inds)
      print(names(co[inds]))
 
      # test_results <- predict(model_fit, newdata = test_data)
      
      test_results = predict(fit.cv, newx = test_data, type = "response", s = "lambda.min")  
      
      # Append the test results to the list
      results_list[[i]] <- c(test_results[[i]], test_results)
    }
  }
  
  # Return the list of results
  return(list(results_list, co))
}
ob_  <- repeated_kfold_cv(Proliferation_6genes, Y, lasso_sample, k=5, repeats=2)



lasso_sample <- function(train_data){
  X_ <- as.matrix(train_data |> select(CCND1, CCNE1, CDKN1A, ESR1, MYC, RB1))
  Y_ <- as.matrix(train_data |>  select(Y))
  fit.cv <- cv.glmnet(X_, Y_)
  return(fit.cv)
  }
  