





library(glmnet)


# standard lasso
lasso01 <- glmnet(as.matrix(X), as.matrix(Y))
plot(lasso01, label = T)
coef(lasso01, s = 0.02)

# lasso w/ cross-validation to determine optimal lambda hyperparameter
lasso02 = cv.glmnet(as.matrix(X), as.matrix(Y), nfolds=5)
co <- coef(lasso02, s = "lambda.min")
inds <- which(co != 0)
variables <- row.names(co)[inds]
variables <- variables[!(variables %in% '(Intercept)')]
show(variables)


# Using LOO cross-validation for testing 
# (going through all observations with LOO as test set)
lasso_loo <- function(X, Y, lambda.min=TRUE){
  X = as.matrix(X)
  Y = as.matrix(Y)
  n = length(Y)
  pred_values <- NULL
  for(i in 1:nrow(X)){
    X_test <- X[i,]
    X_train <- X[-i,]
    Y_train <- Y[-i,]
    
    lasso.cv <- cv.glmnet(X_train, Y_train, nfolds = 5)

    if (lambda.min == TRUE) {
      pred_values[i] = predict(lasso.cv, newx = X_test, type = "response", s = "lambda.min")      
    }
    if (lambda.min == FALSE) {
      pred_values[i] = predict(lasso.cv, newx = X_test, type = "response", s = "lambda.1se")      
    }
  }
  return(pred_values)
}

# Result of lasso with lambda at min:
pred_values_lasso.min <- lasso_loo(X, Y, TRUE)
cat("Correlation for lasso using lambda.min: ")
print(cor(pred_values_lasso.min, Y))
cat("MSE for lasso using lambda.min: ")
print(mean((pred_values_lasso.min - df04$Proliferation.Score)**2))
par(mfrow=c(1,2))
plot(pred_values_lasso.min, df04$Proliferation.Score, main = "Lasso")
abline(lm(df04$Proliferation.Score ~ pred_values_lasso.min), pred_values_lasso.min)

# Result of lasso with lambda at 1 se from min:
pred_values_lasso.1se <- lasso_loo(X, Y, FALSE)
cat("Correlation for lasso using lambda.min: ")
print(cor(pred_values_lasso.1se, PS))
cat("MSE for lasso using lambda.min: ")
print(mean((pred_values_lasso.1se - PS)**2))
plot(pred_values_lasso.1se, PS, main = "Lasso")


