## Various testing of lasso on 6 genes and nodes (depending on input dataframe)
# 1. simple lasso with glmnet
# 2. cross-validation to select features
# 3. LOO-cross validation to tests model
# 4. LOO-cross + 5-fold cross validation to test model
##############################################################################

library(glmnet)


# standard lasso
lasso01 <- glmnet(X, Y)
plot(lasso01, label = T)
coef(lasso01, s = 0.05)

# lasso w/ cross-validation to determine optimal lambda hyper-parameter
lasso02 = cv.glmnet(X, Y, nfolds=5)
co <- coef(lasso02, s = "lambda.min")
inds <- which(co != 0)
variables <- row.names(co)[inds]
variables <- variables[!(variables %in% '(Intercept)')]
show(variables)


pred_ROR <- predict(lasso02, newx = X, type = "response", s = "lambda.min")
cor(pred_ROR, Y)

df$id[df$id == 55 & df$gender == 'm'] <- "60"

df <- mutate(df, id = case_when(
  id == 30 ~ 40, 
  TRUE   ~ id 
))

# Using LOO cross-validation for testing 
# going through all observations with LOO as test set,

lasso_loo <- function(X, Y, lambda.min=TRUE){
  pred_values <- NULL
  for(i in 1:nrow(X)){
    X_train <- X[-i,]
    Y_train <- Y[-i,]
    X_test <- X[i,]
    
    fit <- glmnet(X_train, Y_train)
    
    # print selected predictors
    co <- coef(fit, s=0.05)
    inds <- which(co != 0)
    variables <- row.names(co)[inds]
    variables <- variables[!(variables %in% '(Intercept)')]
    print(variables)
    
    if (lambda.min == TRUE) {
      pred_values[i] = predict(fit, newx = X_test, type = "response", s = 0.05)      
    }
    if (lambda.min == FALSE) {
      pred_values[i] = predict(fit, newx = X_test, type = "response", s = 0.05)      
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



# Using LOO cross-validation for testing 
# going through all observations with LOO as test set,
# but here using 5 fold-cross validation for fitting models (estimating lambda and selecting predictors)

lasso_loo <- function(X, Y, lambda.min=TRUE){
  pred_values <- NULL
  for(i in 1:nrow(X)){
    X_train <- X[-i,]
    Y_train <- Y[-i,]
    X_test <- X[i,]
    
    lasso.cv <- cv.glmnet(X_train, Y_train, nfolds = 5)
    
    # print selected predictors
    co <- coef(lasso.cv, s = "lambda.min")
    inds <- which(co != 0)
    variables <- row.names(co)[inds]
    variables <- variables[!(variables %in% '(Intercept)')]
    print(variables)

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


