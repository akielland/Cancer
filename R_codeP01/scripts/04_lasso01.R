
X <- dplyr::select(df02_LR,  scr_CCND1, scr_CCNE1, scr_CDKN1A, scr_ESR1, scr_MYC, scr_RB1)
Y <- select(df02_LR, sur_ProliferationScore)


# Leave one observation out for testing and going through all observations
lasso_loo <- function(X, Y, lambda.min=TRUE){
  X = as.matrix(X)
  Y = as.matrix(Y)
  n = length(Y)
  pred_values <- NULL
  for(i in 1:nrow(X)){
    X_test <- X[i,]
    X_train <- X[-i,]
    # Y_test <- Y[i,]
    Y_train <- Y[-i,]
    
    lasso.cv <- cv.glmnet(X_train, Y_train)

    if (lambda.min == TRUE) {
      pred_values[i] = predict(lasso.cv, newx = X_test, type = "response", s = "lambda.min")      
    }
    if (lambda.min == FALSE) {
      pred_values[i] = predict(lasso.cv, newx = X_test, type = "response", s = "lambda.1se")      
    }
  }
  return(pred_values)
}

pred_values_lasso.min <- lasso_loo(X, Y, TRUE)
pred_values_lasso.1se <- lasso_loo(X, Y, FALSE)

# Result of lasso with lambda at min:
par(mfrow=c(1,2))
cat("Correlation for lasso using lambda.min: ")
print(cor(pred_values_lasso.min, df02_LR$sur_ProliferationScore))
cat("MSE for lasso using lambda.min: ")
print(mean((pred_values_lasso.min - df02_LR$sur_ProliferationScore)**2))
plot(pred_values_lasso.min, df02_LR$sur_ProliferationScore, main = "Lasso")

# Result of lasso wih lambda at 1 se from min:
cat("Correlation for lasso using lambda.min: ")
print(cor(pred_values_lasso.1se, df02_LR$sur_ProliferationScore))
cat("MSE for lasso using lambda.min: ")
print(mean((pred_values_lasso.1se - df02_LR$sur_ProliferationScore)**2))
plot(pred_values_lasso.1se, df02_LR$sur_ProliferationScore, main = "Lasso")
par(mfrow=c(1,1))

######################################################################

# boostrap
lasso_bootstrap <- function(X, Y, lambda.min=TRUE){
  X = as.matrix(X)
  Y = as.matrix(Y)
  N = length(Y)

  int <- sample(N, size = N*0.632, replace = TRUE)
  X_train = X[int,]
  X_test = X[-int,]
  Y_train = Y[int]
  Y_test = Y[-int]
  
  lasso.cv <- cv.glmnet(X_train, Y_train)
  
  if (lambda.min == TRUE) {
    pred_values = predict(lasso.cv, newx = X_test, type = "response", s = "lambda.min")      
  }
  if (lambda.min == FALSE) {
    pred_values = predict(lasso.cv, newx = X_test, type = "response", s = "lambda.1se")      
  }
  suppressWarnings(out <- cor(pred_values, Y_test))
  return(out)
}

lasso_bootstrap(X, Y, TRUE)

lasso_cor_boot = function(X, Y, n_bootstraps){
  cor_vec <- rep(NA, n_bootstraps)
  for (i in c(1:n_bootstraps)) {
    cor_vec[i] = lasso_bootstrap(X, Y, TRUE)
  }  
  return(cor_vec)
}

cor_vec_boot <- lasso_cor_boot(X,Y,1000)
sum(is.na(cor_vec_boot))/1000
mean(cor_vec_boot, na.rm=TRUE)
var(cor_vec_boot, na.rm=TRUE)
hist(cor_vec_boot, breaks = 50)



