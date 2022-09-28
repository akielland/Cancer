

# Leave one observation out for testing and going through all observations
lasso_loo <- function(X, Y, lambda){
  X = as.matrix(X)
  Y = as.matrix(Y)
  n = length(Y)
  pred_values <- NULL
  for(i in 1:nrow(X)){
    X_test <- X[i,]
    X_train <- X[-i,]
    Y_test <- Y[i,]
    Y_train <- Y[-i,]
    
    lasso.cv <- cv.glmnet(X_train, Y_train)
    if (lambda == one_se) {
      pred_values[i] = predict(lasso.cv, newx = X_test, type = "response", s = "lambda.1se")      
    }
    if (lambda = min) {
      pred_values[i] = predict(lasso.cv, newx = X_test, type = "response", s = "lambda.min")      
    }
  }
  return(pred_values)
}

X <- dplyr::select(df02_LR,  scr_CCND1, scr_CCNE1, scr_CDKN1A, scr_ESR1, scr_MYC, scr_RB1)
Y <- select(df02_LR, sur_ProliferationScore)

pred_values_lasso <- lasso_loo(X, Y, min)

# Result of simple linear model:
par(mfrow=c(1,2))
cat("Correlation for lasso using lambda.min: ")
print(cor(pred_values_lasso, df02_LR$sur_ProliferationScore))
cat("MSE for lasso using lambda.min: ")
print(mean((pred_values_lasso - df02_LR$sur_ProliferationScore)**2))
plot(pred_values_lasso, df02_LR$sur_ProliferationScore, main = "Lasso")

# Result of linear model with interactions:
cat("Correlation for linear model with interactions: ")
print(cor(pred_values_fm02, df02_LR$sur_ProliferationScore))
cat("MSE for linear model with interactions: ")
print(mean((pred_values_fm02 - df02_LR$sur_ProliferationScore)**2))
plot(pred_values_fm02, df02_LR$sur_ProliferationScore, main = "Linear model with interaction terms")
par(mfrow=c(1,1))


pred.lasso.loop = function(X, Y){
  df_ <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(df_) <- c("lambda.min", "Y_test")
  
  for (i in c(1:10)) {
    df_ = rbind(df_, prediction.lasso(X,Y))
  }  
  return(df_)
}

pip.test.l = pred.lasso.loop(X.sim.df, y)

## Results: Lasso prediction

plot(pip.test.l[,1], pip.test.l[,2])
abline(lm(pip.test.l[,2] ~ pip.test.l[,1]), col = "red", lwd=3)
mtext(paste("Correlation:", round(cor(pip.test.l[,2], pip.test.l[,1]), 3)), side=1, line=-1, cex=1.4)

mean((pip.test.l[,1] - pip.test.l[,2])^2)

