#################################################################
##  Bootstrap  Boosting - using original sample as test set ##
#################################################################

## Boosting using bootstrap sample to train 1000 models
## Baselearner:
## - linear
## - spline
## - btree
## Cross-validation (5-fold) used for training/tuning lambda and selecting features
## Test against original data sample
## output: proliferation.score correlation; SEM; Coefficient of genes

library(mboost)


boost_bootstrap_sample <- function(df_train, df_test, method="pearson"){
  # one bootstrap sample
  # output: 
  # 1. correlation between prediction and full input sample  (pearson or spearman should be set)
  # 2. Coefficients of variables 
  # 3. MSE
  n = dim(df_train)[1]
  
  int <- sample.int(n, size = n, replace = TRUE)
  X_train = X[int,]
  Y_train = Y[int]
  
  fit <- gamboost(fm01, data = Proliferation_6genes,
                       baselearner = "bbs", # dfbase=4
                       control = boost_control(mstop = 50))
  
  
  
  co <- as.vector(coef(fit.cv, s = "lambda.min")) # as.vector makes an easier object to work with
  
  pred_fit = predict(fit, Proliferation_6genes)
  
  cor <- suppressWarnings(cor(pred_values, ext.Y, method = method))
  MSE <- mean((ext.Y - pred_values)^2)
  return(list(cor=cor, co=co, MSE=MSE))
}





# linear base learner
fm01 <- Y ~ CCND1 + CCNE1 + CDKN1A + ESR1 + MYC + RB1
boost_m01 = glmboost(fm01, data = Proliferation_6genes)
coef(boost_m01, which = "")
par(mfrow=c(1,2))
plot(boost_m01, off2int=TRUE)
plot(boost_m01)

pred_boost01 = predict(boost_m01, Proliferation_6genes)

cor(pred_boost01, Proliferation_6genes$Y)
plot(pred_boost01, Proliferation_6genes$Y)
abline(lm(Proliferation_6genes$Y ~ pred_boost01))


# smooth - P-spline as base learner
spline01 <- gamboost(fm01, data = Proliferation_6genes,
                     baselearner = "bbs", # dfbase=4
                     control = boost_control(mstop = 50))

coef(spline01, which = "")
par(mfrow=c(2,3))
plot(spline01, off2int=TRUE)

pred_spline01 = predict(spline01, Proliferation_6genes)

cor(pred_spline01, Proliferation_6genes$Y)
par(mfrow=c(1,1))
plot(pred_spline01, Proliferation_6genes$Y)
abline(lm(Proliferation_6genes$Y ~ pred_spline01))

# boosting with stumps
stump01 = gamboost(fm01, data = Proliferation_6genes, 
                   family=Gaussian(), 
                   baselearner='btree',
                   boost_control(mstop = 100))

cols_chosen = unique(stump01$xselect())
colnames(Proliferation_6genes$Y)[cols_chosen+2]

par(mfrow=c(2,3))
plot(stump01)


pred_boost01 = predict(boost_m01, Proliferation_6genes)

cor(pred_boost01, Proliferation_6genes$Y)
plot(pred_boost01, Proliferation_6genes$Y)
abline(lm(Proliferation_6genes$Y ~ pred_boost01))
