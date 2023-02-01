######################################################################
##  mBoost on bootstrap samples - using original sample as test set ##
######################################################################

## Boosting using bootstrap sample to train 1000 models
## Base learners:
## - linear
## - spline
## - btree
## 
## Test against original data sample
## output: proliferation.score correlation; SEM; Coefficient of genes

library(mboost)

# linear base learner
fm01 <- Y ~ CCND1 + CCNE1 + CDKN1A + ESR1 + MYC + RB1
fm05 <- as.formula(paste("Y", paste(genes, collapse="+"), sep=" ~ "))

boost_bootstrap_sample <- function(fm, df_train, df_test, method="pearson"){
  # one bootstrap sample
  # output: 
  # 1. correlation between prediction and full input sample  (pearson or spearman should be set)
  # 2. Coefficients of variables 
  # 3. MSE
  n = dim(df_train)[1]
  
  int <- sample.int(n, size = n, replace = TRUE)
  data_train = df_train[int,]
  
  # fit <- glmboost(fm, data = data_train)
  # fit <- gamboost(fm, data = data_train,
  #                      baselearner = "bbs", # dfbase=4
  #                      control = boost_control(mstop = 50))
  
  fit = gamboost(fm, data = data_train, 
                     family=Gaussian(), 
                     baselearner='btree',
                     boost_control(mstop = 10))
  
  # co <- coef(fit, which = "")
  # print(coef(fit, which = ""))
  cols_chosen = unique(fit$xselect())
  co <- colnames(Proliferation_6genes)[cols_chosen]
  inds <- cols_chosen

  pred_values = predict(fit, df_test)
  
  cor <- suppressWarnings(cor(pred_values, df_test$Y, method = method))
  MSE <- mean((df_test$Y - pred_values)^2)
  return(list(cor=cor, co=co, inds=inds, MSE=MSE))
}


boost_bootstrap_sample(fm01, Proliferation_6genes, Proliferation_6genes)

# 1000 bootstrap fits
boost_boot = function(fm, df_train, df_test, method="pearson", n_bootstraps=1000){
  # run many bootstraps
  # output: - vector with correlations 
  #         - vector with selected features as integers values wrt to X
  #           chronologically added in each bootstrap sample
  #         - vector with SEM
  cor_vec <- rep(NA, n_bootstraps)
  inds_vec <- integer(length = 0)
  MSE_vec <- rep(NA, n_bootstraps)
  for (i in c(1:n_bootstraps)) {
    out <- boost_bootstrap_sample(fm, df_train, df_test, method="pearson")
    cor_vec[i] = as.numeric(out$cor)
    inds_vec <- c(inds_vec, as.integer(out$inds))
    MSE_vec[i] <- as.numeric(out$MSE)
    cat(i," ")
  }  
  return(list(cor_vec=cor_vec, inds_vec=inds_vec, MSE_vec=MSE_vec))
}

set.seed(123)
bb_object_t <- boost_boot(fm01, Proliferation_6genes, Proliferation_6genes, n_bootstraps=2)
bb_object <- boost_boot(fm01, Proliferation_6genes, Proliferation_6genes, n_bootstraps=100)
bb_object <- boost_boot(fm05, dfA03, dfA03, n_bootstraps=1000)


histogram(bb_object$cor_vec, breaks = 99,
          xlab = "Proliferation score", 
          main = "Boosting (trail 1, arm Letro+Ribo)")

# Various objects
save(bb_object, file="bb_object_6genes.RData")
save(bb_object, file="bb_object_nodes01.RData")
save(bb_object, file="bb_object_AllGenes01.RData")
save(bb_object, file="bb_object_ALLGenes_ROR.RData")

# stored object can be loaded to save time
load("bb_object_6Genes.RData")
load("bb_object_nodes01.RData")
load("bb_object_AllGenes01.RData")
load("bb_object_ALLGenes_ROR.RData")







##########################################################
## Earlier code

boost_m01 = glmboost(fm01, data = Proliferation_6genes)
coef(boost_m01, which = "")
par(mfrow=c(1,2))
plot(boost_m01, off2int=TRUE)
plot(boost_m01)

pred_boost01 = predict(boost_m01, Proliferation_6genes)

cor(pred_boost01, Proliferation_6genes$Y)
plot(pred_boost01, Proliferation_6genes$Y)
abline(lm(Proliferation_6genes$Y ~ pred_boost01))

cvm <- cvrisk(boost_m01)
cvm
mstop(cvm)

# smooth - P-spline as base learner
spline01 <- gamboost(fm01, data = Proliferation_6genes,
                     baselearner = "bbs", # dfbase=4
                     control = boost_control(mstop = 40))
spline01$xselect()

coef(spline01, which = "")
par(mfrow=c(2,3))
plot(spline01, off2int=TRUE)

cols_chosen = unique(spline01$xselect())
colnames(Proliferation_6genes)[cols_chosen]


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
colnames(Proliferation_6genes)[cols_chosen+2]

par(mfrow=c(2,3))
plot(stump01)


pred_boost01 = predict(boost_m01, Proliferation_6genes)

cor(pred_boost01, Proliferation_6genes$Y)
plot(pred_boost01, Proliferation_6genes$Y)
abline(lm(Proliferation_6genes$Y ~ pred_boost01))
