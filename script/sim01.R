# Simulating gene expression data
X.sim = ma02[,2:7]
X.sim.df = as.data.frame(X.sim)
y.sd = sd(ma02[,1])


sim01 <- function(X){
  beta = c(0,1,0,1,0,1)
  Y = X %*% beta
  return(Y)
}
Y = sim01(X.sim) 

plot(X, Y)

sim.lm = lm(Y ~ CCND1+CCNE1+CDKN1A+ESR1+MYC+RB1, data = X.sim.df)

plot(Y, sim.lm)

y = predict(sim.lm, newdata = X.sim.df) + rnorm(length(Y), mean = 0, sd = y.sd)

plot(y, Y)

sim.lm.noisy = lm(y ~ CCND1+CCNE1+CDKN1A+ESR1+MYC+RB1, data = X.sim.df)
summary(sim.lm.noisy)


## Lasso prediction
prediction.lasso <- function(X, Y, holdout=5){
  X = as.matrix(X)
  Y = as.matrix(Y)
  n = length(Y)
  sample = sample.int(n=n, size= n - holdout)
  X_train = X[sample,]
  X_test = X[-sample,]
  Y_train = Y[sample,]
  Y_test = Y[-sample,]
  
  lasso.cv <- cv.glmnet(X_train, Y_train)
  #lasso.pred = predict(lasso, newx = as.matrix(X_test), type = "response", s = "lambda.1se")
  lasso.pred = predict(lasso.cv, newx = X_test, type = "response", s = "lambda.min")
  out = data.frame(lasso.pred, Y_test)
  return(out)
}

set.seed(111)

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


## Stump prediction
prediction.stump <- function(df, holdout=5){
  n = length(df[,1])
  sample = sample.int(n=n, size= n - holdout)
  df_train = df[sample,]
  df_test = df[-sample,]
  
  stump = gamboost(y ~ CCND1+CCNE1+CDKN1A+ESR1+MYC+RB1, data = df_train, 
                   family=Gaussian(), 
                   baselearner='btree',
                   boost_control(mstop = 100))
  
  stump.cv = cvrisk(stump)
  stump.pred = predict(stump, newdata = df_test, type = "response")
  
  out = data.frame(stump.pred, df_test[,1])
  return(out)
}

set.seed(111)

pred.stump.loop = function(df){
  df_ <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(df_) <- c("stump.pred", "ProliferationScore")
  
  for (i in c(1:10)) {
    df_ = rbind(df_, prediction.stump(df))
  }  
  return(df_)
}

df.y = as.data.frame(y)
df = cbind(df.y, X.sim.df)
pip.test.stump <- pred.stump.loop(df)


## Results: Stump prediction

plot(pip.test.stump[,1], pip.test.stump[,2])
abline(lm(pip.test.stump[,2] ~ pip.test.stump[,1]), col = "red", lwd=3)
mtext(paste("Correlation:", round(cor(pip.test.stump[,2], pip.test.stump[,1]), 3)), side=1, line=-1, cex=1.4)

mean((pip.test.stump[,1] - pip.test.stump[,2])^2)






# PREDICTION

X  <- subset(df05, select= -c(ProliferationScore, CCNE1, CDKN1A, MYC, RB1))
X  <- subset(df05, select= -c(ProliferationScore))
Y  <- subset(df05, select= ProliferationScore)

par(mfrow=c(2,3), mar=rep(4,4))

for(i in 1:(ncol(X))) {
  x = as.matrix(X[,i])
  y = as.matrix(Y)
  plot(x, y, 
       xlab=colnames(df05)[i+1],
       ylab=colnames(df05)[1],
       pch=20, col='blue',
       main=paste('Plot number',i))
      abline(lm(y ~ x), col = "red", lwd=3)
      mtext(paste("Correlation:", round(cor(x, y), 2)), side=1, line=-1.5, cex=1.2)
      # cat(names[i])
      # print(cor(x,y), method = c("pearson"))
}

## Lasso prediction
prediction.lasso <- function(X, Y, i){
  X = as.matrix(X)
  Y = as.matrix(Y)
  sample = sample.int(n=48, size=45)
  X_train = X[-i,]
  X_test = X[i,]
  Y_train = Y[-i,]
  Y_test = Y[i,]
  
  lasso.cv <- cv.glmnet(X_train, Y_train)
  lasso.pred = predict(lasso.cv, newx = X_test, type = "response", s = "lambda.1se")
  #lasso.pred = predict(lasso.cv, newx = X_test, type = "response", s = "lambda.min")
  out = data.frame(lasso.pred, Y_test)
  #print(coef(lasso.cv, s = "lambda.min"))
  return(out)
}

set.seed(111)

pred.lasso.loop = function(X, Y){
  df_ <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(df_) <- c("lambda.min", "Y_test")
  
  for (i in c(1:48)) {
    df_ = rbind(df_, prediction.lasso(X,Y,i))
  }  
  return(df_)
}

pip.test.l = pred.lasso.loop(X,Y)


## Results: Lasso prediction
par(mfrow=c(1,1))
plot(pip.test.l[,1], pip.test.l[,2])
abline(lm(pip.test.l[,2] ~ pip.test.l[,1]), col = "red", lwd=3)
mtext(paste("Correlation:", round(cor(pip.test.l[,2], pip.test.l[,1]), 3)), side=1, line=-1, cex=1.4)

mean((pip.test.l[,1] - pip.test.l[,2])^2)


## Boost/linear prediction
prediction.boost <- function(X, Y, e){
  df_ = cbind(Y, X)
  sample = sample.int(n=48, size=45)
  df_train = df_[-e,]
  df_test = df_[e,]
  
  boost = gamboost(fm01, data = df_train, 
                   family=Gaussian(), 
                   # baselearner='btree',
                   baselearner="bols",
                   boost_control(mstop = 100))
  
  # stump.cv = cvrisk(boost)
  boost.cv = cvrisk(boost, folds = cv(model.weights(boost), type = "kfold"), papply = mclapply)
  
  boost.pred = predict(boost, newdata = df_test, type = "response")
  
  out = data.frame(boost.pred, df_test[,1])
  return(out)
}

set.seed(111)

pred.boost.loop = function(X, Y){
  df_ <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(df_) <- c("stump.pred", "ProliferationScore")
  
  for (e in c(1:48)) {
    cat(e, ", ")
    df_ = rbind(df_, prediction.boost(X, Y, e))

  }  
  return(df_)
}

pip.test.b <- pred.boost.loop(X, Y)

# Results: Stump prediction
plot(pip.test.b[,1], pip.test.b[,2])
abline(lm(pip.test.b[,2] ~ pip.test.b[,1]), col = "red", lwd=3)
mtext(paste("Correlation:", round(cor(pip.test.b[,2], pip.test.b[,1]), 3)), side=1, line=-1, cex=1.4)

mean((pip.test.b[,1] - pip.test.b[,2])^2)






