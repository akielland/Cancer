# install.packages("glmnet", repos = "https://cran.us.r-project.org")
library(glmnet)

# Importing data
df01 <- readxl::read_xlsx('/Users/anders/Documents/Master/data/SCR_SUR_6genes_noNAN.xlsx')

# Make two dataframes from Letro+Ribo at timepoint SCR and SUR
dfLR <- subset(df01, TrialArmNeo %in% "Letro+Ribo")
dfLR.SCR <- dfLR[dfLR$timepoint %in% "SCR", ]
dfLR.SUR <- dfLR[dfLR$timepoint %in% "SUR", ]

# Make matrix with predictors at timepoint SCR and SUR (mRNA level, also one matrix w/prolif at SCR) and 
# vectors with the response variables (proliferation at SCR and SUR)
df.X01.SCR <- subset(dfLR.SCR, select=-c(UniqueID, timepoint, TrialArmNeo, ProliferationScore))
df.X01.SUR <- subset(dfLR.SUR, select=-c(UniqueID, timepoint, TrialArmNeo, ProliferationScore))
df.X02.SCR <- subset(dfLR.SCR, select=-c(UniqueID, timepoint, TrialArmNeo))

X01.SCR <- as.matrix(df.X01.SCR)
X01.SUR <- as.matrix(df.X01.SUR)
X02.SCR <- as.matrix(df.X02.SCR)
Y.SCR <- as.matrix(dfLR.SCR["ProliferationScore"])
Y.SUR <- as.matrix(dfLR.SUR["ProliferationScore"])

# plot at SCR
par(mfrow=c(2,3), mar=rep(2,4))
names = c('CCND1', 'CCNE1', 'CDKN1A', 'ESR1', 'MYC','RB1') 
for(i in 1:6){
  x = X01.SCR[,i]
  y = Y.SCR
  plot(x, y, main=names[i])
  abline(lm(y ~ x), col = "red", lwd=3)
  #text(paste("Correlation:", round(cor(x, y), 2)), x = -1.5, y = 0.84)
  mtext(paste("Correlation:", round(cor(x, y), 2)), side=1, line=-1.5, cex=1.2)
  print(cor(x,y), method = c("pearson"))
}

mtext("(c)",side=3,line=-1.5, 
      at=par("usr")[1]+0.05*diff(par("usr")[1:2]),
      cex=1.2)
# plot at SUR
for(i in 1:6){
  x = X01.SUR[,i]
  y = Y.SUR
  plot(x, y, main=names[i])
  abline(lm(y ~ x), col = "red", lwd=3)
  mtext(paste("Correlation:", round(cor(x, y), 2)), side=1, line=-1.5, cex=1.2)
  print(cor(x,y), method = c("pearson"))
  }
# plot at genes at SCR against prolif at SUR
for(i in 1:6){
  x = X01.SCR[,i]
  y = Y.SUR
  plot(x, y, main=names[i])
  abline(lm(y ~ x), col = "red", lwd=3)
  mtext(paste("Correlation:", round(cor(x, y), 2)), side=1, line=-1.5, cex=1.2)
  print(cor(x,y), method = c("pearson"))
  }
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1))

# Correlation: proliferation at SCR and SUR
plot(Y.SCR, Y.SUR)
abline(lm(Y.SUR ~ Y.SCR), col='red', lwd=3)
text(paste("Correlation:", round(cor(Y.SCR, Y.SUR), 2)), x = 0.8, y = 0.8)


# Making models
par(mfrow=c(1,2))
m01 <- glmnet(X01.SCR, Y.SUR)
plot(m01, label = T)
print(m01)
coef(m01, s = 0.02)

m02 <- glmnet(X02.SCR, Y.SUR)
plot(m02, label = T)
coef(m02, s = 0.02)

m03 <- glmnet(X01.SCR, Y.SCR)
coef(m03, s = 0.02)
plot(m03, label = T)

m04 <- glmnet(X01.SUR, Y.SUR)
plot(m04, label = T)



# Making models using cross validation
m01.cv <- cv.glmnet(X01.SCR, Y.SUR)
plot(m01.cv)

m02.cv$lambda.min
coef(m02.cv, s = "lambda.1se")

m02.cv <- cv.glmnet(X01.SCR, Y.SCR)
plot(m02.cv)
coef(m02.cv, s = "lambda.min")

m03.cv <- cv.glmnet(X01.SUR, Y.SUR)
plot(m03.cv)

# PREDICTION

# Split in training and test set
pred <- function(X, Y){
  sample = sample(48, 35)
  X_train = X[sample,]
  X_test = X[-sample,]
  Y_train = Y[sample]
  Y_test = Y[-sample]

  cv.m01 <- cv.glmnet(X_train, Y_train)
  plot(cv.m01)
  
  pred.1se = predict(cv.m01, newx = X_test, s = "lambda.1se")
  pred.min = predict(cv.m01, newx = X_test, s = "lambda.min")
  print(pred.1se)
  print(coef(cv.m01, s="lambda.1se"))
  
  par(mfrow=c(1,2))
  plot(Y_test, pred.1se)
  abline(lm(pred.1se ~ Y_test))
  
  plot(Y_test, pred.min)
  abline(lm(pred.min ~ Y_test))
  par(mfrow=c(1,1))
  
  print(mean((pred.1se - Y_test)^2))
  print(mean((pred.min - Y_test)^2))
  cor(Y_test, pred.1se)
  cor(Y_test, pred.min)
}
set.seed(111)
pred(X02.SCR, Y.SUR)



#set.seed(NULL)

> samples <- matrix(NA, nrow=10000, ncol=10)
> for(i in 1:10000) samples[i,] <- sample(infants$weight, 10)
> muhat <- apply(samples, 1, mean)
> hist(muhat,
       breaks=seq(from=2, to=4, by=0.1),
       xlab=expression("Estimated mean birth weight, " * hat(mu)),
       ylab="Relative frequency",
       main="",
       col="blue",
       freq=FALSE)


