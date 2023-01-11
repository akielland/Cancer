## This code uses the boot library to perform bootstrap with 1000 iterations, 
## where each iteration performs linear regression on a random subset of the data 
## and calculates the mean squared error. 
## The bootstrap results, including the mean MSE and 95% confidence interval, 
## are then printed using the boot.ci function.

library(boot)

# Function to fit linear regression model and calculate MSE
linear_regression_mse <- function(fm, data, indices) {
  d <- data[indices,]
  fit <- lm(fm, data=d)
  predict <- predict(fit, newdata = data)
  mse <- mean((data[, "Y"] - predict)^2)
  # print(mse)
  return(predict)
}
print(data1[, "Y"])
# Data
data_in <- Proliferation_6genes
# Formulas
fm01 <- Y ~ CCND1 + CCNE1 + CDKN1A + ESR1 + MYC + RB1
fm02 <- sur_ProliferationScore ~ scr_CCND1*scr_RB1 + scr_CCND1*scr_CDKN1A + scr_CCND1*scr_ESR1 + 
  scr_CCNE1*scr_CDKN1A + scr_CCNE1*scr_RB1 + scr_MYC*scr_ESR1 + scr_MYC*scr_CDKN1A 
fm03 <- sur_ProliferationScore ~ scr_CDKN1A + scr_ESR1 + scr_MYC + 
  scr_ESR1:scr_MYC + scr_CDKN1A:scr_MYC 

# Perform bootstrap with 1000 iterations
b <- boot(data=data_in, statistic=linear_regression_mse, R = 5, fm=fm01)

b$t

# Print the mean MSE and 95% confidence interval
print(boot.ci(b, type = "bca"))




fit <- lm(fm01, data=data1)
dim(data1)
date
prediction <- predict(fit, newdata = data1)





mse <- mean((data$ProliferatonScore - predict)^2)
#print(data$ProliferationScore)
#print(mse)


