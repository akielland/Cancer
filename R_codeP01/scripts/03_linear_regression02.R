## Linear regression with the 6 genes used in he mechanistic model is tested
## Interaction is tested for the gene products which are known to interact

# Bootstrap and Repeated Cross-validation is used to evaluate performance
## Output: mean MSE; 95% confidence interval; histogram/plot

##########################################################################

###############
## Bootstrap ##
###############

## Boot library is used to perform bootstrap with 1000 iterations 
## Fit a linear regression model on the bootstrap sample (a random subset of the data) in each iteration
## Mean squared error is calculated between the prediction of the models and the original data-set in each iteration

library(boot)

# Data: predictor set and output (The output is coded as Y, independent of what it is, from the data clean file)

# Formulas 
# Linear regression:
fm00 <- CCND1 ~ CCND1 + CCNE1
fm01 <- Y ~ CCND1 + CCNE1 + CDKN1A + ESR1 + MYC + RB1
fm02 <- Y ~ CCND1*RB1 + CCND1*CDKN1A + CCND1*ESR1 + CCNE1*CDKN1A + CCNE1*RB1 + MYC*ESR1 + MYC*CDKN1A
fm03 <- Y ~ CDKN1A + ESR1 + MYC + ESR1:MYC + CDKN1A:MYC

# Nodes of mechanistic model:
fm04 <- Y ~ cyclinD1 + cyclinD1Palbo + p21 + cyclinD1p21 + cMyc + cyclinEp21 + Rb1 + ppRb1


# Function to fit linear regression model and calculate MSE on all the data
linear_regression_mse <- function(data, indices, formula) {
  d <- data[indices,]  # allows boot to select sample
  fit <- lm(formula, data = d)
  y_predict <- predict(fit, newdata = data)
  mse <- mean((data$Y - y_predict)^2)
  return(mse)
}

# Make the results reproducible
set.seed(1234)
# Perform bootstrap with 1000 iterations
b <- boot(data=Proliferation_6genes, statistic=linear_regression_mse, R = 1000, formula=fm01)
b <- boot(data=Proliferation_6genes, statistic=linear_regression_mse, R = 1000, formula=fm02)
b <- boot(data=Proliferation_6genes, statistic=linear_regression_mse, R = 1000, formula=fm03)

b <- boot(data=Nodes_Proliferation, statistic=linear_regression_mse, R = 1000, formula=fm04)


# get mean and 95% confidence interval
mean(b$t)
sd(b$t)
# print(boot.ci(b, type = "bca"))    # type = "bca" is throwing an error
boot.ci(b, type="basic")
ci.type <- c("norm","basic", "stud", "perc")
boot.ci(b, type = ci.type)

# plot/histograms of mse for each bootstrap sample
plot(b)
hist(b$t, breaks = 50)


###############################
## Repeated cross-validation ##
###############################

## caret library is used to perform repeated cross_validation 
## 100 k-fold cross-validation is performed (using random splitt)
## Fit a linear regression model on the each k-1 fold k times
## Mean squared error is calculated between the prediction of the models and the k fold (mean of k times)


library(caret)

# set up repeated cross-validation
rep <- 100
folds <- 5
train_control <- trainControl(method = "repeatedcv", 
                              number = folds, repeats = rep, 
                              savePredictions = TRUE)

# train model
rcv <- train(fm01, 
                 data = data_in, 
                 method = "lm", 
                 trControl = train_control)

print(rcv)

# Get RMSE from each k in each iteration
df <- rcv$resample # this get the results (df: resample) of each iteration from the rcv object

# Average RMSE over each iteration
df_RMSE <- df |> 
  mutate(Resample = str_sub(Resample, start = 7L)) |> 
  group_by(Resample) |> 
  summarise(avg = mean(RMSE))

# plot MSE
hist((df_RMSE$avg)**2, breaks = 25)




