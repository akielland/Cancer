## Boot library is used to perform bootstrap with 1000 iterations 
## A linear regression model fitted on the bootstrap sample (a random subset of the data) in each iteration
## Mean squared error is calculated between the prediction of the models and the original data-set in each iteration
## Output: mean MSE; 95% confidence interval; histogram/plot

class(data_in$Y)
data1[, "Y"]

library(boot)

# Data: predictor set and output (The output is coded as Y, independent of what it is, from the data clean file)
data_in <- Proliferation_6genes

# Formulas: linear regressoin
fm00 <- CCND1 ~ CCND1 + CCNE1
fm01 <- Y ~ CCND1 + CCNE1 + CDKN1A + ESR1 + MYC + RB1
fm02 <- Y ~ CCND1*RB1 + CCND1*CDKN1A + CCND1*ESR1 + CCNE1*CDKN1A + CCNE1*RB1 + MYC*ESR1 + MYC*CDKN1A
fm03 <- Y ~ CDKN1A + ESR1 + MYC + ESR1:MYC + CDKN1A:MYC


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
b <- boot(data=data_in, statistic=linear_regression_mse, R = 1000, formula=fm01)
b <- boot(data=data_in, statistic=linear_regression_mse, R = 1000, formula=fm02)
b <- boot(data=data_in, statistic=linear_regression_mse, R = 1000, formula=fm03)

# get mean and 95% confidence interval
mean(b$t)
sd(b$t)
# print(boot.ci(b, type = "bca"))    # type = "bca" is throwing an error
boot.ci(b, type="basic")
ci.type <- c("norm","basic", "stud", "perc")
boot.ci(b, type = ci.type)

# plot/histograms of mse for each bootstrap sample
plot(b)
hist(b$t,breaks = 50)


###############################
## Repeated cross-validation ##
###############################


library(caret)

# set up repeated cross-validation
repeats <- 10
folds <- 5
rcv <- trainControl(method = "repeatedcv", number = folds, repeats = repeats, savePredictions = TRUE)

# train and evaluate your model
results <- train(formula=fm00, data = data_in, method = "lm", trControl = rcv)
results <- train(y = data_in$Y, x = data_in$CCNE1, method = "lm", trControl = rcv)

data_in$Y
data_in$CCNE1

# print the results
print(results)


# load the library
library(caret)
# load the iris dataset
data(iris)
# define training control
train_control <- trainControl(method="repeatedcv", number=10, repeats=3)
# train the model
model <- train(y = data_in[,Y], x = data_in[,CCNE1], trControl=train_control, method="lm")
# summarize results
print(model)

train_control <- trainControl(method="repeatedcv", number=10)
# train the model
model <- train(y = data_in[,Y], x = data_in[,CCNE1], trControl=train_control, method="lm")
# summarize results
print

train.control <- trainControl(method = "repeatedcv", 
                              number = 5, repeats = 100)
# Train the model
model <- train(Y ~ CCND1 + CCNE1 + CDKN1A + ESR1 + MYC + RB1, 
               data = data_in, 
               method = "lm",
               trControl = train.control)
# Summarize the results
print(model)
dt <- model$resample

dt %>% 
  group_by(Resample) %>% 
  summarise(avg = mean(RMSE))

-ends_with("timepoint")




