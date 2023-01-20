###################################
## Testing on the second trail   ##
## Here: lasso                   ##
####################################

# make lasso fit w/cross-validation to determine optimal lambda hyper-parameter
# using all obs from 1. trail to fit model using ROR as dependent variable
# test outcome:
# - PFS_status
# - PFS_months


# Y: Prolif (ca range (-1, 1))
lasso.fit = cv.glmnet(X, Y, nfolds=5)
co <- coef(lasso.fit, s = "lambda.min")
inds <- which(co != 0)
variables <- row.names(co)[inds]
variables <- variables[!(variables %in% '(Intercept)')]
show(variables)
length(variables)


# Y: ROR w/Prolif (count variable, range: 0 - 100)
set.seed(123)
set.seed(0)
lasso.fit = cv.glmnet(X, Y, nfolds=5)
# standardizing features
lasso.fit = cv.glmnet(X_Z, Y, nfolds=5)
# standardizing features and response
lasso.fit = cv.glmnet(X_Z, Y, nfolds=5, standardize.response = TRUE)
# setting setting Y to range (0, 1)
lasso.fit = cv.glmnet(X, Y/100, nfolds=5)

# Get the features
selected_features <- function(){
  co <- coef(lasso.fit, s = "lambda.min")
  inds <- which(co != 0)
  variables <- row.names(co)[inds]
  variables <- variables[!(variables %in% '(Intercept)')]
  print(length(variables))
  return(variables)
}
selected_features()

# using the lasso.fit with lambda at min.SEM to predict ROR from all orignal data
# If want to standardize set: newx = X_Z
pred_ROR <- predict(lasso.fit, newx = X, type = "response", s = "lambda.min")
cor(pred_ROR, Y)
summary(Y)
summary(pred_ROR)

## 2.trail
# using the lasso.fit with lambda at min SEM to predict ROR from the 2.trail
ext.pred_ROR <- predict(lasso.fit, newx = ext.X, type = "response", s = "lambda.min")
# Standardize features
ext.pred_ROR <- predict(lasso.fit, newx = ext.X_Z, type = "response", s = "lambda.min")
summary(ext.pred_ROR)
hist(ext.pred_ROR)

## Correlations
plot(ext.pred_ROR, unlist(ext.PFS_months))
# Pearson
cor(ext.pred_ROR, ext.PFS_months, method = "pearson")
summary(ext.PFS_months)
cor(ext.pred_ROR, ext.PFS_status)
summary(ext.PFS_status)

# Spearman
cor(ext.pred_ROR, ext.PFS_months, method = "spearman")
# p-values
cor.test(ext.pred_ROR, unlist(ext.PFS_months), method = "spearman")
cor.test(as.numeric(as.vector(ext.pred_ROR)), unlist(ext.PFS_months)/100, method = "spearman")




