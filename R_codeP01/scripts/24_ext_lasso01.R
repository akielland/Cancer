###################################
## Testing on the second trail   ##
## Here: lasso                   ##
####################################

# make lasso fit w/cross-validation to determine optimal lambda hyper-parameter
# using all obs from 1. trail
# Y: ROR with Prolif (count variable, range: 0 - 100)
lasso.fit = cv.glmnet(X, Y, nfolds=5)
co <- coef(lasso.fit, s = "lambda.min")
inds <- which(co != 0)

variables <- row.names(co)[inds]
variables <- variables[!(variables %in% '(Intercept)')]
show(variables)
length(variables)

# using the lasso.fit with lambda at min SEM to predict ROR from all data
pred_ROR <- predict(lasso.fit, newx = X, type = "response", s = "lambda.min")
cor(pred_ROR, Y)
summary(Y)
summary(pred_ROR)

## 2.trail
# using the lasso.fit with lambda at min SEM to predict ROR from the 2.trail
ext.pred_ROR <- predict(lasso.fit, newx = ext.X, type = "response", s = "lambda.min")
summary(ext.pred_ROR)
hist(ext.pred_ROR)

cor(ext.pred_ROR, ext.PFS_months)
summary(ext.PFS_months)
cor(ext.pred_ROR, ext.PFS_status)
summary(ext.PFS_status)



# ROR and PFS
# PFS are sensor data


library(pROC)

# extract data
y_labels = df.ROR$ROR_P_GroupSubtypeProliferation

# Using observed ROR score to set probability values
y_pred = df.ROR$Y/100
# Using predicted ROR score to set probability values
y_pred = as.numeric(pred_ROR)/100

# Compute the ROC object and show the curve
ROC_obj <- roc(y_labels, y_pred,
                smoothed = TRUE,
                # arguments for ci
                ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                # arguments for plot
                plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                print.auc=TRUE, show.thres=TRUE)
plot(ROC_obj, print.auc=T)

# Specify SENSITIVITY criteria to meet.
Sn.upper <- 1.0
Sn.lower <- 0.5

# Specify SPECIFICITY criteria to meet.
Sp.upper <- 1.0
Sp.lower <- 0.6

# Extract all coordinate values from the ROC curve.
my.coords <- coords(roc=ROC_obj, x = "all", transpose = FALSE)

# Identify and print all points on the ROC curve that meet the JOINT sensitivity AND specificity criteria.
my.coords[(my.coords$specificity >= Sp.lower & my.coords$specificity <= Sp.upper & 
             my.coords$sensitivity >= Sn.lower & my.coords$sensitivity <= Sn.upper),]
# 0.7 is best
# 0.28 is best
