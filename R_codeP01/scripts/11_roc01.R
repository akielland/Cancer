############
## ROC curve
## - ROR
## 

# Used to estimate cut-off for prediction for PFS 
# PFS are sensor data
# Two types
# - PFS_months
# - PFS_status


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

#plot(ROC_obj, print.auc=T)

# Specify SENSITIVITY criteria to meet.
Sn.upper <- 1.0
Sn.lower <- 0.5

# Specify SPECIFICITY criteria to meet.
Sp.upper <- 1.0
Sp.lower <- 0.5

# Extract all coordinate values from the ROC curve.
my.coords <- coords(roc=ROC_obj, x = "all", transpose = FALSE)

# Identify and print all points on the ROC curve that meet the JOINT sensitivity AND specificity criteria.
my.coords[(my.coords$specificity >= Sp.lower & my.coords$specificity <= Sp.upper & 
             my.coords$sensitivity >= Sn.lower & my.coords$sensitivity <= Sn.upper),]
# 0.7 is best
# 0.28 is best
