############
## ROC curve and Kaplan-Meier plots based on the cut-off value
## - ROR
## 

# Used to estimate cut-off for prediction for PFS 
# PFS are censored data
# Two types
# - PFS_months
# - PFS_status


library(pROC)
library(dplyr)

# install.packages("survival")
# install.packages("survival", dependencies = TRUE)
# install.packages("survminer")
# install.packages("ggsignif")
# install.packages("rstatix")
# install.packages("tidyselect", dependencies = TRUE)

library(survival)
library(survminer)

install.packages("survival", dependencies = TRUE)
library(survival)


# ridge base function returning a trained object
ridge_model <- function(train_data, alpha=0){
  print(dim(train_data))
  X_ <- as.matrix(dplyr::select(train_data, -1))
  Y_ <- as.matrix(dplyr::select(train_data, 1))
  
  fit.cv <- cv.glmnet(X_, Y_, nfolds = 5, alpha=0)
  return(fit.cv)
}

# fit ridge on all data and using 5-fold CV to select the best hyper-parameters
ridge_model_all_trail_data <- ridge_model(ROR_prolif_771genes)

# predict outcome of the trail data using model made on trail data
pred_ROR_trail = predict(ridge_model_all_trail_data, newx = as.matrix(ROR_prolif_771genes)[,-1], type = "response", s = "lambda.min") 


# extract the calssification (high=1, medium/low=0) from the trail data
# y_labels = df.ROR$ROR_P_GroupSubtypeProliferation
y_labels = df.trail_ROR_subtype$ROR_P_GroupSubtypeProliferation

# Using observed ROR score to set probability values
# y_pred = df.ROR$Y/100

# Using predicted ROR score to set probability values
y_pred = as.numeric(pred_ROR_trail)/100

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

# Cut-off values using the prediction on trail data

# 34 0.2626304   0.7894737   0.7272727
# 35 0.2706087   0.8157895   0.7272727
# 36 0.2766501   0.8421053   0.7272727
# 37 0.2903148   0.8421053   0.6363636

# 38 0.2860146   0.8684211   0.6363636
# 39 0.2874070   0.8947368   0.6363636
# 40 0.2884976   0.8947368   0.5454545
# 41 0.2909630   0.9210526   0.5454545


# The code below first divides the data into two groups (high and low) based on the cut-off value for the pred_ROR_cohort variable. 
# Then, it creates two separate Kaplan-Meier plots: 
# - Progression-Free Survival (PFS) and another 
# - Overall Survival (OS). 
# Each plot contains two survival curves, one for the high group and one for the low group, 
# This allows to visually compare the survival probabilities between these two groups over time


ROR_cut_off <- 30

# predict on cohort data
pred_ROR_cohort = predict(ridge_model_all_trail_data, newx = as.matrix(ext.X_Z), type = "response", s = "lambda.min") 

df_pred_PFS <- df04 |> dplyr::select(PFS_months, PFS_status, OS_months, OS_status)
# add the new column to the data frame with the desired name
df_pred_PFS <- df_pred_PFS  |> mutate(pred_ROR_cohort = pred_ROR_cohort[,1])
# df_pred_PFS  |> mutate(test = pred_ROR_cohort)


# Create a new column with the groups based on the cut-off value
df_pred_PFS$group <- ifelse(df_pred_PFS$pred_ROR_cohort >  ROR_cut_off, "high", "low")

# Create the survival objects
PFS_surv <- Surv(df_pred_PFS$PFS_months, df_pred_PFS$PFS_status)
OS_surv <- Surv(df_pred_PFS$OS_months, df_pred_PFS$OS_status)

# Fit the Kaplan-Meier models
PFS_fit <- survfit(PFS_surv ~ group, data = df_pred_PFS)
OS_fit <- survfit(OS_surv ~ group, data = df_pred_PFS)

# Plot the Kaplan-Meier curves
ggsurvplot(PFS_fit,
           data = df_pred_PFS,
           title = "Kaplan-Meier Plot for PFS",
           xlab = "Time (months)",
           ylab = "Progression-Free Survival",
           legend.title = "Group",
           legend.labs = c("High pred_ROR_cohort", "Low pred_ROR_cohort"),
           pval = TRUE,
           pval.method = TRUE)

ggsurvplot(OS_fit,
           data = df_pred_PFS,
           title = "Kaplan-Meier Plot for OS",
           xlab = "Time (months)",
           ylab = "Overall Survival",
           legend.title = "Group",
           legend.labs = c("High pred_ROR_cohort", "Low pred_ROR_cohort"),
           pval = TRUE,
           pval.method = TRUE)

##########
# Plot the Kaplan-Meier curves for PFS
PFS_plot <- ggsurvplot(PFS_fit,
                       data = df_pred_PFS,
                       title = "Kaplan-Meier Plot for PFS",
                       xlab = "Time (months)",
                       ylab = "Progression-Free Survival",
                       legend.title = "Group",
                       pval = TRUE,
                       pval.method = TRUE)

# Modify legend labels for PFS
PFS_plot$plot <- PFS_plot$plot + scale_color_manual(values = c("red", "blue"), labels = c("High pred_ROR_cohort", "Low pred_ROR_cohort"))

# Display PFS plot
print(PFS_plot)

# Plot the Kaplan-Meier curves for OS
OS_plot <- ggsurvplot(OS_fit,
                      data = df_pred_PFS,
                      title = "Kaplan-Meier Plot for OS",
                      xlab = "Time (months)",
                      ylab = "Overall Survival",
                      legend.title = "Group",
                      pval = TRUE,
                      pval.method = TRUE)

# Modify legend labels for OS
OS_plot$plot <- OS_plot$plot + scale_color_manual(values = c("red", "blue"), labels = c("High pred_ROR_cohort", "Low pred_ROR_cohort"))

# Display OS plot
print(OS_plot)




# 0.7 is best
# 0.28 is best (between 0.275 and 0.300 ?)
