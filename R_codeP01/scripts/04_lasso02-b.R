###########################################################
##  Bootstrap Lasso - using original sample as test set ##
###########################################################

## Cross-validation (5-fold) used for training/tuning lambda and selecting features
## output: selected genes; correlation; SEM

## Post Lasso model is made at the bottom

library(glmnet)
library(caret)
library(ggplot2)

# Lasso base function returning a trained object
lasso_sample <- function(train_data){
  X_ <- as.matrix(train_data |> select(-1))
  Y_ <- as.matrix(train_data |>  select(1))
  # NB: NEED TO LOOK INTO the SETTINGS of Lasso
  # fit.cv <- cv.glmnet(X_, Y_, family = "gaussian", alpha = 1, standardize = TRUE, nlambda = 100, nfolds = 5)
  fit.cv <- cv.glmnet(X_, Y_, nfolds = 5)
  return(fit.cv)
}

# Bootstrap
lasso_boot <- function(df_train, df_test, func=lasso_sample, method="pearson", n_bootstraps=8) {
  print(n_bootstraps)
  n <- nrow(df_train)
  cor_vec <- rep(NA, n_bootstraps)
  MSE_vec <- rep(NA, n_bootstraps)
  coef_matrix <- matrix(NA, nrow = n_bootstraps, ncol = ncol(df_train)-1)
  colnames(coef_matrix) <- colnames(df_train[, -1])
  
  for (i in c(1:n_bootstraps)) {
    int <- sample.int(n, size = n, replace = TRUE)
    train_data = df_train[int,]
    
    # Fit the function on the training data and get results
    fit <- func(train_data)
    
    # Extract coefficients of feature
    coef_matrix[i, ] <- coef(fit, s = "lambda.min")[-1]
    
    pred = predict(fit, newx = as.matrix(df_test)[,-1], type = "response", s = "lambda.min") 
    cor_vec[i] <- suppressWarnings(cor(pred, df_test$Y, method = method))
    MSE_vec[i] <- mean((df_test$Y - pred)^2)
    # cor_vec[i]  <- suppressWarnings(cor(pred, df_test[,1], method = method))
    # MSE_vec[i] <- mean((pred - df_test[,1])^2)        
    cat(i, "")
  }
  return(list(cor_vec=cor_vec, coef_matrix=coef_matrix, MSE_vec=MSE_vec))
}

lasso_boot(prolif_6genes, prolif_6genes, func=lasso_sample, method="pearson", n_bootstraps=3)
  
# RUN: lb_obj_6_prolif
set.seed(123)
lb_obj_6_prolif  <- lasso_boot(prolif_6genes, prolif_6genes, lasso_sample, method="pearson", n_bootstraps=1000)
head(lb_obj_6_prolif$coef_matrix)[,1:6]
save(lb_obj_6_prolif, file="/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/lb_obj_6_prolif.RData")
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/lb_obj_6_prolif.RData")

# RUN: lb_obj_6_RORprolif
set.seed(123)
lb_obj_6_RORprolif  <- lasso_boot(RORprolif_6genes, RORprolif_6genes, lasso_sample, method="pearson", n_bootstraps=1000)
head(lb_obj_6_RORprolif$coef_matrix)[,1:6]
save(lb_obj_6_RORprolif, file="/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/lb_obj_6_RORprolif.RData")
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/lb_obj_6_prolif.RData")


# RUN: lb_obj_771_prolif
set.seed(123)
lb_obj_771_prolif  <- lasso_boot(prolif_771genes, prolif_771genes, lasso_sample, method="pearson", n_bootstraps=1000)
head(lb_obj_771_prolif$coef_matrix)[,1:8]
save(lb_obj_771_prolif, file="/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/lb_obj_771_prolif.RData")
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/lb_obj_771_prolif.RData")

# RUN: lb_obj_771_RORprolif
set.seed(123)
lb_obj_771_RORprolif  <- lasso_boot(ROR_prolif_771genes, ROR_prolif_771genes, lasso_sample, method="pearson", n_bootstraps=1000)
head(lb_obj_771_RORprolif$coef_matrix)[,1:8]
save(lb_obj_771_RORprolif, file="/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/lb_obj_771_RORprolif.RData")
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/lb_obj_771_RORprolif.RData")


# RUN: lb_obj_nodes_prolif
set.seed(123)
lb_obj_nodes_prolif  <- lasso_boot(prolif_nodes, prolif_nodes, lasso_sample, method="pearson", n_bootstraps=1000)
head(lb_obj_nodes_prolif$coef_matrix)[,1:8]
save(lb_obj_nodes_prolif, file="/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/lb_obj_nodes_prolif.RData")
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/lb_obj_nodes_prolif.RData")

# RUN: lb_obj_nodes_RORprolif
set.seed(123)
lb_obj_nodes_RORprolif  <- lasso_boot(ROR_prolif_nodes, ROR_prolif_nodes, lasso_sample, method="pearson", n_bootstraps=1000)
head(lb_obj_nodes_RORprolif$coef_matrix)[,1:8]
save(lb_obj_nodes_RORprolif, file="/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/lb_obj_nodes_RORprolif.RData")
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/lb_obj_nodes_RORprolif.RData")




## ANALYSIS
sum(is.na(lasso_k_ob$cor_vec))/length(lasso_k_ob$cor_vec)    # fraction of NA

## Correlations
mean(na.omit(lasso_k_ob$cor_vec))
mean(lasso_k_ob$cor_vec, na.rm=TRUE)
median(lasso_k_ob$cor_vec, na.rm=TRUE)
var(lasso_k_ob$cor_vec, na.rm=TRUE)
sd(na.omit(lasso_k_ob$cor_vec))

# Histogram
correlations_finite <- lasso_k_ob$cor_vec[is.finite(lasso_k_ob$cor_vec)]
cor_df <- data.frame(correlation = correlations_finite)

ggplot(cor_df, aes(x=correlation)) +
  geom_histogram(bins = 30, color = "black", fill = "white") +
  xlab("Correlation") +
  ylab("Frequency") +
  ggtitle("Histogram of Correlation Values")

## MSE
mean(lasso_k_ob$MSE_vec)
sd(lasso_k_ob$MSE_vec)
# Histogram
MSE_df <- data.frame(MSE = lasso_k_ob$MSE_vec)

ggplot(MSE_df, aes(x=MSE)) +
  geom_histogram(bins = 30, color = "black", fill = "white") +
  xlab("MSE") +
  ylab("Frequency") +
  ggtitle("Histogram of MSE Values")


## Most prevalent Features
# Order the features based on their selection frequency
frequency <- data.frame(Feature = colnames(lasso_k_ob$coef_matrix), Frequency = colSums(lasso_k_ob$coef_matrix != 0) / (repeats * folds))
frequency <- frequency[order(frequency$Frequency, decreasing = TRUE),]
frequency[1:5, 1:2]

# Bar plot of the selection frequency of the features
n_best <- 50
ggplot(frequency[1:n_best ,], aes(x = Frequency, y = reorder(Feature, Frequency))) +
  geom_bar(stat = "identity") +
  xlab("Selection Frequency") +
  ylab("Features") +
  ggtitle("Selection Frequency of Features") +
  theme(axis.text.y = element_text(angle = 0, hjust = 0))


## Extracting best features
# calculate the number of features to keep
perc_best <- 0.1
num_features_to_keep <- round(perc_best * ncol(lasso_k_ob$coef_matrix))
# count the frequency of each feature in coef_matrix
counts <- colSums(lasso_k_ob$coef_matrix != 0)
# sort the features based on their frequency
sorted_features <- names(sort(counts, decreasing = TRUE))
sorted_features[1:num_features_to_keep]

## Extracting best features for POST lasso
# extract the top features
perc_best <- 0.1 # the %/100 best fraction 
top_features_with_index = which(colSums(lasso_k_ob$coef_matrix != 0) >= perc_best * repeats)
top_feature_names = colnames(lasso_k_ob$coef_matrix)[top_features_with_index]
# make linear formula of top features
response <- "Y"
formula_string <- paste(response, "~", paste(top_features, collapse = " + "))
formula <- as.formula(formula_string)
# Make data.frame containing only top features
post_lasso_df <- cbind(Proliferation_ALLgenes[, 1], Proliferation_ALLgenes[, top_features_with_index + 1])
colnames(post_lasso_df)[1] <- "Y"




##############################
# older code


## Summery result of lasso bootstrap object
# Correlation
cor_vec <- as.numeric(lb_object[[1]])
sum(is.na(cor_vec))/length(cor_vec)    # fraction of NA
mean(cor_vec, na.rm=TRUE)
median(cor_vec, na.rm=TRUE)
var(cor_vec, na.rm=TRUE)
# par(mfrow=c(1,1))
hist(cor_vec, breaks = 100)

# SEM
SEM_vec <- as.numeric(lb_object[[3]])
mean(SEM_vec)
mean(sqrt(SEM_vec))
hist(SEM_vec, breaks=100)

## Analyzing selected features (genes/nodes...) in lb_object
# count the presence of the individual features for all bootstrap models
covariates_n <- 771  # If ALL genes
# covariates_n <- 6  # If 6 genes
# covariates_n <- 8    # If nodes
vector_1 <- c(1: covariates_n)
vector_2 <- lb_object[[2]] - 1 # shift numbers so intercept becomes 0 and first covariate is 1
covariates_count <- rowSums(outer(vector_1, vector_2, "=="))
# put gene names on the covariates_count vector
covariates_w_names <- setNames(covariates_count, colnames(X))
show(covariates_w_names)

# extract covariates selected at least once
covariates_w_names_ind <- which(covariates_w_names!=0)
covariates_w_names <- covariates_w_names[covariates_w_names_ind]
show(covariates_w_names)
max(covariates_w_names)
ind <- which(covariates_w_names == max(covariates_w_names))
covariates_w_names[ind]
features_ordered <- covariates_w_names[order(covariates_w_names, decreasing = TRUE)]
show(features_ordered)
features_ordered <- as.data.frame(features_ordered)

# extract features selected most times (e.g. > 100)
times_selected <- 100
features_ordered.topp_list <- features_ordered |> filter(features_ordered > times_selected )
cat("Numbers of genes selected and percentages of the 771 genes selected: ")
dim(features_ordered.topp_list)[1]
dim(features_ordered.topp_list)[1]/771 *100
features_ordered.topp_list <- rownames_to_column(features_ordered.topp_list)

# Bar diagram on features ordered by amount of times selected
ggplot(features_ordered.topp_list, aes(y = rowname, x = features_ordered)) +
  geom_col() +
  scale_y_discrete(limits = features_ordered.topp_list$rowname)





