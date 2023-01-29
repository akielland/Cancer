## Repeated cross-validation for testing the Lasso model
##
## Test against the fold
## Cross-validation (5-fold) used for training/tuning lambda and selecting features
## output: selected genes; proliferation.score correlation; SEM

## Post Lasso model is made at the bottom

library(glmnet)
library(caret)
library(ggplot2)

###########################################################
##  Repeated cross validation ##
###########################################################

# Custom function for repeated k-fold cross validation
repeated_kfold_cv <- function(data, func=lasso_sample, folds=5, repeats=100) {
  
  correlations <- numeric(repeats * folds)
  print(length(correlations))
  MSE_vec <- rep(NA, repeats * folds)
  coef_matrix <- matrix(NA, nrow = repeats*folds, ncol = ncol(data[,-1]))
  colnames(coef_matrix) <- colnames(data[, -1])
  coef_matrix_row_index <- 1
  
  # Repeat the cross-validation process
  for (i in 1:repeats) {
    # Create the folds for evaluating the performance
    kf <- caret::createFolds(data[,1], k = folds, list = TRUE, returnTrain = TRUE)
    # Loop through the folds
    for (j in 1:folds) {
      # Get the training and testing data
      train_data <- data[kf[[j]],]
      test_data <- data[-kf[[j]],]
      
      # Fit the function on the training data and get results
      fit <- func(train_data)
      
      coef_matrix[coef_matrix_row_index, ] <- coef(fit, s = "lambda.min")[-1]
      
      pred = predict(fit, newx = as.matrix(test_data)[,-1], type = "response", s = "lambda.min")  
      correlations[coef_matrix_row_index]  <- suppressWarnings(cor(pred, test_data[,1]))
      MSE_vec[coef_matrix_row_index] <- mean((pred - test_data[,1])^2)
      
      cat(coef_matrix_row_index, "")
      coef_matrix_row_index <- coef_matrix_row_index + 1
    }
  }
  return(list(correlations=correlations, coef_matrix=coef_matrix, MSE_vec=MSE_vec))
}

# Set repeats and folds of the cross-validations
repeats = 100
folds = 5
# RUN
lasso_k_ob  <- repeated_kfold_cv(Proliferation_ALLgenes, lasso_sample, folds, repeats)

head(lasso_k_ob$coef_matrix)[,1:10]
mean(na.omit(lasso_k_ob$correlations))

# Order the features based on their selection frequency
frequency <- data.frame(Feature = colnames(lasso_k_ob$coef_matrix), Frequency = colSums(lasso_k_ob$coef_matrix != 0) / (repeats * folds))
frequency <- frequency[order(frequency$Frequency, decreasing = TRUE),]
frequency[1:3, 1:2]

# Bar plot of the selection frequency of the features
ggplot(frequency[1:50,], aes(x = Frequency, y = reorder(Feature, Frequency))) +
  geom_bar(stat = "identity") +
  xlab("Selection Frequency") +
  ylab("Features") +
  ggtitle("Selection Frequency of Features") +
  theme(axis.text.y = element_text(angle = 0, hjust = 0))

# Histogram over correlations
correlations_finite <- lasso_k_ob$correlations[is.finite(lasso_k_ob$correlations)]
MSE_df <- data.frame(correlation = correlations_finite)

ggplot(corr_df, aes(x=correlation)) +
  geom_histogram(bins = 30, color = "black", fill = "white") +
  xlab("Correlation") +
  ylab("Frequency") +
  ggtitle("Histogram of Correlation Values")

# Histogram over MSE
MSE_finite <- lasso_k_ob$MSE_vec[is.finite(lasso_k_ob$MSE_vec)]
corr_df <- data.frame(correlation = MSE_finite)

ggplot(corr_df, aes(x=MSE_vec)) +
  geom_histogram(bins = 30, color = "black", fill = "white") +
  xlab("Correlation") +
  ylab("Frequency") +
  ggtitle("Histogram of Correlation Values")

## Extracting best features
# calculate the number of features to keep
num_features_to_keep <- round(0.1 * ncol(lasso_k_ob$coef_matrix))

# count the frequency of each feature in coef_matrix
counts <- colSums(lasso_k_ob$coef_matrix != 0)

# sort the features based on their frequency
sorted_features <- names(sort(counts, decreasing = TRUE))


# extract the top 10% features
top_features_with_index = which(colSums(lasso_k_ob$coef_matrix != 0) >= 0.1 * repeats)
top_feature_names = colnames(lasso_k_ob$coef_matrix)[top_features_with_index]

# use these top features for post lasso
response <- "Y"
formula_string <- paste(response, "~", paste(top_features, collapse = " + "))
formula <- as.formula(formula_string)

post_lasso_df <- cbind(Proliferation_ALLgenes[, 1], Proliferation_ALLgenes[, top_features_with_index + 1])
colnames(post_lasso_df)[1] <- "Y"

##############################
# older code


histogram(lb_object$cor_vec, breaks = 99,
          xlab = "Proliferation score", 
          main = "Lasso, (trail 1, arm Letro+Ribo)")


# Various objects
save(lb_object, file="lb_object_6genes.RData")
save(lb_object, file="lb_object_nodes01.RData")
save(lb_object, file="lb_object_AllGenes01.RData")
save(lb_object, file="lb_object_ALLGenes_ROR.RData")

# stored object can be loaded to save time
load("lb_object_6Genes.RData")
load("lb_object_nodes01.RData")
load("lb_object_AllGenes01.RData")
load("lb_object_ALLGenes_ROR.RData")

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





