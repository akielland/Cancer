## Lasso using bootstrap sample to train 1000 models
## Cross-validation (5-fold) used for training/tuning lambda and selecting features
## Test against original data sample and leave-out bootstrap sample
## output: selected genes; proliferation.score correlation; SEM

## Post Lasso model is made at the bottom

library(glmnet)

###########################################################
##  Bootstrap  Lasso - using original sample as test set ##
###########################################################


lasso_bootstrap_sample <- function(X, Y, ext.X, ext.Y, lambda.min=TRUE){
  # one bootstrap sample
  # output: 
  # 1. correlation between prediction and full input sample
  # 2. indices of selected variables
  n = length(Y)
  
  int <- sample.int(n, size = n, replace = TRUE)
  X_train = X[int,]
  Y_train = Y[int]
  
  fit.cv <- cv.glmnet(X_train, Y_train, nfolds = 5)
  
  features <- coef(fit.cv, s = "lambda.min")
  inds <- which(features != 0)
  # If want to print selected features
  # inds <- inds[-1] # dropping the first covariates 
  # variables <- row.names(co)[inds]
  # variables <- variables[!(variables %in% '(Intercept)')];
  
  if (lambda.min == TRUE) {
    pred_values = predict(fit.cv, newx = ext.X, type = "response", s = "lambda.min")      
  }
  if (lambda.min == FALSE) {
    pred_values = predict(fit.cv, newx = ext.X, type = "response", s = "lambda.1se")      
  }
  cor <- suppressWarnings(cor(pred_values, ext.Y, method = "spearman"))
  SEM <- mean((ext.Y - pred_values)^2)
  return(list(cor, inds, SEM))
}

lasso_bootstrap_sample(X, Y, X, Y, TRUE)


# 1000 bootstrap fits
lasso_cor_boot = function(X, Y, ext.X, ext.Y, n_bootstraps=100){
  # run many bootstraps
  # output: - vector with correlations 
  #         - vector with selected features as integers values wrt to X
  #           chronologically added in each bootstrap sample
  #         - vector with SEM
  cor_vec <- rep(NA, n_bootstraps)
  inds_vec <- integer(length = 0)
  SEM_vec <- rep(NA, n_bootstraps)
  for (i in c(1:n_bootstraps)) {
    out <- lasso_bootstrap_sample(X, Y, ext.X, ext.Y, TRUE)
    cor_vec[i] = as.numeric(out[1])
    inds_vec <- c(inds_vec, as.integer(out[[2]]))
    SEM_vec[i] <- as.numeric(out[3])
  }  
  return(list(cor_vec, inds_vec, SEM_vec))
  #return(list(cor_vec, inds_vec))
}

# making list object with 
# 1. correlation  
# 2. integer values of features selected in each bootstrap model
# 3. SEM
set.seed(123)
lb_object <- lasso_cor_boot(X_Z, Y, ext.X_Z, ext.Y = ext.PFS_months, 100)

lb_object <- lasso_cor_boot(X_Z, order(Y), ext.X_Z, ext.Y = order(ext.OS_months, decreasing=TRUE), 100)

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




################################################
## Bellow is code I don't use so much anymore ##
################################################

# count the number of times the selected covariates are present in the models
covariates_count_selected <- rowSums(outer(c(1:max(covariates_w_names)), covariates_w_names, "=="))
names(covariates_count_selected) <- c(1:max(covariates_w_names))
show(covariates_count_selected)



regression_bootstrap_test_genes <- function(df, formula, boot_fraction=0.632){
  N = nrow(df)
  
  int <- sample(N, size = N*boot_fraction, replace = TRUE)
  train <- df[int,]
  test <- df[-int,]
  
  lm <- lm(formula, data = train)
  pred_values <- predict(lm, newdata = test)
  out <- cor(pred_values, test$Proliferation.Score)
  return(out)
}

regression_cor_boot_test_genes = function(df, formula, boot_fraction=0.632, n_bootstraps){
  cor_vec <- rep(NA, n_bootstraps)
  for (i in c(1:n_bootstraps)) {
    cor_vec[i] <- regression_bootstrap_test_genes(df, formula)
  }  
  return(cor_vec)
}

regression_test_genes <- regression_cor_boot_test_genes(df04, fm_test_genes, boot_fraction=0.632, 1000)

mean(regression_test_genes, na.rm=TRUE)
var(regression_test_genes, na.rm=TRUE)
par(mfrow=c(1,1))
hist(regression_test_genes, breaks = 50, xlim = c(-0.5, 0.8))



#######

output_minus1 <- output[which(output!=1)]
hist(output_minus1)
test = rowSums(outer(c(1:max(output_minus1)), output_minus1, "=="))
names(test) <- c(1:length(test))
test



##############################################################
##  Bootstrap  Lasso - using out-of-boot sample as test set ##
##############################################################



lasso_bootstrap_sample <- function(X, Y, lambda.min=TRUE){
  # one bootstrap sample
  # output: 
  # 1. correlation between prediction and true of the leave-out of the bootstrap sample
  # 2. indices of selected variables
  X = as.matrix(X)
  Y = as.matrix(Y)
  N = length(Y)

  int <- sample(N, size = N*0.632, replace = FALSE)
  X_train = X[int,]
  X_test = X[-int,]
  Y_train = Y[int]
  Y_test = Y[-int]
  
  lasso.cv <- cv.glmnet(X_train, Y_train, nfolds = 5)
  
  covariates <- coef(lasso.cv, s = "lambda.min")
  inds <- which(covariates != 0)
  # inds <- inds[-1] # dropping the first covariates 
  # variables <- row.names(co)[inds]
  # variables <- variables[!(variables %in% '(Intercept)')];
    
  if (lambda.min == TRUE) {
    pred_values = predict(lasso.cv, newx = X_test, type = "response", s = "lambda.min")      
  }
  if (lambda.min == FALSE) {
    pred_values = predict(lasso.cv, newx = X_test, type = "response", s = "lambda.1se")      
  }
  cor <- suppressWarnings(cor(pred_values, Y_test))
  return(list(cor, inds))
  }

# lasso_bootstrap_sample(X, Y, TRUE)

lasso_cor_boot = function(X, Y, n_bootstraps){
  # run many bootstraps
  # output: - vector with correlations 
  #         - vector with selected covariates as integers values wrt to X
  #           chronologically add in each bootstrap sample
  cor_vec <- rep(NA, n_bootstraps)
  inds_vec <- integer(length = 0)
  
  for (i in c(1:n_bootstraps)) {
    out <- lasso_bootstrap_sample(X, Y, TRUE)
    cor_vec[i] = as.numeric(out[1])
    inds_vec <- c(inds_vec, as.integer(out[[2]]))
  }  
  return(list(cor_vec, inds_vec))
}

# making list object with correlation and integer values of covariates 
# from the different lasso bootstrap models

lb_object <- lasso_cor_boot(X,Y,1000)
lb_object <- lasso_cor_boot(X,order(Y),100)

save(lb_object, file="lb_object_AllGenes01.RData")
load("lb_object_AllGenes01.RData")
save(lb_object, file="lb_object_nodes01.RData")
load("lb_object_nodes01.RData")


# Summery result of lasso bootstrap object
cor_vec <- as.numeric(lb_object[[1]])
sum(is.na(cor_vec))/length(cor_vec)    # fraction of NA
mean(cor_vec, na.rm=TRUE)
var(cor_vec, na.rm=TRUE)
par(mfrow=c(1,1))
hist(cor_vec, breaks = 100)
cor(df08$proliferation, df08$ProliferationScore)

# Analyzing selected covariates (genes/nodes...) in lb_object
# count the presence of the individual covariates for all bootstrap models
covariates_n <- 771  # IF genes
covariates_n <- 8    # IF nodes
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
covariates_w_names[order(covariates_w_names, decreasing = TRUE)]

# count the number of times the selected covariates are present in the models
covariates_count_selected <- rowSums(outer(c(1:max(covariates_w_names)), covariates_w_names, "=="))
names(covariates_count_selected) <- c(1:max(covariates_w_names))
show(covariates_count_selected)

# extract covariates selected at certain amount of times
genes_of_interest <- function(vector, times_selected, above=TRUE){
  # return: names of genes selected a distinct numbers of times
  if (max(vector) < times_selected){
    return(print("argument times_selected has to high value"))
  }
  if (above == TRUE){inds <- which(vector >= times_selected)}
  else{inds <- which(vector == times_selected)}
  
  variable_names <- names(inds)
  return(variable_names)
}
  
test_genes = genes_of_interest(covariates_w_names, 2, TRUE)
show(test_genes)

fm_test_genes = as.formula(paste("Proliferation.Score", "~", paste(test_genes, collapse = "+")))
lm(fm_test_genes, data = df04)

regression_bootstrap_test_genes <- function(df, formula, boot_fraction=0.632){
  N = nrow(df)
  
  int <- sample(N, size = N*boot_fraction, replace = TRUE)
  train <- df[int,]
  test <- df[-int,]
  
  lm <- lm(formula, data = train)
  pred_values <- predict(lm, newdata = test)
  out <- cor(pred_values, test$Proliferation.Score)
  return(out)
}

regression_cor_boot_test_genes = function(df, formula, boot_fraction=0.632, n_bootstraps){
  cor_vec <- rep(NA, n_bootstraps)
  for (i in c(1:n_bootstraps)) {
    cor_vec[i] <- regression_bootstrap_test_genes(df, formula)
  }  
  return(cor_vec)
}

regression_test_genes <- regression_cor_boot_test_genes(df04, fm_test_genes, boot_fraction=0.632, 1000)

mean(regression_test_genes, na.rm=TRUE)
var(regression_test_genes, na.rm=TRUE)
par(mfrow=c(1,1))
hist(regression_test_genes, breaks = 50, xlim = c(-0.5, 0.8))


#########

output_minus1 <- output[which(output!=1)]
hist(output_minus1)
test = rowSums(outer(c(1:max(output_minus1)), output_minus1, "=="))
names(test) <- c(1:length(test))
test



