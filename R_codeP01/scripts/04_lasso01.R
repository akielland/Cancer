# standard lasso
lasso01 <- glmnet(as.matrix(X), as.matrix(Y))
plot(lasso01, label = T)
coef(lasso01, s = 1)

# lasso w/ cross-validation to determine optimal lambda hyperparameter
lasso02 = cv.glmnet(as.matrix(X), as.matrix(Y), nfolds=5)
co <- coef(lasso02, s = "lambda.min")
inds <- which(co!=0)
variables <- row.names(co)[inds]
variables <- variables[!(variables %in% '(Intercept)')];

# Leave one observation out for testing and going through all observations
lasso_loo <- function(X, Y, lambda.min=TRUE){
  X = as.matrix(X)
  Y = as.matrix(Y)
  n = length(Y)
  pred_values <- NULL
  for(i in 1:nrow(X)){
    X_test <- X[i,]
    X_train <- X[-i,]
    Y_train <- Y[-i,]
    
    lasso.cv <- cv.glmnet(X_train, Y_train)

    if (lambda.min == TRUE) {
      pred_values[i] = predict(lasso.cv, newx = X_test, type = "response", s = "lambda.min")      
    }
    if (lambda.min == FALSE) {
      pred_values[i] = predict(lasso.cv, newx = X_test, type = "response", s = "lambda.1se")      
    }
  }
  return(pred_values)
}

# Result of lasso with lambda at min:
pred_values_lasso.min <- lasso_loo(X, Y, TRUE)
par(mfrow=c(1,2))
cat("Correlation for lasso using lambda.min: ")
print(cor(pred_values_lasso.min, df04$Proliferation.Score))
cat("MSE for lasso using lambda.min: ")
print(mean((pred_values_lasso.min - df04$Proliferation.Score)**2))
plot(pred_values_lasso.min, df04$Proliferation.Score, main = "Lasso")
abline(lm(df04$Proliferation.Score ~ pred_values_lasso.min), pred_values_lasso.min)

# Result of lasso with lambda at 1 se from min:
pred_values_lasso.1se <- lasso_loo(X, Y, FALSE)
cat("Correlation for lasso using lambda.min: ")
print(cor(pred_values_lasso.1se, PS))
cat("MSE for lasso using lambda.min: ")
print(mean((pred_values_lasso.1se - PS)**2))
plot(pred_values_lasso.1se, PS, main = "Lasso")

#################
##  bootstrap  ##
#################

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
  inds <- inds[-1]
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
  # output: vector with correlations
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

# Summery result of lasso bootstrap object
cor_vec <- as.numeric(lb_object[[1]])
sum(is.na(cor_vec))/length(cor_vec)    # fraction of NA
mean(cor_vec, na.rm=TRUE)
var(cor_vec, na.rm=TRUE)
par(mfrow=c(1,1))
hist(cor_vec, breaks = 50)

# Analyzing selected covariates/genes in lb_object
# count the presence of the individual covariates for all bootstrap models
vector_1 <- c(1:771)
vector_2 <- lb_object[[2]]
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







output_minus1 <- output[which(output!=1)]
hist(output_minus1)
test = rowSums(outer(c(1:max(output_minus1)), output_minus1, "=="))
names(test) <- c(1:length(test))
test



