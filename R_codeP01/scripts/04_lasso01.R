
# X and Y from the 6 gene data
X <- dplyr::select(df02_LR,  scr_CCND1, scr_CCNE1, scr_CDKN1A, scr_ESR1, scr_MYC, scr_RB1)
Y <- select(df02_LR, sur_ProliferationScore)

# X from the full data set
X <- select(df04, -"Proliferation.Score") |> 
  select(-"Ki67")

# Y's from the full data set
Y <- select(df04, "Proliferation.Score")
Y <- select(df04, "Ki67")
Y[is.na(Y)] = 0

lasso01 <- glmnet(as.matrix(X), as.matrix(Y))
plot(lasso01, label = T)
co <- coef(lasso01, s = 1)

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

pred_values_lasso.min <- lasso_loo(X, Y, TRUE)
pred_values_lasso.1se <- lasso_loo(X, Y, FALSE)

PS <- df04$Proliferation.Score

# Result of lasso with lambda at min:
par(mfrow=c(1,2))
cat("Correlation for lasso using lambda.min: ")
print(cor(pred_values_lasso.min, PS))
cat("MSE for lasso using lambda.min: ")
print(mean((pred_values_lasso.min - PS)**2))
plot(pred_values_lasso.min, PS, main = "Lasso")
abline(lm(PS ~ pred_values_lasso.min), pred_values_lasso.min)

# Result of lasso with lambda at 1 se from min:
cat("Correlation for lasso using lambda.min: ")
print(cor(pred_values_lasso.1se, PS))
cat("MSE for lasso using lambda.min: ")
print(mean((pred_values_lasso.1se - PS)**2))
plot(pred_values_lasso.1se, PS, main = "Lasso")

######################################################################
##  bootstrap  ##
#################

lasso_bootstrap_sample <- function(X, Y, lambda.min=TRUE){
  # one bootstrap sample
  # output: correlation between prediction and true of the leave-out of the bootstrap sample
  X = as.matrix(X)
  Y = as.matrix(Y)
  N = length(Y)

  int <- sample(N, size = N*0.632, replace = TRUE)
  X_train = X[int,]
  X_test = X[-int,]
  Y_train = Y[int]
  Y_test = Y[-int]
  
  lasso.cv <- cv.glmnet(X_train, Y_train, nfolds = 5)
  
  co <- coef(lasso.cv, s = "lambda.min")
  inds <- which(co!=0)
  variables <- row.names(co)[inds]
  variables <- variables[!(variables %in% '(Intercept)')];
    
  if (lambda.min == TRUE) {
    pred_values = predict(lasso.cv, newx = X_test, type = "response", s = "lambda.min")      
  }
  if (lambda.min == FALSE) {
    pred_values = predict(lasso.cv, newx = X_test, type = "response", s = "lambda.1se")      
  }
  cor <- suppressWarnings(cor(pred_values, Y_test))
  return(c(cor, variables))
}

lasso_bootstrap_sample(X, Y, TRUE)

lasso_cor_boot = function(X, Y, n_bootstraps){
  # run many boostraps
  # output: vector with correlations
  cor_vec <- rep(NA, n_bootstraps)
  #var_vec <- data_frame(NA, n_bootstraps)
  
  for (i in c(1:n_bootstraps)) {
    cor_var <- lasso_bootstrap_sample(X, Y, TRUE)
    cor_vec[i] = cor_var[1]
    #var_vec[i] = as.list(cor_var[-1])
  }  
  return(cor_vec)
}

cor_vec_boot <- lasso_cor_boot(X,Y,1000)
sum(is.na(cor_vec_boot))/1000
mean(cor_vec_boot, na.rm=TRUE)
var(cor_vec_boot, na.rm=TRUE)
par(mfrow=c(1,1))
hist(cor_vec_boot, breaks = 50)



