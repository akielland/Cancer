##########################################################################
## Repeated cross-validation for testing the Stacking - Ensemble models 
## Here also with interaction terms
##########################################################################
##
## Test against the fold
## Cross-validation (5-fold) used for training/tuning lambda and selecting features
## output: selected genes (); proliferation.score correlation; SEM

# function to calculate new feature set from level 0 learners
level_0_model <- function(named_list, train_data, alpha0){
  pred_L0 <- data.frame(matrix(nrow = nrow(train_data), ncol = 0))
  fits <- list()
  
  for (key in names(named_list)) {
    char_matrix <- named_list[[key]]
    X_ <- as.matrix(train_data[c(char_matrix)])
    fit_ <- cv.glmnet(X_, train_data$Y, alpha=alpha0, nfolds=5)
    # here prediction is made using same data as used in learning the model... ?
    pred_ = predict(fit_, newx = X_, type = "response", s = "lambda.min")
    colnames(pred_) <- key
    pred_L0 <- cbind(pred_L0, pred_)
    fits[[key]] <- fit_
  }
  return(list(pred_L0=pred_L0, fits=fits))
}
# t <-level_0_model(char_list, prolif_771genes, )
# class(t)

# function to calculate predictions of the level 0 models on some test data
# The only different from the one above is that it takes fitted model as argument instead of fitting them in the function
level_0_test <- function(named_list, test_data, fits){
  pred_test_L0 <- data.frame(matrix(nrow = nrow(test_data), ncol = 0))
  
  for (key in names(named_list)) {
    char_matrix <- named_list[[key]]
    X_ <- as.matrix(test_data[c(char_matrix)])
    fit_ <- fits[[key]]
    test_pred = predict(fit_, newx = X_, type = "response", s = "lambda.min")
    colnames(test_pred) <- key
    pred_test_L0[key] <- test_pred
  }
  return(pred_test_L0)
}
# t1 <- level_0_test(char_list, prolif_771genes, t$fit)


# Function: repeated k-fold cross validation
stacking_rep_cv <- function(df_data, char_list, alpha0=0, alpha1=0.5, folds=5, repeats=1, interactions=FALSE, method="pearson") {
  n_models <- repeats * folds
  print(n_models)
  
  if (interactions){
    binom_coef <- choose(length(char_list), 2)
    coef_matrix <- matrix(NA, nrow = n_models, ncol = length(char_list) + binom_coef)
  } else {
    coef_matrix <- matrix(NA, nrow = n_models, ncol = length(char_list))
  }
  cor_vec <- rep(NA, n_models)
  MSE_vec <- rep(NA, n_models)
  
  count <- 1
  # Repeat the cross-validation process
  for (i in 1:repeats) {
    # Create the folds for evaluating the performance
    kf_1 <- caret::createFolds(df_data[,1], k = folds, list = TRUE, returnTrain = TRUE)
    # Loop through the folds
    for (j in 1:folds) {
      # Get the training and testing data
      train_data <- df_data[kf_1[[j]],]
      test_data <- df_data[-kf_1[[j]],]

      # # Split into training and prediction set
      # ind <- sample.int(2, nrow(train_data), replace=TRUE, prob=c(0.8, 0.2))
      # train0 <- train_data[ind==1,]
      # pred0 <- train_data[ind==2,]
      
      # Partition feature space; run lasso and stack in-sample predictions
      pred_and_fit <- level_0_model(char_list, train_data, alpha0) # returns a dictionary with predictions and models
      pred_L0 <- pred_and_fit$pred_L0
      
      # here make interaction matrix
      if (interactions){
        input_matrix <- model.matrix(~ .^2 - 1, data = pred_L0)
      } else {
        input_matrix = pred_L0
      }
      
      fits <- pred_and_fit$fits   # collect the fitted models for each signature gene sets 
      # fit lasso meta-model
      # meta_model <- cv.glmnet(as.matrix(pred_L0), train_data$Y, nfolds=5)
      meta_model <- cv.glmnet(as.matrix(input_matrix), train_data$Y, nfolds=5, alpha=alpha1)
   
      ####################################################################
      ## test part
      ####################################################################
      # test data prediction from level 0 
      # level_0_test returns the predictions from the the model fits of the individual sign.gene.sets
      pred_test_L0 <- level_0_test(char_list, test_data, fits)
      # here make interaction matrix
      if (interactions){
        input_matrix <- model.matrix(~ .^2 - 1, data = pred_test_L0)
      } else {
        input_matrix = pred_test_L0
      }
      # predict meta_model with test-data
      test_pred <- predict(meta_model, newx = as.matrix(input_matrix), type = "response", s = "lambda.min")
     
      coef_matrix[count, ] <- coef(meta_model, s = "lambda.min")[-1]
      cor_vec[count]  <- suppressWarnings(cor(test_pred, test_data[,1], method=method))
      MSE_vec[count] <- mean((test_pred - test_data[,1])^2)
    
      if (!count %% 10 ) cat(count, "")
      count <- count + 1
    }
  }
  colnames(coef_matrix) <- sub(":", "*", colnames(input_matrix)) # replace ':' with '*' in column names
  
  return(list(cor_vec=cor_vec, coef_matrix=coef_matrix, MSE_vec=MSE_vec))
}

set.seed(1)
t3 <- stacking_rep_cv(ROR_prolif_771genes, char_list, alpha0=0, alpha1=0.5, folds=5, repeats=20, interactions=FALSE, method="pearson")
t3 <- stacking_rep_cv(ROR_prolif_771genes, char_list, alpha0=0, alpha1=0.5, folds=5, repeats=20, interactions=TRUE, method="pearson")
mean(t3$cor_vec)

# Funtion to calculate selection freq
freq_non_zero_non_na <- function(matrix) {
  col_freq <- colSums(!is.na(matrix) & matrix != 0) / nrow(matrix) * 100
  result <- data.frame(column_name = names(col_freq), percentage = col_freq)
  result <- data.frame(percentage = col_freq)
  return(result)
}


# Ridge
# RUN: r_stack_c_771_RORprolif_6
set.seed(123)
r_stack_c_771_RORprolif_6 <- stacking_rep_cv(ROR_prolif_771genes, char_list_6, alpha0=0, alpha1=0, folds=5, repeats=200, interactions=FALSE, method="pearson")
head(r_stack_c_771_RORprolif$coef_matrix)
save(r_stack_c_771_RORprolif, file="/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/r_stack_c_771_RORprolif.RData")
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/r_stack_c_771_RORprolif.RData")
mean(r_stack_c_771_RORprolif$cor_vec, na.rm=T)
sd(r_stack_c_771_RORprolif$cor_vec)
freq_non_zero_non_na(r_stack_c_771_RORprolif$coef_matrix)

mean(r_stack_c_771_RORprolif_6$cor_vec, na.rm=T)
r_stack_c_771_RORprolif_6$cor_vec



# RUN: r_stack_771_RORprolif_5
set.seed(123)
r_stack_771_RORprolif_5 <- stacking_rep_cv(ROR_prolif_771genes,char_list_5, alpha0=0, alpha1=0, folds=5, repeats=200, interactions=FALSE, method="pearson")
head(r_stack_771_RORprolif_5$coef_matrix)
save(r_stack_771_RORprolif_5, file="/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/r_stack_771_RORprolif_5.RData")
mean(r_stack_771_RORprolif_5$cor_vec, na.rm=T)
sd(r_stack_771_RORprolif_5$cor_vec, na.rm=T)

# RUN: r_stack_771_RORprolif_10
set.seed(123)
r_stack_771_RORprolif_10 <- stacking_rep_cv(ROR_prolif_771genes,char_list_10, alpha0=0, alpha1=0, folds=5, repeats=200, interactions=FALSE, method="pearson")
head(r_stack_771_RORprolif_10$coef_matrix)
save(r_stack_771_RORprolif_10, file="/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/r_stack_771_RORprolif_10.RData")
mean(r_stack_771_RORprolif_10$cor_vec, na.rm=T)
sd(r_stack_771_RORprolif_10$cor_vec, na.rm=T)

# Lasso
# RUN: l_stack_c_771_RORprolif
set.seed(123)
l_stack_c_771_RORprolif <- stacking_rep_cv(ROR_prolif_771genes,char_list, alpha0=0, alpha1=1, folds=5, repeats=200, interactions=FALSE, method="pearson")
head(l_stack_c_771_RORprolif$coef_matrix)
save(l_stack_c_771_RORprolif, file="/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/l_stack_c_771_RORprolif.RData")
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/l_stack_c_771_RORprolif.RData")
mean(l_stack_c_771_RORprolif$cor_vec, na.rm=T)
freq_non_zero_non_na(l_stack_c_771_RORprolif$coef_matrix)

# RUN: l_stack_c_771_RORprolif_6
set.seed(123)
l_stack_c_771_RORprolif_6 <- stacking_rep_cv(ROR_prolif_771genes, char_list_6, alpha0=0, alpha1=1, folds=5, repeats=200, interactions=FALSE, method="pearson")
head(l_stack_c_771_RORprolif_6$coef_matrix)
save(l_stack_c_771_RORprolif_6, file="/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/l_stack_c_771_RORprolif_6.RData")
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/l_stack_c_771_RORprolif_6.RData")
mean(l_stack_c_771_RORprolif_6$cor_vec, na.rm=T)
sd(l_stack_c_771_RORprolif_6$cor_vec, na.rm=T)
freq_non_zero_non_na(l_stack_c_771_RORprolif_6$coef_matrix)


# RUN: l_stack_c_771_RORprolif_10
set.seed(123)
l_stack_c_771_RORprolif_10 <- stacking_rep_cv(ROR_prolif_771genes, char_list_10, alpha0=0, alpha1=1, folds=5, repeats=200, interactions=FALSE, method="pearson")
head(l_stack_c_771_RORprolif_10$coef_matrix)
save(l_stack_c_771_RORprolif_10, file="/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/l_stack_c_771_RORprolif_10.RData")
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/l_stack_c_771_RORprolif_10.RData")
mean(l_stack_c_771_RORprolif_10$cor_vec, na.rm=T)
sd(l_stack_c_771_RORprolif_10$cor_vec, na.rm=T)
freq_non_zero_non_na(l_stack_c_771_RORprolif_10$coef_matrix)







# RUN: l_stack_c_interact_771_RORprolif
set.seed(123)
l_stack_c_interact_771_RORprolif <- stacking_rep_cv(ROR_prolif_771genes,char_list, alpha0=0, alpha1=1, folds=5, repeats=200, interactions=TRUE, method="pearson")
head(l_stack_c_interact_771_RORprolif$coef_matrix)
save(l_stack_c_interact_771_RORprolif, file="/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/l_stack_c_interact_771_RORprolif.RData")
mean(l_stack_c_interact_771_RORprolif$cor_vec, na.rm=T)

# Elastic
# RUN: en_stack_c_771_RORprolif
set.seed(123)
en_stack_c_771_RORprolif <- stacking_rep_cv(ROR_prolif_771genes, char_list, alpha0=0, alpha1=0.5, folds=5, repeats=200, interactions=FALSE, method="pearson")
head(en_stack_c_771_RORprolif$coef_matrix)
save(en_stack_c_771_RORprolif, file="/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/en_stack_c_771_RORprolif.RData")
mean(en_stack_c_771_RORprolif$cor_vec, na.rm=T)

# RUN: en_stack_c_interact_771_RORprolif
set.seed(123)
en_stack_c_interact_771_RORprolif <- stacking_rep_cv(ROR_prolif_771genes,char_list, alpha0=0, alpha1=0.5, folds=5, repeats=200, interactions=TRUE, method="pearson")
head(en_stack_c_interact_771_RORprolif$coef_matrix)
save(en_stack_c_interact_771_RORprolif, file="/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/en_stack_c_interact_771_RORprolif.RData")
mean(en_stack_c_interact_771_RORprolif$cor_vec, na.rm=T)





