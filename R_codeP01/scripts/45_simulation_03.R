########################
## Last synergistic. ##
#######################

library(glmnet)

# Function that creates a synthetic dataset (matrix) with interaction effects
syntetic_01 <- function(n, interact, p_gr, p_non_gr, i, scale=TRUE){
  # Create predictors for groups
  X <- matrix(rnorm(n * p_gr), n, p_gr)
  colnames(X) <- paste0("X", 1:p_gr)
  # Add non-gr. predictors
  Z <- matrix(rnorm(n * p_non_gr), nrow = n)
  X <- cbind(X, Z)
  colnames(X)[(p_gr + 1):(p_gr + p_non_gr)] <- paste0("Z", 1:p_non_gr)
  # X_df <- as.data.frame(X)

  
  # Define character lists for each group of features
  char_group1 <- colnames(X)[1:5]
  char_group2 <- colnames(X)[6:10]
  char_group3 <- colnames(X)[11:15]
  char_group4 <- colnames(X)[16:20]
  
  # Create matrices w/ groups of predictors -> predictors also divided into groups
  X_group1 <- X[, char_group1]
  X_group2 <- X[, char_group2]
  X_group3 <- X[, char_group3]
  X_group4 <- X[, char_group4]
  
  char_list_sim <- list(gr1 = char_group1,
                        gr2 = char_group2,
                        gr3 = char_group3,
                        gr4 = char_group4)
  
  # Create true betas
  true_betas <- runif(p_gr, -1, 1)
  # Add non effective betas = 0
  betas <- c(true_betas, rep(0, p_non_gr))
  
  # divide betas into groups
  beta1 <- matrix(true_betas[1:5], ncol=1)
  beta2 <- matrix(true_betas[6:10], ncol=1)
  beta3 <- matrix(true_betas[11:15], ncol=1)
  beta4 <- matrix(true_betas[16:20], ncol=1)
  

  if (interact == 1){
    # Add interaction effects (1): gr1 x gr2.  -> 1 interaction term between the 4 groups: gives 6 simulations
    interaction_effect <- (X_group1 %*% beta1) * (X_group2 %*% beta2) * 1
  } else if (interact == 2) {
    # Add interaction effects (2): gr1 x gr2 + gr1 x gr4
    interaction_effect <- (X_group1 %*% beta1) * (X_group2 %*% beta2) * 1 + (X_group1 %*% beta1) * (X_group4 %*% beta4) * 1
  } else if (interact == 3){
    # Add interaction effects (3): gr1 x gr2 + gr1 x gr4 + gr2 x gr4
    interaction_effect <- (X_group1 %*% beta1) * (X_group2 %*% beta2) * 1 + (X_group1 %*% beta1) * (X_group4 %*% beta4) * 1 + (X_group2 %*% beta2) * (X_group4 %*% beta4) * 1
  }else {
    abort("Error: interact argument must be 1, 2, or 3.")
  }
    

  # use model to create response y
  y <- X %*% betas + interaction_effect + rnorm(n)
  colnames(y) <- "y"
  
  if (scale==TRUE)
    X <- scale(X)
    y <- scale(y)
  
  df <- as.data.frame(cbind(y,X))  
  return(list(df.sim=df, char_list=char_list_sim))
}



# Function: repeated k-fold cross validation: 
simulation_rc <- function(df_data, alpha, groups, folds=5, repeats=2, method="pearson") {
  n_models <- repeats * folds
  print(n_models)
  
  cor_vec <- rep(NA, n_models)
  MSE_vec <- rep(NA, n_models)
  # coef_matrix <- matrix(NA, nrow = n_models, ncol = ncol(df_data[,-1]))
  # colnames(coef_matrix) <- colnames(df_data[,-1])
  
  coef_matrix <- matrix(NA, nrow = n_models, ncol = 6)
  
  row_index <- 1
  
  # Repeat the cross-validation process
  for (i in 1:repeats) {
    # Create the folds for evaluating the performance
    kf <- caret::createFolds(df_data[,1], k = folds, list = TRUE, returnTrain = TRUE)
    # Loop through the folds
    for (j in 1:folds) {
      # Get the training and testing data
      train_data <- df_data[kf[[j]],]
      test_data <- df_data[-kf[[j]],]
      
      # # Fit the function on the training data and get results
      out <- synergistic_adaptive_deviance02(train_data, test_data[,-1], groups, alpha, adaptive=FALSE, lambda_seq = NULL, nfolds = 5, tol = 1e-6, max_iters = 100)
      
      model <- out$model
      # print(out$beta_interaction)
      # print(out$beta_main)
      # coef_matrix[row_index, ] <- out$coef[-1]   # remove intercept
      
      
      
      #pred = predict(model, newx = as.matrix(test_data)[,-1], type = "response", s = "lambda.min")  
      pred <- out$ypred
      cor_vec[row_index]  <- suppressWarnings(cor(pred, test_data[,1], method=method))
      MSE_vec[row_index] <- mean((pred - test_data[,1])^2)
      
      if(!row_index %% 10) cat(row_index, "")
      row_index <- row_index + 1
    }
    cor_mean <- mean(cor_vec)
    MSE_mean <- mean(MSE_vec)
  }
  return(list(cor_mean=cor_mean, MSE_mean=MSE_mean, coef_matrix=coef_matrix))
}






# Function to run many simulations
run_sim_02 <- function(n_simulations, interact, alpha, adaptive, tol){

  #coef_matrix <- matrix(NA, nrow = n_simulations, ncol = 6)
  correlation_mean <- rep(NA, n_simulations)
  MSE_mean <- rep(NA, n_simulations)
  
  row_index <- 1
  for (i in c(1:n_simulations)){
    # i = i+12345
    set.seed(i) 
    
    data <- syntetic_01(n, interact, p_gr, p_non_gr, i)
    
    results <- simulation_rc(df_data=data$df.sim, alpha, data$char_list)
    
    # print(out$beta_interaction)
    #coef_matrix[row_index, ] <- out$beta_interaction
    cat(row_index, "")
    row_index <- row_index + 1
    
    correlation_mean[i] <- results$cor_mean
    MSE_mean[i] <- results$MSE_mean
  }
  # colnames(coef_matrix) <- c("1x2", "1x3", "1x4", "2x3", "2x4", "3x4")
  return(list(cor_mean = correlation_mean, MSE_mean = MSE_mean))
}




n <- 50
n <- 100
n <- 500
n <- 1000
p_gr <- 20
p_non_gr <- 100
interact <- 1
interact <- 2
interact <- 3

n_simulations <- 50

t1 <- run_sim_02(n_simulations, interact, alpha = 1, adaptive = FALSE, tol=1e-4)
t2 <- run_sim_02(n_simulations, interact, alpha = 0.5, adaptive = FALSE, tol=1e-6)

mean(t2$cor_mean, na.rm=T)
sd(t2$cor_mean, na.rm=T)

mean(t2$MSE_mean, na.rm=T)
sd(t2$MSE_mean, na.rm=T)


# Count the number of non-zero and non-NA values in each column
count_non_zero_non_na <- apply(t1, 2, function(x) sum(!is.na(x) & x != 0))
print(count_non_zero_non_na*2)



save(t1, file="/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/syn_sim_rc_500n_interact3.RData")
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/syn_sim_rc_500n_interact3.RData")



synergistic_adaptive_deviance02 <- function(df, X_test, char_list, alpha = 0.5, adaptive=FALSE, lambda_seq = NULL, nfolds = 5, tol = 1e-4, max_iters = 100) {
  # Extract predictor matrix and response vector; assuming response variable in first column
  X <- as.matrix(df[, -1])
  X <- scale(X)
  X_test <- scale(X_test)
  y <- df[, 1]
  
  # Create a list of groups, each group containing the column name and indices of df corresponding to the feature names in char_list
  groups <- lapply(char_list, function(char_group) {
    sapply(char_group, function(feature_name) which(colnames(df) == feature_name) - 1)
  })
  
  # function caculating the interaction terms
  interaction_terms <- function(X, X_test, betas, group1, group2) {

    interaction_X <- (X[, group1, drop = FALSE] %*% betas[group1]) * (X[, group2, drop = FALSE] %*% betas[group2])
    interaction_X_test <- (X_test[, group1, drop = FALSE] %*% betas[group1]) * (X_test[, group2, drop = FALSE] %*% betas[group2])
    return(list(interaction_X=interaction_X, interaction_X_test=interaction_X_test))
  }
  
  # learn betas without interaction terms to get starting beta values, using ridge
  fit_beta0 <- cv.glmnet(X, y, alpha = 0, lambda = lambda_seq, nfolds = nfolds)
  best_lambda <- fit_beta0$lambda.min
  beta_main <- coef(fit_beta0, s = best_lambda)[-1]  # exclude intercept
  
  # calculate interaction terms for all combinations of groups
  n_groups <- length(groups)
  interaction_indices <- t(combn(n_groups, 2))
  
  # initialize interaction matrix
  n_combinations <- ncol(combn(n_groups, 2))
  interactions <- matrix(0, nrow(X), n_combinations)
  interactions_test <- matrix(0, nrow(X_test), n_combinations)
  
  # Calculate interaction terms using the current betas
  for (i in seq_len(nrow(interaction_indices))) {
    group1 <- groups[[interaction_indices[i, 1]]]
    group2 <- groups[[interaction_indices[i, 2]]]

    tmp <- interaction_terms(X, X_test, beta_main, group1, group2)
  
    interactions[, i] <- tmp$interaction_X
    interactions_test[, i] <- tmp$interaction_X_test
  }
  
  # make feature matrix with original features and interactions
  X_with_interactions <- cbind(X, interactions)
  X_test_with_interactions <- cbind(X_test, interactions_test)
  
  # initialize variables for convergence checking
  converged <- FALSE
  iter <- 0
  obj_diff <- NULL
  betas <- rep(1, ncol(X_with_interactions))
  dev_glmnet <- obj_dev <- Inf
  
  
  # Learn betas outside the interaction term until convergence by deviance
  while (!converged && iter < max_iters) {
    iter <- iter + 1
    if(!iter %% 10) cat("iter =", iter, "; dev_glmnet =", dev_glmnet, "; obj_dev =", obj_dev, "\n")
    # Calculate interaction terms using the current betas
    for (i in seq_len(nrow(interaction_indices))) {
      group1 <- groups[[interaction_indices[i, 1]]]
      group2 <- groups[[interaction_indices[i, 2]]]
      tmp <- interaction_terms(X, X_test, beta_main, group1, group2)
      
      interactions[, i] <- tmp$interaction_X
      interactions_test[, i] <- tmp$interaction_X_test
    }
    
    # Remove interactions terms = 0
    interactions_sd = apply(interactions, 2, sd)
    interactions_scaled = interactions[, interactions_sd != 0]
    interactions_test_scaled = interactions_test[, interactions_sd != 0]
    
    # Make feature matrix with original features and interactions
    X_with_interactions <- cbind(X, interactions_scaled)
    X_test_with_interactions <- cbind(X_test, interactions_test_scaled)
    
    # Fit model using adaptive elastic-net
    if (adaptive==TRUE) {
      pf <- rep(1, ncol(X_with_interactions))
      pf[1:ncol(X)] <- betas[1:ncol(X)] 
      pf[pf == 0] <- 1e-2
      pf <- 1 / pf
      if (!is.null(ncol(interactions_scaled)) && ncol(interactions_scaled) > 0) {
        pf[(ncol(X) + 1):ncol(X_with_interactions)] <- 0
      }
      fit_with_interactions <- cv.glmnet(X_with_interactions, y, alpha = alpha, lambda = lambda_seq, nfolds = nfolds, penalty.factor = pf)
      # Fit model with just elastiv net
    } else {
      fit_with_interactions <- cv.glmnet(X_with_interactions, y, alpha = alpha, lambda = lambda_seq, nfolds = nfolds)
    }
    best_lambda <- fit_with_interactions$lambda.min
    
    # Update betas using the new beta values
    new_betas <- as.vector(coef(fit_with_interactions))[-1]  # exclude intercept
    new_beta_main <- new_betas[1:ncol(X)]
    
    best_lambda_idx <- which(fit_with_interactions$lambda == best_lambda)
    
    # calculates the deviance for the current model using the best lambda value. 
    # The deviance is calculated as (1 - dev.ratio) * nulldev since we can only extract dev.ratio from the model object
    # dev.ratio is the deviance ratio at the best lambda: (null deviance - model deviance) / null deviance
    # nulldev is the null deviance (the deviance of a model with no features -> just the intercept).
    dev_glmnet_curr <-  (1 - fit_with_interactions$glmnet$dev.ratio[best_lambda_idx]) * fit_with_interactions$glmnet$nulldev
    obj_dev <- abs(dev_glmnet_curr - dev_glmnet)/dev_glmnet_curr
    if ( obj_dev < tol ) {  # Only check convergence for the betas without interaction terms
      converged <- TRUE
    } else {
      dev_glmnet <- dev_glmnet_curr
      betas <- new_betas  # Update betas only if not converged
      beta_main <- new_beta_main
    }
    
    obj_diff <- c(obj_diff, dev_glmnet)
  }
  
  beta_interaction <- new_betas[ncol(X) + 1:ncol(interactions)]

  ypred <- X_test_with_interactions %*% matrix(new_betas, ncol = 1) + as.vector(coef(fit_with_interactions))[1]
  
  # make named vectors for beta_main and beta_interaction
  named_beta_main <- setNames(beta_main, colnames(df)[-1])
  
  # make named vector for beta_interaction
  interaction_names <- apply(interaction_indices, 1, function(idx) {
    paste(names(char_list)[idx[1]], "*", names(char_list)[idx[2]])
  })
  named_beta_interaction <- setNames(beta_interaction, interaction_names)
  
  print("beta_interaction = ")
  print(beta_interaction)
  
  return(list(model=fit_with_interactions, beta_main = named_beta_main, beta_interaction=beta_interaction, ypred=ypred))
  # return(list(obj_diff=obj_diff, beta_main = named_beta_main, beta_interaction = named_beta_interaction, iterations = iter, best_lambda = best_lambda))
}

#result <- synergistic(ROR_prolif_771genes, char_list, alpha = 0.1, lambda_seq = NULL, nfolds = 5, tol = 1e-5, max_iters = 100)

t1 <- run_sim_02(n_simulations, interact, alpha = 1, adaptive = FALSE, tol=1e-4)









char_list_sim <- list(gr1 = char_group1,
                      gr2 = char_group2,
                      gr3 = char_group3,
                      gr4 = char_group4)




# Replace NA values with 0
mat_no_na <- replace(t1, is.na(t1), 0)
mat_no_na
# Calculate column means
colMeans(mat_no_na)
colSds(mat_no_na)

# Calculate the mean of non-NA values in each column
mean_non_na <- apply(t1, 2, function(x) mean(x, na.rm = TRUE))
print("Mean of non-NA values:")
print(mean_non_na)

# Calculate the proportion of NA values in each column
mean_na <- apply(t1, 2, function(x) mean(is.na(x)))
mean_na


print("Proportion of NA values:")
print(mean_na/n_simulations)

# Count the number of numeric values in each column
count_numeric <- apply(t2, 2, function(x) sum(sapply(x, is.numeric)))
print(count_numeric)
# Count the number of non-zero and non-NA values in each column
count_non_zero_non_na <- apply(t2, 2, function(x) sum(!is.na(x) & x != 0))
print(2*(count_non_zero_non_na))




# Run the elastic net interaction function
result <- syn_elastic_net_interaction(X_df, y, char_list, alpha = 0.01, lambda_seq = NULL, nfolds = 5, tol = 1e-8, max_iters = 100)
# Print the results
print(result)




## Bootstrapping a large data set in order to capture > 1 interactions
set.seed(123)
std_df <- scale(ROR_prolif_771genes)
dim(std_df)

boot_idx <- sample(1:dim(std_df)[1], 1000, replace = TRUE)

df_boot <- std_df[boot_idx, ]

result <- synergistic(df_boot, char_list, alpha = 0.1, lambda_seq = NULL, nfolds = 5, tol = 1e-6, max_iters = 100)




X_df <- matrix(rnorm(n * p), n, p)
colnames(X_df) <- paste0("X", 1:p)
X_df <- as.data.frame(X_df)
X_df <- cbind(X_df, matrix(rnorm(n*100), nrow = n))


# Create true betas
true_betas <- runif(p, -1, 1)
true_betas <- c(true_betas, rep(0,100))

beta1 <- matrix(true_betas[1:5], ncol=1)
beta2 <- matrix(true_betas[6:10], ncol=1)

beta4 <- matrix(true_betas[16:20], ncol=1)

# Define character lists for each group of features
char_group1 <- colnames(X_df)[1:5]
char_group2 <- colnames(X_df)[6:10]
char_group3 <- colnames(X_df)[11:15]
char_group4 <- colnames(X_df)[16:20]

X_df <- as.matrix(X_df)

X_group1 <- X_df[, char_group1]
X_group2 <- X_df[, char_group2]
X_group4 <- X_df[, char_group4]


# Add interaction effects

# 1 interaction term between the 4 groups: gives 6 simulations
interaction_effect <- (X_group1 %*% beta1) * (X_group2 %*% beta2) * 1 
interaction_effect <- (X_group1 %*% beta1) * (X_group3 %*% beta3) * 0.5 
interaction_effect <- (X_group1 %*% beta1) * (X_group4 %*% beta4) * 0.5 
interaction_effect <- (X_group2 %*% beta2) * (X_group3 %*% beta3) * 1 
interaction_effect <- (X_group2 %*% beta2) * (X_group4 %*% beta4) * 0.5 
interaction_effect <- (X_group3 %*% beta3) * (X_group4 %*% beta4) * 1 

# All possible interaction term between the 4 groups: gives 1 simulation
interaction_effect <- 
  (X_group1 %*% beta1) * (X_group2 %*% beta2) * 1 
+ (X_group1 %*% beta1) * (X_group3 %*% beta3) * 0.5 
+ (X_group1 %*% beta1) * (X_group4 %*% beta4) * 0.5 
+ (X_group2 %*% beta2) * (X_group3 %*% beta3) * 1 
+ (X_group2 %*% beta2) * (X_group4 %*% beta4) * 0.5 
+ (X_group3 %*% beta3) * (X_group4 %*% beta4) * 1 



X_mat <- as.matrix(X_df)
y <- X_mat %*% true_betas + interaction_effect + rnorm(n)

char_list <- list(group1 = char_group1,
                  group2 = char_group2,
                  group3 = char_group3,
                  group4 = char_group4)









X_df <- matrix(rnorm(n * p), n, p)
colnames(X_df) <- paste0("X", 1:p)

X_df <- cbind(X_df, matrix(rnorm(n*100, sd=0.5), nrow = n))

# Create true betas
true_betas <- runif(p, -1, 1)
true_betas <- c(true_betas, rep(0,100))

beta1 <- matrix(true_betas[1:5], ncol=1)
beta2 <- matrix(true_betas[6:10], ncol=1)

# Define character lists for each group of features
char_group1 <- colnames(X_df)[1:5]
char_group2 <- colnames(X_df)[6:10]
char_group3 <- colnames(X_df)[11:15]
char_group4 <- colnames(X_df)[16:20]

X_group1 = X_df[, 1:5]
X_group2 = X_df[, char_group2]

# Add interaction effects
interaction_effect <- (X_group1 %*% beta1) * (X_group2 %*% beta2) * 1
X_mat <- as.matrix(X_df)
y <- X_mat %*% true_betas + interaction_effect + rnorm(n)

char_list <- list(group1 = char_group1,
                  group2 = char_group2,
                  group3 = char_group3,
                  group4 = char_group4)

# Run the elastic net interaction function
result <- elastic_net_interaction(X_df, y, char_list, alpha = 0.01, lambda_seq = NULL, nfolds = 5, tol = 1e-8, max_iters = 100)
# Print the results
print(result)


# Create a synthetic dataset with just random values
# X_df <- as.data.frame(matrix(rnorm(n * p), n, p))
# colnames(X_df) <- paste0("X", 1:p)
# y <- matrix(rnorm(n), n, 1)
# Run the elastic net interaction function
# result <- elastic_net_interaction(X_df, y, char_list, alpha = 0.5, lambda_seq = NULL, nfolds = 5, tol = 1e-4, max_iters = 1000)
# Print the results
# print(result)