
library(sparsegl)
set.seed(123)

elastic_net_interaction_sparsegl <- function(x_df, y, groups, lambda1, lambda2, tol = 1e-4, max_iters = 100) {
  
  x <- as.matrix(x_df)
  
  groups_indices <- lapply(groups, function(group) {
    sapply(group, function(feature_name) which(colnames(x_df) == feature_name))
  })
  
  interaction_terms <- function(x, betas, group1, group2) {
    (x[, group1, drop = FALSE] %*% betas[group1]) * (x[, group2, drop = FALSE] %*% betas[group2])
  }
  
  n_groups <- length(groups_indices)
  interaction_indices <- t(combn(n_groups, 2))
  
  n_combinations <- ncol(combn(n_groups, 2))
  interactions <- matrix(0, nrow(x), n_combinations)
  
  betas <- rep(1, ncol(x))
  
  converged <- FALSE
  iter <- 0
  obj_diff <- NULL
  
  while (!converged && iter < max_iters) {
    iter <- iter + 1
    
    for (i in seq_len(nrow(interaction_indices))) {
      group1 <- groups_indices[[interaction_indices[i, 1]]]
      group2 <- groups_indices[[interaction_indices[i, 2]]]
      interactions[, i] <- interaction_terms(x, betas, group1, group2)
    }
    
    x_with_interactions <- cbind(x, interactions)
    
    fit_with_interactions <- sparsegl(x_with_interactions, y, group = unlist(groups_indices), lambda1 = lambda1, lambda2 = lambda2)
    
    new_betas <- fit_with_interactions$beta
    new_beta_main <- new_betas[1:ncol(x)]
    
    obj_dev <- sqrt(sum((betas - new_betas)^2)) / sqrt(sum(new_betas^2))
    
    if (obj_dev < tol) {
      converged <- TRUE
    } else {
      betas <- new_betas
    }
    
    obj_diff <- c(obj_diff, obj_dev)
  }
  
  beta_interaction <- betas[(ncol(x) + 1):length(betas)]
  
  named_beta_main <- setNames(new_beta_main, colnames(x_df))
  
  interaction_names <- apply(interaction_indices, 1, function(idx) {
    paste(names(groups)[idx[1]], "*", names(groups)[idx[2]])
  })
  named_beta_interaction <- setNames(beta_interaction, interaction_names)
  
  return(list(obj_diff = obj_diff, beta_main = named_beta_main, beta_interaction = named_beta_interaction, iterations = iter))
}




# Make sure to install and load the required package
# install.packages("sparsegl")
library(sparsegl)

sparse_group_lasso_interaction <- function(x_df, y, groups, lambda_seq = NULL, nfolds = 5, tol = 1e-4, max_iters = 100) {
  x <- as.matrix(x_df)
  
  # Create a list of groups, each group containing the column indices of x corresponding to the feature names in char_list
  groups <- lapply(char_list, function(char_group) {
    sapply(char_group, function(feature_name) which(colnames(x_df) == feature_name))
  })
  
  # Define interaction terms
  interaction_terms <- function(x, betas, group1, group2) {
    (x[, group1, drop = FALSE] %*% betas[group1]) * (x[, group2, drop = FALSE] %*% betas[group2])
  }
  
  # learn betas without interaction terms to get starting beta values
  fit_beta0 <- sparsegl(x, y, groups = unlist(groups), lambda = lambda_seq, nfolds = nfolds)
  best_lambda <- fit_beta0$lambda.min
  beta_main <- fit_beta0$beta
  
  # calculate interaction terms for all combinations of groups
  n_groups <- length(groups)
  interaction_indices <- t(combn(n_groups, 2))
  
  # initialize interaction matrix
  n_combinations <- ncol(combn(n_groups, 2))
  interactions <- matrix(0, nrow(x), n_combinations)
  
  # Calculate interaction terms using the current betas
  for (i in seq_len(nrow(interaction_indices))) {
    group1 <- groups[[interaction_indices[i, 1]]]
    group2 <- groups[[interaction_indices[i, 2]]]
    interactions[, i] <- interaction_terms(x, beta_main, group1, group2)
  }
  
  # make feature matrix with original features and interactions
  x_with_interactions <- cbind(x, interactions)
  
  # initialize variables for convergence checking
  converged <- FALSE
  iter <- 0
  obj_diff <- NULL
  betas <- rep(1, ncol(x_with_interactions))
  dev_glmnet <- obj_dev <- Inf
  
  # Learn betas outside the interaction term until convergence
  while (!converged && iter < max_iters) {
    iter <- iter + 1
    if(!iter %% 10) cat("iter =", iter, "; dev_glmnet =", dev_glmnet, "; obj_dev =", obj_dev, "\n")
    # Calculate interaction terms using the current betas
    for (i in seq_len(nrow(interaction_indices))) {
      group1 <- groups[[interaction_indices[i, 1]]]
      group2 <- groups[[interaction_indices[i, 2]]]
      interactions[, i] <- interaction_terms(x, beta_main, group1, group2)
    }
    
    # Make feature matrix with original features and interactions
    x_with_interactions <- cbind(x, interactions)
    
    

    library(sparsegl)
    
    elastic_net_interaction_sparsegl <- function(x_df, y, groups, lambda_seq = NULL, nfolds = 5, tol = 1e-4, max_iters = 100) {
      
      x <- as.matrix(x_df)
      
      # Create a list of groups, each group containing the column indices of x corresponding to the feature names in char_list
      groups <- lapply(char_list, function(char_group) {
        sapply(char_group, function(feature_name) which(colnames(x_df) == feature_name))
      })
      
      # Define interaction terms
      interaction_terms <- function(x, betas, group1, group2) {
        (x[, group1, drop = FALSE] %*% betas[group1]) * (x[, group2, drop = FALSE] %*% betas[group2])
      }
      
      # learn betas without interaction terms to get starting beta values
      fit_beta0 <- sparsegl(x, y, group = unlist(groups), lambda = lambda_seq)
      best_lambda <- fit_beta0$lambda.min
      beta_main <- fit_beta0$beta
      
      # calculate interaction terms for all combinations of groups
      n_groups <- length(groups)
      interaction_indices <- t(combn(n_groups, 2))
      
      # initialize interaction matrix
      n_combinations <- ncol(combn(n_groups, 2))
      interactions <- matrix(0, nrow(x), n_combinations)
      
      # Calculate interaction terms using the current betas
      for (i in seq_len(nrow(interaction_indices))) {
        group1 <- groups[[interaction_indices[i, 1]]]
        group2 <- groups[[interaction_indices[i, 2]]]
        interactions[, i] <- interaction_terms(x, beta_main, group1, group2)
      }
      
      # make feature matrix with original features and interactions
      x_with_interactions <- cbind(x, interactions)
      
      # initialize variables for convergence checking
      converged <- FALSE
      iter <- 0
      obj_diff <- NULL
      betas <- rep(1, ncol(x_with_interactions))
      dev_sparsegl <- obj_dev <- Inf
      
      # Learn betas outside the interaction term until convergence
      while (!converged && iter < max_iters) {
        iter <- iter + 1
        if(!iter %% 10) cat("iter =", iter, "; dev_sparsegl =", dev_sparsegl, "; obj_dev =", obj_dev, "\n")
        # Calculate interaction terms using the current betas
        for (i in seq_len(nrow(interaction_indices))) {
          group1 <- groups[[interaction_indices[i, 1]]]
          group2 <- groups[[interaction_indices[i, 2]]]
          interactions[, i] <- interaction_terms(x, beta_main, group1, group2)
        }
        
        # Make feature matrix with original features and interactions
        x_with_interactions <- cbind(x, interactions)
        
        # Fit model
        fit_with_interactions <- sparsegl(x_with_interactions, y, group = unlist(groups), lambda = lambda_seq)
        best_lambda <- fit_with_interactions$lambda.min
        
        # Update betas using the new beta values
        new_betas <- fit_with_interactions$beta
        
        new_beta_main <- new_betas[1:ncol(x)]
        
        best_lambda_idx <- which(fit_with_interactions$lambda == best_lambda)
        dev_sparsegl_curr <-  (1 -
                                 



# Create a synthetic dataset with some interaction effect
set.seed(42)
n <- 100
p <- 20

X_df <- matrix(rnorm(n * p), n, p)
colnames(X_df) <- paste0("X", 1:p)

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

X_group1 = X_df[, char_group1]
X_group2 = X_df[, char_group2]
X_group4 = X_df[, char_group4]

# Add interaction effects
interaction_effect <- (X_group1 %*% beta1) * (X_group2 %*% beta2) * 1 + (X_group1 %*% beta1) * (X_group4 %*% beta4) * 0.5
X_mat <- as.matrix(X_df)
y <- X_mat %*% true_betas + interaction_effect + rnorm(n)

char_list <- list(group1 = char_group1,
                  group2 = char_group2,
                  group3 = char_group3,
                  group4 = char_group4)

x_df <- scale(X_df)
y <- scale(y)
# Run the elastic net interaction function
set.seed(123)
result <- elastic_net_interaction(x_df, y, char_list, alpha = 0.05, lambda_seq = NULL, nfolds = 5, tol = 1e-4, max_iters = 300)
# Print the results
print(result)
plot(result$obj_diff, type = "l", lty = 1)



