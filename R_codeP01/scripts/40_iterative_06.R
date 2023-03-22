library(glmnet)
set.seed(123)

elastic_net_interaction <- function(x_df, y, groups, alpha = 0.5, lambda_seq = NULL, nfolds = 5, tol = 1e-4, max_iters = 1000) {
  # x_df: data frame of features (n x p)
  # y: response vector (n x 1)
  # char_list: list of character vectors representing groups of features
  # alpha: elastic-net mixing parameter (0 <= alpha <= 1)
  # lambda: regularization parameter
  # tol: tolerance for convergence
  # max_iters: maximum number of iterations
  
  # Convert data frame to matrix
  x <- as.matrix(x_df)
  
  # If no lambda sequence is provided, use cv.glmnet's default
  if (is.null(lambda_seq)) {
    lambda_seq <- NULL
  } else {
    lambda_seq <- as.vector(lambda_seq)
  }
  
  # Create a list of groups, each group containing the column indices of x corresponding to the feature names in char_list
  groups <- lapply(char_list, function(char_group) {
    sapply(char_group, function(feature_name) which(colnames(x_df) == feature_name))
  })
  
  # Define interaction terms
  interaction_terms <- function(x, betas, group1, group2) {
    x[, group1, drop = FALSE] %*% betas[group1] * x[, group2, drop = FALSE] %*% betas[group2]
  }
  
  # learn betas without interaction terms to get starting beta values, using ridge
  fit_beta0 <- cv.glmnet(x, y, alpha = 0, lambda = lambda_seq, nfolds = nfolds)
  best_lambda <- fit_beta0$lambda.min
  beta_main <- coef(fit_beta0, s = best_lambda)[-1]  # exclude intercept
  
  # Calculate interaction terms for all combinations of groups
  n_groups <- length(groups)
  interaction_indices <- t(combn(n_groups, 2))
  
  # Initialize interaction matrix
  n_combinations <- ncol(combn(n_groups, 2))
  interactions <- matrix(0, nrow(x), n_combinations)
  
  # Calculate interaction terms using the current betas
  for (i in seq_len(nrow(interaction_indices))) {
    group1 <- groups[[interaction_indices[i, 1]]]
    group2 <- groups[[interaction_indices[i, 2]]]
    interactions[, i] <- interaction_terms(x, beta_main, group1, group2)
  }
  
  # Make feature matrix with original features and interactions
  x_with_interactions <- cbind(x, interactions)
  
  # Initialize variables for convergence check
  converged <- FALSE
  iter <- 0
  obj_diff <- NULL
  betas <- rep(0, ncol(x_with_interactions))
  
  # Learn betas outside the interaction term until convergence
  while (!converged && iter < max_iters) {
    iter <- iter + 1
    
    # Calculate interaction terms using the current betas
    for (i in seq_len(nrow(interaction_indices))) {
      group1 <- groups[[interaction_indices[i, 1]]]
      group2 <- groups[[interaction_indices[i, 2]]]
      interactions[, i] <- interaction_terms(x, beta_main, group1, group2)
    }
    
    # Make feature matrix with original features and interactions
    x_with_interactions <- cbind(x, interactions)
    
    # Fit model
    fit_with_interactions <- cv.glmnet(x_with_interactions, y, alpha = alpha, lambda = lambda_seq, nfolds = nfolds)
    best_lambda <- fit_with_interactions$lambda.min
    
    # Update betas using the new beta values
    new_betas <- as.vector(coef(fit_with_interactions))[-1]  # exclude intercept
    
    # Check for convergence
    new_beta_main <- new_betas[1:ncol(x)]
    
    obj_diff <- c(obj_diff, mean(abs(new_betas - betas)))
    # print(mean(abs(new_betas - betas)))
    
    if (max(abs(new_beta_main - beta_main)) < tol ) {  # Only check convergence for the betas without interaction terms
      converged <- TRUE
    } else {
      beta_main <- new_beta_main
      betas <- new_betas  # Update betas only if not converged
    }
  
  }
  # betas <- new_betas
  beta_interaction <- betas[(ncol(x) + 1):length(betas)]
  
  # Create named vectors for beta_main and beta_interaction
  named_beta_main <- setNames(beta_main, colnames(x_df))
  
  # Create named vector for beta_interaction
  interaction_names <- apply(interaction_indices, 1, function(idx) {
    paste(names(char_list)[idx[1]], "*", names(char_list)[idx[2]])
  })
  named_beta_interaction <- setNames(beta_interaction, interaction_names)
  
  return(list(obj_diff=obj_diff, beta_main = named_beta_main, beta_interaction = named_beta_interaction, iterations = iter, best_lambda = best_lambda))
}

# Create a synthetic dataset with some interaction effect
set.seed(42)
n <- 100
p <- 20

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