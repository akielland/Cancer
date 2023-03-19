library(glmnet)


elastic_net_interaction <- function(x_df, y, char_list, alpha = 0.5, lambda = 0.1, n_iters = 3) {
  # x_df: data frame of features (n x p)
  # y: response vector (n x 1)
  # char_list: list of character vectors representing groups of features
  # alpha: elastic-net mixing parameter (0 <= alpha <= 1)
  # lambda: regularization parameter
  # n_iters: number of iterations for learning betas and gamma
  
  # Convert data frame to matrix
  x <- as.matrix(x_df)
  
  # Create a list of groups, each group containing the column indices of x corresponding to the feature names in char_list
  groups <- lapply(char_list, function(char_group) {
    sapply(char_group, function(feature_name) which(colnames(x_df) == feature_name))
  })
  
  # Define interaction terms
  interaction_terms <- function(x, betas, group1, group2) {
    x[, group1, drop = FALSE] %*% betas[group1] * x[, group2, drop = FALSE] %*% betas[group2]
  }
  
  # Learn betas without interaction terms to get starting beta values, using ridge
  fit_beta0 <- glmnet(x, y, alpha = 0, lambda = lambda)
  beta0 <- coef(fit_beta0)[-1]  # exclude intercept
  
  # Calculate interaction terms for all combinations of groups
  n_groups <- length(groups)
  interaction_indices <- t(combn(n_groups, 2))
  
  # Initialized interactions matrix with the number of rows equal to the number of observations in x 
  # and the number of columns equal to the number of combinations of groups
  # interactions <- matrix(0, nrow(x), ncol(interaction_indices))
  # Initialize interaction matrix
  n_combinations <- ncol(combn(n_groups, 2))
  interactions <- matrix(0, nrow(x), n_combinations)
  
  
  # Learn betas outside the interaction term
  for (iter in 1:n_iters) {
    # Calculate interaction terms using the current betas
    # updated calculation of interaction terms to consider all combinations of group
    for (i in seq_len(nrow(interaction_indices))) {
      group1 <- groups[[interaction_indices[i, 1]]]
      group2 <- groups[[interaction_indices[i, 2]]]
      interactions[, i] <- interaction_terms(x, beta0, group1, group2)
    }
    
    # Make feature matrix with original features and interactions
    x_with_interactions <- cbind(x, interactions)
    
    # Fit model
    fit_with_interactions <- glmnet(x_with_interactions, y, alpha = alpha, lambda = lambda)    
    # Update betas using the new beta values
    beta0 <- coef(fit_with_interactions)[-1]  # exclude intercept
  }
  
  betas <- beta0[1:ncol(x)]
  beta_interaction <- beta0[(ncol(x) + 1):length(beta0)]
  
  return(list(betas = betas, beta_interaction = beta_interaction))
}
# Run with test with syntetic data

# Create a synthetic dataset
set.seed(123)
n <- 100
p <- 20

X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("V", 1:p)
y <- rnorm(n)

# Define feature groups as character lists
char_group1 <- c("V1", "V2", "V3", "V4", "V5")
char_group2 <- c("V6", "V7", "V8", "V9", "V10")
char_group3 <- c("V11", "V12", "V13", "V14", "V15")
char_group4 <- c("V16", "V17", "V18", "V19", "V20")

# Create a named list with feature groups
char_list <- list(group1 = char_group1,
                  group2 = char_group2,
                  group3 = char_group3,
                  group4 = char_group4)

# Convert the feature matrix X to a data frame
X_df <- as.data.frame(X)

# Run the elastic_net_interaction function
result <- elastic_net_interaction(X_df, y, char_list, alpha = 0.5, lambda = 0.1, n_iters = 3)

# Print betas and interaction betas
print(result$betas)
print(result$beta_interaction)

