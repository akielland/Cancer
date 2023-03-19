# Load required libraries
library(glmnet)

# Custom elastic-net function with interaction terms
elastic_net_interaction <- function(x, y, groups, alpha = 0.5, lambda = 0.1, n_iters = 3) {
  # x: feature matrix (n x p)
  # y: response vector (n x 1)
  # groups: list of vectors containing column indices of x for each group
  # alpha: elastic-net mixing parameter (0 <= alpha <= 1)
  # lambda: regularization parameter
  # n_iters: number of iterations for learning betas and gamma
  
  # Define interaction terms
  interaction_terms <- function(x, betas, group1, group2) {
    x[, group1, drop = FALSE] %*% betas[group1] * x[, group2, drop = FALSE] %*% betas[group2]
  }
  
  # learn betas without interaction terms to get starting beta values, using ridge
  fit_beta0 <- glmnet(x, y, alpha = 0, lambda = lambda)
  beta0 <- coef(fit_beta0)[-1]  # exclude intercept
  
  # starting values for the interaction term
  interactions = interaction_terms(x, beta0, groups[[1]], groups[[2]])
  
  # Learn betas outside the interaction term
  for (iter in 1:n_iters) {
    # make feature matrix with original features and interactions
    x_with_interactions <- cbind(x, interactions) # correct the function name
    # fit model
    fit_with_interactions <- glmnet(x_with_interactions, y, alpha = alpha, lambda = lambda)    
    # update interaction terms using the new beta values
    betas_and_interaction_betas <- coef(fit_with_interactions)[-1] # exclude intercept
    interactions = interaction_terms(x, betas_and_interaction_betas[1:ncol(x)], groups[[1]], groups[[2]])
  }
  
  betas <- betas_and_interaction_betas[1:ncol(x)] # correct the index range
  beta_interaction <- betas_and_interaction_betas[(ncol(x) + 1):length(betas_and_interaction_betas)]
  
  return(list(betas = betas, beta_interaction = beta_interaction))
}


# Example usage
n <- 100
p1 <- 3
p2 <- 2
p3 <- 4
x1 <- matrix(rnorm(n * p1), n, p1)
x2 <- matrix(rnorm(n * p2), n, p2)
x3 <- matrix(rnorm(n * p3), n, p3)
x <- cbind(x1, x2, x3)
y <- rnorm(n)

groups <- list(1:p1, (p1 + 1):(p1 + p2))
result <- elastic_net_interaction(x, y, groups)

print(result)
