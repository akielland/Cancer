# Load required libraries
library(glmnet)

# Custom elastic-net function with interaction terms
elastic_net_interaction <- function(x, y, groups, alpha = 0.5, lambda = 0.1, n_iters = 100) {
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
  
  # Initialize betas
  betas <- rep(0, ncol(x))
  beta_interaction <- 0
  
  for (iter in 1:n_iters) {
    # Learn betas outside the interaction term
    x_no_interaction <- x[, -c(groups[[1]], groups[[2]])]
    fit_no_interaction <- glmnet(x_no_interaction, y, alpha = alpha, lambda = lambda)
    betas[-c(groups[[1]], groups[[2]])] <- coef(fit_no_interaction)[-1]  # exclude intercept
    
    # Update the interaction term
    interactions <- interaction_terms(x, betas, groups[[1]], groups[[2]])
    
    # Learn beta_interaction
    fit_interaction <- glmnet(interactions, y - x_no_interaction %*% betas[-c(groups[[1]], groups[[2]])], alpha = alpha, lambda = lambda)
    beta_interaction <- coef(fit_interaction)[-1]  # exclude intercept
  }
  
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
