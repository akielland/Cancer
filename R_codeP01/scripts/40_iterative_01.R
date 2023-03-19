# Load required libraries
library(glmnet)

# Custom elastic-net function with interaction terms
elastic_net_interaction <- function(x, y, groups, alpha = 0.5, lambda = 0.1, n_iters = 100) {
  # x: feature matrix (n x p)
  # y: response vector (n x 1)
  # groups: list of vectors containing column indices of x for each group
  # alpha: elastic-net mixing parameter (0 <= alpha <= 1)
  # lambda: regularization parameter
  # n_iters: number of iterations for learning betas and gammas
  
  # Define interaction terms
  interaction_terms <- function(x, groups, betas) {
    interactions <- 1
    for (i in seq_along(groups)) {
      interactions <- interactions * (x[, groups[[i]], drop = FALSE] %*% betas[[i]])
    }
    return(interactions)
  }
  
  # Initialize betas and gamma
  betas <- lapply(groups, function(g) rep(0, length(g)))
  gamma <- 0
  
  for (iter in 1:n_iters) {
    # Update betas
    for (i in seq_along(groups)) {
      # Update the betas for each group by training an elastic-net model
      # Exclude the group to be updated from the interaction terms
      other_groups <- groups[-i]
      other_betas <- betas[-i]
      interactions <- interaction_terms(x, other_groups, other_betas)
      
      x_group <- x[, groups[[i]], drop = FALSE] * interactions
      fit <- glmnet(x_group, y, alpha = alpha, lambda = lambda)
      
      # Update betas for the current group
      betas[[i]] <- coef(fit)[-1]  # exclude intercept
    }
    
    # Update gamma
    interactions <- interaction_terms(x, groups, betas)
    fit_gamma <- glmnet(interactions, y, alpha = alpha, lambda = lambda)
    gamma <- coef(fit_gamma)[-1]  # exclude intercept
  }
  
  return(list(betas = betas, gamma = gamma))
}

# Example usage
n <- 100
p1 <- 3
p2 <- 2
x1 <- matrix(rnorm(n * p1), n, p1)
x2 <- matrix(rnorm(n * p2), n, p2)
x <- cbind(x1, x2)
y <- rnorm(n)

groups <- list(1:p1, (p1 + 1):(p1 + p2))
result <- elastic_net_interaction(x, y, groups)

print(result)
