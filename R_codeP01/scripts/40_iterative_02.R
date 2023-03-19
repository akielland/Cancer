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
  interaction_terms <- function(x, groups) {
    n <- ncol(x)
    interactions <- matrix(0, nrow(x), choose(length(groups), 2))
    counter <- 1
    for (i in 1:(length(groups) - 1)) {
      for (j in (i + 1):length(groups)) {
        interactions[, counter] <- rowSums(outer(x[, groups[[i]], drop = FALSE], x[, groups[[j]], drop = FALSE], "*"))
        counter <- counter + 1
      }
    }
    return(interactions)
  }
  
  # Calculate the interaction terms
  interactions <- interaction_terms(x, groups)
  
  # Create a combined matrix with original features and interaction terms
  XX <- cbind(x, interactions)
  
  # Initialize betas
  betas <- rep(0, ncol(XX))
  
  # Fit the elastic-net model
  for (iter in 1:n_iters) {
    fit <- glmnet(XX, y, alpha = alpha, lambda = lambda)
    betas <- coef(fit)[-1]  # exclude intercept
    print(betas)
  }
  
  return(betas)
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

groups <- list(1:p1, (p1 + 1):(p1 + p2), (p1 + p2 + 1):(p1 + p2 + p3))
betas <- elastic_net_interaction(x, y, groups)

print(betas)
