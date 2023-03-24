library(glmnet)
set.seed(123)

elastic_net_interaction <- function(x_df, y, groups, alpha = 0.5, lambda_seq = NULL, nfolds = 5, tol = 1e-4, max_iters = 1000) {

  x <- as.matrix(x_df)
  
  # Create a list of groups, each group containing the column indices of x corresponding to the feature names in char_list
  groups <- lapply(char_list, function(char_group) {
    sapply(char_group, function(feature_name) which(colnames(x_df) == feature_name))
  })
  
  # Define interaction terms
  interaction_terms <- function(x, betas, group1, group2) {
    (x[, group1, drop = FALSE] %*% betas[group1]) * (x[, group2, drop = FALSE] %*% betas[group2])
  }
  
  # learn betas without interaction terms to get starting beta values, using ridge
  fit_beta0 <- cv.glmnet(x, y, alpha = 0, lambda = lambda_seq, nfolds = nfolds)
  best_lambda <- fit_beta0$lambda.min
  beta_main <- coef(fit_beta0, s = best_lambda)[-1]  # exclude intercept
  
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
    
    # Fit model
    pf <- betas
    pf[pf == 0] <- 1e-2
    pf <- 1 / pf
    fit_with_interactions <- cv.glmnet(x_with_interactions, y, alpha = alpha, lambda = lambda_seq, nfolds = nfolds, penalty.factor = pf)
    best_lambda <- fit_with_interactions$lambda.min
    
    # Update betas using the new beta values
    new_betas <- as.vector(coef(fit_with_interactions))[-1]  # exclude intercept
    

    new_beta_main <- new_betas[1:ncol(x)]
    
    best_lambda_idx <- which(fit_with_interactions$lambda == best_lambda)
    dev_glmnet_curr <-  (1 - fit_with_interactions$glmnet$dev.ratio[best_lambda_idx]) * fit_with_interactions$glmnet$nulldev
    obj_dev <- abs(dev_glmnet_curr - dev_glmnet)/dev_glmnet_curr
    
    # obj_dev <- abs(new_beta_main - beta_main)
    # obj_dev <- max(obj_dev)
    # if (obj_dev == 0) {
    #   obj_dev <- 1
    # } else {
    #   obj_dev <- mean(abs((new_beta_main - beta_main))[new_beta_main != 0 | beta_main != 0])
    # }
    
    if ( obj_dev < tol ) {  # Only check convergence for the betas without interaction terms
      converged <- TRUE
    } else {
      dev_glmnet <- dev_glmnet_curr
      beta_main <- new_beta_main
      betas <- new_betas  # Update betas only if not converged
    }
    
    obj_diff <- c(obj_diff, dev_glmnet)
  }
  # betas <- new_betas
  beta_interaction <- betas[(ncol(x) + 1):length(betas)]
  
  # make named vectors for beta_main and beta_interaction
  named_beta_main <- setNames(beta_main, colnames(x_df))
  
  # make named vector for beta_interaction
  interaction_names <- apply(interaction_indices, 1, function(idx) {
    paste(names(char_list)[idx[1]], "*", names(char_list)[idx[2]])
  })
  named_beta_interaction <- setNames(beta_interaction, interaction_names)
  
  return(list(obj_diff=obj_diff, beta_main = named_beta_main, beta_interaction = named_beta_interaction, iterations = iter, best_lambda = best_lambda))
  # return(list(beta_main = beta_main, beta_interaction = beta_interaction, iterations = iter, best_lambda = best_lambda))
}

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
boot_idx <- sample(1:n, 1000, replace = TRUE)

## Bootstrapping a large data set in order to capture > 1 interactions
set.seed(123)
x_df <- x_df[boot_idx, ]
y <- y[boot_idx]
  
# Run the elastic net interaction function
set.seed(123)
result <- elastic_net_interaction(x_df, y, char_list, alpha = 0.5, lambda_seq = NULL, nfolds = 5, tol = 1e-3, max_iters = 300)
# Print the results
print(result)
plot(result$obj_diff, type = "l", lty = 1)

