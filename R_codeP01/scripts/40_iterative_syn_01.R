

synergistic <- function(df, char_list, alpha=0.5, nfolds=5, tol=1e-4, max_iters=100) {
  # Extract predictor matrix and response vector; assuming response variable in first column
  X <- as.matrix(df[, -1])
  X <- scale(X)
  y <- df[, 1]
  
  # Create a list of groups, each group containing the column name and indices of df corresponding to the feature names in char_list
  groups <- lapply(char_list, function(char_group) {
    sapply(char_group, function(feature_name) which(colnames(df) == feature_name) - 1)
  })

  # function caculating the interaction terms
  interaction_terms <- function(X, betas, group1, group2) {
    (X[, group1, drop = FALSE] %*% betas[group1]) * (X[, group2, drop = FALSE] %*% betas[group2])
  }
  
  # learn betas without interaction terms to get starting beta values, using ridge
  fit_beta0 <- cv.glmnet(X, y, alpha = 0, nfolds = nfolds)
  best_lambda <- fit_beta0$lambda.min
  beta_main <- coef(fit_beta0, s = best_lambda)[-1]  # exclude intercept
  
  # calculate interaction terms for all combinations of groups
  n_groups <- length(groups)
  interaction_indices <- t(combn(n_groups, 2))
  
  # initialize interaction matrix
  n_combinations <- ncol(combn(n_groups, 2))
  interactions <- matrix(0, nrow(X), n_combinations)
  
  # Calculate interaction terms using the current betas
  for (i in seq_len(nrow(interaction_indices))) {
    group1 <- groups[[interaction_indices[i, 1]]]
    group2 <- groups[[interaction_indices[i, 2]]]
    interactions[, i] <- interaction_terms(X, beta_main, group1, group2)
  }
  
  # make feature matrix with original features and interactions
  X_with_interactions <- cbind(X, interactions)
  
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
      interactions[, i] <- interaction_terms(X, beta_main, group1, group2)
    }
    
    # Make feature matrix with original features and interactions
    X_with_interactions <- cbind(X, interactions)
    
    # Fit model using adaptive elastic-net
    pf <- betas
    pf[pf == 0] <- 1e-2
    pf <- 1 / pf
    pf[(ncol(X) + 1):length(betas)] <- 0
    fit_with_interactions <- cv.glmnet(X_with_interactions, y, alpha=alpha, nfolds=nfolds, penalty.factor = pf)
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
      beta_main <- new_beta_main
      betas <- new_betas  # Update betas only if not converged
    }
    
    obj_diff <- c(obj_diff, dev_glmnet)
  }
  # betas <- new_betas
  beta_interaction <- betas[(ncol(X) + 1):length(betas)]
  
  # make named vectors for beta_main and beta_interaction
  named_beta_main <- setNames(beta_main, colnames(df)[-1])
  
  # make named vector for beta_interaction
  interaction_names <- apply(interaction_indices, 1, function(idx) {
    paste(names(char_list)[idx[1]], "*", names(char_list)[idx[2]])
  })
  named_beta_interaction <- setNames(beta_interaction, interaction_names)
  
  return(list(model=fit_with_interactions, obj_diff=obj_diff, beta_main = named_beta_main, beta_interaction = named_beta_interaction, iterations = iter, best_lambda = best_lambda))
 }

result <- synergistic(ROR_prolif_771genes, char_list, alpha=0.001, nfolds = 5, tol=1e-5, max_iters=100)


# Function: repeated k-fold cross validation: 
synergistic_rep_cv <- function(df_data, char_list, alpha, folds=5, repeats=1, tol=1e-5, max_iters=100, method="pearson") {
  n_models <- repeats * folds
  print(n_models)
  
  cor_vec <- rep(NA, n_models)
  MSE_vec <- rep(NA, n_models)
  coef_matrix <- matrix(NA, nrow = n_models, ncol = ncol(df_data[,-1]))
  colnames(coef_matrix) <- colnames(df_data[,-1])
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
      out <- synergistic(df_data, char_list, alpha, nfolds, tol, max_iters)

      model <- out$model
      coef_matrix[row_index, ] <- out$coef[-1]   # remove intercept
      
      pred = predict(model, newx = as.matrix(test_data)[,-1], type = "response", s = "lambda.min")  
      cor_vec[row_index]  <- suppressWarnings(cor(pred, test_data[,1], method=method))
      MSE_vec[row_index] <- mean((pred - test_data[,1])^2)
      
      if(!row_index %% 10) cat(row_index, "")
      row_index <- row_index + 1
    }
  }
  return(list(cor_vec=cor_vec, MSE_vec=MSE_vec, coef_matrix=coef_matrix))
}


char_list <- list(imm_inf_ = char_immune_inf,
                  prolif_ = char_prolif,
                  ER_sing_ = char_ER_signaling,
                  anti_pres_ = char_antigen_present,
                  angiogen_ = char_angiogenesis
                  )

char_list <- list(ER_sing_ = char_ER_signaling,
                  prolif_ = char_prolif
                  )


set.seed(123)
t1 <- synergistic(ROR_prolif_771genes, char_list, alpha=0.1, nfolds=5, repeats=5, tol=1e-6, max_iters=100)



dat <- list(ROR_prolif_771genes=ROR_prolif_771genes, char_list=char_list, synergistic=synergistic)
save(dat, file = "dat.RData")

res## Bootstrapping a large data set in order to capture > 1 interactions
set.seed(123)
std_df <- scale(ROR_prolif_771genes)
dim(std_df)

boot_idx <- sample(1:dim(std_df)[1], 1000, replace = TRUE)

df_boot <- std_df[boot_idx, ]

result <- synergistic(df_boot, char_list, alpha = 0.01, lambda_seq = NULL, nfolds = 5, tol = 1e-6, max_iters = 100)






# Create a synthetic dataset with some interaction effect
set.seed(42)
n <- 100
p <- 20

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
interaction_effect <- (X_group1 %*% beta1) * (X_group2 %*% beta2) * 1 + (X_group1 %*% beta1) * (X_group4 %*% beta4) * 0.5
interaction_effect <- (X_group1 %*% beta1) * (X_group2 %*% beta2) * 1 
X_mat <- as.matrix(X_df)
y <- X_mat %*% true_betas + interaction_effect + rnorm(n)

char_list <- list(group1 = char_group1,
                  group2 = char_group2,
                  group3 = char_group3,
                  group4 = char_group4)


# Run the elastic net interaction function
set.seed(123)
result <- synergistic(X_df, char_list, alpha = 0.05, lambda_seq = NULL, nfolds = 5, tol = 1e-2, max_iters = 300)
# Print the results
print(result)
plot(result$obj_diff, type = "l", lty = 1)

