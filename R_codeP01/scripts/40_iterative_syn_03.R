

synergistic <- function(df, char_list, alpha = 0.5, lambda_seq = NULL, nfolds = 5, tol = 1e-4, max_iters = 100) {
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
  fit_beta0 <- cv.glmnet(X, y, alpha = 0, lambda = lambda_seq, nfolds = nfolds)
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
    
    # Remove interactions terms = 0
    interactions_sd = apply(interactions, 2, sd)
    interactions_scaled = interactions[, interactions_sd != 0]
    
    # Make feature matrix with original features and interactions
    X_with_interactions <- cbind(X, interactions_scaled)
    
    # Fit model using adaptive elastic-net
    pf <- rep(1, ncol(X_with_interactions))
    pf[1:ncol(X)] <- betas[1:ncol(X)] 
    pf[pf == 0] <- 1e-2
    pf <- 1 / pf
    if (ncol(interactions_scaled) > 0) {
      pf[(ncol(X) + 1):ncol(X_with_interactions)] <- 0
    }
    fit_with_interactions <- cv.glmnet(X_with_interactions, y, alpha = alpha, lambda = lambda_seq, nfolds = nfolds, penalty.factor = pf)
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

  beta_interaction <- betas[ncol(X) + 1:ncol(interactions)]
  
  # make named vectors for beta_main and beta_interaction
  named_beta_main <- setNames(beta_main, colnames(df)[-1])
  
  # make named vector for beta_interaction
  interaction_names <- apply(interaction_indices, 1, function(idx) {
    paste(names(char_list)[idx[1]], "*", names(char_list)[idx[2]])
  })
  named_beta_interaction <- setNames(beta_interaction, interaction_names)
  
  
  return(list(model=fit_with_interactions, coef=beta_interaction))
  # return(list(obj_diff=obj_diff, beta_main = named_beta_main, beta_interaction = named_beta_interaction, iterations = iter, best_lambda = best_lambda))
 }

#result <- synergistic(ROR_prolif_771genes, char_list, alpha = 0.1, lambda_seq = NULL, nfolds = 5, tol = 1e-5, max_iters = 100)



# Function: repeated k-fold cross validation: 
syn_sym_rep_cv <- function(df_data=sim.data, groups=char_list, folds=5, repeats=2, method="pearson") {
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
      out <- synergistic(train_data, char_list, alpha = 0.1, lambda_seq = NULL, nfolds = 5, tol = 1e-5, max_iters = 100)
 
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

t <- syn_sym_rep_cv()
