## Bootstrapping a large data set in order to capture > 1 interactions
set.seed(123)
std_df <- scale(ROR_prolif_771genes)
dim(std_df)

boot_idx <- sample(1:dim(std_df)[1], 1000, replace = TRUE)

df_boot <- std_df[boot_idx, ]

result <- synergistic(df_boot, char_list, alpha = 0.1, lambda_seq = NULL, nfolds = 5, tol = 1e-6, max_iters = 100)


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