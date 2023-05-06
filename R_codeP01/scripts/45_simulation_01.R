library(glmnet)

# Function that creates a synthetic dataset (matrix) with interaction effects
syntetic_01 <- function(n, interact, p_gr, p_non_gr, i, scale=TRUE){
  # Create predictors for groups
  X <- matrix(rnorm(n * p_gr), n, p_gr)
  colnames(X) <- paste0("X", 1:p_gr)
  # Add non-gr. predictors
  Z <- matrix(rnorm(n * p_non_gr), nrow = n)
  X <- cbind(X, Z)
  colnames(X)[(p_gr + 1):(p_gr + p_non_gr)] <- paste0("Z", 1:p_non_gr)
  # X_df <- as.data.frame(X)

  
  # Define character lists for each group of features
  char_group1 <- colnames(X)[1:5]
  char_group2 <- colnames(X)[6:10]
  char_group3 <- colnames(X)[11:15]
  char_group4 <- colnames(X)[16:20]
  
  # Create matrices w/ groups of predictors -> predictors also divided into groups
  X_group1 <- X[, char_group1]
  X_group2 <- X[, char_group2]
  X_group3 <- X[, char_group3]
  X_group4 <- X[, char_group4]
  
  char_list_sim <- list(gr1 = char_group1,
                        gr2 = char_group2,
                        gr3 = char_group3,
                        gr4 = char_group4)
  
  # Create true betas
  true_betas <- runif(p_gr, -1, 1)
  # Add non effective betas = 0
  betas <- c(true_betas, rep(0, p_non_gr))
  
  # divide betas into groups
  beta1 <- matrix(true_betas[1:5], ncol=1)
  beta2 <- matrix(true_betas[6:10], ncol=1)
  beta3 <- matrix(true_betas[11:15], ncol=1)
  beta4 <- matrix(true_betas[16:20], ncol=1)
  

  if (interact == 1){
    # Add interaction effects (1): gr1 x gr2.  -> 1 interaction term between the 4 groups: gives 6 simulations
    interaction_effect <- (X_group1 %*% beta1) * (X_group2 %*% beta2) * 1
  } else if (interact == 2) {
    # Add interaction effects (2): gr1 x gr2 + gr1 x gr4
    interaction_effect <- (X_group1 %*% beta1) * (X_group2 %*% beta2) * 1 + (X_group1 %*% beta1) * (X_group4 %*% beta4) * 1
  } else if (interact == 3){
    # Add interaction effects (3): gr1 x gr2 + gr1 x gr4 + gr2 x x gr4
    interaction_effect <- (X_group1 %*% beta1) * (X_group2 %*% beta2) * 1 + (X_group1 %*% beta1) * (X_group4 %*% beta4) * 1 + (X_group2 %*% beta2) * (X_group4 %*% beta4) * 1
  }else {
    abort("Error: interact argument must be 1, 2, or 3.")
  }
    

  # use model to create response y
  y <- X %*% betas + interaction_effect + rnorm(n)
  colnames(y) <- "y"
  
  if (scale==TRUE)
    X <- scale(X)
    y <- scale(y)
  
  df <- as.data.frame(cbind(y,X))  
  return(list(df.sim=df, char_list=char_list_sim))
}
  
n <- 50
n <- 100
n <- 500
n <- 1000
p_gr <- 20
p_non_gr <- 100
interact <- 1

char_list_sim <- list(gr1 = char_group1,
                      gr2 = char_group2,
                      gr3 = char_group3,
                      gr4 = char_group4)


sim.data <- syntetic_01(n, interact, p_gr, p_non_gr, i=1, scale=FALSE)
head(sim.data)
typeof(sim.data); class(sim.data)
hist(sim.data$y)

# Function to run many simulations
run_sim <- function(n_simulations, interact, alpha, adaptive, tol){

  coef_matrix <- matrix(NA, nrow = n_simulations, ncol = 6)
  row_index <- 1
  for (i in c(1:n_simulations)){
    # i = i+12345
    set.seed(i) 
    data <- syntetic_01(n, interact, p_gr, p_non_gr, i)
    out <- synergistic_adaptive_deviance(data$df.sim, data$char_list, alpha, adaptive, lambda_seq = NULL, nfolds = 5, tol = 1e-6, max_iters = 100)
    
    # print(out$beta_interaction)
    coef_matrix[row_index, ] <- out$beta_interaction
    cat(row_index, "")
    row_index <- row_index + 1
    
  }
  colnames(coef_matrix) <- c("1x2", "1x3", "1x4", "2x3", "2x4", "3x4")
  return(coef_matrix)
}


n <- 50
n <- 100
n <- 500
n <- 1000
p_gr <- 20
p_non_gr <- 80
interact <- 1
interact <- 2
interact <- 3


n_simulations <- 100
set.seed(123)
t1 <- run_sim(n_simulations, interact, alpha = 1, adaptive = FALSE, tol=1e-4)

n_simulations <- 50
set.seed(1000)
t1 <- run_sim(n_simulations, interact, alpha = 0.5, adaptive = FALSE, tol=1e-6)

# Count the number of non-zero and non-NA values in each column
count_non_zero_non_na <- apply(t1, 2, function(x) sum(!is.na(x) & x != 0))
print(count_non_zero_non_na*2)


# Replace NA values with 0
mat_no_na <- replace(t1, is.na(t1), 0)
mat_no_na
# Calculate column means
colMeans(mat_no_na)
colSds(mat_no_na)

# Calculate the mean of non-NA values in each column
mean_non_na <- apply(t1, 2, function(x) mean(x, na.rm = TRUE))
print("Mean of non-NA values:")
print(mean_non_na)

# Calculate the proportion of NA values in each column
mean_na <- apply(t1, 2, function(x) mean(is.na(x)))
mean_na


print("Proportion of NA values:")
print(mean_na/n_simulations)

# Count the number of numeric values in each column
count_numeric <- apply(t2, 2, function(x) sum(sapply(x, is.numeric)))
print(count_numeric)
# Count the number of non-zero and non-NA values in each column
count_non_zero_non_na <- apply(t2, 2, function(x) sum(!is.na(x) & x != 0))
print(2*(count_non_zero_non_na))




# Run the elastic net interaction function
result <- syn_elastic_net_interaction(X_df, y, char_list, alpha = 0.01, lambda_seq = NULL, nfolds = 5, tol = 1e-8, max_iters = 100)
# Print the results
print(result)




## Bootstrapping a large data set in order to capture > 1 interactions
set.seed(123)
std_df <- scale(ROR_prolif_771genes)
dim(std_df)

boot_idx <- sample(1:dim(std_df)[1], 1000, replace = TRUE)

df_boot <- std_df[boot_idx, ]

result <- synergistic(df_boot, char_list, alpha = 0.1, lambda_seq = NULL, nfolds = 5, tol = 1e-6, max_iters = 100)




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

# 1 interaction term between the 4 groups: gives 6 simulations
interaction_effect <- (X_group1 %*% beta1) * (X_group2 %*% beta2) * 1 
interaction_effect <- (X_group1 %*% beta1) * (X_group3 %*% beta3) * 0.5 
interaction_effect <- (X_group1 %*% beta1) * (X_group4 %*% beta4) * 0.5 
interaction_effect <- (X_group2 %*% beta2) * (X_group3 %*% beta3) * 1 
interaction_effect <- (X_group2 %*% beta2) * (X_group4 %*% beta4) * 0.5 
interaction_effect <- (X_group3 %*% beta3) * (X_group4 %*% beta4) * 1 

# All possible interaction term between the 4 groups: gives 1 simulation
interaction_effect <- 
  (X_group1 %*% beta1) * (X_group2 %*% beta2) * 1 
+ (X_group1 %*% beta1) * (X_group3 %*% beta3) * 0.5 
+ (X_group1 %*% beta1) * (X_group4 %*% beta4) * 0.5 
+ (X_group2 %*% beta2) * (X_group3 %*% beta3) * 1 
+ (X_group2 %*% beta2) * (X_group4 %*% beta4) * 0.5 
+ (X_group3 %*% beta3) * (X_group4 %*% beta4) * 1 



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