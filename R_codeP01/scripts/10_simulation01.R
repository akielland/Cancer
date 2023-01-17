




set.seed(123)
n <- 100 # number of samples
p <- 20 # number of genetic features
X_ <- matrix(rnorm(n*p), nrow=n, ncol=p) # generate random data for genetic features

# assign effect to half of the genes
p_half <- p/2
beta <- c(rnorm(p_half), rep(0, p_half)) # generate random coefficients for the first half of the genes and 0 for the other half

Y_ <- X_ %*% beta + rnorm(n) # generate outcome variable Y










set.seed(123)
n <- 100 # number of samples
p <- 20 # number of genetic features

# Correlated genes
p_correlated <- round(p * 0.3) # 30% of the genes
X_correlated <- matrix(rnorm(n*p_correlated), nrow=n, ncol=p_correlated) # generate correlated genes data
X_correlated <- X_correlated + rnorm(n, mean = 0, sd = 0.1) # add small random noise

# Non-correlated genes
p_non_correlated <- p - p_correlated
X_non_correlated <- matrix(rnorm(n*p_non_correlated), n, p_non_correlated) # generate non-correlated genes data

# Merge correlated and non-correlated genes
X <- cbind(X_correlated, X_non_correlated)

beta <- c(rnorm(p_half), rep(0, p_half)) # generate random coefficients for the first half of the genes and 0 for the other half
Y <- X %*% beta + rnorm(n) # generate outcome variable Y