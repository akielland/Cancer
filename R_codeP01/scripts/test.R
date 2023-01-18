library(glmnet)

# load data
data(mtcars)
x = model.matrix(mpg ~ ., mtcars)[,-1]
y = mtcars$mpg

# specify a range of alpha values to search
alphas = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)

# specify lambda values
lambdas = 10^seq(3, -3, length.out = 100)

# run elastic net with cross-validation and alpha search
fit = cv.glmnet(x, y, alpha = alphas, lambda = lambdas, nfolds = 5)

# find the row with the minimum cross-validated error
idx = which.min(fit$cvm)

# extract the corresponding alpha and lambda values
best_alpha = fit$alpha[idx]
best_lambda = fit$lambda[idx]

# print the best alpha and lambda values
print(paste("Best alpha:", best_alpha))
print(paste("Best lambda:", best_lambda))

____________
min_error = min(fit$cvm)
idx = which(fit$cvm == min_error)
best_alpha = fit$alpha[idx]
best_lambda = fit$lambda[idx]


df <- data.frame(col1 = rnorm(100), 
                 col2 = rnorm(100), 
                 col3 = rnorm(100))

# Use ggplot to create a histogram
ggplot(df, aes(x = names(df))) + 
  geom_histogram(aes(y = ..count..), 
                 binwidth = 0.5, 
                 fill = "blue", 
                 color = "black") + 
  ggtitle("Histogram of Column Values")
