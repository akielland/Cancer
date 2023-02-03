# Load required libraries
library(prcomp)
library(ggplot2)
library(caret)

features <- prolif_771genes |> select(-Y)

pca_results <- function(fm, data){
  X_ <- model.matrix(fm, data = data)[,-1]  
  pca_results <- prcomp(X_, scale = TRUE) # create PCA object
  # Get the number of components that explain at least X % of the variance
  n_components <- sum(pca_results$sdev^2 / sum(pca_results$sdev^2) >= 0.1)
  # Extract the first n_components principal components
  df_first_pca <- data.frame(pca_results$x[, 1:n_components])
  return(df_first_pca)
}


immune_inf <- pca_results(fm_immune_inf, data=prolif_771genes)

summary(immune_inf)










# Perform PCA on the data
pca_results <- prcomp(X_, scale = TRUE)

# Get the number of components that explain at least 50% of the variance
n_components <- sum(pca_results$sdev^2 / sum(pca_results$sdev^2) >= 0.1)

# Extract the first n_components principal components and use them as predictors in a linear regression model
prolif_771genes_pca <- data.frame(pca_results$x[, 1:n_components], Y = prolif_771genes$Y)

model_fit <- lm(Y ~ ., data = prolif_771genes_pca)

# Print summary of the model
summary(model_fit)

# Plot the regression results
ggplot(prolif_771genes_pca, aes(x = PC2, y = Y)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("PCA-Regression Results")
