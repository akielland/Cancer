

path <- "/Users/anders/Documents/MASTER/Cancer/model_instances/lk_AllGenesProlif.RData"

report02_out <- function(){
  load(path)  
  cat(sum(is.na(lasso_k_ob$correlations))/length(lasso_k_ob$correlations))
}
  
report02_out()

## ANALYSIS
sum(is.na(lasso_k_ob$correlations))/length(lasso_k_ob$correlations)    # fraction of NA

## Correlations
mean(na.omit(lasso_k_ob$correlations))
mean(lasso_k_ob$correlations, na.rm=TRUE)
median(lasso_k_ob$correlations, na.rm=TRUE)
var(lasso_k_ob$correlations, na.rm=TRUE)
sd(na.omit(lasso_k_ob$correlations))

# Histogram
correlations_finite <- lasso_k_ob$correlations[is.finite(lasso_k_ob$correlations)]
cor_df <- data.frame(correlation = correlations_finite)

ggplot(cor_df, aes(x=correlation)) +
  geom_histogram(bins = 30, color = "black", fill = "white") +
  xlab("Correlation") +
  ylab("Frequency") +
  ggtitle("Histogram of Correlation Values")

## MSE
mean(lasso_k_ob$MSE_vec)
sd(lasso_k_ob$MSE_vec)
# Histogram
MSE_df <- data.frame(MSE = lasso_k_ob$MSE_vec)

ggplot(MSE_df, aes(x=MSE)) +
  geom_histogram(bins = 30, color = "black", fill = "white") +
  xlab("MSE") +
  ylab("Frequency") +
  ggtitle("Histogram of MSE Values")


## Most prevalent Features
# Order the features based on their selection frequency
frequency <- data.frame(Feature = colnames(lasso_k_ob$coef_matrix), Frequency = colSums(lasso_k_ob$coef_matrix != 0) / (repeats * folds))
frequency <- frequency[order(frequency$Frequency, decreasing = TRUE),]
frequency[1:5, 1:2]

# Bar plot of the selection frequency of the features
n_best <- 50
ggplot(frequency[1:n_best ,], aes(x = Frequency, y = reorder(Feature, Frequency))) +
  geom_bar(stat = "identity") +
  xlab("Selection Frequency") +
  ylab("Features") +
  ggtitle("Selection Frequency of Features") +
  theme(axis.text.y = element_text(angle = 0, hjust = 0))


## Extracting best features
# calculate the number of features to keep
perc_best <- 0.1
num_features_to_keep <- round(perc_best * ncol(lasso_k_ob$coef_matrix))
# count the frequency of each feature in coef_matrix
counts <- colSums(lasso_k_ob$coef_matrix != 0)
# sort the features based on their frequency
sorted_features <- names(sort(counts, decreasing = TRUE))
sorted_features[1:num_features_to_keep]

## Extracting best features for POST lasso
# extract the top features
perc_best <- 0.1 # the %/100 best fraction 
top_features_with_index = which(colSums(lasso_k_ob$coef_matrix != 0) >= perc_best * repeats)
top_feature_names = colnames(lasso_k_ob$coef_matrix)[top_features_with_index]
# make linear formula of top features
response <- "Y"
formula_string <- paste(response, "~", paste(top_features, collapse = " + "))
formula <- as.formula(formula_string)
# Make data.frame containing only top features
post_lasso_df <- cbind(Proliferation_ALLgenes[, 1], Proliferation_ALLgenes[, top_features_with_index + 1])
colnames(post_lasso_df)[1] <- "Y"

```

