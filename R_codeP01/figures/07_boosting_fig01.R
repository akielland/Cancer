###################
## LASSO figures ##
###################
library(reshape2)

## IMPORT AND PREPARE DATA

# 1000 bootstrap fits

## 771 genes -> proliferation score - bootstrap
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/xg_b_obj_771_p.RData")
head(xg_b_obj_771_p$coef_matrix)[,1:6]

# Make data frame for correlations
cor_vec_b <- e_b_obj_771_p$cor_vec
cor_vec_b_finite <- cor_vec_b[is.finite(cor_vec_b)]
corr_df_b_p <- data.frame(correlation = cor_vec_b_finite)


# 771 genes -> ROR score - bootstrap
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/xg_b_obj_771_ROR_p.RData")
head(xg_b_obj_771_ROR_p$coef_matrix)[,1:6]
e_b_obj_771_ROR_p$co_matrix <- t(e_b_obj_771_ROR_p$co_matrix)

# Make data frame for correlations
cor_vec_b_ROR <- e_b_obj_771_ROR_p$cor_vec
cor_vec_b_ROR_finite <- cor_vec_b_ROR[is.finite(cor_vec_b_ROR)]
corr_df_b_ROR <- data.frame(correlation = cor_vec_b_ROR_finite)


# 200 repeats 5-folds cross-validations

### 771 genes -> proliferation score - repeated cross-validation
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/xg_c_obj_771_prolif.RData")
head(xg_c_obj_771_prolif$coef_matrix)[,1:6]
e_c_obj_771_prolif$co_matrix <- t(e_c_obj_771_prolif$co_matrix)

# Make data frame for correlations
cor_vec_c <- e_c_obj_771_prolif$cor_vec
cor_vec_c_finite <- cor_vec_c[is.finite(cor_vec_c)]
corr_df_r_p <- data.frame(correlation = cor_vec_c_finite)


# 771 genes -> ROR score - repeated cross-validation
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/xg_c_obj_771_RORprolif.RData")
head(xg_c_obj_771_RORprolif$coef_matrix)[,1:6]
e_c_obj_771_RORprolif$co_matrix <- t(e_c_obj_771_RORprolif$co_matrix)

# Make data frame for correlations
cor_vec_c_ROR <- e_c_obj_771_RORprolif$cor_vec
cor_vec_c_ROR_finite <- cor_vec_c_ROR[is.finite(cor_vec_c_ROR)]
datcorr_df_r_ROR <- data.frame(correlation = cor_vec_c_ROR_finite)




##############
## SLECTION ##
##############


gene_freq_boost <- function(matrix, name, n_models=1000, hight=15){
  # Order features based on their selection frequency in a df
  coef_matrix <- matrix
  frequency <- data.frame(Feature = colnames(coef_matrix), Frequency = colSums(!is.na(coef_matrix)) / (n_models))
  frequency <- frequency[order(frequency$Frequency, decreasing = TRUE),]
  
  # extract the top features
  perc_best50 <- 0.5 # how often they were selected in percentage
  perc_best10 <- 0.10
  perc_best30 <- 0.30
  
  # Calculate top_features selected 50%
  top_features_with_index = which(colSums(!is.na(coef_matrix)) >= perc_best50 * n_models)
  top_feature_names = colnames(coef_matrix)[top_features_with_index]
  cat("\n")
  cat("Features selected 50% or more times:", "\n")
  if (length(top_feature_names) == 0) {
    cat("None selected that many times \n")
  } else {
    cat(top_feature_names, "\n")
  }
  cat("\n")
  
  num_features_to_keep <- 20 
  # count the frequency of each feature in coef_matrix
  counts <- colSums(!is.na(coef_matrix))
  # sort the features based on their frequency
  sorted_features <- names(sort(counts, decreasing = TRUE))
  cat("Top 20 featrues:", "\n")
  print(sorted_features[1:num_features_to_keep])
  
  # Number of genes that were selected more than 10% of the times
  more_than_10_percent <- nrow(frequency[frequency$Frequency >= perc_best10, ])
  cat("\n")
  cat("Number of genes selected more than 10% of the times:", more_than_10_percent, "\n")
  
  # Number of genes that were selected more than 30% of the times
  more_than_30_percent <- nrow(frequency[frequency$Frequency >= perc_best30, ])
  cat("\n")
  cat("Number of genes selected more than 30% of the times:", more_than_30_percent, "\n")
  
  # Create a horizontal histogram of the features that are selected more than 10% of the times
  h <- ggplot(frequency[frequency$Frequency >= perc_best10,], aes(x = Frequency, y = reorder(Feature, Frequency))) +
    geom_bar(stat = "identity") +
    xlab("Selection Frequency") +
    ylab("Genes") +
    theme(axis.text.y = element_text(angle = 0, hjust = 0))

  ggsave(paste0("figures/", name, ".pdf"), plot = h, width = 10, height = hight)
  print(h)
}

gene_freq_boost(xg_b_obj_771_p$coef_matrix, "boost_90perc_b_p", hight = 12.2)
gene_freq_boost(xg_b_obj_771_ROR_p$coef_matrix, "boost_90perc_b_ROR", hight = 10.1)
gene_freq_boost(xg_c_obj_771_prolif$coef_matrix, "boost_90perc_rc_p", hight = 7.8)
gene_freq_boost(xg_c_obj_771_RORprolif$coef_matrix, "boost_90perc_rc_ROR", hight = 3.7)


######################################################
# Selecting top 20 coefficients by size in a BoxPlot #
######################################################
# No sizes for stumps


#################
## PERFORMANCE ##
#################

# Histogram correlations
correlations_finite <- cor_vec[is.finite(cor_vec)]
cor_df <- data.frame(correlation = correlations_finite)

h1 <- ggplot(cor_df, aes(x=correlation)) +
  geom_histogram(bins = 30, color = "black", fill = "white") +
  xlab("Correlation") +
  ylab("Frequency") +
  ggtitle("Histogram of Correlation Values")
# show(h1)

# Histogram MSE
MSE_df <- data.frame(MSE = MSE_vec)

h2 <- ggplot(MSE_df, aes(x=MSE)) +
  geom_histogram(bins = 30, color = "black", fill = "white") +
  xlab("MSE") +
  ylab("Frequency") +
  ggtitle("Histogram of MSE Values")
# show(h2)

# Display the histograms side by side
grid.arrange(h1, h2, ncol=2)




