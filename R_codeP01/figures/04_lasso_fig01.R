###################
## LASSO figures ##
###################
library(reshape2)


## IMPORT AND PREPARE DATA

# 1000 bootstrap fits


## 771 genes -> proliferation score - bootstrap
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/lb_obj_771_prolif.RData")
head(lb_obj_771_prolif$coef_matrix)[,1:6]

# Make data_frame for coefficients
co_matrix <- r_b_obj_771_p$co_matrix
# Convert matrix to a data_frame and add row names as a new column (gene names)
cof_df_b_p <- as.data.frame(co_matrix)
cof_df_b_p$gene <- rownames(cof_df_b_p)

# Make data frame for correlations
cor_vec_b <- r_b_obj_771_p$cor_vec
cor_vec_b_finite <- cor_vec_b[is.finite(cor_vec_b)]
corr_df_b_p <- data.frame(correlation = cor_vec_b_finite)


# 771 genes -> ROR score - bootstrap
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/lb_obj_771_RORprolif.RData")
head(lb_obj_771_RORprolif$coef_matrix)[,1:6]

# Make data_frame for coefficients
co_matrix <- r_b_obj_771_ROR_p$co_matrix
# Convert matrix to a data_frame and add row names as a new column (gene names)
cof_df_b_ROR <- as.data.frame(co_matrix)
cof_df_b_ROR$gene <- rownames(cof_df_b_ROR)

# Make data frame for correlations
cor_vec_b_ROR <- r_b_obj_771_ROR_p$cor_vec
cor_vec_b_ROR_finite <- cor_vec_b_ROR[is.finite(cor_vec_b_ROR)]
corr_df_b_ROR <- data.frame(correlation = cor_vec_b_ROR_finite)


# 200 repeats 5-folds cross-validations

### 771 genes -> proliferation score - repeated cross-validation
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/lc_obj_771_prolif.RData")
head(lc_obj_771_prolif$coef_matrix)[,1:6]

# Make data_frame for coefficients
#NB: need to transpose the coef matrix for it to match with the bootstrap data
cof_df_r_p_transposed <- t(r_c_771_prolif$coef_matrix)
co_matrix <- cof_df_r_p_transposed
# co_matrix <- r_c_771_prolif$co_matrix
# Convert matrix to a data_frame and add row names as a new column (gene names)
cof_df_r_p <- as.data.frame(co_matrix)
cof_df_r_p$gene <- rownames(cof_df_r_p)

# Make data frame for correlations
cor_vec_c <- r_c_771_prolif$cor_vec
cor_vec_c_finite <- cor_vec_c[is.finite(cor_vec_c)]
corr_df_r_p <- data.frame(correlation = cor_vec_c_finite)


# 771 genes -> ROR score - repeated cross-validation
# RUN: r_c_771_RORprolif <- ridge_rep_cv(ROR_prolif_771genes, func=ridge_sample, folds, repeats, method="pearson")
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/lc_obj_771_RORprolif.RData")
head(lc_obj_771_RORprolif$coef_matrix)[,1:6]

#NB: need to transpose the coef matrix for it to match with the bootstrap data
cof_df_r_p_transposed <- t(r_c_771_RORprolif$coef_matrix)
co_matrix <- cof_df_r_p_transposed
# Convert matrix to a data_frame and add row names as a new column (gene names)
cof_df_r_ROR <- as.data.frame(co_matrix)
cof_df_r_ROR$gene <- rownames(cof_df_r_ROR)

# Make data frame for correlations
cor_vec_c_ROR <- r_c_771_RORprolif$cor_vec
cor_vec_c_ROR_finite <- cor_vec_c_ROR[is.finite(cor_vec_c_ROR)]
datcorr_df_r_ROR <- data.frame(correlation = cor_vec_c_ROR_finite)




##############
## SLECTION ##
##############


gene_freq <- function(object, name, n_models=1000){
  # Order features based on their selection frequency in a df
  coef_matrix <- object$coef_matrix
  frequency <- data.frame(Feature = colnames(coef_matrix), Frequency = colSums(coef_matrix != 0) / (n_models))
  frequency <- frequency[order(frequency$Frequency, decreasing = TRUE),]
  
  # extract the top features
  perc_best50 <- 0.5 # how often they were selected in percentage
  perc_best10 <- 0.1
  top_features_with_index = which(colSums(coef_matrix != 0) >= perc_best50 * n_models)
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
  counts <- colSums(coef_matrix != 0)
  # sort the features based on their frequency
  sorted_features <- names(sort(counts, decreasing = TRUE))
  cat("Top 20 featrues:", "\n")
  print(sorted_features[1:num_features_to_keep])
  
  # Number of genes that were selected more than 10% of the times
  more_than_10_percent <- nrow(frequency[frequency$Frequency >= perc_best10, ])
  cat("\n")
  cat("Number of genes selected more than 10% of the times:", more_than_10_percent, "\n")
  
  
  # Create a horizontal histogram of the features that are selected more than 10% of the times
  h <- ggplot(frequency[frequency$Frequency >= perc_best10,], aes(x = Frequency, y = reorder(Feature, Frequency))) +
    geom_bar(stat = "identity") +
    xlab("Selection Frequency") +
    ylab("Genes") +
    theme(axis.text.y = element_text(angle = 0, hjust = 0))
  
  ggsave(paste0("figures/", name, ".pdf"), plot = h, width = 10, height = 15)
  print(h)
}

gene_freq(lb_obj_771_prolif, "lasso_90perc_b_p")
gene_freq(lb_obj_771_RORprolif, "lasso_90perc_b_ROR")
gene_freq(lc_obj_771_prolif, "lasso_90perc_r_p")
gene_freq(lc_obj_771_RORprolif, "lasso_90perc_r_ROR")


######################################################
# Selecting top 20 coefficients by size in a BoxPlot #
######################################################


top20_boxplot <- function(coef_matrix, name) {
  # Calculate the frequency of each feature being selected (non-zero coefficients)
  frequency <- colSums(coef_matrix != 0)
  
  # Find the top 20 features based on frequency
  top20_features <- names(sort(frequency, decreasing = TRUE))[1:20]
  
  # Filter the coef_matrix to keep only the top 20 features
  top20_coef_matrix <- coef_matrix[, top20_features]
  
  # Convert the filtered matrix to a long format data frame
  top20_coef_long <- melt(top20_coef_matrix, varnames = "model", value.name = "coefficient")
  
  # Add the gene names to the data frame
  top20_coef_long$gene <- rep(top20_features, each = nrow(top20_coef_matrix))
  
  # Remove zero coefficient values from the data frame
  top20_coef_long <- top20_coef_long[top20_coef_long$coefficient != 0, ]
  
  # Update the gene column to be a factor with levels sorted in alphabetical order
  top20_coef_long$gene <- factor(top20_coef_long$gene, levels = sort(unique(top20_coef_long$gene), decreasing = FALSE))
  
  # Create a box plot of the top 20 features with gene names on the y-axis and coefficient values on the x-axis
  box_plot <- ggplot(top20_coef_long, aes(y = gene, x = coefficient, group = gene)) +
    geom_boxplot() +
    scale_y_discrete(limits = rev(levels(top20_coef_long$gene))) +
    theme_minimal() +
    labs(y = "Gene", x = "Coefficient values") +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14))
  
  ggsave(paste0("figures/", name, ".pdf"), plot = box_plot, width = 10, height = 8)
  print(box_plot)
}

top20_boxplot(lb_obj_771_prolif$coef_matrix, "lasso_top20_b_p")
top20_boxplot(lb_obj_771_RORprolif$coef_matrix, "lasso_top20_b_ROR")
top20_boxplot(lc_obj_771_prolif$coef_matrix, "lasso_top20_r_p")
top20_boxplot(lc_obj_771_RORprolif$coef_matrix, "lasso_top20_r_ROR")


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




