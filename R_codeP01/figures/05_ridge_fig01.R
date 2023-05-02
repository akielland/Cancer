###################################
# Ridge Figures
###################################
library(dplyr)
library(tidyr)
library(xtable)
library(reshape2)       

# objects / dataframes
# cor_vec=cor_vec, co_matrix=co_matrix, MSE_vec = MSE_vec))

## IMPORT AND PREPARE DATA

# 1000 bootstrap fits

## 771 genes -> proliferation score - bootstrap
# RUN: r_b_obj_771_p <- ridge_boot(X, Y, X, Y)
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/r_b_obj_771_p.RData")
head(r_b_obj_771_p$co_matrix)[,1:6]

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
# RUN: r_b_obj_771_ROR_p <- ridge_boot(X, Y, X, Y)
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/r_b_obj_771_ROR_p.RData")
head(r_b_obj_771_ROR_p$co_matrix)[,1:6]

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
# RUN: r_c_771_prolif <- ridge_rep_cv(prolif_771genes, func=ridge_sample, folds, repeats, method="pearson")
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/r_c_771_prolif.RData")
head(r_c_771_prolif$coef_matrix)[,1:6]

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
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/r_c_771_RORprolif.RData")
head(r_c_771_RORprolif$coef_matrix)[,1:6]

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



##########################
## COEFFICIENT ANALYSIS ##
##########################


## Coefficient size distributions 

# Calculate the average coefficient for each gene one the model ridge_b_p
# (code below should be fixed...!!!!!!!!!!!!!!!)
cof_df_b_p$avg_coefficient <- rowMeans(cof_df_b_p[, 1:ncol(co_matrix)])
cof_df_b_p_tests$avg_coefficient <- rowMeans(cof_df_b_p[, 1:ncol(cof_df_b_p)])

# Create a histogram of the average coefficients with a specified range using ggplot2
histogram_plot <- ggplot(co_df, aes(x = avg_coefficient)) +
  geom_histogram(binwidth = 0.001, fill = "blue", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(x = "Mean Coefficient", y = "Frequency") +
  scale_x_continuous(limits = c(-0.015, 0.015))

# Display the histogram plot
print(histogram_plot)

# Save the plot as a PDF file
ggsave("figures/ridge_histogram_plot.pdf", plot = histogram_plot, width = 10, height = 8)

######################################################
# Selecting top 20 coefficients by size in a BoxPlot #
######################################################

top20 <- function(co_df, name){
  # Melt the data frame to a long format
  co_df_long <- melt(co_df, id.vars = "gene", variable.name = "model", value.name = "coefficient")
  
  # Calculate the absolute values of coefficients and find the top 20 genes with largest coefficients
  co_df_long$abs_coefficient <- abs(co_df_long$coefficient)
  top_genes <- co_df_long %>% group_by(gene) %>% summarise(mean_abs_coefficient = mean(abs_coefficient)) %>%
    top_n(20, mean_abs_coefficient) %>% select(gene)
  
  # Filter the coefficients of the top 20 genes
  co_df_long_top_genes <- co_df_long %>% filter(gene %in% top_genes$gene)
  
  # Update the gene column to be a factor with levels sorted in reverse alphabetical order
  co_df_long_top_genes$gene <- factor(co_df_long_top_genes$gene, levels = sort(unique(co_df_long_top_genes$gene), decreasing = TRUE))
  
  # Create the horizontal box plot and store it in a variable
  box_plot <- ggplot(co_df_long_top_genes, aes(x = gene, y = coefficient)) +
    geom_boxplot() +
    coord_flip() +
    theme_minimal() +
    labs(x = "Gene", y = "Coefficient values") +
    theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14))
  
  # Save plot with custom name
  ggsave(paste0("figures/", name, ".pdf"), plot = box_plot, width = 10, height = 8)
  print(box_plot)
}

top20(cof_df_b_p, "ridge_top20_b_p")
top20(cof_df_b_ROR, "ridge_top20_b_ROR")
top20(cof_df_r_p, "ridge_top20_r_p") 
top20(cof_df_r_ROR, "ridge_top20_r_ROR")


##############################################################
# Table of the top 20 genes for each model: with mean and SD #
##############################################################

# Function to process the data frame and calculate mean and SD for each gene
process_df <- function(df, model_name) {
  df_long <- melt(df, id.vars = "gene", variable.name = "model", value.name = "coefficient")
  df_summary <- df_long %>% group_by(gene) %>% 
    summarise(mean = signif(mean(coefficient), 3), sd = signif(sd(coefficient), 3), .groups = 'drop') %>% 
    mutate(model = model_name) %>% 
    slice_max(order_by = abs(mean), n = 20) %>%  # Select the top 20 genes sorted by the absolute value of the mean coefficients
    arrange(gene)  # Arrange the rows alphabetically by gene names
  return(df_summary)
}


# Process each data frame
cof_summary_b_p <- process_df(cof_df_b_p, "Bootstrap_P")
cof_summary_b_ROR <- process_df(cof_df_b_ROR, "Bootstrap_ROR")
cof_summary_r_p <- process_df(cof_df_r_p, "RepeatedCV_P")
cof_summary_r_ROR <- process_df(cof_df_r_ROR, "RepeatedCV_ROR")

# Combine the top 20 genes from each model class into a single data frame
top20_all_models <- bind_rows(cof_summary_b_p, cof_summary_b_ROR, cof_summary_r_p, cof_summary_r_ROR)

# Reshape the data frame to a wide format an rearrange the columns
top20_wide <- top20_all_models %>%
  pivot_wider(id_cols = gene, names_from = model, values_from = c(mean, sd), names_sep = "_") %>%
  arrange(gene) %>%  # Arrange the rows alphabetically by gene names
  # Reorder the columns:
  select(gene,
         mean_Bootstrap_P, sd_Bootstrap_P,
         mean_Bootstrap_ROR, sd_Bootstrap_ROR,
         mean_RepeatedCV_P, sd_RepeatedCV_P,
         mean_RepeatedCV_ROR, sd_RepeatedCV_ROR) 

# Print the combined table
print(top20_wide, n=40)

# Replace 'your_dataframe' with the actual name of your data frame
ridge_feature_df_table <- top20_wide

# Replace NAs with "-", dont work I think
# ridge_feature_df_table[is.na(ridge_feature_df_table)] <- "-"

# Convert numeric columns to character and replace NAs with ""
ridge_feature_df_table[sapply(ridge_feature_df_table, is.numeric)] <- lapply(ridge_feature_df_table[sapply(ridge_feature_df_table, is.numeric)], function(x) { ifelse(is.na(x), "", as.character(x)) })

# Merge the first two columns and the last two columns
ridge_feature_df_table <- ridge_feature_df_table %>%
  mutate(Bootstrap_P = paste(mean_Bootstrap_P, sd_Bootstrap_P, sep = " / "),
         Bootstrap_ROR = paste(mean_Bootstrap_ROR, sd_Bootstrap_ROR, sep = " / "),
         RepeatedCV_P = paste(mean_RepeatedCV_P, sd_RepeatedCV_P, sep = " / "),
         RepeatedCV_ROR = paste(mean_RepeatedCV_ROR, sd_RepeatedCV_ROR, sep = " / ")) %>%
  select(gene, Bootstrap_P, Bootstrap_ROR, RepeatedCV_P, RepeatedCV_ROR)

# Convert the data frame to LaTeX
latex_code <- xtable(ridge_feature_df_table)

# Print the LaTeX code
print(latex_code, type = "latex", include.rownames = FALSE)



#########################
# Performance Analysis ##
#########################

# Set custom theme with larger font size
theme_custom <- theme(plot.title = element_text(size = 14, face = "bold"),
                      axis.title = element_text(size = 12),
                      axis.text = element_text(size = 12))


# Create individual ggplot histograms
p1 <- ggplot(data1, aes(x = correlation)) + geom_histogram(bins = 30, fill = "steelblue", alpha = 0.8) +
  labs(title = "a)") + theme_custom + ylab("Count") + theme(axis.title.x = element_blank())

p2 <- ggplot(data2, aes(x = correlation)) + geom_histogram(bins = 30, fill = "steelblue", alpha = 0.8) +
  labs(title = "b)") + theme_custom + theme(axis.title.x = element_blank(), axis.title.y = element_blank())

p3 <- ggplot(data3, aes(x = correlation)) + geom_histogram(bins = 30, fill = "steelblue", alpha = 0.8) +
  labs(title = "c)") + theme_custom + xlab("correlation") + ylab("Count")

p4 <- ggplot(data4, aes(x = correlation)) + geom_histogram(bins = 30, fill = "steelblue", alpha = 0.8) +
  labs(title = "d)") + theme_custom + xlab("correlation") + theme(axis.title.y = element_blank())

# Combine the histograms in a 2x2 grid
grid_histograms <- grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)

# Save the grid as a PDF
ggsave("figures/ridge_histogram.pdf", grid_histograms, width = 11, height = 8.5)
