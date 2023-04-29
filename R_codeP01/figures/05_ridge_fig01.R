###################################
# Ridge Figures
###################################

# objects / dataframes
# cor_vec=cor_vec, co_matrix=co_matrix, MSE_vec = MSE_vec))

## IMPORT AND PREPARE DATA

# 1000 bootstrap fits
# RUN: r_b_obj_771_p <- ridge_boot(X, Y, X, Y)
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/r_b_obj_771_p.RData")
head(r_b_obj_771_p$co_matrix)[,1:6]

# RUN: r_b_obj_771_ROR_p <- ridge_boot(X, Y, X, Y)
head(r_b_obj_771_ROR_p$co_matrix)[,1:6]
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/r_b_obj_771_ROR_p.RData")

# 200 repeats 5-folds cross-validations

# RUN: r_c_771_prolif <- ridge_rep_cv(prolif_771genes, func=ridge_sample, folds, repeats, method="pearson")
head(r_c_771_prolif$coef_matrix)[,1:6]
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/r_c_771_prolif.RData")

# RUN: r_c_771_RORprolif <- ridge_rep_cv(ROR_prolif_771genes, func=ridge_sample, folds, repeats, method="pearson")
head(r_c_771_RORprolif$coef_matrix)[,1:6]
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/r_c_771_RORprolif.RData")



# 771 genes -> proliferation score - bootstrap
# RUN: r_b_obj_771_p <- ridge_boot(X, Y, X, Y)
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/r_b_obj_771_p.RData")
head(r_b_obj_771_p$co_matrix)[,1:6]

# Assuming the coefficient matrix is stored in r_b_obj_771_p$co_matrix
co_matrix <- r_b_obj_771_p$co_matrix
# Convert the matrix to a data frame and add row names as a new column (gene names)
co_df <- as.data.frame(co_matrix)
co_df$gene <- rownames(co_df)

# Make data frame of correlations
cor_vec_b <- r_b_obj_771_p$cor_vec
cor_vec_b_finite <- cor_vec_b[is.finite(cor_vec_b)]
data1 <- data.frame(correlation = cor_vec_b_finite)



# 771 genes -> ROR-proliferation score - bootstrap
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/r_b_obj_771_ROR_p.RData")
cor_vec_b_ROR <- r_b_obj_771_ROR_p$cor_vec
vec_b_ROR <- xg_b_obj_771_ROR_p$cor_vec

cor_vec_b_ROR_finite <- cor_vec_b_ROR[is.finite(cor_vec_b_ROR)]
data2 <- data.frame(correlation = cor_vec_b_ROR_finite)


### 771 genes -> proliferation score - repeated cross-validation
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/r_c_obj_771_prolif.RData")
cor_vec_c <- r_c_obj_771_prolif$cor_vec

cor_vec_c_finite <- cor_vec_c[is.finite(cor_vec_c)]
data3 <- data.frame(correlation = cor_vec_c_finite)


# 771 genes -> ROR-proliferation score - repeated cross-validation
# (ridge)
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/r_c_771_RORprolif.RData")
cor_vec_c_ROR <- r_c_771_RORprolif$cor_vec

cor_vec_c_ROR_finite <- cor_vec_c_ROR[is.finite(cor_vec_c_ROR)]
data4 <- data.frame(correlation = cor_vec_c_ROR_finite)




###########
## COEFFICIENT ANALYSIS
############


## Coefficient size distributions

# Calculate the average coefficient for each gene
co_df$avg_coefficient <- rowMeans(co_df[, 1:ncol(co_matrix)])

# Create a histogram of the average coefficients with a specified range using ggplot2
histogram_plot <- ggplot(co_df, aes(x = avg_coefficient)) +
  geom_histogram(binwidth = 0.001, fill = "blue", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(x = "Average Coefficient", y = "Frequency") +
  scale_x_continuous(limits = c(-0.015, 0.015))

# Display the histogram plot
print(histogram_plot)

# Save the plot as a PDF file
ggsave("figures/ridge_histogram_plot.pdf", plot = histogram_plot, width = 10, height = 8)


## Selecting top 20 coefficients by size

# Melt the data frame to a long format
co_df_long <- melt(co_df, id.vars = "gene", variable.name = "ridge_run", value.name = "coefficient")

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
  # labs(title = "Horizontal Box Plot of Top 20 Coefficients", x = "Gene", y = "Coefficient values") +
  labs(x = "Gene", y = "Coefficient values") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot as a PDF file
ggsave("figures/ridge_boxplot_top20.pdf", plot = box_plot, width = 10, height = 8)


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
