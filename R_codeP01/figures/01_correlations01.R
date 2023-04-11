# Install and load extrafont package
install.packages("extrafont")
library(extrafont)

# Load required libraries
library(ggplot2)
library(gridExtra)



# 771 genes -> proliferation score - bootstrap
# (ridge)
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/r_b_obj_771_p.RData")
cor_vec_b <- r_b_obj_771_p$cor_vec
# (lasso)
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/lb_obj_771_prolif.RData")
cor_vec_b <- lb_obj_771_prolif

cor_vec_b_finite <- cor_vec_b[is.finite(cor_vec_b)]
data1 <- data.frame(correlation = cor_vec_b_finite)


# 771 genes -> ROR-proliferation score - bootstrap
# (ridge)
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/r_b_obj_771_ROR_p.RData")
cor_vec_b_ROR <- r_b_obj_771_ROR_p$cor_vec
# (lasso)
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/lb_obj_771_RORprolif.RData")
cor_vec_b_ROR <- lb_obj_771_RORprolif

cor_vec_b_ROR_finite <- cor_vec_b_ROR[is.finite(cor_vec_b_ROR)]
data2 <- data.frame(correlation = cor_vec_b_ROR_finite)


### 771 genes -> proliferation score - repeated cross-validation
# (ridge) 
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/r_c_obj_771_prolif.RData")
cor_vec_c <- r_c_obj_771_prolif$cor_vec
# (lasso)
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/lc_obj_771_prolif.RData")
cor_vec_c <- lc_obj_771_prolif

cor_vec_c_finite <- cor_vec_c[is.finite(cor_vec_c)]
data3 <- data.frame(correlation = cor_vec_c_finite)


# 771 genes -> ROR-proliferation score (repeated cross-validation)
# (ridge)
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/r_c_771_RORprolif.RData")
cor_vec_c_ROR <- r_c_771_RORprolif$cor_vec
# (lasso)
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/lc_obj_771_RORprolif.RData")
cor_vec_c_ROR <- report02_out(lc_obj_771_RORprolif)

cor_vec_c_ROR_finite <- cor_vec_c_ROR[is.finite(cor_vec_c_ROR)]
data4 <- data.frame(correlation = cor_vec_c_ROR_finite)






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
ggsave("figures/lasso_histogram.pdf", grid_histograms, width = 11, height = 8.5)
ggsave("figures/elastic_histogram.pdf", grid_histograms, width = 11, height = 8.5)




# Histogram correlations
# Create individual ggplot histograms
# Create individual ggplot histograms
p1 <- ggplot(data1, aes(x = correlation)) +  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.8) + labs(title = "Histogram 1") + theme_custom
p2 <- ggplot(data2, aes(x = correlation)) + geom_histogram(bins = 30, fill = "steelblue", alpha = 0.8) + labs(title = "Histogram 2") + theme_custom
p3 <- ggplot(data3, aes(x = correlation)) + geom_histogram(bins = 30, fill = "steelblue", alpha = 0.8) + labs(title = "Histogram 3") + theme_custom
p4 <- ggplot(data4, aes(x = correlation)) + geom_histogram(bins = 30, fill = "steelblue", alpha = 0.8) + labs(title = "Histogram 4") + theme_custom

p1 <- ggplot(data1, aes(x = correlation, y = ..density..)) + geom_histogram(aes(y = ..density..), bins = 30, fill = "steelblue", alpha = 0.8) +
  labs(title = "Histogram 1") + theme_custom + ylab("Frequency")

p2 <- ggplot(data2, aes(x = correlation, y = ..density..)) + geom_histogram(aes(y = ..density..), bins = 30, fill = "steelblue", alpha = 0.8) +
  labs(title = "Histogram 2") + theme_custom

p3 <- ggplot(data3, aes(x = correlation, y = ..density..)) + geom_histogram(aes(y = ..density..), bins = 30, fill = "steelblue", alpha = 0.8) +
  labs(title = "Histogram 3") + theme_custom + xlab("X-axis label") + ylab("Frequency")

p4 <- ggplot(data4, aes(x = correlation, y = ..density..)) + geom_histogram(aes(y = ..density..), bins = 30, fill = "steelblue", alpha = 0.8) +
  labs(title = "Histogram 4") + theme_custom + xlab("X-axis label")


correlations_finite1 <- cor_vec1[is.finite(cor_vec1)]
cor_df_1 <- data.frame(correlation = correlations_finite1)

h1 <- ggplot(cor_df_1, aes(x=correlation)) +
  geom_histogram(bins = 30, color = "black", fill = "white") +
  xlab("Correlation") +
  ylab("Frequency") +
  ggtitle("Histogram of Correlation Values")


# Import fonts available in your system
font_import()

# Use the same font family as your LaTeX document, e.g., 'Times New Roman'
theme_set(theme_bw(base_family = "Times New Roman"))


ggplot(...) + theme_bw(base_size = 12, base_family = "Times New Roman")
ggplot(...) + theme_bw(base_size = 12, base_family = "Times New Roman") +
  theme(line = element_line(size = 1))

# Example of a custom color palette
my_colors <- c("blue", "red", "green")
# Use custom colors in ggplot
ggplot(...) +
  geom_line(aes(color = factor(group))) +
  scale_color_manual(values = my_colors)

