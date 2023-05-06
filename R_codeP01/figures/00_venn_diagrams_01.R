





install.packages("VennDiagram")

install.packages("futile.options")
install.packages("lambda.r")
install.packages("eulerr")

install.packages("ggvenn")


library(eulerr)

# Create sets with named elements


ridge_b_P <- c(RPS6KA5, PSMB7 , PRKACA POLR2A PMS2, PIK3CA PFDN2 PCNA, OAZ1 NSD1 , MTOR MAPK1, HDAC2 GATA3, G6PD ERCC1 , CETN2 BAX, ADD1, SKP1)


lasso_b_P <- c("LEFTY2", "GATA3", "CACNA1H", "EFNA3", "HOXA9", "CAMK2B", "BMPR1B", "NSD1", "CA12", "HOXA7", "JAG1", "APOE", "PLA2G2A", "TAPBP", "S100A7", "CALML5", "HDAC2", "CHIT1", "CBLC", "FGF13")
lasso_b_ROR <- c( "CHIT1"   "LEFTY2"  "CA12"    "CACNA1H" "PMS2"    "E2F5"    "FGF13" "HOXA7"   "FZD9"    "ACTR3B"  "EFNA3"   "APOE"    "PROM1"   "BBOX1", "CETN2"   "ITGB1"   "HDAC2"   "IFT140"  "RELN"    "ACVR1B")
lasso_rc_P <- c("GATA3"   "LEFTY2"  "CACNA1H" "HOXA9"   "CA12"    "BMPR1B"  "EFNA3" "NSD1"    "HOXA7"   "CAMK2B"  "APOE"    "THY1"    "JAG1"    "CBLC", "CHIT1"   "TAPBP"   "CALML5"  "RPS6KB1" "PIM1"    "EGLN2")
lasso_rc_ROR <- c("CHIT1"   "CA12"    "LEFTY2"  "PMS2"    "E2F5"    "ITGB1"   "CACNA1H", "APOE"    "HOXA7"   "HDAC2"   "THY1"    "EFNA3", "PROM1", "EGF", "FGF13"   "FZD9",  "IFT140",  "WDR77",  "CCL4",  "BBOX1")

elastic_b_P <- c("LEFTY2", "GATA3", "CACNA1H", "EFNA3", "CAMK2B", "NSD1", "BMPR1B", "HOXA9", "CA12", "APOE", "JAG1", "PLA2G2A", "HOXA7", "FGF13", "TAPBP", "FAM198B", "HDAC2", "CALML5", "EYA2", "S100A7")




boosting_b_P <- c("EFNA3", "BMPR1B", "CHIT1", "CACNA1H", "CAMK2B", "LEFTY2", "CA12", "HOXA7", "CALML5", "EGLN2", "DKK1", "CDCA7L", "CKB", "NSD1", "OLFML2B", "SFRP4", "CD84", "KIT", "ZFYVE9", "GATA3")




# Create a named list of sets
venn_data_b_P <- list(Lasso = lasso_v, `Elastic Net` = elastic_v, Boosting = boosting_v)
venn_data_b_ROR <- list(Lasso = lasso_v, `Elastic Net` = elastic_v, Boosting = boosting_v)
venn_data_rc_P <- list(Lasso = lasso_v, `Elastic Net` = elastic_v, Boosting = boosting_v)
venn_data_rc_ROR <- list(Lasso = lasso_v, `Elastic Net` = elastic_v, Boosting = boosting_v)
# Fit the Euler diagram
fit <- euler(venn_data)

# Plot the Euler diagram
plot(fit)
# Extract the intersections
sets

pdf("figures/venn_diagram_b_p.pdf", width = 10, height = 8)
plot(fit)
dev.off()




##############################
# 
##############################

# Create a matrix with 0's
intersection_matrix <- matrix(0, nrow = length(unique(unlist(venn_data))), ncol = length(venn_data))

# Set row names as unique gene names
rownames(intersection_matrix) <- unique(unlist(venn_data))

# Set column names as set names
colnames(intersection_matrix) <- names(venn_data)

# Fill the matrix with 1's for the presence of genes in each set
for (i in seq_along(venn_data)) {
  intersection_matrix[rownames(intersection_matrix) %in% venn_data[[i]], i] <- 1
}

# Convert the matrix to a data frame
intersection_df <- apply(intersection_matrix, 2, function(x) ifelse(x == 0, " ", "X"))
intersection_df <- as.data.frame(intersection_df, stringsAsFactors = FALSE)


# Calculate the selection frequency for each gene (based on the number of "X" in each row)
intersection_df$selection_frequency <- apply(intersection_df, 1, function(x) sum(x == "X"))

# Sort the data frame by selection frequency (descending) and gene names (ascending)
intersection_df_sorted <- intersection_df[order(-intersection_df$selection_frequency, rownames(intersection_df)), ]

# Remove the selection_frequency column
intersection_df_sorted$selection_frequency <- NULL

# Print the sorted data frame
print(intersection_df_sorted)

# Convert the sorted data frame to a LaTeX table
latex_table <- xtable(intersection_df_sorted, caption = "Gene presence in Lasso, ElasticNet, and Boosting sets")

# Set the table positioning to 'H' (place the table exactly here)
print(latex_table, floating.environment = "table", table.placement = "!h", sanitize.text.function = identity)




plot_eulerr_with_labels <- function(fit) {
  plot(fit)
  
  # Extract intersection matrix
  intersections <- as.data.frame(fit$intersections)
  
  # Loop through each intersection
  for (i in 1:nrow(intersections)) {
    row <- intersections[i, ]
    genes <- row[[1]]
    
    # If the intersection has genes, display them
    if (length(genes) > 0) {
      gene_labels <- paste(genes, collapse = ", ")
      text(row$x, row$y, labels = gene_labels, cex = 0.8, font = 2)
    }
  }
}

# Fit the Euler diagram
fit <- euler(venn_data)

# Plot the Euler diagram with gene names
plot_eulerr_with_labels(fit)


plot_eulerr_with_labels <- function(fit) {
  plot(fit)
  
  # Extract intersection matrix
  intersections <- as.data.frame(fit$intersections)
  
  # Loop through each intersection
  for (i in 1:nrow(intersections)) {
    row <- intersections[i, ]
    genes <- row[["label"]]
    
    # If the intersection has genes, display them
    if (length(genes) > 0) {
      gene_labels <- paste(genes, collapse = ", ")
      text(row$x, row$y, labels = gene_labels, cex = 0.8, font = 2)
    }
  }
}

# Fit the Euler diagram
fit <- euler(venn_data)

# Plot the Euler diagram with gene names
plot_eulerr_with_labels(fit)

################################

library(ggvenn)


# Create a named list of sets
venn_data <- list(Lasso = lasso_v, ElasticNet = elastic_v, Boosting = boosting_v)

# Create a dataframe for ggvenn
venn_df <- data.frame(
  Gene = unlist(venn_data),
  Lasso = unlist(lapply(venn_data, `%in%`, x = lasso_v)),
  ElasticNet = unlist(lapply(venn_data, `%in%`, x = elastic_v)),
  Boosting = unlist(lapply(venn_data, `%in%`, x = boosting_v))
)

# Plot the Venn diagram
ggvenn(venn_df, columns = c("Lasso", "ElasticNet", "Boosting")) +
  labs(title = "Venn diagram with gene names") +
  theme_bw()

# Get the ggvenn data
venn_results <- ggvenn::ggvenn_data(venn_df, columns = c("Lasso", "ElasticNet", "Boosting"))

# Print the gene names for each intersection
print(venn_results)


################################
# CODE bellow not in use

library(venn)

# Calculate intersections
venn_intersection <- venn(venn_data)
# Print the gene names for each intersection
for (intersection in names(venn_intersection)) {
  cat(intersection, ":\n")
  print(venn_intersection[[intersection]])
  cat("\n")
}



library(VennDiagram)

# Corrected variables
lasso_v <- c("LEFTY2", "GATA3", "CACNA1H", "EFNA3", "HOXA9", "CAMK2B", "BMPR1B", "NSD1", "CA12", "HOXA7", "JAG1",
             "APOE", "PLA2G2A", "TAPBP", "S100A7", "CALML5", "HDAC2", "CHIT1", "CBLC", "FGF13")
elastic_v <- c("LEFTY2", "GATA3", "CACNA1H", "EFNA3", "CAMK2B", "NSD1", "BMPR1B", "HOXA9", "CA12", "APOE", "JAG1",
               "PLA2G2A", "HOXA7", "FGF13", "TAPBP", "FAM198B", "HDAC2", "CALML5", "EYA2", "S100A7")
boosting_v <- c("EFNA3", "BMPR1B", "CHIT1", "CACNA1H", "CAMK2B", "LEFTY2", "CA12", "HOXA7", "CALML5", "EGLN2",
                "DKK1", "CDCA7L", "CKB", "NSD1", "OLFML2B", "SFRP4", "CD84", "KIT", "ZFYVE9", "GATA3")

venn_data <- list(Set1 = lasso_v, Set2 = elastic_v, Set3 = boosting_v)


# Plot the Venn diagram
venn.plot <- venn.diagram(
  x = venn_data,
  filename = NULL,
  category.names = c("Lasso", "Elastic Net", "Boosting"),
  fill = c("red", "green", "blue"),
  alpha = 0.50,
  cex = 2,
  show.plot.labels = TRUE
)

# Use grid.newpage() to avoid plotting on top of a previous plot
grid.newpage()
grid.draw(venn.plot)

# Modified function
venn_diagram <- function(venn_data, name) {
  venn.diagram(
    x = venn_data,
    filename = paste0("figures/", name, ".pdf"),
    category.names = c("Lasso", "Elastic Net", "Boosting"),
    fill = c("red", "green", "blue"),
    alpha = 0.50,
    cex = 2,
    show.plot.labels = TRUE
  )
}

venn_diagram(venn_data, "venn_diagram_b_p")
