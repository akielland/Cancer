# Install and load extrafont package# 
install.packages("extrafont")
library(extrafont)

# Load required libraries
library(ggplot2)
library(gridExtra)

data1 <- data.frame(value = prolif_771genes$Y)
data2 <- data.frame(value = ROR_prolif_771genes$Y)


# Set custom theme with larger font size
theme_custom <- theme(plot.title = element_text(size = 14, face = "bold"),
                      axis.title = element_text(size = 12),
                      axis.text = element_text(size = 12))


h1 <- ggplot(data1, aes(x = value)) +  
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.8) +
  labs(title = "a)") +
  theme_custom +
  ylab("Count") +
  xlab("proliferation score") +
  scale_y_continuous(limits = c(0, 8))  # Set the y-axis range

h2 <- ggplot(data2, aes(x = value)) + 
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.8) +
  labs(title = "ab") +
  theme_custom +
  xlab("ROR score") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),      # Remove y-axis text
        axis.ticks.y = element_blank()) +   # Remove y-axis ticks
  scale_y_continuous(limits = c(0, 8))  # Set the y-axis range

# Combine the histograms in a 2x2 grid
grid_histograms <- grid.arrange(h1, h2, ncol = 2, nrow = 1)


# Save the grid as a PDF
ggsave("figures/response_histogram.pdf", grid_histograms, width = 11, height = 8.5)


