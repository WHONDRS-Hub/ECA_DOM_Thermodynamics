# Load required libraries
library(WGCNA)
library(ggplot2)
library(reshape2)

moduleEigengenes = read.csv('Module_Eigengenes_11-24.csv', row.names = 1)
moduleEigengenes <- moduleEigengenes[, colnames(moduleEigengenes) != "MEgrey"]

# 2. Perform PCA
pca <- prcomp(moduleEigengenes, scale. = FALSE) # No scaling
explained_variance <- summary(pca)$importance[2,] * 100
pca_scores <- as.data.frame(pca$x) # PCA scores for samples
pca_loadings <- as.data.frame(pca$rotation) # Loadings for modules

# Add treatment information for plotting
pca_scores$treatment <- case_when(grepl("W", row.names(pca_scores)) ~ "Wet",grepl("D", row.names(pca_scores)) ~ "Dry")


ggplot(pca_scores, aes(x = PC1, y = PC2, color = treatment)) +
  geom_point(size = 3) +
  geom_segment(data = pca_loadings, aes(x = 0, y = 0, xend = PC1 * 5, yend = PC2 * 5), 
               arrow = arrow(length = unit(0.2, "cm")), color = "darkgrey") +
  geom_text(data = pca_loadings, aes(x = PC1 * 5, y = PC2 * 5, label = rownames(pca_loadings)),
            hjust = 0.5, vjust = 1.5, size = 4, color = "darkgrey") +
  theme_bw() +
  labs(title = "PCA of Module Eigengenes",
       x = paste0("PC1: ", round(explained_variance[1], 1), "%"),
       y = paste0("PC2: ", round(explained_variance[2], 1), "%")) + 
  scale_color_manual(values = c("darkorange", "lightblue"))

# 4. Identify the Principal Component with Maximum Treatment Separation
# Calculate treatment means for each PC
mean_diff <- colMeans(subset(pca_scores, treatment == "Wet")[, -ncol(pca_scores)]) - 
  colMeans(subset(pca_scores, treatment == "Dry")[, -ncol(pca_scores)])

# Find the PC with the highest absolute mean difference
max_diff_pc <- names(which.max(abs(mean_diff)))

# 5. Correlate Modules with the Principal Component
correlations <- sapply(colnames(moduleEigengenes), function(mod) {
  cor(moduleEigengenes[[mod]], pca_scores[[max_diff_pc]], method = "pearson")
})

# Convert correlations to a data frame for visualization
cor_df <- data.frame(Module = names(correlations), Correlation = correlations)
cor_df <- cor_df[order(-abs(cor_df$Correlation)), ]

# 6. Visualize Correlation Strengths
ggplot(cor_df, aes(x = reorder(Module, -abs(Correlation)), y = Correlation)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_bw() +
  labs(title = paste("Module Correlations with", max_diff_pc),
       x = "Module", y = "Correlation Strength")

# 7. Output Results
cat("The PC showing the greatest separation between Wet and Dry is:", max_diff_pc, "\n")
cat("Correlations of modules with", max_diff_pc, ":\n")
print(cor_df)