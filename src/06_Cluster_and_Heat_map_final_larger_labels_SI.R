# Load required libraries
rm(list = ls(all = TRUE))

library(tidyverse)
library(ggplot2)
library(pheatmap)
library(viridis)
library(dendextend)
library(dplyr)
library(RColorBrewer)
library(tibble)
library(grid)
library(colorspace)
library(corrplot)
library(Hmisc)
library(stringr)
library(gridExtra)
library(ggcorrplot)
library(factoextra)
library(PMCMRplus)
library(multcompView)
library(dunn.test)
library(ggpubr)
library(patchwork)
library(MASS)

# ==== Defining paths and working directories ======
figure_path <- "Figures/"
dir.create(figure_path, showWarnings = FALSE, recursive = TRUE)

# ===== Publication-ready plotting settings =====
cluster_colors <- c("1" = "#74C0FC", "2" = "#FFD43B", "3" = "#FF8FAB")

pub_theme_s3 <- theme_bw(base_size = 22) +
  theme(
    aspect.ratio = 1,
    legend.position = "none",
    axis.text = element_text(size = 24, color = "black"),
    axis.title = element_text(size = 28, color = "black", face = "bold"),
    plot.title = element_text(size = 30, color = "black", face = "bold", hjust = 0),
    strip.text = element_text(size = 24, face = "bold"),
    plot.margin = margin(18, 18, 28, 18)
  )

# ====== Read in data ======
data <- read.csv("Data/Medians_of_Median_molecular_properties_per_site_and_treatment_unique_formulas.csv")

row.names(data) <- paste0(data$site, "_", data$Treatment)

sample_data <- read_csv(
  "v4_CM_SSS_Data_Package_1/Sample_Data/v3_CM_SSS_Sediment_Sample_Data_Summary.csv",
  comment = "#",
  na = c("N/A", -9999)
) %>%
  slice(-(1:11)) %>%
  mutate_at(vars(-Sample_Name, -Field_Name, -IGSN, -Material), as.numeric) %>%
  dplyr::select(
    site = Sample_Name,
    Mean_ATP_picomoles_per_g,
    Mean_Specific_Surface_Area_m2_per_g,
    C_percent_per_mg = "01395_C_percent_per_mg",
    N_percent_per_mg = "01397_N_percent_per_mg",
    Percent_Tot_Sand,
    Mean_Fe_mg_per_kg,
    Mean_Gravimetric_Moisture_g_per_g = "Mean_62948_Gravimetric_Moisture_g_per_g"
  )

sample_data$site <- gsub("CM", "EC", sample_data$site)
sample_data$site <- gsub("_Sediment", "", sample_data$site)

field_metadata <- read.csv("EC_Data_Package/EC_Field_Metadata.csv") %>%
  dplyr::select(site = Parent_ID, State)

ecoregions <- read.csv("Data/EC_Ecoregions.csv") %>%
  dplyr::select(site = "Parent_ID", Ecoregion = Name) %>%
  mutate(Ecoregion = str_trim(str_extract(Ecoregion, "(?<=- ).*")))

field_metadata <- merge(field_metadata, ecoregions, by = "site")

# ==== Set up data =====
dom_data <- data %>%
  dplyr::filter(site != "EC_023") %>%
  dplyr::filter(!(site %in% c("EC_012", "EC_011", "EC_053", "EC_057", "EC_052"))) %>%
  dplyr::select(site, Treatment, Median_delGcoxPerCmol, Median_Lambda)

sites <- as.data.frame(unique(dom_data$site))
names(sites) <- "site"

explanatory_data <- merge(sample_data, field_metadata, by = "site")
explanatory_data <- merge(explanatory_data, sites, by = "site")

field_sample_data <- merge(sample_data, sites, by = "site") %>%
  dplyr::select(-Mean_Gravimetric_Moisture_g_per_g)

site_metadata <- explanatory_data %>%
  dplyr::select(-Mean_Gravimetric_Moisture_g_per_g)

# ===== Calculations ======
treatment_effects <- dom_data %>%
  filter(Treatment %in% c("Wet", "Dry")) %>%
  dplyr::select(site, Treatment, starts_with("Median_")) %>%
  pivot_wider(
    id_cols = site,
    names_from = Treatment,
    values_from = starts_with("Median_")
  ) %>%
  mutate(
    effect_delGcoxPerCmol = Median_delGcoxPerCmol_Wet / Median_delGcoxPerCmol_Dry,
    effect_Lambda = Median_Lambda_Wet / Median_Lambda_Dry
  ) %>%
  dplyr::select(site, starts_with("effect_"))

treatment_df <- treatment_effects

effect_matrix <- treatment_effects %>%
  column_to_rownames("site") %>%
  as.matrix()

effect_matrix_scaled <- scale(effect_matrix)

# ===== Cube root transforming data before correlations =====
cube_root <- function(x) sign(x) * abs(x)^(1 / 3)

cube_field_means <- field_sample_data %>%
  filter(Mean_Fe_mg_per_kg < 10) %>%
  mutate(across(where(is.numeric), cube_root)) %>%
  rename_with(where(is.numeric), .fn = ~ paste0("cube_", .x)) %>%
  dplyr::select(-contains("per_L"))

cube_treatment_effects <- treatment_effects %>%
  mutate(across(where(is.numeric), cube_root)) %>%
  rename_with(where(is.numeric), .fn = ~ paste0("cube_", .x))

dat_field <- merge(cube_treatment_effects, cube_field_means, by = "site") %>%
  column_to_rownames("site")

# ==== Correlation matrix =====
cor_matrix <- cor(dat_field, use = "complete.obs", method = "pearson")
print(cor_matrix)

corrplot(
  cor_matrix,
  method = "color",
  type = "upper",
  order = "hclust",
  tl.cex = 1.2,
  tl.col = "black",
  tl.srt = 45
)

ggcorrplot(
  cor_matrix,
  hc.order = TRUE,
  type = "lower",
  lab = TRUE,
  lab_size = 5,
  method = "circle",
  colors = c("blue", "white", "red"),
  title = "Pearson Correlation Matrix",
  ggtheme = theme_bw(base_size = 18)
)

# ==== Creating color objects for heat map =====
site_annotation <- site_metadata %>%
  dplyr::select(
    site,
    State,
    Ecoregion,
    Mean_ATP_picomoles_per_g,
    Mean_Specific_Surface_Area_m2_per_g,
    C_percent_per_mg,
    N_percent_per_mg,
    Percent_Tot_Sand,
    Mean_Fe_mg_per_kg
  ) %>%
  column_to_rownames("site")

site_annotation <- site_annotation[rownames(effect_matrix), ]

categorical_colors <- list(
  Ecoregion = setNames(
    qualitative_hcl(
      n = length(unique(site_annotation$Ecoregion)),
      palette = "Dark 3"
    ),
    unique(site_annotation$Ecoregion)
  ),
  State = setNames(
    lighten(
      qualitative_hcl(length(unique(site_annotation$State)), palette = "Dynamic"),
      0.1
    ),
    unique(site_annotation$State)
  )
)

numerical_colors <- list(
  Mean_ATP_picomoles_per_g = colorRampPalette(c("#FFFFFF", "#E6E1F9", "#B197FC", "#845EF7"))(100),
  Mean_Specific_Surface_Area_m2_per_g = colorRampPalette(c("#FFFFFF", "#D0EBFF", "#74C0FC", "#1971C2"))(100),
  C_percent_per_mg = colorRampPalette(c("#FFFFFF", "#E9ECEF", "#ADB5BD", "#495057"))(100),
  N_percent_per_mg = colorRampPalette(c("#FFFFFF", "#D3F9D8", "#69DB7C", "#2B8A3E"))(100),
  Percent_Tot_Sand = colorRampPalette(c("#FFFFFF", "#FFF3BF", "#FFD43B", "#F08C00"))(100),
  Mean_Fe_mg_per_kg = colorRampPalette(c("#FFFFFF", "#FFDEEB", "#FF8FAB", "#E03131"))(100),
  Mean_Gravimetric_Moisture_g_per_g = colorRampPalette(c("#8D6E63", "#D7CCC8", "#81D4FA", "#0277BD"))(100)
)

# ===== Clustering based on the effect ======
site_dist <- dist(effect_matrix_scaled, method = "euclidean")
site_hclust <- hclust(site_dist, method = "ward.D2")

# fviz_nbclust(effect_matrix_scaled, hcut, method = "wss") +
#   labs(title = "Elbow Method")

site_clusters <- cutree(site_hclust, k = 3)
site_annotation$ResponseCluster <- factor(site_clusters)

categorical_colors$ResponseCluster <- setNames(
  cluster_colors,
  levels(factor(site_clusters))
)

annotation_colors <- c(categorical_colors, numerical_colors)

# ===== Visualize clusters ======
site_cluster_data <- data.frame(
  site = names(site_clusters),
  cluster = site_clusters
) %>%
  left_join(site_metadata, by = "site") %>%
  mutate(cluster = factor(cluster))

dom_cluster_data <- data.frame(
  site = names(site_clusters),
  cluster = site_clusters
) %>%
  left_join(treatment_df, by = "site") %>%
  mutate(cluster = factor(cluster))

# === Figure S3 boxplots: field variables within clusters ====

p44 <- ggplot(site_cluster_data, aes(x = factor(cluster), fill = Ecoregion)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = categorical_colors$Ecoregion) +
  theme_bw(base_size = 22) +
  theme(
    axis.text = element_text(size = 24, color = "black"),
    axis.title = element_text(size = 28, color = "black", face = "bold"),
    legend.text = element_text(size = 20, color = "black"),
    legend.title = element_text(size = 22, color = "black", face = "bold"),
    plot.title = element_text(size = 30, face = "bold")
  ) +
  labs(
    title = "A",
    x = "Cluster",
    y = "Proportion"
  )

p1 <- ggplot(site_cluster_data, aes(x = cluster, y = Mean_Specific_Surface_Area_m2_per_g)) +
  geom_boxplot(aes(fill = cluster), linewidth = 1.1, outlier.size = 2.8) +
  scale_fill_manual(values = cluster_colors) +
  pub_theme_s3 +
  labs(
    title = "A",
    x = "Cluster",
    y = expression(Specific~Surface~Area~(m^2~g^{-1})),
    fill = "Cluster"
  )

p2 <- ggplot(site_cluster_data, aes(x = cluster, y = C_percent_per_mg)) +
  geom_boxplot(aes(fill = cluster), linewidth = 1.1, outlier.size = 2.8) +
  scale_fill_manual(values = cluster_colors) +
  pub_theme_s3 +
  labs(
    title = "B",
    x = "Cluster",
    y = "Carbon (%)",
    fill = "Cluster"
  )

p3 <- ggplot(site_cluster_data, aes(x = cluster, y = Mean_ATP_picomoles_per_g)) +
  geom_boxplot(aes(fill = cluster), linewidth = 1.1, outlier.size = 2.8) +
  scale_fill_manual(values = cluster_colors) +
  pub_theme_s3 +
  labs(
    title = "C",
    x = "Cluster",
    y = expression(ATP~(picomoles~g^{-1})),
    fill = "Cluster"
  )

p4 <- ggplot(site_cluster_data, aes(x = cluster, y = Percent_Tot_Sand)) +
  geom_boxplot(aes(fill = cluster), linewidth = 1.1, outlier.size = 2.8) +
  scale_fill_manual(values = cluster_colors) +
  pub_theme_s3 +
  labs(
    title = "D",
    x = "Cluster",
    y = "Total Sand (%)",
    fill = "Cluster"
  )

p5 <- ggplot(site_cluster_data, aes(x = cluster, y = Mean_Fe_mg_per_kg)) +
  geom_boxplot(aes(fill = cluster), linewidth = 1.1, outlier.size = 2.8) +
  scale_fill_manual(values = cluster_colors) +
  pub_theme_s3 +
  labs(
    title = "E",
    x = "Cluster",
    y = expression(Fe~(mg~kg^{-1})),
    fill = "Cluster"
  )

p6 <- ggplot(site_cluster_data, aes(x = cluster, y = N_percent_per_mg)) +
  geom_boxplot(aes(fill = cluster), linewidth = 1.1, outlier.size = 2.8) +
  scale_fill_manual(values = cluster_colors) +
  pub_theme_s3 +
  labs(
    title = "F",
    x = "Cluster",
    y = "Nitrogen (%)",
    fill = "Cluster"
  )

# ===== Larger statistical annotation function for Figure S3 =====
add_kw_stats <- function(
    plot,
    data,
    y_var,
    group_var = "cluster",
    kw_text_size = 8,
    letter_text_size = 9
) {
  y_data <- data[[y_var]]
  
  kw_formula <- as.formula(paste(y_var, "~", group_var))
  kw_result <- kruskal.test(kw_formula, data = data)
  
  p_value_formatted <- ifelse(
    kw_result$p.value < 0.001,
    "p < 0.001",
    paste0("p = ", signif(kw_result$p.value, 3))
  )
  
  y_max <- max(y_data, na.rm = TRUE)
  y_min <- min(y_data, na.rm = TRUE)
  y_range <- y_max - y_min
  
  stat_y <- y_max + 0.22 * y_range
  letter_y <- y_max + 0.07 * y_range
  upper_y <- y_max + 0.32 * y_range
  
  if (kw_result$p.value < 0.05) {
    posthoc_result <- kwAllPairsDunnTest(
      kw_formula,
      data = data,
      p.adjust.method = "bonf"
    )
    
    pvalues <- posthoc_result$p.value
    diag(pvalues) <- 1
    
    letters_result <- multcompLetters(
      pvalues,
      compare = "<",
      threshold = 0.05,
      Letters = letters
    )
    
    letter_df <- data.frame(
      cluster = names(letters_result$Letters),
      letter = unname(letters_result$Letters),
      stringsAsFactors = FALSE
    )
    
    enhanced_plot <- plot +
      annotate(
        "text",
        x = 2,
        y = stat_y,
        label = paste0("K-W H = ", round(as.numeric(kw_result$statistic), 1), ", ", p_value_formatted),
        size = kw_text_size,
        fontface = "bold",
        color = "black"
      ) +
      geom_text(
        data = letter_df,
        aes(x = cluster, y = letter_y, label = letter),
        inherit.aes = FALSE,
        size = letter_text_size,
        fontface = "bold",
        color = "black"
      )
  } else {
    enhanced_plot <- plot +
      annotate(
        "text",
        x = 2,
        y = stat_y,
        label = paste0("K-W H = ", round(as.numeric(kw_result$statistic), 1), ", ", p_value_formatted),
        size = kw_text_size,
        fontface = "bold",
        color = "black"
      )
  }
  
  enhanced_plot <- enhanced_plot +
    coord_cartesian(ylim = c(NA, upper_y), clip = "off")
  
  return(enhanced_plot)
}

p1_enhanced <- add_kw_stats(p1, site_cluster_data, "Mean_Specific_Surface_Area_m2_per_g")
p2_enhanced <- add_kw_stats(p2, site_cluster_data, "C_percent_per_mg")
p3_enhanced <- add_kw_stats(p3, site_cluster_data, "Mean_ATP_picomoles_per_g")
p4_enhanced <- add_kw_stats(p4, site_cluster_data, "Percent_Tot_Sand")
p5_enhanced <- add_kw_stats(p5, site_cluster_data, "Mean_Fe_mg_per_kg")
p6_enhanced <- add_kw_stats(p6, site_cluster_data, "N_percent_per_mg")

grid_arranged_plots <- grid.arrange(
  p1_enhanced, p2_enhanced, p3_enhanced,
  p4_enhanced, p5_enhanced, p6_enhanced,
  ncol = 2
)

ggsave(
  filename = "Figures/FigureS3_Field_variables_within_clusters.pdf",
  plot = grid_arranged_plots,
  width = 14,
  height = 18,
  units = "in",
  dpi = 300,
  device = cairo_pdf
)

ggsave(
  filename = "Figures/FigureS3_Field_variables_within_clusters.png",
  plot = grid_arranged_plots,
  width = 14,
  height = 18,
  units = "in",
  dpi = 600
)

# ==== Heat map plot Field data =====
pdf("Figures/Figure2_Heatmap_Field_data_and_Clusters.pdf", width = 14, height = 10)
pheatmap(
  effect_matrix_scaled,
  cluster_rows = site_hclust,
  cluster_cols = TRUE,
  clustering_method = "ward.D2",
  color = colorRampPalette(c("#3D5A80", "#FFFFFF", "#E07A5F"))(100),
  breaks = seq(-2.5, 2.5, length.out = 101),
  fontsize_number = 8,
  angle_col = 45,
  annotation_row = site_annotation,
  annotation_colors = annotation_colors,
  cutree_rows = 3,
  annotation_legend = TRUE,
  legend = TRUE,
  fontsize_row = 8,
  cellwidth = 15,
  cellheight = 12,
  border_color = NA
)
dev.off()

png("Figures/Figure2_Heatmap_Field_data_and_Clusters.png", width = 4200, height = 3000, res = 300)
pheatmap(
  effect_matrix_scaled,
  cluster_rows = site_hclust,
  cluster_cols = TRUE,
  clustering_method = "ward.D2",
  color = colorRampPalette(c("#3D5A80", "#FFFFFF", "#E07A5F"))(100),
  breaks = seq(-2.5, 2.5, length.out = 101),
  fontsize_number = 8,
  angle_col = 45,
  annotation_row = site_annotation,
  annotation_colors = annotation_colors,
  cutree_rows = 3,
  annotation_legend = TRUE,
  legend = TRUE,
  fontsize_row = 8,
  cellwidth = 15,
  cellheight = 12,
  border_color = NA
)
dev.off()

# === DOM effects by cluster ====
kruskal_gibbs <- kruskal.test(effect_delGcoxPerCmol ~ cluster, data = dom_cluster_data)
dunn_gibbs <- dunn.test(
  dom_cluster_data$effect_delGcoxPerCmol,
  dom_cluster_data$cluster,
  method = "bonferroni"
)

kruskal_lambda <- kruskal.test(effect_Lambda ~ cluster, data = dom_cluster_data)
dunn_lambda <- dunn.test(
  dom_cluster_data$effect_Lambda,
  dom_cluster_data$cluster,
  method = "bonferroni"
)

effectGibbs <- ggplot(dom_cluster_data, aes(y = effect_delGcoxPerCmol, x = cluster, fill = cluster)) +
  geom_boxplot(linewidth = 1) +
  scale_fill_manual(values = cluster_colors) +
  annotate(
    "text",
    x = 1.2,
    y = max(dom_cluster_data$effect_delGcoxPerCmol, na.rm = TRUE) + 0.005,
    label = "K-W p < 0.001",
    size = 7,
    fontface = "bold"
  ) +
  annotate("text", x = 1, y = max(dom_cluster_data$effect_delGcoxPerCmol, na.rm = TRUE), label = "ac", size = 7, fontface = "bold") +
  annotate("text", x = 2, y = max(dom_cluster_data$effect_delGcoxPerCmol, na.rm = TRUE) + 0.005, label = "ab", size = 7, fontface = "bold") +
  annotate("text", x = 3, y = max(dom_cluster_data$effect_delGcoxPerCmol, na.rm = TRUE) - 0.005, label = "bc", size = 7, fontface = "bold") +
  theme_bw(base_size = 20) +
  theme(
    aspect.ratio = 1,
    legend.position = "bottom",
    axis.text = element_text(size = 18, color = "black"),
    axis.title = element_text(size = 22, color = "black", face = "bold"),
    legend.text = element_text(size = 18, color = "black"),
    legend.title = element_text(size = 20, color = "black", face = "bold"),
    plot.title = element_text(size = 24, color = "black", face = "bold"),
    legend.key.size = unit(1.2, "cm")
  ) +
  labs(
    x = "",
    y = expression(paste(Effect~Size~Delta, G[cox], ~(kJ~Cmol^{-1}))),
    fill = "Cluster",
    title = "B"
  )

effectlambda <- ggplot(dom_cluster_data, aes(y = effect_Lambda, x = cluster, fill = cluster)) +
  geom_boxplot(linewidth = 1) +
  scale_fill_manual(values = cluster_colors) +
  annotate(
    "text",
    x = 1.2,
    y = max(dom_cluster_data$effect_Lambda, na.rm = TRUE) + 0.05,
    label = "K-W p < 0.001",
    size = 7,
    fontface = "bold"
  ) +
  annotate("text", x = 1, y = max(dom_cluster_data$effect_Lambda, na.rm = TRUE) - 0.02, label = "a", size = 7, fontface = "bold") +
  annotate("text", x = 2, y = max(dom_cluster_data$effect_Lambda, na.rm = TRUE) - 0.02, label = "b", size = 7, fontface = "bold") +
  annotate("text", x = 3, y = max(dom_cluster_data$effect_Lambda, na.rm = TRUE) + 0.02, label = "ab", size = 7, fontface = "bold") +
  theme_bw(base_size = 20) +
  theme(
    aspect.ratio = 1,
    legend.position = "bottom",
    axis.text = element_text(size = 18, color = "black"),
    axis.title = element_text(size = 22, color = "black", face = "bold"),
    legend.text = element_text(size = 18, color = "black"),
    legend.title = element_text(size = 20, color = "black", face = "bold"),
    plot.title = element_text(size = 24, color = "black", face = "bold"),
    legend.key.size = unit(1.2, "cm")
  ) +
  labs(
    x = "",
    y = expression(paste(Effect~Size~lambda)),
    fill = "Cluster",
    title = "C"
  )

print("ΔGcox Kruskal-Wallis Test:")
print(kruskal_gibbs)
print("ΔGcox Dunn's Test:")
print(dunn_gibbs)

print("Lambda Kruskal-Wallis Test:")
print(kruskal_lambda)
print("Lambda Dunn's Test:")
print(dunn_lambda)

combined_plot <- effectGibbs + effectlambda +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  "Figures/Figure2_effect_sizes_by_cluster.png",
  combined_plot,
  width = 10,
  height = 5,
  dpi = 300
)

ggsave(
  "Figures/Figure2_effect_sizes_by_cluster.pdf",
  combined_plot,
  width = 10,
  height = 5,
  device = cairo_pdf
)

# ==== Perform a Linear Discriminant Analysis for Field Data ====
site_lda_data <- data.frame(
  site = names(site_clusters),
  cluster = factor(site_clusters),
  stringsAsFactors = FALSE
) %>%
  left_join(site_metadata, by = "site") %>%
  dplyr::select(-Ecoregion) %>%
  na.omit()

cat("Cluster distribution:\n")
print(table(site_lda_data$cluster))
cat("\nNumber of sites per cluster:\n")
print(site_lda_data %>% count(cluster))

lda_model <- lda(
  cluster ~
    Mean_ATP_picomoles_per_g +
    Mean_Specific_Surface_Area_m2_per_g +
    C_percent_per_mg +
    N_percent_per_mg +
    Percent_Tot_Sand +
    Mean_Fe_mg_per_kg,
  data = site_lda_data
)

cat("\n===== LDA RESULTS =====\n")
print(lda_model)

lda_predictions <- predict(lda_model, site_lda_data)

confusion_matrix <- table(
  Predicted = lda_predictions$class,
  Actual = site_lda_data$cluster
)

cat("\n===== CONFUSION MATRIX =====\n")
print(confusion_matrix)

accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
cat("\nOverall Classification Accuracy:", round(accuracy * 100, 1), "%\n")

lda_coefficients <- as.data.frame(lda_model$scaling)
lda_coefficients$Variable <- rownames(lda_coefficients)

cat("\n===== DISCRIMINANT FUNCTION COEFFICIENTS =====\n")
print(lda_coefficients)

lda_plot_data <- data.frame(
  site = site_lda_data$site,
  cluster = site_lda_data$cluster,
  LD1 = lda_predictions$x[, 1],
  LD2 = lda_predictions$x[, 2]
)

lda_stats <- list(
  prop_var = lda_model$svd^2 / sum(lda_model$svd^2),
  accuracy = sum(lda_predictions$class == site_lda_data$cluster) / nrow(site_lda_data)
)

lda_plot <- ggplot(lda_plot_data, aes(x = LD1, y = LD2, color = cluster, fill = cluster)) +
  stat_ellipse(level = 0.95, alpha = 0.2, linewidth = 1) +
  geom_point(size = 3.5, alpha = 0.8) +
  scale_color_manual(values = cluster_colors) +
  scale_fill_manual(values = cluster_colors) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = sprintf("Accuracy: %.1f%%", lda_stats$accuracy * 100),
    hjust = 1.1,
    vjust = 1.1,
    size = 5.5,
    color = "black"
  ) +
  labs(
    x = sprintf("LD1 (%.1f%% of variance)", lda_stats$prop_var[1] * 100),
    y = sprintf("LD2 (%.1f%% of variance)", lda_stats$prop_var[2] * 100),
    color = "Response\nCluster",
    fill = "Response\nCluster",
    title = "D"
  ) +
  theme_bw(base_size = 18) +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 20, color = "black", face = "bold"),
    legend.text = element_text(size = 16, color = "black"),
    legend.title = element_text(size = 18, color = "black", face = "bold"),
    plot.title = element_text(size = 24, face = "bold"),
    legend.key.size = unit(1.2, "cm")
  )

ggsave("Figures/Figure_2_LDA_discriminant_space.png", plot = lda_plot, width = 8, height = 6, dpi = 300)
ggsave("Figures/Figure_2_LDA_discriminant_space.pdf", plot = lda_plot, width = 8, height = 6, dpi = 300)

var_importance <- lda_coefficients %>%
  mutate(
    LD1_abs = abs(LD1),
    LD2_abs = abs(LD2),
    Total_importance = sqrt(LD1^2 + LD2^2)
  ) %>%
  arrange(desc(Total_importance))

importance_plot <- ggplot(var_importance, aes(x = reorder(Variable, Total_importance), y = Total_importance)) +
  geom_col(fill = "steelblue", alpha = 0.7) +
  coord_flip() +
  labs(
    title = "Variable Importance in Cluster Discrimination",
    x = "Environmental Variables",
    y = "Discriminant Importance"
  ) +
  theme_bw(base_size = 18) +
  theme(
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 20, color = "black", face = "bold"),
    plot.title = element_text(size = 22, face = "bold")
  )

ggsave("Figures/LDA_variable_importance.png", plot = importance_plot, width = 8, height = 6, dpi = 300)
ggsave("Figures/LDA_variable_importance.pdf", plot = importance_plot, width = 8, height = 6, dpi = 300)

top_vars <- head(var_importance$Variable, 3)

plot_list <- list()

for (i in 1:3) {
  var_name <- top_vars[i]
  
  p <- ggplot(site_lda_data, aes_string(x = "cluster", y = var_name, fill = "cluster")) +
    geom_boxplot(alpha = 0.7, linewidth = 1) +
    scale_fill_manual(values = cluster_colors) +
    labs(
      title = paste("Cluster Differences:", var_name),
      x = "Response Cluster",
      y = var_name
    ) +
    theme_bw(base_size = 18) +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 16, color = "black"),
      axis.title = element_text(size = 20, color = "black", face = "bold"),
      plot.title = element_text(size = 20, face = "bold")
    )
  
  plot_list[[i]] <- p
}

combined_boxplots <- grid.arrange(grobs = plot_list, ncol = 3)

ggsave(
  "Figures/LDA_top_variables_boxplots.png",
  plot = combined_boxplots,
  width = 12,
  height = 4,
  dpi = 300
)

ggsave(
  "Figures/LDA_top_variables_boxplots.pdf",
  plot = combined_boxplots,
  width = 12,
  height = 4,
  dpi = 300
)

lda_results <- list(
  model = lda_model,
  predictions = lda_predictions,
  accuracy = accuracy,
  confusion_matrix = confusion_matrix,
  variable_importance = var_importance
)

saveRDS(lda_results, "Figures/lda_cluster_analysis.rds")