rm(list=ls(all=T))

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
if(!require(corrplot)) install.packages("corrplot", dependencies = TRUE)
library(corrplot)
library(Hmisc)
library(stringr)

if(!require(factoextra)) install.packages("factoextra", dependencies = TRUE)
library(factoextra)

if(!require(MASS)) install.packages("MASS", dependencies = TRUE)
library(MASS)

if(!require(gridExtra)) install.packages("gridExtra", dependencies = TRUE)
library(gridExtra)

if(!require(ggcorrplot)) install.packages("ggcorrplot", dependencies = TRUE)
library(ggcorrplot)

if(!require(cowplot)) install.packages("cowplot", dependencies = TRUE)
library(cowplot)

# ==== Defining paths and working directories ======
figure_path = "Figures/"

# ====== Read in data ======
data = read.csv("Data/Medians_of_Median_molecular_properties_per_site_and_treatment_unique_formulas.csv")
row.names(data) = paste0(data$site, "_", data$Treatment)

large_data = read.csv("Data/EC_Field_Geospatial_NPP_PET_VGC.csv")
# NOTE: keeping your change
names(large_data)[2] = "code"
names(large_data)[1] = "site"

# ==== Set up data =====
dom_data = data %>%
  dplyr::filter(site != "EC_023") %>%
  dplyr::filter(!(site %in% c("EC_012","EC_011","EC_053","EC_057","EC_052"))) %>%
  dplyr::select(site, Treatment, Median_delGcoxPerCmol, Median_Lambda)

sites = as.data.frame(unique(dom_data$site))
names(sites) = "site"
explanatory_data = merge(sites, large_data, by = "site") %>%
  distinct()

site_metadata = explanatory_data

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
    effect_delGcoxPerCmol = (Median_delGcoxPerCmol_Wet / Median_delGcoxPerCmol_Dry),
    effect_Lambda        = (Median_Lambda_Wet        / Median_Lambda_Dry)
  ) %>%
  dplyr::select(site, starts_with("effect_"))

effect_matrix <- treatment_effects %>%
  column_to_rownames("site") %>%
  as.matrix()

effect_matrix_scaled <- scale(effect_matrix)

# ===== Cube root transforming data before correlations =====
cube_root <- function(x) sign(x) * (abs(x))^(1/3)

cube_field_means = explanatory_data %>%
  mutate(across(where(is.numeric), cube_root)) %>%
  rename_with(where(is.numeric), .fn = ~ paste0("cube_", .x))

cube_treatment_effects = treatment_effects %>%
  mutate(across(where(is.numeric), cube_root)) %>%
  rename_with(where(is.numeric), .fn = ~ paste0("cube_", .x))

dat_field = merge(cube_treatment_effects, cube_field_means, by = "site") %>%
  column_to_rownames("site") %>%
  dplyr::select(-code)


# ==== Correlation matrix =====
cor_matrix <- cor(dat_field, use = "complete.obs", method = "pearson")
print(cor_matrix)

pretty_names <- colnames(cor_matrix) %>%
  stringr::str_remove("^cube_") %>%     # drop cube_ prefix
  stringr::str_remove("^effect_") %>%   # drop effect_ prefix (if present)
  stringr::str_replace_all("_", " ") %>%# underscores -> spaces
  stringr::str_wrap(width = 10)         # wrap long labels onto multiple lines

dimnames(cor_matrix) <- list(pretty_names, pretty_names)

# Plot A 
pA <- ggcorrplot(
  cor_matrix,
  method = "square",
  type = "upper",
  hc.order = TRUE,
  lab = TRUE
) +
  theme_bw() +
  theme(
    plot.title = element_blank(),
    
    axis.text.x = element_text(
      size = 7,
      angle = 45,      # try 60 or 90 if still tight
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(size = 7),
    
    plot.margin = margin(t = 5, r = 5, b = 25, l = 5)
  )

ggsave("Figures/Figure_S3_correlation_plot_A.png", pA, width = 8, height = 6, dpi = 300)
ggsave("Figures/Figure_S3_correlation_plot_A.pdf", pA, width = 8, height = 6, dpi = 300)

# ===== Clustering based on the effect ======
site_dist <- dist(effect_matrix_scaled, method = "euclidean")
site_hclust <- hclust(site_dist, method = "ward.D2")

fviz_nbclust(effect_matrix_scaled, hcut, method = "wss") +
  labs(title = "Elbow Method")

site_clusters <- cutree(site_hclust, k = 3)

# ==== LDA ====
site_lda_data <- data.frame(
  site = names(site_clusters),
  cluster = factor(site_clusters),
  stringsAsFactors = FALSE
) %>%
  left_join(explanatory_data, by = "site") %>%
  na.omit()

cat("Cluster distribution:\n")
print(table(site_lda_data$cluster))
cat("\nNumber of sites per cluster:\n")
print(site_lda_data %>% count(cluster))

lda_model <- lda(cluster ~
                   slope +
                   elevation +
                   AridityWs +
                   Pct_shrub +
                   precipitation +
                   PctFst +
                   PctAg +
                   NPP_mean_kgC_m2_yr +
                   ET_mean_mm_yr+
                   PET_mean_mm_yr,
                 data = site_lda_data)

cat("\n===== LDA RESULTS =====\n")
print(lda_model)

lda_predictions <- predict(lda_model, site_lda_data)

confusion_matrix <- table(Predicted = lda_predictions$class, Actual = site_lda_data$cluster)
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
  LD1 = lda_predictions$x[,1],
  LD2 = lda_predictions$x[,2]
)

lda_stats <- list(
  prop_var  = lda_model$svd^2 / sum(lda_model$svd^2),
  accuracy  = sum(lda_predictions$class == site_lda_data$cluster) / nrow(site_lda_data)
)

# Plot B (LDA)
pB <- ggplot(lda_plot_data, aes(x = LD1, y = LD2, color = cluster, fill = cluster)) +
  stat_ellipse(level = 0.95, alpha = 0.2, size = 0.9) +
  geom_point(size = 2.6, alpha = 0.75) +
  scale_color_manual(values = c("1" = "#74C0FC", "2" = "#FFD43B", "3" = "#FF8FAB")) +
  scale_fill_manual(values  = c("1" = "#74C0FC", "2" = "#FFD43B", "3" = "#FF8FAB")) +
  annotate(
    "text", x = Inf, y = Inf,
    label = sprintf("Accuracy: %.1f%%", lda_stats$accuracy * 100),
    hjust = 1.08, vjust = 1.2, size = 3.2,
    color = "black"
  ) +
  labs(
    x = sprintf("LD1 (%.1f%%)", lda_stats$prop_var[1] * 100),
    y = sprintf("LD2 (%.1f%%)", lda_stats$prop_var[2] * 100),
    color = "Response\nCluster",
    fill  = "Response\nCluster"
  ) +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    legend.position = "right",
    axis.text  = element_text(size = 11),
    axis.title = element_text(size = 12),
    legend.text  = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.key.size = unit(0.7, "cm")
  )

ggsave("Figures/Figure_S3_LDA_discriminant_space.png", plot = pB, width = 8, height = 6, dpi = 300)
ggsave("Figures/Figure_S3_LDA_discriminant_space.pdf", plot = pB, width = 8, height = 6, dpi = 300)

var_importance <- lda_coefficients %>%
  mutate(
    Total_importance = sqrt(LD1^2 + LD2^2)
  ) %>%
  arrange(desc(Total_importance))

# ====== Boxplots =====

vars_to_plot <- site_lda_data %>%
  dplyr::select(where(is.numeric)) %>%
  dplyr::select(-any_of(c("code"))) %>%
  names()

# Ensure cluster ordering is stable
site_lda_data <- site_lda_data %>%
  mutate(cluster = factor(cluster, levels = sort(unique(as.integer(as.character(cluster))))))

# Compute Kruskal–Wallis stats for every variable
kw_table <- purrr::map_dfr(vars_to_plot, function(v) {
  d <- site_lda_data %>%
    dplyr::select(cluster, all_of(v)) %>%
    tidyr::drop_na()
  
  if (n_distinct(d$cluster) < 2 || n_distinct(d[[v]]) < 2) {
    return(tibble::tibble(variable = v, H = NA_real_, p_value = NA_real_, n = nrow(d)))
  }
  
  kt <- kruskal.test(d[[v]] ~ d$cluster)
  
  tibble::tibble(
    variable = v,
    H = as.numeric(kt$statistic),
    p_value = as.numeric(kt$p.value),
    n = nrow(d)
  )
})

print(kw_table)
write.csv(kw_table, file = "Figures/S3_KruskalWallis_stats_all_variables.csv", row.names = FALSE)

plot_list <- vector("list", length(vars_to_plot))

for (i in seq_along(vars_to_plot)) {
  var_name <- vars_to_plot[i]
  stat_row <- kw_table %>% dplyr::filter(variable == var_name)
  
  # One-line annotation
  if (is.na(stat_row$H) || is.na(stat_row$p_value)) {
    lab <- "K-W: not testable"
  } else {
    lab <- sprintf(
      "K-W H = %.2f, p = %s",
      stat_row$H,
      format.pval(stat_row$p_value, digits = 2, eps = 1e-3)
    )
  }
  
  # Wrap long y-axis labels so they don't blow up the layout
  ylab_wrapped <- stringr::str_wrap(var_name, width = 18)
  
  p_box <- ggplot(site_lda_data, aes(x = cluster, y = .data[[var_name]], fill = cluster)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.35, linewidth = 0.25) +
    scale_fill_manual(values = c("1" = "#74C0FC", "2" = "#FFD43B", "3" = "#FF8FAB")) +
    labs(
      x = NULL,
      y = ylab_wrapped
    ) +
    annotate(
      "text",
      x = Inf, y = Inf,
      label = lab,
      hjust = 1.05, vjust = 1.2,
      size = 2.6,
      color = "black"
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_blank(),
      axis.text.x  = element_text(size = 7),
      axis.text.y  = element_text(size = 7),
      axis.title.y = element_text(size = 7),
      axis.title.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )
  
  plot_list[[i]] <- p_box
}

# ==== Export figure =============================================================

AB_panel <- cowplot::plot_grid(
  pA, pB,
  labels = c("A", "B"),
  ncol = 2,
  label_size = 18,
  rel_widths = c(1.25, 1)   # was ~1.05; widen A
)



ncol_C <- 4

C_panel <- cowplot::plot_grid(
  plotlist = plot_list,
  ncol = ncol_C,
  align = "hv"
)


C_panel_labeled <- cowplot::ggdraw(C_panel) +
  cowplot::draw_label("C", x = 0.01, y = 0.99, hjust = 0, vjust = 1,
                      fontface = "bold", size = 18)


final_single_page <- cowplot::plot_grid(
  AB_panel,
  C_panel_labeled,
  ncol = 1,
  rel_heights = c(1, 1.7)
)


ggsave(
  filename = "Figures/Figure_S3_SINGLE_PAGE_COMPOSITE_A_B_C.pdf",
  plot = final_single_page,
  width = 16, height = 20, units = "in", dpi = 300
)


ggsave(
  filename = "Figures/Figure_S3_SINGLE_PAGE_COMPOSITE_A_B_C.png",
  plot = final_single_page,
  width = 16, height = 20, units = "in", dpi = 300
)

# ===== Save results =====
lda_results <- list(
  model = lda_model,
  predictions = lda_predictions,
  accuracy = accuracy,
  confusion_matrix = confusion_matrix,
  variable_importance = var_importance,
  kw_table = kw_table
)

saveRDS(lda_results, "Figures/S3_lda_cluster_analysis.rds")
