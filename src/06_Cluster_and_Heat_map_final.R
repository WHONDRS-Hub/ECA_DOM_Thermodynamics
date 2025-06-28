# Load required libraries
rm(list=ls(all=T))

library(tidyverse)
library(ggplot2)
library(pheatmap)
library(viridis)
library(dendextend) # For enhanced dendrograms
library(dplyr)
library(viridis)
library(RColorBrewer)
library(pheatmap)
library(tibble)
library(grid)
library(colorspace)
if(!require(corrplot)) install.packages("corrplot")
library(corrplot)
library(Hmisc)
library(stringr)
# ==== Defining paths and working directories ======
github = 'C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/ECA_DOM_Thermodynamics/'
data_path = paste0(github,'Data/')
figure_path = paste0(github,'Figures/')

# ====== Read in data ======
data = read.csv(paste0(data_path,'Medians_of Median_molecular_properties_per_site_and_treatment_unique_formulas.csv'))

row.names(data) = paste0(data$site,'_',data$Treatment)

sample_data = read_csv(paste0(github,'v4_CM_SSS_Data_Package/Sample_Data/v3_CM_SSS_Sediment_Sample_Data_Summary.csv'),comment = '#', na = c('N/A', -9999)) %>%
  slice(-(1:11)) %>%  # Remove the first 11 rows
  mutate_at(vars(-Sample_Name, -Field_Name, -IGSN, -Material), as.numeric) %>%  
  dplyr::select(site = Sample_Name, Mean_ATP_picomoles_per_g,
  Mean_Specific_Surface_Area_m2_per_g,
  C_percent_per_mg = '01395_C_percent_per_mg',                
  N_percent_per_mg = '01397_N_percent_per_mg',
  Percent_Tot_Sand, Mean_Fe_mg_per_kg,
  Mean_Gravimetric_Moisture_g_per_g = "Mean_62948_Gravimetric_Moisture_g_per_g")

sample_data$site = gsub('CM','EC',sample_data$site)
sample_data$site = gsub('_Sediment','',sample_data$site)

field_metadata = read.csv(paste0(github,'EC_Data_Package/EC_Field_Metadata.csv')) %>%
  dplyr::select(site = Parent_ID,State)

ecoregions = read.csv('Data/EC_Ecoregions.csv')%>%
  dplyr::select(site = 'Parent_ID', Ecoregion = Name)%>%
  mutate(Ecoregion = str_trim(str_extract(Ecoregion, "(?<=- ).*")))

field_metadata = merge(field_metadata,ecoregions, by = 'site')

# ==== Set up data =====
dom_data = data %>%
  dplyr::filter(site != "EC_023") %>% 
  dplyr::filter(!(site %in% c("EC_012","EC_011","EC_053","EC_057","EC_052"))) %>%
  #Remove EC_011, EC_012, EC_023, EC_052, EC_053, and EC_057 for having too much water added (11 and 12), no mg/kg calculation (23), and being duplicated NEON sites (52, 53, and 57)
  dplyr::select(site,Treatment,Median_delGcoxPerCmol,Median_Lambda)

# Doing this to ensure my explanatory variables are from the same sites as my dom variables
sites = as.data.frame(unique(dom_data$site))
names(sites) = 'site'
explanatory_data = merge(sample_data,field_metadata, by = 'site')
explanatory_data = merge(explanatory_data,sites, by = 'site')

field_sample_data = merge(sample_data, sites, by = 'site')
site_metadata = explanatory_data

# ===== Calculations ======
# Calculate treatment effects size
treatment_effects <- dom_data %>%
  filter(Treatment %in% c("Wet", "Dry")) %>%
  select(site, Treatment, starts_with("Median_")) %>%
  pivot_wider(
    id_cols = site,
    names_from = Treatment,
    values_from = starts_with("Median_")
  ) %>%
  mutate(
    effect_delGcoxPerCmol = (Median_delGcoxPerCmol_Wet / Median_delGcoxPerCmol_Dry),
    effect_Lambda = (Median_Lambda_Wet / Median_Lambda_Dry)
  ) %>%
  select(site, starts_with("effect_"))

treatment_df = treatment_effects
# Prepare the effect matrix for heatmap
effect_matrix <- treatment_effects %>%
  column_to_rownames("site") %>%
  as.matrix()

effect_matrix_scaled <- scale(effect_matrix) # z-score normalization for clustering. Clustering uses Eucledean distance which is sensitive to scale

# ===== Cube root transforming data before correlations =====
cube_root <- function(x) sign(x) * (abs(x))^(1/3)

cube_field_means = field_sample_data %>% 
  filter(Mean_Fe_mg_per_kg < 10) %>%  # remove Fe outlier for analysis site EC_068
  mutate(across(where(is.numeric), cube_root)) %>% # cube root transform data
  rename_with(where(is.numeric), .fn = ~ paste0("cube_", .x)) %>% 
  #column_to_rownames("Sample_Name") %>%
  select(-contains("per_L")) # remove per_L data, analysis ran on _per_kg

cube_treatment_effects = treatment_effects %>%  
  mutate(across(where(is.numeric), cube_root)) %>% # cube root transform data
  rename_with(where(is.numeric), .fn = ~ paste0("cube_", .x))

dat_field = merge(cube_treatment_effects,cube_field_means, by = 'site') %>%
  column_to_rownames('site')

# ==== Creating color objects for heat map =====
# For field data
site_annotation <- site_metadata %>% 
  select(
    site, 
    Ecoregion,
    State,
    Mean_ATP_picomoles_per_g,
    Mean_Specific_Surface_Area_m2_per_g,
    C_percent_per_mg,
    N_percent_per_mg,
    Percent_Tot_Sand,
    Mean_Fe_mg_per_kg,
    Mean_Gravimetric_Moisture_g_per_g
  ) %>%
  column_to_rownames("site")

# Ensure annotations match the order of sites in the heatmap
site_annotation <- site_annotation[rownames(effect_matrix), ]

# Define categorical palettes with improved colors

library(colorspace)

categorical_colors <- list(
  Ecoregion = setNames(
    qualitative_hcl(n = length(unique(site_annotation$Ecoregion)), 
                    palette = "Dark 3"),  # or "Set 2", "Cold", "Warm"
    unique(site_annotation$Ecoregion)),
  # State - use a brighter palette
  State = setNames(
    lighten(qualitative_hcl(length(unique(site_annotation$State)), palette = "Dynamic"), 0.1),
    unique(site_annotation$State)
  )
)

# Create custom gradients with better visibility
numerical_colors <- list(
  # ATP - light to medium purple
  Mean_ATP_picomoles_per_g = colorRampPalette(c("#FFFFFF", "#E6E1F9", "#B197FC", "#845EF7"))(100),
  
  # SSA - light blue to dark blue gradient
  Mean_Specific_Surface_Area_m2_per_g = colorRampPalette(c("#FFFFFF", "#D0EBFF", "#74C0FC", "#1971C2"))(100),
  
  # Carbon content (improved grayscale)
  C_percent_per_mg = colorRampPalette(c("#FFFFFF", "#E9ECEF", "#ADB5BD", "#495057"))(100),
  
  # Nitrogen content (light green to dark green)
  `N_percent_per_mg` = colorRampPalette(c("#FFFFFF", "#D3F9D8", "#69DB7C", "#2B8A3E"))(100),
  
  # Sand content (light yellow to orange)
  Percent_Tot_Sand = colorRampPalette(c("#FFFFFF", "#FFF3BF", "#FFD43B", "#F08C00"))(100),
  
  # Iron content (light pink to red)
  Mean_Fe_mg_per_kg = colorRampPalette(c("#FFFFFF", "#FFDEEB", "#FF8FAB", "#E03131"))(100),
  
  Mean_Gravimetric_Moisture_g_per_g = colorRampPalette(c("#8D6E63", "#D7CCC8", "#81D4FA", "#0277BD"))(100)
)


# ===== Clustering based on the effect ======
# Perform hierarchical clustering to group similar sites
site_dist <- dist(effect_matrix_scaled, method = "euclidean")
site_hclust <- hclust(site_dist, method = "ward.D2")

# Determine the number of clusters visually using the elbow method
library(factoextra)
# Elbow method
fviz_nbclust(effect_matrix_scaled, hcut, method = "wss") +
  labs(title = "Elbow Method")

# Based on the plot it seems the k = 3 is a good number of clusters
# Cut the dendrogram to create site clusters
site_clusters <- cutree(site_hclust, k = 3)
site_annotation$ResponseCluster <- factor(site_clusters)

# Add ResponseCluster to categorical_colors
categorical_colors$ResponseCluster <- setNames(
  c("#74C0FC", "#FFD43B", "#FF8FAB"),  # Light blue, yellow, pink
  levels(factor(site_clusters))
)

# Combine color lists
annotation_colors <- c(categorical_colors, numerical_colors)
# ===== Visualize clusters ======
site_cluster_data <- data.frame(
  site = names(site_clusters),
  cluster = site_clusters
) %>%
  left_join(
    site_metadata,
    by = "site"
  )

dom_cluster_data <- data.frame(
  site = names(site_clusters),
  cluster = site_clusters
) %>%
  left_join(
    treatment_df,
    by = "site"
  )


# === Stacked ploys and Boxplots clusters ====
# Visualize cluster composition by Region
p44 = ggplot(site_cluster_data, aes(x = factor(cluster), fill = Ecoregion)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = categorical_colors$Ecoregion) +
  theme_bw() +
  labs(title = "A",
       x = "Cluster",
       y = "Proportion")

p45 = ggplot(site_cluster_data, aes(x = factor(cluster), fill = State)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = categorical_colors$State) +
  theme_bw() +
  labs(title = 'B',
       x = "Cluster",
       y = "Proportion")
library(gridExtra)
# Arrange all plots in a grid - adjust layout as needed
arranged_plots <- grid.arrange(
  p44, p45,
  ncol = 1
)


# Save as PDF and PNG
ggsave(
  filename = "Figures/FigureS2_Ecoregion_and_State_clusters.pdf",
  plot = arranged_plots,
  width = 10,
  height = 14,
  units = "in",
  dpi = 300
)

ggsave(
  filename = "Figures/FigureS2_Ecoregion_and_State_clusters.png",
  plot = arranged_plots,
  width = 10,
  height = 14,
  units = "in",
  dpi = 300
)

# Visualize numerical variables by clusters with higher Pearson correlation coefficients to the effects using field data 

# SSA
p1 <- ggplot(site_cluster_data, aes(x = factor(cluster), y = Mean_Specific_Surface_Area_m2_per_g)) +
  geom_boxplot(aes(fill = factor(cluster))) +
  scale_fill_manual(values = c("1" = "#74C0FC", "2" = "#FFD43B", "3" = "#FF8FAB")) +
  theme_bw() +theme(legend.position="none")+ theme(aspect.ratio=1)+
  labs(x = "Cluster",
       y = "Specific Surface Area",
       fill = "Cluster")

# Carbon content
p2 <- ggplot(site_cluster_data, aes(x = factor(cluster), y = C_percent_per_mg)) +
  geom_boxplot(aes(fill = factor(cluster))) +
  scale_fill_manual(values = c("1" = "#74C0FC", "2" = "#FFD43B", "3" = "#FF8FAB")) +
  theme_bw() +theme(legend.position="none")+theme(aspect.ratio=1)+
  labs(x = "Cluster",
       y = "Carbon (%)",
       fill = "Cluster")

# ATP content
p3 <- ggplot(site_cluster_data, aes(x = factor(cluster), y = Mean_ATP_picomoles_per_g)) +
  geom_boxplot(aes(fill = factor(cluster))) +
  scale_fill_manual(values = c("1" = "#74C0FC", "2" = "#FFD43B", "3" = "#FF8FAB")) +
  theme_bw() +theme(legend.position="none")+theme(aspect.ratio=1)+
  labs(x = "Cluster",
       y = "ATP (picomoles/g)",
       fill = "Cluster")

# Sand content
p4 <- ggplot(site_cluster_data, aes(x = factor(cluster), y = Percent_Tot_Sand)) +
  geom_boxplot(aes(fill = factor(cluster))) +
  scale_fill_manual(values = c("1" = "#74C0FC", "2" = "#FFD43B", "3" = "#FF8FAB")) +
  theme_bw() + theme(legend.position="none")+theme(aspect.ratio=1)+
  labs(x = "Cluster",
       y = "Total Sand (%)",
       fill = "Cluster")

p5 = ggplot(site_cluster_data, aes(x = factor(cluster), y = Mean_Fe_mg_per_kg)) +
  geom_boxplot(aes(fill = factor(cluster))) +
  scale_fill_manual(values = c("1" = "#74C0FC", "2" = "#FFD43B", "3" = "#FF8FAB")) +
  theme_bw() + theme(legend.position="none")+theme(aspect.ratio=1)+
  labs(x = "Cluster",
       y = "Fe (mg/kg)",
       fill = "Cluster")

p6 <- ggplot(site_cluster_data, aes(x = factor(cluster), y = N_percent_per_mg)) +
  geom_boxplot(aes(fill = factor(cluster))) +
  scale_fill_manual(values = c("1" = "#74C0FC", "2" = "#FFD43B", "3" = "#FF8FAB")) +
  theme_bw() +theme(legend.position="none")+theme(aspect.ratio=1)+
  labs(x = "Cluster",
       y = "Nitrogen (%)",
       fill = "Cluster")

p7 <- ggplot(site_cluster_data, aes(x = factor(cluster), y = Mean_Gravimetric_Moisture_g_per_g)) +
  geom_boxplot(aes(fill = factor(cluster))) +
  scale_fill_manual(values = c("1" = "#74C0FC", "2" = "#FFD43B", "3" = "#FF8FAB")) +
  theme_bw() + 
  theme(legend.position="none") + 
  theme(aspect.ratio=1) +
  labs( x = "Cluster",
       y = "Moisture (g/g)",
       fill = "Cluster")

# Doing some stats before exporting 
# if(!require(PMCMRplus)) install.packages("PMCMRplus")
# if(!require(multcompView)) install.packages("multcompView")
# if(!require(gridExtra)) install.packages("gridExtra")
library(PMCMRplus)
library(multcompView)
library(ggplot2)


# Function to add Kruskal-Wallis statistics and post-hoc letters to a boxplot
add_kw_stats <- function(plot, data, y_var, group_var = "cluster") {
  # Extract y variable from the plot
  y_data <- data[[y_var]]
  
  # Get y-axis label from the plot
  y_lab <- ggplot_build(plot)$layout$panel_params[[1]]$y$name
  
  # Perform Kruskal-Wallis test
  kw_formula <- as.formula(paste(y_var, "~", group_var))
  kw_result <- kruskal.test(kw_formula, data = data)
  
  # Format p-value
  p_value_formatted <- ifelse(kw_result$p.value < 0.001, "p < 0.001***", 
                              ifelse(kw_result$p.value < 0.01, paste0("p = ", round(kw_result$p.value, 3), "**"),
                                     ifelse(kw_result$p.value < 0.05, paste0("p = ", round(kw_result$p.value, 3), "*"),
                                            paste0("p = ", round(kw_result$p.value, 3), " "))))
  
  # Perform post-hoc test if KW is significant
  if(kw_result$p.value < 0.05) {
    posthoc_result <- kwAllPairsDunnTest(kw_formula, data = data, p.adjust.method = "bonf")
    
    # Extract p-values from the post-hoc test
    pvalues <- posthoc_result$p.value
    diag(pvalues) <- 1
    
    # Generate letters for the groups
    letters <- multcompLetters(pvalues, compare = "<", threshold = 0.05, Letters = letters)
    
    # Create a data frame for the letters
    letter_df <- data.frame(
      cluster = names(letters$Letters), 
      letter = unname(letters$Letters),
      stringsAsFactors = FALSE
    )
    
    # Get max y-value for plotting
    y_max <- max(y_data, na.rm = TRUE)
    
    # Add stats and letters to plot
    enhanced_plot <- plot +
      # Add Kruskal-Wallis result
      annotate("text", x = 2, y = y_max * 1.15, 
               label = paste("H =", round(kw_result$statistic, 1), ",", p_value_formatted),
               size = 3) +
      # Add letters for significant differences
      geom_text(data = letter_df, 
                aes(x = cluster, y = y_max * 0.9, label = letter),
                size = 4)
  } else {
    # If not significant, just add KW stats
    y_max <- max(y_data, na.rm = TRUE)
    enhanced_plot <- plot +
      annotate("text", x = 2, y = y_max * 1.15, 
               label = paste("K-W H =", round(kw_result$statistic, 1), ",", p_value_formatted),
               size = 3)
  }
  
  # Make sure the plot extends high enough to show annotations
  enhanced_plot <- enhanced_plot + 
    coord_cartesian(ylim = c(NA, y_max * 1.2))
  
  return(enhanced_plot)
}

# Apply the function to each plot with appropriate variables
p1_enhanced <- add_kw_stats(p1, site_cluster_data, "Mean_Specific_Surface_Area_m2_per_g")
p2_enhanced <- add_kw_stats(p2, site_cluster_data, "C_percent_per_mg")
p3_enhanced <- add_kw_stats(p3, site_cluster_data, "Mean_ATP_picomoles_per_g")
p4_enhanced <- add_kw_stats(p4, site_cluster_data, "Percent_Tot_Sand")
p5_enhanced <- add_kw_stats(p5, site_cluster_data, "Mean_Fe_mg_per_kg")
p6_enhanced <- add_kw_stats(p6, site_cluster_data, "N_percent_per_mg")
p7_enhanced <- add_kw_stats(p7, site_cluster_data, "Mean_Gravimetric_Moisture_g_per_g")

# Arrange all plots in a grid - adjust layout as needed
grid_arranged_plots <- grid.arrange(
  p1_enhanced, p2_enhanced, p3_enhanced, 
  p4_enhanced, p5_enhanced, p6_enhanced,p7_enhanced,
  ncol = 2
)


# Save as PDF and PNG
ggsave(
  filename = "Figures/FigureS3_Field_variables_within_clusters.pdf",
  plot = grid_arranged_plots,
  width = 10,
  height = 14,
  units = "in",
  dpi = 300
)

ggsave(
  filename = "Figures/FigureS3_Field_variables_within_clusters.png",
  plot = grid_arranged_plots,
  width = 10,
  height = 14,
  units = "in",
  dpi = 300
)


# ==== Heat map plot Field data =====
#Create main heatmap with all annotations and legends
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

# Also create a PNG version
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
library(dunn.test)  # For Dunn's test
library(ggpubr)     # For statistical annotations

# Perform tests for ΔGcox
kruskal_gibbs <- kruskal.test(effect_delGcoxPerCmol ~ cluster, data = dom_cluster_data)
dunn_gibbs <- dunn.test(dom_cluster_data$effect_delGcoxPerCmol, 
                        dom_cluster_data$cluster, 
                        method = "bonferroni")

# Perform tests for Lambda  
kruskal_lambda <- kruskal.test(effect_Lambda ~ cluster, data = dom_cluster_data)
dunn_lambda <- dunn.test(dom_cluster_data$effect_Lambda, 
                         dom_cluster_data$cluster, 
                         method = "bonferroni")


effectGibbs <- ggplot(dom_cluster_data, aes(y = effect_delGcoxPerCmol, x = cluster, fill = factor(cluster))) +
  geom_boxplot() +
  scale_fill_manual(values = c("1" = "#74C0FC", "2" = "#FFD43B", "3" = "#FF8FAB")) +
  annotate("text", x = 1.2, y = max(dom_cluster_data$effect_delGcoxPerCmol, na.rm = TRUE) + 0.005,
           label = paste("K-W p < 0.001"),
           size = 5) +
  annotate("text", x = 1, y = max(dom_cluster_data$effect_delGcoxPerCmol, na.rm = TRUE),
           label = 'ac', size = 5) +
  annotate("text", x = 2, y = max(dom_cluster_data$effect_delGcoxPerCmol, na.rm = TRUE)+0.005 ,
           label = 'ab', size = 5) +
  annotate("text", x = 3, y = max(dom_cluster_data$effect_delGcoxPerCmol, na.rm = TRUE)-0.005 ,
           label = 'bc', size = 5) +
  
  theme_bw() + theme(aspect.ratio = 1, legend.position = "bottom",
                     axis.text = element_text(size = 14, color = "black"),      # Axis numbers
                     axis.title = element_text(size = 16, color = "black"),     # Axis titles  
                     legend.text = element_text(size = 14, color = "black"),    # Legend text
                     legend.title = element_text(size = 16, color = "black"),   # Legend title
                     plot.title = element_text(size = 18, color = "black"),     # Plot title
                     legend.key.size = unit(1.2, "cm")                         # Larger legend symbols
  ) +
  labs(x = '', y = expression(paste(Effect~Size~Delta, G[cox], ~(kJ~Cmol^{-1}))),
       fill = 'Cluster', title = 'B')



effectlambda <- ggplot(dom_cluster_data, aes(y = effect_Lambda, x = cluster, fill = factor(cluster))) +
  geom_boxplot() +
  scale_fill_manual(values = c("1" = "#74C0FC", "2" = "#FFD43B", "3" = "#FF8FAB")) +
  # Add Kruskal-Wallis result
  annotate("text", x = 1.2, y = max(dom_cluster_data$effect_Lambda, na.rm = TRUE) + 0.05,
           label = paste("K-W p < 0.001"),
           size = 5) +
  # Add letters above each boxplot
  annotate("text", x = 1, y = max(dom_cluster_data$effect_Lambda, na.rm = TRUE)- 0.02,
           label = 'a', size = 5) +
  annotate("text", x = 2, y = max(dom_cluster_data$effect_Lambda, na.rm = TRUE) -0.02,
           label = 'b', size = 5) +
  annotate("text", x = 3, y = max(dom_cluster_data$effect_Lambda, na.rm = TRUE) + 0.02,
           label = 'ab', size = 5) +
  theme_bw() + 
  theme(aspect.ratio = 1, legend.position = "bottom",
        axis.text = element_text(size = 14, color = "black"),      # Axis numbers
        axis.title = element_text(size = 16, color = "black"),     # Axis titles  
        legend.text = element_text(size = 14, color = "black"),    # Legend text
        legend.title = element_text(size = 16, color = "black"),   # Legend title
        plot.title = element_text(size = 18, color = "black"),     # Plot title
        legend.key.size = unit(1.2, "cm")                         # Larger legend symbols
  )  +
  labs(x = '', 
       y = expression(paste(Effect~Size~lambda)),
       fill = 'Cluster', 
       title = 'C')

# Print test results to console
print("ΔGcox Kruskal-Wallis Test:")
print(kruskal_gibbs)
print("ΔGcox Dunn's Test:")
print(dunn_gibbs)

print("Lambda Kruskal-Wallis Test:")
print(kruskal_lambda)  
print("Lambda Dunn's Test:")
print(dunn_lambda)

library(patchwork)  # For combining plots
library(ggplot2)

# Combine the plots using patchwork
combined_plot <- effectGibbs + effectlambda + 
  plot_layout(guides = "collect") &  # Collect legends together
  theme(legend.position = "bottom")

# Save the combined plot
ggsave(
  "Figures/Figure2_effect_sizes_by_cluster.png",
  combined_plot,
  width = 10,
  height = 5,
  dpi = 300
)

# For a PDF version
ggsave(
  "Figures/Figure2_effect_sizes_by_cluster.pdf",
  combined_plot,
  width = 10,
  height = 5,
  device = cairo_pdf
)

# ==== Perform a linear Discriminant analysis for Field Data ====
library(MASS)
# library(dplyr)
library(ggplot2)
library(gridExtra)

# ===== Prepare data for LDA =====
# Create site-level dataset combining clusters with environmental variables
site_lda_data <- data.frame(
  site = names(site_clusters),
  cluster = factor(site_clusters),  # Make sure it's a factor
  stringsAsFactors = FALSE
) %>%
  left_join(site_metadata, by = "site") %>%
  # Remove any rows with missing values
  na.omit()

# Check your data
cat("Cluster distribution:\n")
print(table(site_lda_data$cluster))
cat("\nNumber of sites per cluster:\n")
print(site_lda_data %>% count(cluster))

# ===== Perform LDA =====
# LDA to distinguish clusters based on environmental variables
lda_model <- lda(cluster ~ 
                   Mean_ATP_picomoles_per_g +
                   Mean_Specific_Surface_Area_m2_per_g +
                   C_percent_per_mg +
                   N_percent_per_mg +
                   Percent_Tot_Sand +
                   Mean_Fe_mg_per_kg +
                   Mean_Gravimetric_Moisture_g_per_g,
                 data = site_lda_data)

# Print LDA results
cat("\n===== LDA RESULTS =====\n")
print(lda_model)

# ===== Assess classification accuracy =====
# Predict cluster membership
lda_predictions <- predict(lda_model, site_lda_data)

# Create confusion matrix
confusion_matrix <- table(Predicted = lda_predictions$class, 
                          Actual = site_lda_data$cluster)
cat("\n===== CONFUSION MATRIX =====\n")
print(confusion_matrix)

# Calculate accuracy
accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
cat("\nOverall Classification Accuracy:", round(accuracy * 100, 1), "%\n")

# ===== Variable importance in discrimination =====
# Extract coefficients for each discriminant function
lda_coefficients <- as.data.frame(lda_model$scaling)
lda_coefficients$Variable <- rownames(lda_coefficients)

cat("\n===== DISCRIMINANT FUNCTION COEFFICIENTS =====\n")
print(lda_coefficients)

# ===== Visualizations =====
# 1. LDA plot (discriminant space)
lda_plot_data <- data.frame(
  site = site_lda_data$site,
  cluster = site_lda_data$cluster,
  LD1 = lda_predictions$x[,1],
  LD2 = lda_predictions$x[,2]
)

lda_stats <- list(
  # Proportion of variance explained by each LD
  prop_var = lda_model$svd^2 / sum(lda_model$svd^2),
  # Classification accuracy (if you have predictions)
  accuracy = sum(lda_predictions$class == site_lda_data$cluster) / nrow(site_lda_data)
)

# Enhanced plot with ellipses and internal stats
p1 <- ggplot(lda_plot_data, aes(x = LD1, y = LD2, color = cluster, fill = cluster)) +
  # Add confidence ellipses (95% confidence)
  stat_ellipse(level = 0.95, alpha = 0.2, size = 1) +
  # Original points
  geom_point(size = 3, alpha = 0.7) +
  # Colors
  scale_color_manual(values = c("1" = "#74C0FC", "2" = "#FFD43B", "3" = "#FF8FAB")) +
  scale_fill_manual(values = c("1" = "#74C0FC", "2" = "#FFD43B", "3" = "#FF8FAB")) +
  # Add statistics text inside plot
  annotate("text", x = Inf, y = Inf, 
           label = sprintf("Accuracy: %.1f%%",lda_stats$accuracy * 100),
           hjust = 1.1, vjust = 1.1, size = 3.5, 
           color = "black", fontface = "plain") +
  # Labels
  labs( x = sprintf("LD1 (%.1f%% of variance)", lda_stats$prop_var[1] * 100),
          y = sprintf("LD2 (%.1f%% of variance)", lda_stats$prop_var[2] * 100),
       color = "Response\nCluster",
       fill = "Response\nCluster",
       title = 'D') +
  theme_bw() +
  theme(legend.position = "right",
        axis.text = element_text(size = 14),        # Axis numbers
        axis.title = element_text(size = 16),       # Axis titles  
        legend.text = element_text(size = 14),      # Legend text
        legend.title = element_text(size = 16),     # Legend title
        plot.title = element_text(size = 18),       # Plot title
        legend.key.size = unit(1.2, "cm")          # Larger legend symbols
  )

ggsave("Figures/Figure_2_LDA_discriminant_space.png", plot = p1, width = 8, height = 6, dpi = 300)

ggsave("Figures/Figure_2_LDA_discriminant_space.pdf", plot = p1, width = 8, height = 6, dpi = 300)
# 2. Variable importance plot
# Calculate the magnitude of coefficients for each variable
var_importance <- lda_coefficients %>%
  mutate(
    LD1_abs = abs(LD1),
    LD2_abs = abs(LD2),
    Total_importance = sqrt(LD1^2 + LD2^2)
  ) %>%
  arrange(desc(Total_importance))

p2 <- ggplot(var_importance, aes(x = reorder(Variable, Total_importance), y = Total_importance)) +
  geom_col(fill = "steelblue", alpha = 0.7) +
  coord_flip() +
  labs(title = "Variable Importance in Cluster Discrimination",
       x = "Environmental Variables",
       y = "Discriminant Importance") +
  theme_bw()

ggsave("Figures/LDA_variable_importance.png", plot = p2, width = 8, height = 6, dpi = 300)

ggsave("Figures/LDA_variable_importance.pdf", plot = p2, width = 8, height = 6, dpi = 300)
# 3. Boxplots of top discriminating variables by cluster
top_vars <- head(var_importance$Variable, 3)

# Create boxplots for top 3 variables
plot_list <- list()
for(i in 1:3) {
  var_name <- top_vars[i]
  
  p <- ggplot(site_lda_data, aes_string(x = "cluster", y = var_name, fill = "cluster")) +
    geom_boxplot(alpha = 0.7) +
    scale_fill_manual(values = c("1" = "#74C0FC", "2" = "#FFD43B", "3" = "#FF8FAB")) +
    labs(title = paste("Cluster Differences:", var_name),
         x = "Response Cluster",
         y = var_name) +
    theme_bw() +
    theme(legend.position = "none")
  
  plot_list[[i]] <- p
}

combined_boxplots <- grid.arrange(grobs = plot_list, ncol = 3)

# Save combined boxplots
ggsave("Figures/LDA_top_variables_boxplots.png", plot = combined_boxplots, 
       width = 12, height = 4, dpi = 300)

ggsave("Figures/LDA_top_variables_boxplots.pdf", plot = combined_boxplots, 
       width = 12, height = 4, dpi = 300)

# ===== Save results =====
# Save LDA results and predictions
lda_results <- list(
  model = lda_model,
  predictions = lda_predictions,
  accuracy = accuracy,
  confusion_matrix = confusion_matrix,
  variable_importance = var_importance)

saveRDS(lda_results, "Figures/lda_cluster_analysis.rds")
