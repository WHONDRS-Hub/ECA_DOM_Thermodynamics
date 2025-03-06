# Load required libraries
rm(list=ls(all=T))

library(tidyverse)
library(pheatmap)
library(viridis)
library(dendextend) # For enhanced dendrograms
# Load required libraries
library(dplyr)
library(viridis)
library(RColorBrewer)
library(pheatmap)
library(tibble)
library(grid)
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
  Percent_Tot_Sand, Mean_Fe_mg_per_kg)

sample_data$site = gsub('CM','EC',sample_data$site)
sample_data$site = gsub('_Sediment','',sample_data$site)


field_metadata = read.csv(paste0(github,'EC_Data_Package/EC_Field_Metadata.csv')) %>%
  dplyr::select(site = Parent_ID,State)


# ==== Set up data =====
dom_data = data %>%
  dplyr::filter(site != "EC_023") %>% 
  dplyr::filter(!(site %in% c("EC_012","EC_011","EC_053","EC_057","EC_052"))) %>%
  #Remove EC_011, EC_012, EC_023, EC_052, EC_053, and EC_057 for having too much water added (11 and 12), no mg/kg calculation (23), and being duplicated NEON sites (52, 53, and 57)
  dplyr::select(site,Treatment,Median_delGcoxPerCmol,Median_Lambda)
 
sites = as.data.frame(unique(dom_data$site))
names(sites) = 'site'
explanatory_data = merge(sample_data,field_metadata, by = 'site')
explanatory_data = merge(explanatory_data,sites, by = 'site')

site_metadata = explanatory_data
# ===== Calculations ======
# Load required libraries
library(tidyverse)
library(pheatmap)
library(viridis)
library(dendextend)

# Step 1: Calculate treatment effects as before
treatment_effects <- dom_data %>%
  filter(Treatment %in% c("Wet", "Dry")) %>%
  select(site, Treatment, starts_with("Median_")) %>%
  pivot_wider(
    id_cols = site,
    names_from = Treatment,
    values_from = starts_with("Median_")
  ) %>%
  mutate(
    effect_delGcoxPerCmol = log(Median_delGcoxPerCmol_Wet / Median_delGcoxPerCmol_Dry),
    effect_Lambda = log(Median_Lambda_Wet / Median_Lambda_Dry)
  ) %>%
  select(site, starts_with("effect_"))

# Step 2: Prepare the effect matrix for heatmap
effect_matrix <- treatment_effects %>%
  column_to_rownames("site") %>%
  as.matrix()

effect_matrix_scaled <- scale(effect_matrix)

# Step 3: Create ecoregion mapping based on full state names
ecoregion_map <- list(
  "Northeast" = c("Maine", "New Hampshire", "Vermont", "Massachusetts", "Rhode Island", 
                  "Connecticut", "New York", "New Jersey", "Pennsylvania", "Delaware", "Maryland"),
  "Southeast" = c("Virginia", "North Carolina", "South Carolina", "Georgia", "Florida", 
                  "Alabama", "Mississippi", "Tennessee", "Kentucky", "West Virginia"),
  "Midwest" = c("Ohio", "Indiana", "Illinois", "Michigan", "Wisconsin", 
                "Minnesota", "Iowa", "Missouri"),
  "Great Plains" = c("North Dakota", "South Dakota", "Nebraska", "Kansas", "Oklahoma", "Texas"),
  "Southwest" = c("Arizona", "New Mexico", "Nevada", "Utah"),
  "Rocky Mountains" = c("Montana", "Idaho", "Wyoming", "Colorado"),
  "Pacific Northwest" = c("Washington", "Oregon"),
  "Pacific" = c("California", "Hawaii"),
  "Alaska" = c("Alaska")
)

# Function to map state to ecoregion
state_to_ecoregion <- function(state) {
  for (region in names(ecoregion_map)) {
    if (state %in% ecoregion_map[[region]]) {
      return(region)
    }
  }
  return("Other")  # Default for any unmapped states
}

# Load required libraries
library(pheatmap)
library(viridis)
library(colorspace)
library(dplyr)

# Add Region to metadata and keep numerical variables as-is
site_annotation <- site_metadata %>% 
  mutate(Ecoregion = sapply(State, state_to_ecoregion)) %>%
  select(
    site, 
    Region = Ecoregion,
    State,
    Mean_ATP_picomoles_per_g,
    Mean_Specific_Surface_Area_m2_per_g,
    C_percent_per_mg,
    N_percent_per_mg,
    Percent_Tot_Sand,
    Mean_Fe_mg_per_kg
  ) %>%
  column_to_rownames("site")

# Ensure annotations match the order of sites in the heatmap
site_annotation <- site_annotation[rownames(effect_matrix), ]

# Define categorical palettes with improved colors
categorical_colors <- list(
  # Region - use viridis with alpha transparency to lighten
  Region = setNames(
    lighten(viridis(length(unique(site_annotation$Region))), 0.2),
    unique(site_annotation$Region)
  ),
  # State - use a brighter palette
  State = setNames(
    lighten(qualitative_hcl(length(unique(site_annotation$State)), palette = "Pastel 1"), 0.1),
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
  Mean_Fe_mg_per_kg = colorRampPalette(c("#FFFFFF", "#FFDEEB", "#FF8FAB", "#E03131"))(100)
)

# Perform hierarchical clustering to group similar sites
site_dist <- dist(effect_matrix_scaled, method = "euclidean")
site_hclust <- hclust(site_dist, method = "ward.D2")

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

# Print all column names to debug
cat("Column names in site_annotation:\n")
print(colnames(site_annotation))
cat("\nKeys in annotation_colors:\n")
print(names(annotation_colors))

# MUCH SIMPLER APPROACH: Create one high-quality plot with all legends
# 1. Create main heatmap with all annotations and legends
pdf("complete_heatmap_with_legends_2.pdf", width = 14, height = 10)
pheatmap(
  effect_matrix_scaled,
  cluster_rows = site_hclust,
  cluster_cols = TRUE,
  clustering_method = "ward.D2",
  color = colorRampPalette(c("#3D5A80", "#FFFFFF", "#E07A5F"))(100),
  breaks = seq(-2, 2, length.out = 101),
  display_numbers = matrix(
    ifelse(abs(effect_matrix_scaled) > 1.5, "*", ""), 
    nrow = nrow(effect_matrix_scaled)
  ),
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
png("complete_heatmap_with_legends.png", width = 4200, height = 3000, res = 300)
pheatmap(
  effect_matrix_scaled,
  cluster_rows = site_hclust,
  cluster_cols = TRUE,
  clustering_method = "ward.D2",
  color = colorRampPalette(c("#3D5A80", "#FFFFFF", "#E07A5F"))(100),
  breaks = seq(-2, 2, length.out = 101),
  display_numbers = matrix(
    ifelse(abs(effect_matrix_scaled) > 1.5, "*", ""), 
    nrow = nrow(effect_matrix_scaled)
  ),
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

# 2. Alternative approach: Create main heatmap only
pdf("main_heatmap_only.pdf", width = 11, height = 10)
pheatmap(
  effect_matrix_scaled,
  cluster_rows = site_hclust,
  cluster_cols = TRUE,
  clustering_method = "ward.D2",
  color = colorRampPalette(c("#3D5A80", "#FFFFFF", "#E07A5F"))(100),
  breaks = seq(-2, 2, length.out = 101),
  display_numbers = matrix(
    ifelse(abs(effect_matrix_scaled) > 1.5, "*", ""), 
    nrow = nrow(effect_matrix_scaled)
  ),
  fontsize_number = 8,
  angle_col = 45,
  annotation_row = site_annotation,
  annotation_colors = annotation_colors,
  cutree_rows = 3,
  annotation_legend = FALSE,
  legend = FALSE,
  fontsize_row = 8,
  cellwidth = 15,
  cellheight = 12,
  border_color = NA
)
dev.off()

# Create a simple legend file for manually combining later
pdf("main_effect_legend.pdf", width = 6, height = 2)
# Empty plot
plot(c(0,1), c(0,1), type = "n", axes = FALSE, xlab = "", ylab = "")

# Add the legend manually
legend("center", 
       legend = c("Strong Negative", "Negative", "Neutral", "Positive", "Strong Positive"),
       fill = colorRampPalette(c("#3D5A80", "#98C1D9", "#FFFFFF", "#F9C784", "#E07A5F"))(5),
       horiz = TRUE,
       title = "Effect Size",
       cex = 1.2,
       bty = "n")
dev.off()




# Step 7: Analyze site clusters with the original metadata
site_cluster_data <- data.frame(
  site = names(site_clusters),
  cluster = site_clusters
) %>%
  left_join(
    site_metadata %>% mutate(Ecoregion = sapply(State, state_to_ecoregion)),
    by = "site"
  )

# Visualize cluster composition by Region
ggplot(site_cluster_data, aes(x = factor(cluster), fill = Ecoregion)) +
  geom_bar(position = "fill") +
  theme_bw() +
  labs(title = "Composition of Response Clusters by Ecoregion",
       x = "Cluster",
       y = "Proportion")

# Visualize numerical variables by cluster - creating more informative boxplots
# Moisture content
p1 <- ggplot(site_cluster_data, aes(x = factor(cluster), y = Mean_Specific_Surface_Area_m2_per_g)) +
  geom_boxplot(aes(fill = factor(cluster))) +
  scale_fill_manual(values = c("1" = "#74C0FC", "2" = "#FFD43B", "3" = "#FF8FAB")) +
  theme_bw() +
  labs(x = "Cluster",
       y = "Specific Surface Area",
       fill = "Cluster")

# Carbon content
p2 <- ggplot(site_cluster_data, aes(x = factor(cluster), y = C_percent_per_mg)) +
  geom_boxplot(aes(fill = factor(cluster))) +
  scale_fill_manual(values = c("1" = "#74C0FC", "2" = "#FFD43B", "3" = "#FF8FAB")) +
  theme_bw() +
  labs(title = "Carbon Content by Cluster",
       x = "Cluster",
       y = "Carbon (%)",
       fill = "Cluster")

# ATP content
p3 <- ggplot(site_cluster_data, aes(x = factor(cluster), y = Mean_ATP_picomoles_per_g)) +
  geom_boxplot(aes(fill = factor(cluster))) +
  scale_fill_manual(values = c("1" = "#74C0FC", "2" = "#FFD43B", "3" = "#FF8FAB")) +
  theme_bw() +
  labs(title = "ATP Content by Cluster",
       x = "Cluster",
       y = "ATP (picomoles/g)",
       fill = "Cluster")

# Sand content
p4 <- ggplot(site_cluster_data, aes(x = factor(cluster), y = Percent_Tot_Sand)) +
  geom_boxplot(aes(fill = factor(cluster))) +
  scale_fill_manual(values = c("1" = "#74C0FC", "2" = "#FFD43B", "3" = "#FF8FAB")) +
  theme_bw() +
  labs(title = "Sand Content by Cluster",
       x = "Cluster",
       y = "Total Sand (%)",
       fill = "Cluster")

# Combine plots with gridExtra if available
# install.packages("gridExtra")
 library(gridExtra)
grid.arrange(p1, p2, p3, p4, ncol = 2)

# Step 8: Analyze average treatment effect by cluster
effect_by_cluster <- treatment_effects %>%
  mutate(cluster = site_clusters[site]) %>%
  group_by(cluster) %>%
  summarize(across(starts_with("effect_"), mean)) %>%
  pivot_longer(cols = -cluster, names_to = "variable", values_to = "avg_effect") %>%
  # Clean up variable names for plotting
  mutate(variable = str_replace(variable, "effect_", ""))

# Plot average effects by cluster
ggplot(effect_by_cluster, aes(x = variable, y = avg_effect, fill = factor(cluster))) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("1" = "#74C0FC", "2" = "#FFD43B", "3" = "#FF8FAB")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Average Treatment Effects by Response Cluster",
       x = "DOM Property",
       y = "Average Log Response Ratio",
       fill = "Cluster")

# ===== Identify if clusters differ in their site characteristics =====
# Compare site characteristics across clusters statistically
cluster_stats <- site_cluster_data %>%
  group_by(cluster) %>%
  summarize(
    mean_sand = mean(Percent_Tot_Sand, na.rm = TRUE),
    mean_carbon = mean(`C_percent_per_mg`, na.rm = TRUE),
    mean_nitrogen = mean(`N_percent_per_mg`, na.rm = TRUE)
  )

# ANOVA to test for significant differences
sand_anova <- aov(Percent_Tot_Sand ~ factor(cluster), data = site_cluster_data)
summary(sand_anova)

library(vegan)
site_characteristics <- site_cluster_data %>%
  select(Percent_Tot_Sand, `C_percent_per_mg`, `N_percent_per_mg`) %>%
  scale()

adonis2(site_characteristics ~ factor(site_cluster_data$cluster), method = 'euclidean')

# ===== Identify if clusters differ in their site characteristics =====
library(factoextra)  # For PCA visualization
library(FactoMineR)  # For PCA
missing_values <- colSums(is.na(site_cluster_data))
print("Missing values per column:")
print(missing_values)

#  remove rows with NA values in numeric columns
site_data_clean <- site_cluster_data %>%
  filter(!is.na(Mean_ATP_picomoles_per_g) & 
           !is.na(Mean_Specific_Surface_Area_m2_per_g) &
           !is.na(C_percent_per_mg) &
           !is.na(N_percent_per_mg) &
           !is.na(Percent_Tot_Sand) &
           !is.na(Mean_Fe_mg_per_kg))

# Convert cluster to factor if not already
site_data_clean$cluster <- as.factor(site_data_clean$cluster)

# Select only numerical environmental variables for multivariate analysis
env_vars <- site_data_clean %>%
  select(Mean_ATP_picomoles_per_g, Mean_Specific_Surface_Area_m2_per_g,
         C_percent_per_mg, N_percent_per_mg, Percent_Tot_Sand, Mean_Fe_mg_per_kg)

# Get categorical variables
cat_vars <- site_data_clean %>%
  select(State, Ecoregion)

# 1. Principal Component Analysis (PCA)
# Standardize data (important for variables with different scales)
pca_result <- PCA(env_vars, scale.unit = TRUE, graph = FALSE)

# Visualize PCA with clusters
fviz_pca_biplot(pca_result, 
                habillage = site_data_clean$cluster,  # Color by cluster
                geom = "point",
                pointsize = 2,
                palette = c("#74C0FC", "#FFD43B", "#FF8FAB"),                addEllipses = TRUE,
                ellipse.type = "confidence",
                ggtheme = theme_bw(),
                title = " ") +
  theme(legend.title = element_text(face = "bold"))

# 2. PERMANOVA - test if clusters differ significantly in environmental space
perm_result <- adonis2(env_vars ~ cluster, data = site_data_clean, 
                       method = "euclidean", permutations = 999)
print("PERMANOVA Results:")
print(perm_result)
