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
  Percent_Tot_Sand, Mean_Fe_mg_per_kg)

sample_data$site = gsub('CM','EC',sample_data$site)
sample_data$site = gsub('_Sediment','',sample_data$site)


field_metadata = read.csv(paste0(github,'EC_Data_Package/EC_Field_Metadata.csv')) %>%
  dplyr::select(site = Parent_ID,State)

ecoregions = read.csv('Data/EC_Ecoregions.csv')%>%
  dplyr::select(site = 'Parent_ID', Ecoregion = Name)%>%
  mutate(Ecoregion = str_trim(str_extract(Ecoregion, "(?<=- ).*")))

field_metadata = merge(field_metadata,ecoregions, by = 'site')

bulk_medians = read.csv('Data/Bulk_medians_and_effect_size.csv')%>%
  mutate(site = gsub('_all','',Sample_Name)) %>%
  dplyr::select(-Sample_Name)

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
bulk_medians = merge(bulk_medians,sites, by = 'site')
field_sample_data = merge(sample_data, sites, by = 'site')
site_metadata = explanatory_data
bulk_metadata = merge(field_metadata,bulk_medians, by = 'site')
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

# Prepare the effect matrix for heatmap
effect_matrix <- treatment_effects %>%
  column_to_rownames("site") %>%
  as.matrix()

effect_matrix_scaled <- scale(effect_matrix) # z-score normalization for clustering. Clustering uses Eucledean distance which is sensitive to scale

# ===== Cube root transforming data before correlations =====
cube_root <- function(x) sign(x) * (abs(x))^(1/3)

# Fe outlier not included in analysis - remove from DF
cube_bulk_medians = bulk_medians %>% 
  mutate(across(where(is.numeric), cube_root)) %>% # cube root transform data
  rename_with(where(is.numeric), .fn = ~ paste0("cube_", .x)) %>% 
  #column_to_rownames("Sample_Name") %>%
  filter(cube_Effect_Size_Fe_mg_per_kg > -1) %>%  # remove Fe outlier for analysis site EC_032
  select(-contains("per_L")) %>% # remove per_L data, analysis ran on _per_kg
  relocate(cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H, .before = cube_Effect_Size_SpC_microsiemens_per_cm)

cube_field_means = field_sample_data %>% 
  filter(Mean_Fe_mg_per_kg < 10) %>%  # remove Fe outlier for analysis site EC_068
  mutate(across(where(is.numeric), cube_root)) %>% # cube root transform data
  rename_with(where(is.numeric), .fn = ~ paste0("cube_", .x)) %>% 
  #column_to_rownames("Sample_Name") %>%
  select(-contains("per_L")) # remove per_L data, analysis ran on _per_kg

cube_treatment_effects = treatment_effects %>%  
  mutate(across(where(is.numeric), cube_root)) %>% # cube root transform data
  rename_with(where(is.numeric), .fn = ~ paste0("cube_", .x))

dat_bulk = merge(cube_treatment_effects,cube_bulk_medians, by = 'site') %>%
  column_to_rownames('site')

dat_field = merge(cube_treatment_effects,cube_field_means, by = 'site') %>%
  column_to_rownames('site')

# ====== Correlation matrices for lab and field data =====
bulk_correlation <- cor(dat_bulk, method = "pearson")
bulk_significance <- rcorr(as.matrix(dat_bulk), type = "pearson")
# Plot
pdf("Figures/FigureS2_Pearson_Correlation_Treatment_Effects_vs_bulk_medians.pdf", width = 14, height = 14)  

par(mfrow=c(1,1))
corrplot(bulk_correlation, method = "circle", type = "upper", 
         tl.col = "black", tl.srt = 45, 
         tl.cex = 0.8)

dev.off()

png("Figures/FigureS2_Pearson_Correlation_Treatment_Effects_vs_bulk_medians.png", width = 14, height = 14, units = "in", res = 300)  

par(mfrow=c(1,1))
corrplot(bulk_correlation, method = "circle", type = "upper", 
         tl.col = "black", tl.srt = 45, 
         tl.cex = 0.8)

dev.off()

field_correlation <- cor(dat_field, method = "pearson")
field_significance <- rcorr(as.matrix(dat_field), type = "pearson")

# Plot
pdf("Figures/FigureS3_Pearson_Correlation_Treatment_Effects_vs_field_means.pdf", width = 14, height = 14)  

par(mfrow=c(1,1))
corrplot(field_correlation, method = "circle", type = "upper", 
         tl.col = "black", tl.srt = 45, 
         tl.cex = 1)

dev.off()

png("Figures/FigureS3_Pearson_Correlation_Treatment_Effects_vs_field_means.png", width = 14, height = 14, units = "in", res = 300)  

par(mfrow=c(1,1))
corrplot(field_correlation, method = "circle", type = "upper", 
         tl.col = "black", tl.srt = 45, 
         tl.cex = 1)

dev.off()

# Export in a csv all the correlations
# Calculate correlations and p-values between treatment effects and bulk variables

# Create dataframe for all bulk correlations
bulk_all_corr <- data.frame(
  treatment_effect = character(),
  bulk_variable = character(),
  correlation = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Extract treatment effect names and bulk variable names
treatment_effects <- colnames(dat_bulk)[1:2]
bulk_variables <- colnames(dat_bulk)[3:ncol(dat_bulk)]

# Fill dataframe with all correlations and p-values
for (effect in treatment_effects) {
  for (var in bulk_variables) {
    effect_idx <- which(rownames(bulk_significance$r) == effect)
    var_idx <- which(colnames(bulk_significance$r) == var)
    
    if (length(effect_idx) > 0 && length(var_idx) > 0) {
      corr_value <- bulk_significance$r[effect_idx, var_idx]
      p_value <- bulk_significance$P[effect_idx, var_idx]
      
      bulk_all_corr <- rbind(bulk_all_corr, data.frame(
        treatment_effect = effect,
        bulk_variable = var,
        correlation = corr_value,
        p_value = p_value,
        stringsAsFactors = FALSE
      ))
    }
  }
}

# Sort by treatment effect and absolute correlation value
bulk_all_corr <- bulk_all_corr[order(bulk_all_corr$treatment_effect, -abs(bulk_all_corr$correlation)), ]

# Export all bulk correlations to CSV
write.csv(bulk_all_corr, file = "Figures/SI_all_bulk_correlations.csv", row.names = FALSE)

# Field data
# Create dataframe for all field correlations
field_all_corr <- data.frame(
  treatment_effect = character(),
  field_variable = character(),
  correlation = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Extract field variable names
field_variables <- colnames(dat_field)[3:ncol(dat_field)]

# Fill dataframe with all correlations and p-values
for (effect in treatment_effects) {
  for (var in field_variables) {
    effect_idx <- which(rownames(field_significance$r) == effect)
    var_idx <- which(colnames(field_significance$r) == var)
    
    if (length(effect_idx) > 0 && length(var_idx) > 0) {
      corr_value <- field_significance$r[effect_idx, var_idx]
      p_value <- field_significance$P[effect_idx, var_idx]
      
      field_all_corr <- rbind(field_all_corr, data.frame(
        treatment_effect = effect,
        field_variable = var,
        correlation = corr_value,
        p_value = p_value,
        stringsAsFactors = FALSE
      ))
    }
  }
}

# Sort by treatment effect and absolute correlation value
field_all_corr <- field_all_corr[order(field_all_corr$treatment_effect, -abs(field_all_corr$correlation)), ]

# Export all field correlations to CSV
write.csv(field_all_corr, file = "Figures/SI_all_field_correlations.csv", row.names = FALSE)

# === Filtering variables with higher correlations =====
high_corr_vars <- list()
high_corr_vars_field <- list()

# Loop through each treatment effect variable
for (effect_var in colnames(cube_treatment_effects)[2:3]) {
  # For bulk data
  if(effect_var %in% rownames(bulk_correlation)) {
    correlations_bulk <- bulk_correlation[effect_var, ]
    strong_correlations_bulk <- abs(correlations_bulk) > 0.20
    strong_var_names_bulk <- names(correlations_bulk)[strong_correlations_bulk]
    
    if (length(strong_var_names_bulk) > 0) {
      high_corr_vars[[effect_var]] <- data.frame(
        treatment_effect = effect_var,  # Add this column
        env_variable = strong_var_names_bulk,
        correlation = correlations_bulk[strong_correlations_bulk],
        stringsAsFactors = FALSE
      )
      high_corr_vars[[effect_var]] <- high_corr_vars[[effect_var]][
        order(-abs(high_corr_vars[[effect_var]]$correlation)), ]
    } else {
      high_corr_vars[[effect_var]] <- data.frame(
        treatment_effect = effect_var,  # Add this column
        env_variable = "None > 0.2",
        correlation = NA,
        stringsAsFactors = FALSE
      )
    }
  }
  
  # For field data
  if(effect_var %in% rownames(field_correlation)) {
    correlations_field <- field_correlation[effect_var, ]
    strong_correlations_field <- abs(correlations_field) > 0.2
    strong_var_names_field <- names(correlations_field)[strong_correlations_field]
    
    if (length(strong_var_names_field) > 0) {
      high_corr_vars_field[[effect_var]] <- data.frame(
        treatment_effect = effect_var,  # Add this column
        env_variable = strong_var_names_field,
        correlation = correlations_field[strong_correlations_field],
        stringsAsFactors = FALSE
      )
      high_corr_vars_field[[effect_var]] <- high_corr_vars_field[[effect_var]][
        order(-abs(high_corr_vars_field[[effect_var]]$correlation)), ]  # Fixed this line
    } else {
      high_corr_vars_field[[effect_var]] <- data.frame(
        treatment_effect = effect_var,  # Add this column
        env_variable = "None > 0.2",
        correlation = NA,
        stringsAsFactors = FALSE
      )
    }
  }
}

# Combine results
high_corr_combined <- do.call(rbind, high_corr_vars)

# Add p-values to the combined dataframe - Fixed
high_corr_combined$p_value <- NA
for(i in 1:nrow(high_corr_combined)) {
  if(high_corr_combined$env_variable[i] != "None > 0.2") {
    effect_var <- high_corr_combined$treatment_effect[i]
    env_var <- high_corr_combined$env_variable[i]
    
    # Check if both variables exist in the p-value matrix
    if(effect_var %in% rownames(bulk_significance$P) && 
       env_var %in% colnames(bulk_significance$P)) {
      
      effect_idx <- which(rownames(bulk_significance$P) == effect_var)
      env_idx <- which(colnames(bulk_significance$P) == env_var)
      
      # Check that we found valid indices
      if(length(effect_idx) > 0 && length(env_idx) > 0) {
        high_corr_combined$p_value[i] <- bulk_significance$P[effect_idx, env_idx]
      }
    }
  }
}


# Export to CSV
write.csv(high_corr_combined, file = "Figures/SI_higher_correlation_bulk_variables.csv", row.names = FALSE)

high_corr_combined_field <- do.call(rbind, high_corr_vars_field)

high_corr_combined_field$p_value <- NA

for(i in 1:nrow(high_corr_combined_field)) {
  if(high_corr_combined_field$env_variable[i] != "None > 0.2") {
    effect_var <- high_corr_combined_field$treatment_effect[i]
    env_var <- high_corr_combined_field$env_variable[i]
    
    # Check if both variables exist in the p-value matrix
    if(effect_var %in% rownames(field_significance$P) && 
       env_var %in% colnames(field_significance$P)) {
      
      effect_idx <- which(rownames(field_significance$P) == effect_var)
      env_idx <- which(colnames(field_significance$P) == env_var)
      
      # Check that we found valid indices
      if(length(effect_idx) > 0 && length(env_idx) > 0) {
        high_corr_combined_field$p_value[i] <- field_significance$P[effect_idx, env_idx]
      }
    }
  }
}

# Export to CSV
write.csv(high_corr_combined_field, file = "Figures/SI_higher_correlation_field_variables.csv", row.names = FALSE)


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
    Mean_Fe_mg_per_kg
  ) %>%
  column_to_rownames("site")

# Ensure annotations match the order of sites in the heatmap
site_annotation <- site_annotation[rownames(effect_matrix), ]

# Define categorical palettes with improved colors
categorical_colors <- list(
  Ecoregion = setNames(
    lighten(viridis(length(unique(site_annotation$Ecoregion))), 0.2),
    unique(site_annotation$Ecoregion)
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

# For bulk data
bulk_annotation <- bulk_metadata %>% 
  select(
    site, 
    Ecoregion,
    State,
    Median_ATP_picomoles_per_g,
    Mean_Specific_Surface_Area_m2_per_g,
    Median_X01395_C_percent_per_mg,
    Median_X01397_N_percent_per_mg,
    Percent_Tot_Sand,
    Median_Fe_mg_per_kg,
    median_Dry_Final_Gravimetric,
    Median_pH
  ) %>%
  column_to_rownames("site")

# Ensure annotations match the order of sites in the heatmap
bulk_annotation <- bulk_annotation[rownames(effect_matrix), ]

# Define categorical palettes with improved colors
categorical_colors_bulk <- list(
  Ecoregion = setNames(
    lighten(viridis(length(unique(site_annotation$Ecoregion))), 0.2),
    unique(site_annotation$Ecoregion)
  ),
  # State - use a brighter palette
  State = setNames(
    lighten(qualitative_hcl(length(unique(site_annotation$State)), palette = "Pastel 1"), 0.1),
    unique(site_annotation$State)
  )
)

# Create custom gradients with better visibility
numerical_colors_bulk <- list(
  # ATP - light to medium purple
  Median_ATP_picomoles_per_g = colorRampPalette(c("#FFFFFF", "#E6E1F9", "#B197FC", "#845EF7"))(100),
  
  # SSA - light blue to dark blue gradient
  Mean_Specific_Surface_Area_m2_per_g = colorRampPalette(c("#FFFFFF", "#D0EBFF", "#74C0FC", "#1971C2"))(100),
  
  # Carbon content (improved grayscale)
  Median_X01395_C_percent_per_mg = colorRampPalette(c("#FFFFFF", "#E9ECEF", "#ADB5BD", "#495057"))(100),
  
  # Nitrogen content (light green to dark green)
  Median_X01397_N_percent_per_mg = colorRampPalette(c("#FFFFFF", "#D3F9D8", "#69DB7C", "#2B8A3E"))(100),
  
  # Sand content (light yellow to orange)
  Percent_Tot_Sand = colorRampPalette(c("#FFFFFF", "#FFF3BF", "#FFD43B", "#F08C00"))(100),
  
  # Iron content (light pink to red)
  Median_Fe_mg_per_kg = colorRampPalette(c("#FFFFFF", "#FFDEEB", "#FF8FAB", "#E03131"))(100),
  # Gravimetric Moisture - light teal to dark teal
  median_Dry_Final_Gravimetric = colorRampPalette(c("#FFFFFF", "#C3FAE8", "#63E6BE", "#0B7285"))(100),
  
  # pH - light peach to burgundy
  Median_pH = colorRampPalette(c("#FFFFFF", "#FFE8CC", "#FFA94D", "#9C4221"))(100)
)

# For effect size
effect_annotation <- bulk_metadata %>% 
  select(
    site, 
    Ecoregion,
    State,
    Effect_Size_ATP_picomoles_per_g,
    Mean_Specific_Surface_Area_m2_per_g,
    Effect_Size_C_percent_per_mg,
    Effect_Size_N_percent_per_mg,
    Percent_Tot_Sand,
    Effect_Size_Fe_mg_per_kg,
    median_Dry_Final_Gravimetric,
    Effect_Size_pH
  ) %>%
  column_to_rownames("site")

# Ensure annotations match the order of sites in the heatmap
effect_annotation <- effect_annotation[rownames(effect_matrix), ]

# Define categorical palettes with improved colors
categorical_colors_effect <- list(
  Ecoregion = setNames(
    lighten(viridis(length(unique(site_annotation$Ecoregion))), 0.2),
    unique(site_annotation$Ecoregion)
  ),
  # State - use a brighter palette
  State = setNames(
    lighten(qualitative_hcl(length(unique(site_annotation$State)), palette = "Pastel 1"), 0.1),
    unique(site_annotation$State)
  )
)

# Create custom gradients with better visibility
numerical_colors_effect <- list(
  # ATP - light to medium purple
  Effect_Size_ATP_picomoles_per_g = colorRampPalette(c("#FFFFFF", "#E6E1F9", "#B197FC", "#845EF7"))(100),
  
  # SSA - light blue to dark blue gradient
  Mean_Specific_Surface_Area_m2_per_g = colorRampPalette(c("#FFFFFF", "#D0EBFF", "#74C0FC", "#1971C2"))(100),
  
  # Carbon content (improved grayscale)
  Effect_Size_C_percent_per_mg = colorRampPalette(c("#FFFFFF", "#E9ECEF", "#ADB5BD", "#495057"))(100),
  
  # Nitrogen content (light green to dark green)
  Effect_Size_N_percent_per_mg = colorRampPalette(c("#FFFFFF", "#D3F9D8", "#69DB7C", "#2B8A3E"))(100),
  
  # Sand content (light yellow to orange)
  Percent_Tot_Sand = colorRampPalette(c("#FFFFFF", "#FFF3BF", "#FFD43B", "#F08C00"))(100),
  
  # Iron content (light pink to red)
  Effect_Size_Fe_mg_per_kg = colorRampPalette(c("#FFFFFF", "#FFDEEB", "#FF8FAB", "#E03131"))(100),
  # Gravimetric Moisture - light teal to dark teal
  median_Dry_Final_Gravimetric = colorRampPalette(c("#FFFFFF", "#C3FAE8", "#63E6BE", "#0B7285"))(100),
  
  # pH - light peach to burgundy
  Effect_Size_pH = colorRampPalette(c("#FFFFFF", "#FFE8CC", "#FFA94D", "#9C4221"))(100)
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
bulk_colors <- c(categorical_colors_bulk, numerical_colors_bulk)
effect_colors <- c(categorical_colors_effect, numerical_colors_effect)
# ===== Visualize clusters ======
site_cluster_data <- data.frame(
  site = names(site_clusters),
  cluster = site_clusters
) %>%
  left_join(
    site_metadata,
    by = "site"
  )

# Visualize cluster composition by Region
ggplot(site_cluster_data, aes(x = factor(cluster), fill = Ecoregion)) +
  geom_bar(position = "fill") +
  theme_bw() +
  labs(title = "Composition of Response Clusters by Ecoregion",
       x = "Cluster",
       y = "Proportion")

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
  labs(title = " ",
       x = "Cluster",
       y = "Carbon (%)",
       fill = "Cluster")

# ATP content
p3 <- ggplot(site_cluster_data, aes(x = factor(cluster), y = Mean_ATP_picomoles_per_g)) +
  geom_boxplot(aes(fill = factor(cluster))) +
  scale_fill_manual(values = c("1" = "#74C0FC", "2" = "#FFD43B", "3" = "#FF8FAB")) +
  theme_bw() +theme(legend.position="none")+theme(aspect.ratio=1)+
  labs(title = " ",
       x = "Cluster",
       y = "ATP (picomoles/g)",
       fill = "Cluster")

# Sand content
p4 <- ggplot(site_cluster_data, aes(x = factor(cluster), y = Percent_Tot_Sand)) +
  geom_boxplot(aes(fill = factor(cluster))) +
  scale_fill_manual(values = c("1" = "#74C0FC", "2" = "#FFD43B", "3" = "#FF8FAB")) +
  theme_bw() + theme(legend.position="none")+theme(aspect.ratio=1)+
  labs(title = " ",
       x = "Cluster",
       y = "Total Sand (%)",
       fill = "Cluster")

p5 = ggplot(site_cluster_data, aes(x = factor(cluster), y = Mean_Fe_mg_per_kg)) +
  geom_boxplot(aes(fill = factor(cluster))) +
  scale_fill_manual(values = c("1" = "#74C0FC", "2" = "#FFD43B", "3" = "#FF8FAB")) +
  theme_bw() + theme(legend.position="none")+theme(aspect.ratio=1)+
  labs(title = " ",
       x = "Cluster",
       y = "Fe (mg/kg)",
       fill = "Cluster")

p6 <- ggplot(site_cluster_data, aes(x = factor(cluster), y = N_percent_per_mg)) +
  geom_boxplot(aes(fill = factor(cluster))) +
  scale_fill_manual(values = c("1" = "#74C0FC", "2" = "#FFD43B", "3" = "#FF8FAB")) +
  theme_bw() +theme(legend.position="none")+theme(aspect.ratio=1)+
  labs(title = " ",
       x = "Cluster",
       y = "Nitrogen (%)",
       fill = "Cluster")


# Doing some stats before exporting 
# if(!require(PMCMRplus)) install.packages("PMCMRplus")
# if(!require(multcompView)) install.packages("multcompView")
# if(!require(gridExtra)) install.packages("gridExtra")
library(PMCMRplus)
library(multcompView)
library(ggplot2)
library(gridExtra)

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

# Arrange all plots in a grid - adjust layout as needed
grid_arranged_plots <- grid.arrange(
  p1_enhanced, p2_enhanced, p3_enhanced, 
  p4_enhanced, p5_enhanced, p6_enhanced,
  ncol = 2
)


# Save as PDF and PNG
ggsave(
  filename = "Figures/FigureS_Field_variables_within_clusters.pdf",
  plot = grid_arranged_plots,
  width = 10,
  height = 14,
  units = "in",
  dpi = 300
)

ggsave(
  filename = "Figures/FigureS_Field_variables_within_clusters.png",
  plot = grid_arranged_plots,
  width = 10,
  height = 14,
  units = "in",
  dpi = 300
)

# Experiment data
site_cluster_data_bulk <- data.frame(
  site = names(site_clusters),
  cluster = site_clusters
) %>%
  left_join(
    bulk_metadata,
    by = "site"
  )

# SSA
p1 <- ggplot(site_cluster_data_bulk, aes(x = factor(cluster), y = Mean_Specific_Surface_Area_m2_per_g)) +
  geom_boxplot(aes(fill = factor(cluster))) +
  scale_fill_manual(values = c("1" = "#74C0FC", "2" = "#FFD43B", "3" = "#FF8FAB")) +
  theme_bw() +theme(legend.position="none")+ theme(aspect.ratio=1)+
  labs(x = "Cluster",
       y = "Specific Surface Area",
       fill = "Cluster")

# Carbon content
p2 <- ggplot(site_cluster_data_bulk, aes(x = factor(cluster), y = Median_X01395_C_percent_per_mg)) +
  geom_boxplot(aes(fill = factor(cluster))) +
  scale_fill_manual(values = c("1" = "#74C0FC", "2" = "#FFD43B", "3" = "#FF8FAB")) +
  theme_bw() +theme(legend.position="none")+theme(aspect.ratio=1)+
  labs(title = " ",
       x = "Cluster",
       y = "Carbon (%)",
       fill = "Cluster")

# ATP content
p3 <- ggplot(site_cluster_data_bulk, aes(x = factor(cluster), y = Median_ATP_picomoles_per_g)) +
  geom_boxplot(aes(fill = factor(cluster))) +
  scale_fill_manual(values = c("1" = "#74C0FC", "2" = "#FFD43B", "3" = "#FF8FAB")) +
  theme_bw() +theme(legend.position="none")+theme(aspect.ratio=1)+
  labs(title = " ",
       x = "Cluster",
       y = "ATP (picomoles/g)",
       fill = "Cluster")

# Sand content
p4 <- ggplot(site_cluster_data_bulk, aes(x = factor(cluster), y = Percent_Tot_Sand)) +
  geom_boxplot(aes(fill = factor(cluster))) +
  scale_fill_manual(values = c("1" = "#74C0FC", "2" = "#FFD43B", "3" = "#FF8FAB")) +
  theme_bw() + theme(legend.position="none")+theme(aspect.ratio=1)+
  labs(title = " ",
       x = "Cluster",
       y = "Total Sand (%)",
       fill = "Cluster")

p5 = ggplot(site_cluster_data_bulk, aes(x = factor(cluster), y = Median_Fe_mg_per_kg)) +
  geom_boxplot(aes(fill = factor(cluster))) +
  scale_fill_manual(values = c("1" = "#74C0FC", "2" = "#FFD43B", "3" = "#FF8FAB")) +
  theme_bw() + theme(legend.position="none")+theme(aspect.ratio=1)+
  labs(title = " ",
       x = "Cluster",
       y = "Fe (mg/kg)",
       fill = "Cluster")

p6 <- ggplot(site_cluster_data_bulk, aes(x = factor(cluster), y = Median_X01397_N_percent_per_mg)) +
  geom_boxplot(aes(fill = factor(cluster))) +
  scale_fill_manual(values = c("1" = "#74C0FC", "2" = "#FFD43B", "3" = "#FF8FAB")) +
  theme_bw() +theme(legend.position="none")+theme(aspect.ratio=1)+
  labs(title = " ",
       x = "Cluster",
       y = "Nitrogen (%)",
       fill = "Cluster")


# Doing some stats before exporting 
# Apply the function to each plot with appropriate variables
p1_enhanced <- add_kw_stats(p1, site_cluster_data_bulk, "Mean_Specific_Surface_Area_m2_per_g")
p2_enhanced <- add_kw_stats(p2, site_cluster_data_bulk, "Median_X01395_C_percent_per_mg")
p3_enhanced <- add_kw_stats(p3, site_cluster_data_bulk, "Median_ATP_picomoles_per_g")
p4_enhanced <- add_kw_stats(p4, site_cluster_data_bulk, "Percent_Tot_Sand")
p5_enhanced <- add_kw_stats(p5, site_cluster_data_bulk, "Median_Fe_mg_per_kg")
p6_enhanced <- add_kw_stats(p6, site_cluster_data_bulk, "Median_X01397_N_percent_per_mg")

# Arrange all plots in a grid - adjust layout as needed
grid_arranged_plots <- grid.arrange(
  p1_enhanced, p2_enhanced, p3_enhanced, 
  p4_enhanced, p5_enhanced, p6_enhanced,
  ncol = 2
)


# Save as PDF and PNG
ggsave(
  filename = "Figures/FigureS_Bulk_variables_within_clusters.pdf",
  plot = grid_arranged_plots,
  width = 10,
  height = 14,
  units = "in",
  dpi = 300
)

ggsave(
  filename = "Figures/FigureS_Bulk_variables_within_clusters.png",
  plot = grid_arranged_plots,
  width = 10,
  height = 14,
  units = "in",
  dpi = 300
)

# Effect size
site_cluster_data_effect <- data.frame(
  site = names(site_clusters),
  cluster = site_clusters
) %>%
  left_join(
    bulk_metadata,
    by = "site"
  )

# SSA
p1 <- ggplot(site_cluster_data_effect, aes(x = factor(cluster), y = Mean_Specific_Surface_Area_m2_per_g)) +
  geom_boxplot(aes(fill = factor(cluster))) +
  scale_fill_manual(values = c("1" = "#74C0FC", "2" = "#FFD43B", "3" = "#FF8FAB")) +
  theme_bw() +theme(legend.position="none")+ theme(aspect.ratio=1)+
  labs(x = "Cluster",
       y = "Specific Surface Area",
       fill = "Cluster")

# Carbon content
p2 <- ggplot(site_cluster_data_effect, aes(x = factor(cluster), y = Effect_Size_C_percent_per_mg)) +
  geom_boxplot(aes(fill = factor(cluster))) +
  scale_fill_manual(values = c("1" = "#74C0FC", "2" = "#FFD43B", "3" = "#FF8FAB")) +
  theme_bw() +theme(legend.position="none")+theme(aspect.ratio=1)+
  labs(title = " ",
       x = "Cluster",
       y = "Effect_Size_Carbon (%)",
       fill = "Cluster")

# ATP content
p3 <- ggplot(site_cluster_data_effect, aes(x = factor(cluster), y = Effect_Size_ATP_picomoles_per_g)) +
  geom_boxplot(aes(fill = factor(cluster))) +
  scale_fill_manual(values = c("1" = "#74C0FC", "2" = "#FFD43B", "3" = "#FF8FAB")) +
  theme_bw() +theme(legend.position="none")+theme(aspect.ratio=1)+
  labs(title = " ",
       x = "Cluster",
       y = "Effect_Size_ATP (picomoles/g)",
       fill = "Cluster")

# Sand content
p4 <- ggplot(site_cluster_data_effect, aes(x = factor(cluster), y = Percent_Tot_Sand)) +
  geom_boxplot(aes(fill = factor(cluster))) +
  scale_fill_manual(values = c("1" = "#74C0FC", "2" = "#FFD43B", "3" = "#FF8FAB")) +
  theme_bw() + theme(legend.position="none")+theme(aspect.ratio=1)+
  labs(title = " ",
       x = "Cluster",
       y = "Total Sand (%)",
       fill = "Cluster")

p5 = ggplot(site_cluster_data_effect, aes(x = factor(cluster), y = Effect_Size_Fe_mg_per_kg)) +
  geom_boxplot(aes(fill = factor(cluster))) +
  scale_fill_manual(values = c("1" = "#74C0FC", "2" = "#FFD43B", "3" = "#FF8FAB")) +
  theme_bw() + theme(legend.position="none")+theme(aspect.ratio=1)+
  labs(title = " ",
       x = "Cluster",
       y = "Effect_Size_Fe (mg/kg)",
       fill = "Cluster")

p6 <- ggplot(site_cluster_data_effect, aes(x = factor(cluster), y = Effect_Size_N_percent_per_mg)) +
  geom_boxplot(aes(fill = factor(cluster))) +
  scale_fill_manual(values = c("1" = "#74C0FC", "2" = "#FFD43B", "3" = "#FF8FAB")) +
  theme_bw() +theme(legend.position="none")+theme(aspect.ratio=1)+
  labs(title = " ",
       x = "Cluster",
       y = "Effect_Size_Nitrogen (%)",
       fill = "Cluster")


# Doing some stats before exporting 
# Apply the function to each plot with appropriate variables
p1_enhanced <- add_kw_stats(p1, site_cluster_data_effect, "Mean_Specific_Surface_Area_m2_per_g")
p2_enhanced <- add_kw_stats(p2, site_cluster_data_effect, "Effect_Size_C_percent_per_mg")
p3_enhanced <- add_kw_stats(p3, site_cluster_data_effect, "Effect_Size_ATP_picomoles_per_g")
p4_enhanced <- add_kw_stats(p4, site_cluster_data_effect, "Percent_Tot_Sand")
p5_enhanced <- add_kw_stats(p5, site_cluster_data_effect, "Effect_Size_Fe_mg_per_kg")
p6_enhanced <- add_kw_stats(p6, site_cluster_data_effect, "Effect_Size_N_percent_per_mg")

# Arrange all plots in a grid - adjust layout as needed
grid_arranged_plots <- grid.arrange(
  p1_enhanced, p2_enhanced, p3_enhanced, 
  p4_enhanced, p5_enhanced, p6_enhanced,
  ncol = 2
)


# Save as PDF and PNG
ggsave(
  filename = "Figures/FigureS_effect_variables_within_clusters.pdf",
  plot = grid_arranged_plots,
  width = 10,
  height = 14,
  units = "in",
  dpi = 300
)

ggsave(
  filename = "Figures/FigureS_effect_variables_within_clusters.png",
  plot = grid_arranged_plots,
  width = 10,
  height = 14,
  units = "in",
  dpi = 300
)


# ==== Heat map plot Field data =====
#Create main heatmap with all annotations and legends
pdf("Figures/Heatmap_Field_data_and_Clustes.pdf", width = 14, height = 10)
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
png("Figures/Heatmap_Field_data_and_Clustes.png", width = 4200, height = 3000, res = 300)
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


# ==== Heat map plot Bulk data =====
#Create main heatmap with all annotations and legends
pdf("Figures/Heatmap_Bulk_data_and_Clustes.pdf", width = 14, height = 10)
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
  annotation_row = bulk_annotation,
  annotation_colors = bulk_colors,
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
png("Figures/Heatmap_Bulk_data_and_Clustes.png", width = 4200, height = 3000, res = 300)
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
  annotation_row = bulk_annotation,
  annotation_colors = bulk_colors,
  cutree_rows = 3,
  annotation_legend = TRUE,
  legend = TRUE,
  fontsize_row = 8,
  cellwidth = 15,
  cellheight = 12,
  border_color = NA
)
dev.off()

# ==== Heat map plotEffect data =====
#Create main heatmap with all annotations and legends
pdf("Figures/Heatmap_Effect_data_and_Clustes.pdf", width = 14, height = 10)
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
  annotation_row = effect_annotation,
  annotation_colors = effect_colors,
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
png("Figures/Heatmap_Effect_data_and_Clustes.png", width = 4200, height = 3000, res = 300)
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
  annotation_row = effect_annotation,
  annotation_colors = effect_colors,
  cutree_rows = 3,
  annotation_legend = TRUE,
  legend = TRUE,
  fontsize_row = 8,
  cellwidth = 15,
  cellheight = 12,
  border_color = NA
)
dev.off()

# ===== Identify if clusters differ in their site characteristics using field data =====
# Seems that field data is the one that provides the most variability within clusters so only looking at field data

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

# Principal Component Analysis (PCA)
# Standardize data (important for variables with different scales)
pca_result <- PCA(env_vars, scale.unit = TRUE, graph = FALSE)

# Calculate PERMANOVA
perm_result <- adonis2(env_vars ~ cluster, data = site_data_clean,
                       method = "euclidean", permutations = 999)

# Format PERMANOVA results for display
perm_f_value <- round(perm_result$F[1], 2)
perm_r2 <- round(perm_result$R2[1], 3)
perm_p_value <- perm_result$`Pr(>F)`[1]

# Format p-value with appropriate significance stars
if(perm_p_value < 0.001) {
  p_text <- "p < 0.001 ***"
} else if(perm_p_value < 0.01) {
  p_text <- paste0("p = ", round(perm_p_value, 3), " **")
} else if(perm_p_value < 0.05) {
  p_text <- paste0("p = ", round(perm_p_value, 3), " *")
} else {
  p_text <- paste0("p = ", round(perm_p_value, 3), " ")
}

# Create annotation text
perm_text <- paste0("PERMANOVA: ",  p_text)

# Get eigenvalues to calculate variance explained
eigenvalues <- pca_result$eig
pc1_var <- round(eigenvalues[1, 2], 1)  # % variance explained by PC1
pc2_var <- round(eigenvalues[2, 2], 1)  # % variance explained by PC2

# Create enhanced PCA plot
pca_plot <- fviz_pca_biplot(pca_result,
                            habillage = site_data_clean$cluster,  # Color by cluster
                            geom = "point",
                            pointsize = 3,
                            palette = c("#74C0FC", "#FFD43B", "#FF8FAB"),
                            addEllipses = TRUE,
                            ellipse.type = "confidence",
                            ggtheme = theme_bw(),
                            title = " ") +
  # Add axes labels with variance explained
  xlab(paste0("PC1 (", pc1_var, "% explained var.)")) +
  ylab(paste0("PC2 (", pc2_var, "% explained var.)")) +
  # Add PERMANOVA results
  annotate("text", x = -2 ,
           y = -2,
           label = perm_text,
           hjust = 0,
           size = 3.5,
           fontface = "bold") +
  # Improve legend
  theme(legend.title = element_text(face = "bold"),
        legend.position = "right",
        plot.title = element_text(face = "bold", size = 14),
        panel.grid.minor = element_blank())

# Print plot to screen
print(pca_plot)

# Save as PDF
ggsave(
  filename = "Figures/FigureS_PCA_cluster_field_data.pdf",
  plot = pca_plot,
  width = 10,
  height = 8,
  units = "in",
  dpi = 300
)

# Save as PNG
ggsave(
  filename = "Figures/FigureS_PCA_cluster_field_data.png",
  plot = pca_plot,
  width = 10,
  height = 8,
  units = "in",
  dpi = 300
)

var_contrib <- fviz_contrib(pca_result, choice = "var", axes = 1:2, top = 6)
print(var_contrib)

# Loading scores of variables on principal components
print(pca_result$var$coord)

cluster_profiles <- site_data_clean %>%
  group_by(cluster) %>%
  summarise(dplyr::across(c(Mean_ATP_picomoles_per_g, Mean_Specific_Surface_Area_m2_per_g,
                     C_percent_per_mg, N_percent_per_mg, 
                     Percent_Tot_Sand, Mean_Fe_mg_per_kg), 
                   list(mean = ~mean(., na.rm = TRUE),
                        sd = ~sd(., na.rm = TRUE))))

print(cluster_profiles)
