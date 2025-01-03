# ==== Loading libraries =========
rm(list=ls(all=T))

library(stringr); library(devtools);  library("plyr")
library("readr");  library(readxl);library(crayon); library(vegan)
# Load in necessary libraries first
library(reshape2)
library(ggpubr) # For to combine plots
library(dplyr) # For reorganization
library(stringr) # For string manipulation
# ==== Defining paths and working directories ======
github = 'C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/ECA_DOM_Thermodynamics/'
# ====== Read in data ======
# Processed ICR Data
data = read.csv(list.files(path = github, pattern = "*unique_formulas_Data.csv", full.names = T),row.names = 1)
mol = read.csv(list.files(path = github, pattern = "*cal_pts_Mol.csv"), row.names = 1)
# Fixing colnames 
colnames(data) = gsub('SIR.','SIR-',colnames(data))

sample_data = read_csv(paste0(github,'EC_Data_Package/Sample_Data/EC_Sediment_Sample_Data_Summary.csv'),comment = '#', na = c('N/A', -9999)) %>%
  slice(-(1:11))%>%
  mutate_at(vars(-Sample_Name,-Field_Name,-IGSN,-Material), as.numeric)

effect_size = read_csv(paste0(github,'EC_Data_Package/Sample_Data/EC_Sediment_Effect_Size.csv'),comment = '#', na = c('N/A', -9999)) %>%
  slice(-(1:11))%>%
  dplyr::select('Sample_Name',"Effect_Size_Respiration_Rate_mg_DO_per_L_per_H","Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H","Effect_Size_Initial_Gravimetric_Moisture_g_per_g","Effect_Size_Final_Gravimetric_Moisture_g_per_g","Effect_Size_Extractable_NPOC_mg_per_kg","Effect_Size_Extractable_TN_mg_per_L")%>%
  mutate_at(vars(-Sample_Name), as.numeric)


effect_size$site = effect_size$Sample_Name
effect_size$site = gsub('_all','',effect_size$site)

factors = data.frame(Sample_Name = sample_data$Sample_Name, site = sample_data$Sample_Name, Treatment = sample_data$Sample_Name)
factors$site = str_extract(factors$site, "EC_0[0-9]{2}|EC_([A-Za-z0-9]+)")
factors$Treatment = str_extract(factors$Treatment, "W|D|Blk")

factors$Treatment = gsub('W','Wet', factors$Treatment)
factors$Treatment = gsub('D','Dry', factors$Treatment)

sample_data = merge(factors,sample_data, by = 'Sample_Name')

site_order = read.csv('Site_order.csv')
# ========= Plots ======
# Select the mol variables of interest
mol2 = mol %>% dplyr::select(AI_Mod,DBE_1,NOSC,delGcoxPerCmol,delGd,lamO2)

# Calculating weighed avg metris and thermodynamics per sample
df.merge = merge(data,mol2, by = 'row.names')

df.stats = as.data.frame(matrix(NA, nrow = ncol(data), ncol = 7))
colnames(df.stats) = c('Sample_Name','Weighted_Avg_AI_mod','Weighted_Avg_NOSC','Weighted_Avg_DBE','Weighted_Avg_delGcoxPerCmol','Weighted_Avg_delGcoxPerCompmol','Weighted_Avg_Lambda')

for (i in 2:(ncol(data)+1)){
  df.stats$Sample_Name[i-1] = colnames(df.merge[i])
  df.stats$Weighted_Avg_delGcoxPerCmol[i-1] = weighted.mean(df.merge$delGcoxPerCmol[which(df.merge[, i] > 0)],df.merge[,i][which(df.merge[, i] > 0)])
  df.stats$Weighted_Avg_delGcoxPerCompmol[i-1] = weighted.mean(df.merge$delGd[which(df.merge[, i] > 0)],df.merge[,i][which(df.merge[, i] > 0)])
  df.stats$Weighted_Avg_Lambda[i-1] = weighted.mean(df.merge$lamO2[which(df.merge[, i] > 0)],df.merge[,i][which(df.merge[, i] > 0)])
  df.stats$Weighted_Avg_AI_mod[i-1] = weighted.mean(df.merge$AI_Mod[which(df.merge[, i] > 0)],df.merge[,i][which(df.merge[, i] > 0)])
  df.stats$Weighted_Avg_NOSC[i-1] = weighted.mean(df.merge$NOSC[which(df.merge[, i] > 0)],df.merge[,i][which(df.merge[, i] > 0)])
  df.stats$Weighted_Avg_DBE[i-1] = weighted.mean(df.merge$DBE_1[which(df.merge[, i] > 0)],df.merge[,i][which(df.merge[, i] > 0)])
}

# Adding factors manually to df.stats to not complicate merging
df.stats$site = str_extract(df.stats$Sample_Name, "EC_0[0-9]{2}|EC_([A-Za-z0-9]+)")
df.stats$Treatment = str_extract(df.stats$Sample_Name, "W|D|Blk")

df.stats$Treatment = gsub('W','Wet', df.stats$Treatment)
df.stats$Treatment = gsub('D','Dry', df.stats$Treatment)

write.csv(df.stats,'Weighted_averages_for_molecular_properties_per_sample.csv', row.names = F)

wilcox_test <- function(df, variable, group_var) {
  results <- data.frame()
  unique_sites <- unique(df[[group_var]])
  for (site in unique_sites) {
    site_data <- df[df[[group_var]] == site, ]
    wet_data <- site_data[site_data$Treatment == 'Wet', variable]
    dry_data <- site_data[site_data$Treatment == 'Dry', variable]
    
    if (length(wet_data) > 0 & length(dry_data) > 0) {
      if (length(wet_data) > 1 & length(dry_data) > 1) {
        test_result <- wilcox.test(wet_data, dry_data, exact = FALSE)
        results <- rbind(results, data.frame(Site = site, Variable = variable, 
                                             W = test_result$statistic, 
                                             p.value = test_result$p.value))
      } else {
        results <- rbind(results, data.frame(Site = site, Variable = variable, 
                                             W = NA, 
                                             p.value = NA))
      }
    } else {
      results <- rbind(results, data.frame(Site = site, Variable = variable, 
                                           W = NA, 
                                           p.value = NA))
    }
  }
  return(results)
}

variables_to_test <- c('Weighted_Avg_AI_mod', 'Weighted_Avg_NOSC', 'Weighted_Avg_DBE', 'Weighted_Avg_delGcoxPerCmol','Weighted_Avg_delGcoxPerCompmol', 'Weighted_Avg_Lambda')
wilcox_results <- data.frame()

# Perform Wilcoxon tests for each variable and combine the results
for (variable in variables_to_test) {
  wilcox_results <- rbind(wilcox_results, wilcox_test(df.stats, variable, 'site'))
}

#Filter for significant p-values (e.g., p < 0.05)
significant_results <- wilcox_results %>% filter(p.value < 0.05)

# Count how many sites have significant differences across multiple mol values
significant_count <- significant_results %>%
  group_by(Site) %>%
  summarise(Num_Significant_Variables = n())

df_filtered <- df.stats %>% filter(site %in% significant_count$Site)


ggplot(df_filtered, aes(x = site, y = Gibbs, fill = Treatment)) +
  geom_boxplot() +
  labs(x = "Location", y = "Intensity weighted avg Gibbs per C mol Unique peaks") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(df_filtered, aes(x = site, y = Lambda, fill = Treatment)) +
  geom_boxplot() +
  labs(x = "Location", y = "Intensity weighted avg Lambda Unique peaks") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



ggplot(df.stats, aes(x = site, y = Gibbs, fill = Treatment)) +
  geom_boxplot() +
  labs(x = "Location", y = "Intensity weighted avg Gibbs per C mol Unique peaks") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(df.stats, aes(x = site, y = Lambda, fill = Treatment)) +
  geom_boxplot() +
  labs(x = "Location", y = "Intensity weighted avg Lambda Unique peaks") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ===== Bar plots with average values for mol characteristics for each site ====

df_melted <- df.stats %>%
  pivot_longer(
    cols = starts_with("Weighted_Avg"),  # Choose columns to pivot (all those starting with "Weighted_Avg")
    names_to = "Property",               # New column for the names of the properties
    values_to = "Value"                  # New column for the values of the properties
  )

error_data <- df_melted %>%
  group_by(site, Treatment,Property) %>%
  summarise(
    mean = mean(Value),
    sd = sd(Value),
    se = sd / sqrt(n()),  # Standard error of the mean
    .groups = "drop"
  )


ggplot(error_data, aes(x = site, y = mean, fill = Treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black") +  # Bar plot
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),  # Error bars based on SE
    position = position_dodge(width = 0.9),  # Match dodge width of bars
    width = 0.2  # Error bar width
  ) +
  labs(
    x = "Site",
    y = "Average Molecular Property\nwithin a site and treatment",
    fill = "Treatment"
  ) +
  scale_fill_manual(values = c("Wet" = "lightblue", "Dry" = "darkorange")) +  # Custom colors
  theme_bw() +  # Minimal theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top") +  # Rotate x-axis labels for readability
  facet_wrap(~ Property, scales = "free_y")  # Facet by Property with free y scale

# Bar plots for sites in group 1

sites_in_order_group_1 <- site_order %>%
  filter(order_group == 1) %>%  # Only keep sites in order_group 1
  pull(site)  # Extract the site names


sites_in_order_group_2 <- site_order %>%
  filter(order_group == 2) %>%  # Only keep sites in order_group 1
  pull(site)  # Extract the site names

filtered_data <- error_data %>%
  filter(site %in% sites_in_order_group_1)  

ggplot(filtered_data, aes(x = site, y = mean, fill = Treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black") +  # Bar plot
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),  # Error bars based on SE
    position = position_dodge(width = 0.9),  # Match dodge width of bars
    width = 0.2  # Error bar width
  ) +
  labs(
    x = "Site",
    y = "Average Molecular Property\nwithin a site and treatment",
    fill = "Treatment"
  ) +
  scale_fill_manual(values = c("Wet" = "lightblue", "Dry" = "darkorange")) +  # Custom colors
  theme_bw() +  # Minimal theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top") +  # Rotate x-axis labels for readability
  facet_wrap(~ Property, scales = "free_y")  # Facet by Property with free y scale

ggplot(df_melted, aes(x = site, y = Value, fill = Treatment)) +
  geom_boxplot() +
  labs(
    x = "Site",
    y = "Average Molecular Property\nwithin a site and treatment",
    fill = "Treatment"
  ) +
  scale_fill_manual(values = c("Wet" = "lightblue", "Dry" = "darkorange")) +  # Custom colors
  theme_bw() +  # Minimal theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top") +  # Rotate x-axis labels for readability
  facet_wrap(~ Property, scales = "free_y")  # Facet by Property with free y scale

filtered_data <- df_melted %>%
  filter(site %in% sites_in_order_group_1)  

ggplot(filtered_data, aes(x = Treatment, y = Value, fill = Treatment)) +
  geom_boxplot() +
  labs(
    x = "Site",
    y = "Average Molecular Property\nwithin a site and treatment",
    fill = "Treatment"
  ) +
  scale_fill_manual(values = c("Wet" = "lightblue", "Dry" = "darkorange")) +  # Custom colors
  theme_bw() +  # Minimal theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top") +  # Rotate x-axis labels for readability
  facet_wrap(~ Property, scales = "free_y")  # Facet by Property with free y scale

filtered_data <- sample_data %>%
  filter(site %in% sites_in_order_group_1)  

filtered_data2 <- sample_data %>%
  filter(site %in% sites_in_order_group_2)  

effect_size_ratio2 <- filtered_data2 %>%
  group_by(site, Treatment) %>%
  summarise(
    median_moisture = median(Median_62948_Final_Gravimetric_Moisture_g_per_g, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  spread(key = Treatment, value = median_moisture) %>%  # Spread the data for Wet and Dry
  mutate(
    effect_size = Wet / Dry  # Calculate the ratio of Wet to Dry median moisture
  )

effect_size_ratio1 <- filtered_data %>%
  group_by(site, Treatment) %>%
  summarise(
    median_moisture = median(Median_62948_Final_Gravimetric_Moisture_g_per_g, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  spread(key = Treatment, value = median_moisture) %>%  # Spread the data for Wet and Dry
  mutate(
    effect_size = Wet - Dry  # Calculate the ratio of Wet to Dry median moisture
  )

# ===== Calculate a median of these metrics per site and treatment ====
df.medians = df.stats %>%
  group_by(site,Treatment) %>%
  summarise(
    Median_Weighted_Avg_AI_mod = median(Weighted_Avg_AI_mod, na.rm = TRUE),
    Median_Weighted_Avg_NOSC = median(Weighted_Avg_NOSC, na.rm = TRUE),
    Median_Weighted_Avg_DBE = median(Weighted_Avg_DBE, na.rm = TRUE),
    Median_Weighted_Avg_delGcoxPerCmol = median(Weighted_Avg_delGcoxPerCmol, na.rm = TRUE),
    Median_Weighted_Avg_delGcoxPerCompmol = median(Weighted_Avg_delGcoxPerCmol, na.rm = TRUE),
    Median_Weighted_Avg_Lambda = median(Weighted_Avg_Lambda, na.rm = TRUE)
  )

write.csv(df.medians, 'Medians_of_Weighted_averages_for_molecular_properties_per_site_and_treatment.csv', row.names = F)

# ==== PCA of median mol properties =====
library(ggfortify)  # For PCA visualization
pca_data <- df.medians %>%
  select(starts_with("Median_")) %>%
  mutate(across(everything(), as.numeric))

pca_data = pca_data[,2:7]

# Step 2: Perform PCA
pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)

pca_scores <- as.data.frame(pca_result$x)
pca_scores <- cbind(pca_scores, df.medians[, c("site", "Treatment")])
# Extract PCA loadings and scale them to fit the plot
loadings <- as.data.frame(pca_result$rotation[, 1:2])  # Loadings for PC1 and PC2
loadings$Variable <- rownames(loadings)  # Add variable names
loadings <- loadings %>%
  mutate(PC1 = PC1 * max(abs(pca_scores$PC1)) * 0.3,  # Scale loadings to match PCA plot
         PC2 = PC2 * max(abs(pca_scores$PC2)) * 0.3)

library(Polychrome)

# Generate a palette for 40 sites
color_palette <- createPalette(50, c("#000000", "#FFFFFF"), M = 500)
names(color_palette) <- unique(pca_scores$site)  # Assign site names to colors

# PCA Plot with Custom Colors
ggplot(pca_scores, aes(x = PC1, y = PC2, color = site, shape = Treatment)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_segment(
    data = loadings, 
    aes(x = 0, y = 0, xend = PC1, yend = PC2), 
    arrow = arrow(length = unit(0.2, "cm")), 
    inherit.aes = FALSE,
    color = "black", 
    alpha = 0.8
  ) +
  geom_text(
    data = loadings, 
    aes(x = PC1, y = PC2, label = Variable),
    inherit.aes = FALSE,
    size = 4, vjust = 1, hjust = 1, color = "black"
  ) +
  scale_color_manual(values = color_palette) +  # Use the Polychrome palette
  labs(
    x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "% Variance Explained)"),
    y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "% Variance Explained)"),
    color = "Site",
    shape = "Treatment"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

# ==== Bar plots with medians =====
df_melted <- df.stats %>%
  pivot_longer(
    cols = starts_with("Weighted_Avg"),  # Choose columns to pivot (all those starting with "Weighted_Avg")
    names_to = "Property",               # New column for the names of the properties
    values_to = "Value"                  # New column for the values of the properties
  )

# Step 2: Summarize the data by site, treatment, and property using median and IQR
error_data <- df_melted %>%
  group_by(site, Treatment, Property) %>%
  summarise(
    median_value = median(Value),           # Median value of the property
    IQR_value = IQR(Value),                 # Interquartile range (IQR) as error
    .groups = "drop"
  )

# Step 3: Plot the data with IQR as error bars
ggplot(error_data, aes(x = site, y = median_value, fill = Treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black") +  # Bar plot
  geom_errorbar(
    aes(ymin = median_value - IQR_value / 2, ymax = median_value + IQR_value / 2),  # Error bars based on IQR
    position = position_dodge(width = 0.9),  # Match dodge width of bars
    width = 0.2  # Error bar width
  ) +
  labs(
    x = "Site",
    y = "Median Molecular Property\nwithin a site and treatment",
    fill = "Treatment"
  ) +
  scale_fill_manual(values = c("Wet" = "lightblue", "Dry" = "darkorange")) +  # Custom colors
  theme_bw() +  # Minimal theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top") +  # Rotate x-axis labels for readability
  facet_wrap(~ Property, scales = "free_y")  # Facet by Property with free y scale

# ===== Rations and histograms =====
df_medians_melted <- df.medians %>%
  pivot_longer(
    cols = starts_with("Median_"),  
    names_to = "Property",               # New column for the names of the properties
    values_to = "Value"                  # New column for the values of the properties
  )

pivot_data <- df_medians_melted %>%
  pivot_wider(
    names_from = Treatment,  # Spread out the Treatment variable into separate columns for Wet and Dry
    values_from = Value      # Values from the Value column
  )

# Step 2: Filter out rows with NA in Wet or Dry columns and calculate ratio
ratio_data <- pivot_data %>%
  filter(!is.na(Wet) & !is.na(Dry)) %>%  # Keep only rows where both Wet and Dry values are present
  mutate(
    ratio_wet_dry = Wet / Dry  # Calculate Wet/Dry ratio for each property
  ) %>%
  select(site, Property, ratio_wet_dry)  # Keep necessary columns (site, Property, ratio)



ggplot(ratio_data, aes(x = ratio_wet_dry)) +
  geom_histogram(binwidth = 0.05, color = "black", fill = "lightblue", alpha = 0.7) +  # Smaller binwidth for more bins
  labs(
    x = "Wet/Dry Ratio",
    y = "Frequency",
    title = "Histogram of Wet/Dry Ratios for Each Property"
  ) +
  facet_wrap(~ Property, scales = "free") +  # Free scales for both x and y axes
  theme_bw() +  # Minimal theme
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 12),
    legend.position = "top"
  )
# ==== ADD ====
df.all = merge(sample_data,df.medians, by = 'Sample_Name')%>%
  mutate(
    Median_Respiration_Rate_mg_DO_per_kg_per_H = ifelse(        Median_Respiration_Rate_mg_DO_per_kg_per_H == 0,
                                                                -min(abs(Median_Respiration_Rate_mg_DO_per_kg_per_H[
                                                                  Median_Respiration_Rate_mg_DO_per_kg_per_H != 0 & 
                                                                    !is.na(Median_Respiration_Rate_mg_DO_per_kg_per_H)
                                                                ])) / 2, 
                                                                Median_Respiration_Rate_mg_DO_per_kg_per_H
    ),                                       
    Cubic_Root_Respiration = abs(Median_Respiration_Rate_mg_DO_per_kg_per_H)^(1/3)
  )%>%
  filter(!is.na(Cubic_Root_Respiration))

ggplot(df.all, aes(x = Median_Extractable_NPOC_mg_per_kg, y = Gibbs_median, color = Treatment)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(x = "Extractable_NPOC_mg_per_kg", y = "Gibbs Median") +
  theme_bw() +
  theme(legend.position = 'top')

ggplot(df.all, aes(x = Median_Extractable_NPOC_mg_per_kg, y = Cubic_Root_Respiration, color = Treatment)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(x = "Extractable_NPOC_mg_per_L", y = "Cubic Root Respiration") +
  theme_bw() +
  theme(legend.position = 'top')


ggplot(df.all, aes(x = (1/Median_Extractable_NPOC_mg_per_kg), y = Cubic_Root_Respiration, color = Treatment)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(x = "1/Extractable_NPOC_mg_per_L", y = "Cubic Root Respiration") +
  theme_bw() +
  theme(legend.position = 'top')
########################################
# ===== Some plots with the different sites ======
# Calculate effect size as the ratio of W/D median

df.medians <- df.medians %>%
  mutate(Treatment = ifelse(str_detect(Sample_Name, "-W"), 'Wet', 'Dry'))

# Split the data into Wet and Dry
df_wet <- df.medians %>% filter(Treatment == 'Wet')
df_dry <- df.medians %>% filter(Treatment == 'Dry')

# Rename columns to indicate treatment type
df_wet <- df_wet %>%
  rename_with(~ paste0(., "_Wet"), -c(site, Sample_Name, Treatment))

df_dry <- df_dry %>%
  rename_with(~ paste0(., "_Dry"), -c(site, Sample_Name, Treatment))

# Join the Wet and Dry dataframes by site
df_joined <- left_join(df_wet, df_dry, by = "site")

# Calculate the ratios
df_ratios <- df_joined %>%
  mutate(
    AI_mod_ratio = AI_mod_median_Wet / AI_mod_median_Dry,
    NOSC_ratio = NOSC_median_Wet / NOSC_median_Dry,
    DBE_ratio = DBE_median_Wet / DBE_median_Dry,
    Gibbs_ratio = Gibbs_median_Wet / Gibbs_median_Dry,
    Lambda_ratio = Lambda_median_Wet / Lambda_median_Dry
  ) %>%
  select(site, AI_mod_ratio, NOSC_ratio, DBE_ratio, Gibbs_ratio, Lambda_ratio)

#############################
# ===== Scatter plots ====
ggplot(df.all, aes(x = AI_mod_median, y = Cubic_Root_Respiration, color = Treatment)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(x = "AI_mod Median", y = "Cubic Root Respiration") +
  theme_bw() +
  theme(legend.position = 'top')

ggplot(df.all, aes(x = NOSC_median, y = Cubic_Root_Respiration, color = Treatment)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(x = "NOSC Median", y = "Cubic Root Respiration") +
  theme_bw() +
  theme(legend.position = 'top')

ggplot(df.all, aes(x = DBE_median, y = Cubic_Root_Respiration, color = Treatment)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(x = "DBE Median", y = "Cubic Root Respiration") +
  theme_bw() +
  theme(legend.position = 'top')

ggplot(df.all, aes(x = Gibbs_median, y = Cubic_Root_Respiration, color = Treatment)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(x = "Gibbs Median", y = "Cubic Root Respiration") +
  theme_bw() +
  theme(legend.position = 'top')

ggplot(df.all, aes(x = Lambda_median, y = Cubic_Root_Respiration, color = Treatment)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(x = "Lambda Median", y = "Cubic Root Respiration") +
  theme_bw() +
  theme(legend.position = 'top')

ggplot(df.all, aes(x = Median_62948_Final_Gravimetric_Moisture_g_per_g, y = Gibbs_median, color = Treatment)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(x = "Final Gravimetric Moisture (g per g)", y = "Gibbs Median") +
  theme_bw() +
  theme(legend.position = 'top')

ggplot(df.all, aes(x = Median_62948_Final_Gravimetric_Moisture_g_per_g, y = Lambda_median, color = Treatment)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(x = "Final Gravimetric Moisture (g per g)", y = "Lambda Median") +
  theme_bw() +
  theme(legend.position = 'top')


ggplot(df.all, aes(x = Median_Extractable_NPOC_mg_per_kg, y = Gibbs_median, color = Treatment)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(x = "Extractable_NPOC_mg_per_kg", y = "Gibbs Median") +
  theme_bw() +
  theme(legend.position = 'top')

ggplot(df.all, aes(x = Median_Extractable_NPOC_mg_per_kg, y = Lambda_median, color = Treatment)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(x = "Extractable_NPOC_mg_per_kg", y = "Lambda Median") +
  theme_bw() +
  theme(legend.position = 'top')

ggplot(df.all, aes(x = Median_Extractable_NPOC_mg_per_L, y = Gibbs_median, color = Treatment)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(x = "Extractable_NPOC_mg_per_L", y = "Gibbs Median") +
  theme_bw() +
  theme(legend.position = 'top')

ggplot(df.all, aes(x = Median_Extractable_NPOC_mg_per_L, y = Lambda_median, color = Treatment)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(x = "Extractable_NPOC_mg_per_L", y = "Lambda Median") +
  theme_bw() +
  theme(legend.position = 'top')
##############################
# Calculate median per Location and Treatment
median_data <- df %>%
  group_by(Location, Treatment) %>%
  summarise_at(median)

# Calculate ratio between W and D
median_data <- median_data %>%
  mutate(ratio_W_to_D = W / D)

# Plot histogram of the ratios
ggplot(median_data, aes(x = ratio_W_to_D)) +
  geom_histogram(binwidth = 0.1, fill = "skyblue", color = "black") +
  labs(x = "Ratio of W to D", y = "Frequency") +
  theme_minimal()

cols_of_interest <- 2:6  # Columns of interest
plots <- list()  # Initialize a list to store plots

for (col in cols_of_interest) {
  plot <- ggplot(df, aes_string(x = colnames(df)[col], y = "log_rate", color = "Treatment")) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +  # Add linear regression line
    facet_wrap(~ Treatment) +  # Separate plots by Treatment
    labs(x = colnames(df)[col], y = "Log Rate") +  # Label axes
    theme_minimal()
  
  plots[[as.character(col)]] <- plot  # Store each plot in the list
}

# Composite plot with all individual plots
composite_plot <- cowplot::plot_grid(plotlist = plots, align = "hv")

df2 = subset(df, df$log_rate<= 1)
for (col in cols_of_interest) {
  plot <- ggplot(df2, aes_string(x = colnames(df)[col], y = "log_rate", color = "Treatment")) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +  # Add linear regression line
    facet_wrap(~ Treatment) +  # Separate plots by Treatment
    labs(x = colnames(df2)[col], y = "Log Rate") +  # Label axes
    theme_minimal()
  
  plots[[as.character(col)]] <- plot  # Store each plot in the list
}

# Composite plot with all individual plots
composite_plot <- cowplot::plot_grid(plotlist = plots, align = "hv")


for (col in cols_of_interest) {
  plot <- ggplot(df2, aes_string(x = colnames(df)[col], y = 'Respiration_Rate_mg_DO_per_L_per_H', color = "Treatment")) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +  # Add linear regression line
    facet_wrap(~ Treatment) +  # Separate plots by Treatment
    labs(x = colnames(df2)[col], y = "Rate") +  # Label axes
    theme_minimal()
  
  plots[[as.character(col)]] <- plot  # Store each plot in the list
}

# Composite plot with all individual plots
composite_plot <- cowplot::plot_grid(plotlist = plots, align = "hv")
