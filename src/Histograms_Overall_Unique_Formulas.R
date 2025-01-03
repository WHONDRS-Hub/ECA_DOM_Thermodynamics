# ==== Loading libraries =========
rm(list=ls(all=T))

library(stringr); library(devtools);  library("plyr")
library("readr"); library(tidyr); library(readxl);library(crayon); library(vegan)
# Load in necessary libraries first
library(reshape2)
library(ggpubr) # For to combine plots
library(dplyr) # For reorganization
library(stringr) # For string manipulation
# ==== Defining paths and working directories ======
github = 'C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/ECA_DOM_Thermodynamics/'

# ====== Read in and clean up data ======
# Processed ICR Data

data = read.csv(list.files(path = github, pattern = "*unique_formulas_Data.csv", full.names = T),row.names = 1)
mol = read.csv(list.files(path = github, pattern = "*cal_pts_Mol.csv"), row.names = 1)
# Fixing colnames 
colnames(data) = gsub('SIR.','SIR-',colnames(data))

site_order = read.csv('Site_order.csv')

# ========= Plots ======

# Renaming bs1_classes
mol$bs1_class[grep(";", mol$bs1_class)] = "Multi-class"

# Adding mass
mol$Mass = as.numeric(as.character(row.names(mol)))

# Creating factors sheet
factors = data.frame(Samples = colnames(data), Location = colnames(data), Treatment = colnames(data))
factors$Location = str_extract(factors$Location, "EC_0[0-9]{2}|EC_([A-Za-z0-9]+)")
factors$Treatment = str_extract(factors$Treatment, "W|D|Blk")

sample_info <- data.frame(
  sample = colnames(data),
  site = str_extract(colnames(data), "EC_[A-Z0-9]+"),
  treatment = case_when(grepl("W", colnames(data)) ~ "Wet",
                        grepl("D", colnames(data)) ~ "Dry",TRUE~"Samples"))
 
# === Density plots of OM properties ====

# === Function to process and prepare the data for Wet and Dry treatments ===
process_treatment_data <- function(data, mol, treatment, property_cols) {
  
  # Subset data based on Wet or Dry treatment
  treatment_data <- data %>% select(matches(treatment))
  treatment_mol <- mol %>% select(all_of(property_cols))
  
  # Merge treatment data with molecular properties
  df_merge <- merge(treatment_data,treatment_mol, by = 'row.names', all = FALSE)
  
  # Clean data by removing rows with all zeros (or NAs in some cases)
  df_cleaned <- df_merge %>%
    pivot_longer(cols = starts_with("EC"), names_to = "Sample", values_to = "Intensity") %>%
    filter(Intensity > 0) %>%
    pivot_longer(cols = c(property_cols), names_to = "Property", values_to = "Value") %>%
    mutate(Treatment = treatment)
  
  return(df_cleaned)
}

# === Process both Wet and Dry treatments for OM properties ===
om_properties <- c("AI_Mod", "DBE_1", "NOSC")

# Process Wet and Dry treatment data for OM properties
wet_om_data <- process_treatment_data(data, mol, "W", om_properties)
dry_om_data <- process_treatment_data(data, mol, "D", om_properties)

# Combine Wet and Dry data
om_combined_data <- bind_rows(wet_om_data, dry_om_data)

# === Process both Wet and Dry treatments for Thermodynamic properties ===
thermo_properties <- c("delGcoxPerCmol", "lamO2", "delGd")

# Process Wet and Dry treatment data for Thermodynamic properties
wet_thermo_data <- process_treatment_data(data, mol, "W", thermo_properties)
dry_thermo_data <- process_treatment_data(data, mol, "D", thermo_properties)

# Combine Wet and Dry thermodynamic data
thermo_combined_data <- bind_rows(wet_thermo_data, dry_thermo_data)

# === Combine both OM and Thermodynamic data into one dataframe ===
final_combined_data <- bind_rows(om_combined_data, thermo_combined_data)

final_combined_data$Treatment = gsub('W','Wet',final_combined_data$Treatment)
final_combined_data$Treatment = gsub('D','Dry',final_combined_data$Treatment)

metabolite_count_data <- data %>%
  summarise(across(everything(), ~sum(. != 0))) %>%  # Summing non-zero values
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "MetaboliteCount")

# Merge metabolite count with treatment information (Wet vs Dry)
metabolite_count_data <- metabolite_count_data %>%
  left_join(sample_info, by = c("Sample" = "sample"))

# === Combine Metabolite Count with other data for final plot ===
final_combined_data_with_metabolite <- bind_rows(final_combined_data, 
                                                 data.frame(Property = "Metabolite Count",
                                                            Value = metabolite_count_data$MetaboliteCount,
                                                            Sample = metabolite_count_data$Sample,
                                                            Treatment = metabolite_count_data$treatment))

# === Plot Histograms for each Property ===
ggplot(final_combined_data_with_metabolite, aes(x = Value, fill = Treatment)) + 
 # geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  geom_histogram(aes(fill = Treatment), alpha = 1, color = 'black')+
  facet_wrap(~ Property, scales = "free") +
  scale_fill_manual(values = c("darkorange","lightblue")) +
  theme_bw() +
  labs(x = " ", y = "Frequency", fill = "Treatment") +
  theme(legend.position = "top") 
  #ggtitle("Histograms of OM and Thermodynamic Properties by Treatment")

# ==== Box plots for sites in group 1 (higher counts on wet) ===
# Step 1: Filter the sites in order_group == 1 from the site_order dataframe
sites_in_order_group_1 <- site_order %>%
  filter(order_group == 1) %>%  # Only keep sites in order_group 1
  pull(site)  # Extract the site names

final_combined_data_with_metabolite = final_combined_data_with_metabolite %>%
  mutate(Site = str_extract(final_combined_data_with_metabolite$Sample, "EC_0[0-9]{2}|EC_([A-Za-z0-9]+)"))
# Step 2: Filter final_combined_data_with_metabolite to include only sites in order_group 1
filtered_data <- final_combined_data_with_metabolite %>%
  filter(Site %in% sites_in_order_group_1)  # Keep only the samples (sites) in order_group 1

# Step 3: Plot histograms for the filtered data
ggplot(filtered_data, aes(x = Value,  fill = Treatment)) + 
  # geom_histogram(aes(fill = Treatment), alpha = 1, color = 'black') +
  geom_density()+
  facet_wrap(~ Property, scales = "free") +
  scale_fill_manual(values = c("darkorange", "lightblue")) +
  theme_bw() +
  labs(x = " ", y = "Frequency", fill = "Treatment") +
  theme(legend.position = "top") 

 # ==== PCA =========
# Transpose the data to have samples as rows and molecular formulas as columns
data_transposed <- t(data)

# Perform PCA
pca_result <- prcomp(data_transposed, scale. = TRUE)

# Get explained variance for each PC
explained_variance <- summary(pca_result)$importance[2,] * 100  # % variance explained

# Convert PCA results to a data frame for ggplot
pca_data_for_plot <- as.data.frame(pca_result$x)
pca_data_for_plot$Treatment <- case_when(grepl("W", row.names(pca_data_for_plot)) ~ "Wet",grepl("D", row.names(pca_data_for_plot)) ~ "Dry")

# Create a PCA plot using ggplot2

ggplot(pca_data_for_plot, aes(x = PC1, y = PC2, color = Treatment)) +
  geom_point(size = 3) +  # Plot the PCA points
  stat_ellipse(aes(group = Treatment), level = 0.95, linetype = "solid") +  # Add ellipses for Wet/Dry groups
  labs(title = "PCA of Molecular Formula (Wet vs Dry)",
       x = paste0("PC1: ", round(explained_variance[1], 1), "%"),
       y = paste0("PC2: ", round(explained_variance[2], 1), "%")) + 
  scale_color_manual(values = c("darkorange", "lightblue")) +  # Color by Wet and Dry treatment
  theme_bw() +
  theme(legend.title = element_blank())  # Clean up the legend title

# ===== NMDS ===============

library(vegan)
dist_matrix <- vegdist(data, method = "jaccard")  # Adjust if needed (e.g., "jaccard", "euclidean")

nmds_result <- metaMDS(dist_matrix, k = 2, trymax = 100)

# Step 3: Extract NMDS coordinates
nmds_data <- as.data.frame(scores(nmds_result))  # Extract NMDS scores (coordinates)
nmds_data$Sample <- rownames(nmds_data)  # Add sample names
nmds_data$Treatment <- factors$Treatment  # Assuming you have a column in `factors` for Wet/Dry treatment

# Step 4: Plot NMDS with ggplot
ggplot(nmds_data, aes(x = MDS1, y = MDS2, color = Treatment)) +
  geom_point(size = 3) +  # Plot the NMDS points
  stat_ellipse(aes(group = Treatment), level = 0.95, linetype = "solid") +  # Add ellipses for Wet/Dry groups
  labs(title = "NMDS of Molecular Data (Wet vs Dry)",
       x = paste0("MDS1: ", round(nmds_result$stress * 100, 1), "% Stress"),
       y = paste0("MDS2")) + 
  scale_color_manual(values = c("Wet" = "darkorange", "Dry" = "lightblue")) +  # Color by Wet and Dry treatment
  theme_bw() +
  theme(legend.title = element_blank())  # Clean up the legend title

