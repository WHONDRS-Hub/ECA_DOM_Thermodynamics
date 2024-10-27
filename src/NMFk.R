rm(list=ls(all=T))
# ==== Install and load necessary libraries =====
#install.packages("NMF")
#install.packages("cluster")

library(NMF)
library(cluster)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(readr)
library(caret)

# ==== Defining paths and working directories ======
input = 'C:/Users/gara009/OneDrive - PNNL/Documents - Core Richland and Sequim Lab-Field Team/Data Generation and Files/ECA/FTICR/03_ProcessedData/CoreMS/EC_Data_Processed_FTICR/'

github = 'C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/ECA_DOM_Thermodynamics/'
out_plots = paste0(github,'CoreMS/Plots/')
out_data = paste0(github,'CoreMS/Data/')

# ====== Read in data ======
# Processed ICR Data
data = read.csv(list.files(pattern = "*Processed_CoreMS_Clean_Data.csv", path = input, full.names = T), row.names = 1)
mol = read.csv(list.files(pattern = "*Processed_CoreMS_Clean_Mol.csv", path = input, full.names = T), row.names = 1)

#Matching names for data
colnames(data) = gsub('SIR.','INC-',colnames(data))

resp = read_csv(paste0(github,'EC_Data_Package/Sample_Data/EC_Sediment_SpC_pH_Temp_Respiration.csv'), comment = '#', na = c('N/A', -9999)) %>%
  slice(-(1:11))%>%
  select(Sample_Name,Respiration_Rate_mg_DO_per_kg_per_H) %>%
  mutate(Respiration_Rate_mg_DO_per_kg_per_H = as.numeric(Respiration_Rate_mg_DO_per_kg_per_H),
         Respiration_Rate_mg_DO_per_kg_per_H = ifelse(is.na(Respiration_Rate_mg_DO_per_kg_per_H) | Respiration_Rate_mg_DO_per_kg_per_H == 0, -min(abs(Respiration_Rate_mg_DO_per_kg_per_H[Respiration_Rate_mg_DO_per_kg_per_H != 0 & !is.na(Respiration_Rate_mg_DO_per_kg_per_H)])) / 2, 
                                                      Respiration_Rate_mg_DO_per_kg_per_H), # Replacing zeros and NA with half teh lowest respiration
         treatment = case_when(
           str_detect(Sample_Name, "D") ~ "Dry",
           str_detect(Sample_Name, "W") ~ "Wet",
           TRUE ~ NA_character_  # If neither condition is met
         ))
resp$Respiration_Rate_mg_DO_per_kg_per_H = abs(resp$Respiration_Rate_mg_DO_per_kg_per_H)
# ==== Normalize data =====
# Check the proportion of zeros
zero_proportion <- mean(fticr_data == 0)
print(paste("Proportion of zeros in the dataset:", zero_proportion))


fticr_data_t <- t(data)
fticr_data_t <- as.data.frame(fticr_data_t)

fticr_data_t$Sample_Name <- rownames(fticr_data_t)
combined_data <- merge(fticr_data_t, resp, by = "Sample_Name")

intensity_cols <- grep("^[0-9]+", colnames(combined_data))
X <- combined_data[, intensity_cols]
Y_treatment <- as.factor(combined_data$treatment)
Y_respiration <- combined_data$respiration

# Min-Max Scaling
min_max_scale <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}
X_min_max_scaled <- apply(X, 2, min_max_scale)

# Check for NA values
any(is.na(X_min_max_scaled)) # Should be FALSE after the previous round of cleaning

# Check for rows that are entirely NA or zero (i.e., need to drop these)
rows_with_issues <- apply(X_min_max_scaled, 1, function(row) {
  all(is.na(row)) || all(row == 0)
})
sum(rows_with_issues)  # Number of problematic rows

# Ensure there are no such rows by removing them
X_clean <- X_min_max_scaled[!rows_with_issues, ]

# Verify clean data
any(is.na(X_clean))       # Should be FALSE
any(apply(X_clean, 1, function(row) all(is.na(row))))   # Should be FALSE
any(apply(X_clean, 1, function(row) all(row == 0)))     # Should be FALSE



# ==== Apply NMF =====

# Set a random seed for reproducibility
set.seed(123)

# Determine the optimal number of components
estim.r <- nmf(X_clean, 2:6, nrun = 30, method = "brunet", seed = 1234)
plot(estim.r)

# Apply NMF with the selected rank (e.g., k = 3)
nmf_model <- nmf(X_filtered, 3, nrun = 30, method = "brunet", seed = 1234)
summary(nmf_model)

# ==== Extract basis and coefficient matrices ====

# Extract basis (W) and coefficient (H) matrices
W_matrix <- basis(nmf_model)
H_matrix <- coef(nmf_model)

# Transpose H_matrix if needed
H_matrix <- t(H_matrix)

# ===== Apply k-means clustering to coefficient matrix ====

# Apply k-means clustering
num_clusters <- 3
set.seed(123)
kmeans_result <- kmeans(H_matrix, centers = num_clusters, nstart = 25)

# Add cluster assignments to metadata
combined_data$cluster <- as.factor(kmeans_result$cluster)

# Visualize clustering results
pca_result <- prcomp(H_matrix, scale. = TRUE)
ggplot(data.frame(pca_result$x), aes(PC1, PC2, color = combined_data$cluster)) + 
  geom_point(alpha = 0.7) + 
  theme_minimal() + 
  labs(title = "NMF Clustering Results", x = "PC1", y = "PC2") +
  scale_color_discrete(name = "Cluster")

# Analyze clusters with treatment and respiration rates
table(combined_data$treatment, combined_data$cluster)

combined_data %>%
  group_by(cluster, treatment) %>%
  summarise(mean_respiration = mean(respiration))
