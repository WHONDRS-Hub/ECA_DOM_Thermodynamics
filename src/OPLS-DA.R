# ==== Loading libraries =========
rm(list=ls(all=T))

#install.packages("BiocManager")
#BiocManager::install("ropls")

library(ropls)
library(tidyverse);library(ggplot2); library(caret)

# ==== Defining paths and working directories ======
input = 'C:/Users/gara009/OneDrive - PNNL/Documents - Core Richland and Sequim Lab-Field Team/Data Generation and Files/ECA/FTICR/03_ProcessedData/CoreMS/EC_Data_Processed_FTICR/'

github = 'C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/ECA_DOM_Thermodynamics/'
out_plots = paste0(github,'CoreMS/Plots/')
out_data = paste0(github,'CoreMS/Data/')

# ====== Read in data ======
# Processed ICR Data
data = read.csv(list.files(pattern = "*Processed_CoreMS_Clean_Data.csv", path = input, full.names = T), row.names = 1)
#mol = read.csv(list.files(pattern = "*Processed_CoreMS_Clean_Mol.csv", path = input, full.names = T), row.names = 1)

# Step 3: Normalize the data
# Apply log1p() transformation, which handles zeros by computing log(1 + intensity)
data_log <- log1p(data)  # log1p(x) is log(x + 1), handles zeros without creating artificial values

# Calculate the variance of each feature
feature_variances <- apply(data_log, 1, var)

# Set a threshold for variance (e.g., keep features with variance greater than a specified value)
threshold <- 1e-4  # Adjust this threshold based on your data

# Remove features with more than a certain percentage of zeros (e.g., >90%)
zero_threshold <- 0.9  # Adjust this as needed
non_zero_fraction <- rowMeans(data_log == 0)

# Combine conditions: keep features with sufficient variance and less than 90% zeros
data_filtered <- data_log[feature_variances > threshold & non_zero_fraction < zero_threshold, ]

# First step remove singletons, i.e. formula that only appear one time 
# Try first the filtering and the transformation last
# Think about variance partitioning more
# Try WGCNA co-correlates each of the intensities of the MF and has an algorithm for identifying clustering sort of like k-mean clustering but no ML. You might end with some clusters that have similar features

# Check final filtered data
dim(data_filtered)  # Should have some columns remaining

data_t <- t(data_filtered)

# Step 7: Extract sample names and treatment information
sample_names <- colnames(data)
treatment <- ifelse(grepl("SIR.D", sample_names), "Dry", "Wet")

# Step 9: Apply OPLS-DA for treatment separation
opls_model_treatment <- opls(data_t, treatment, predI = 1, orthoI = 1, permI = 100)

# Results
#OPLS-DA
# 555 samples x 3312 variables and 1 response
# standard scaling of predictors and response(s)
# R2X(cum) R2Y(cum) Q2(cum) RMSEE pre ort pR2Y  pQ2
# Total    0.185    0.202  -0.235 0.448   1   1 0.17 0.94

# ======= Other way ======

# Step 2: Log-transform the data to handle skewness and improve normality
data_log <- log(data + 1)  # Adding 1 to avoid log(0)

# Step 3: Filter features with low variance and high sparsity
feature_variances <- apply(data_log, 1, var)
zero_threshold <- 0.9  # Allow maximum 90% zeros
non_zero_fraction <- rowMeans(data_log == 0)

# Filter features: Keep those with variance > threshold and < 90% zeros
threshold <- 1e-4  # Adjust this threshold based on your data
data_filtered <- data_log[feature_variances > threshold & non_zero_fraction < zero_threshold, ]

# Check if data_filtered has columns left
if (ncol(data_filtered) == 0) {
  stop("No features left after filtering. Adjust your thresholds.")
}

# Step 4: Handle NAs by imputing with the mean (or median) of each feature
data_clean <- apply(data_filtered, 2, function(x) {
  ifelse(is.na(x), mean(x, na.rm = TRUE), x)
})

# Convert back to matrix if necessary
data_clean <- as.data.frame(data_clean)

# Step 5: Transpose the data for OPLS-DA
data_t <- t(data_clean)

# Ensure it’s numeric
data_t <- as.data.frame(apply(data_t, 2, as.numeric))

# Step 6: Define treatment variable based on sample names
treatment <- ifelse(grepl("SIR.D", sample_names), "Dry", "Wet")

# Step 7: Check for NAs in data_t
if (sum(is.na(data_t)) > 0) {
  stop("Data still contains NAs. Check your imputation step.")
}

# Step 8: Run OPLS-DA with reduced complexity (e.g., fewer components)
opls_model_treatment <- opls(data_t, treatment, predI = 1, orthoI = 0, permI = 100)  # No orthogonal components

# Step 9: Visualize the results
plot(opls_model_treatment, typeVc = "x-score", parAsColFcVn = treatment)

# Step 10: Evaluate the model's performance
summary(opls_model_treatment)

# Additional Step: Consider running PCA for exploratory analysis
pca_result <- prcomp(data_t, center = TRUE, scale. = TRUE)
plot(pca_result$x, col = treatment, pch = 19, main = "PCA of FTICR-MS Data")
# Result
# PLS-DA
# 555 samples x 3312 variables and 1 response
# standard scaling of predictors and response(s)
# R2X(cum) R2Y(cum) Q2(cum) RMSEE pre ort pR2Y pQ2
# Total    0.152   0.0248 -0.0393 0.495   1   0 0.99 0.5
# Warning message:
#   Single component model: only 'overview' and 'permutation' (in case of single response (O)PLS(-DA)) plots available


# ====== OPLS-DA for respiration =====
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

data_log <- log1p(data)  # log1p(x) is log(x + 1), handles zeros without creating artificial values

# Calculate the variance of each feature
feature_variances <- apply(data_log, 1, var)

# Set a threshold for variance (e.g., keep features with variance greater than a specified value)
threshold <- 1e-4  # Adjust this threshold based on your data

# Remove features with more than a certain percentage of zeros (e.g., >90%)
zero_threshold <- 0.9  # Adjust this as needed
non_zero_fraction <- rowMeans(data_log == 0)

# Combine conditions: keep features with sufficient variance and less than 90% zeros
data_filtered <- data_log[feature_variances > threshold & non_zero_fraction < zero_threshold, ]

# Check final filtered data
dim(data_filtered)  # Should have some columns remaining

data_t <- t(data_filtered)
data_t = as.data.frame(data_t)
# Merge ICR data with respiration data
data_t$Sample_Name <- rownames(data_t)
combined_data <- merge(data_t, resp, by = "Sample_Name")
Y_respiration <- combined_data$Respiration_Rate_mg_DO_per_kg_per_H
#Identify the columns corresponding to intensities
intensity_cols <- grep("^[0-9]+", colnames(combined_data))

# Extract predictor matrix X and response vectors Y
X <- combined_data[, intensity_cols]
opls_model_respiration <- opls(X, Y_respiration, predI = 1, orthoI = 1, permI = 100)

#Results
OPLS
554 samples x 3312 variables and 1 response
standard scaling of predictors and response(s)
R2X(cum) R2Y(cum) Q2(cum) RMSEE pre ort pR2Y  pQ2
Total    0.188    0.175  0.0264   582   1   1  0.4 0.01

# ===== Format data ====

# Transpose ICR data
data_t <- t(data) 
# Convert the transposed data to a data frame
data_t <- as.data.frame(data_t)

# Merge ICR data with respiration data
data_t$Sample_Name <- rownames(data_t)
combined_data <- merge(data_t, resp, by = "Sample_Name")

# ===== Prepare data for OPLS-DA =====
# Separate predictors from the response variable

# Identify the columns corresponding to intensities
intensity_cols <- grep("^[0-9]+", colnames(combined_data))

# Extract predictor matrix X and response vectors Y
X <- combined_data[, intensity_cols]
Y_treatment <- as.factor(combined_data$treatment)
Y_respiration <- combined_data$Respiration_Rate_mg_DO_per_kg_per_H

# ===== Running OPLS-DA for H1 =====
# Remove features with near zero variance
nzv <- nearZeroVar(X_normalized)
X_filtered <- X_normalized[, -nzv] # Maybe do both times

# Normalize the intensity data
X_normalized <- scale(X, center = TRUE, scale = TRUE) # z-score transformation
# This might not be ideal because z-score transformations might have been used for real concentrations
# we are assuming that each of the intensity values are link of counts but not really. It is just the concepts are related to each other. 
# Normalize as relative abundance/total sum transformation. sum constrained data the longer an ICR cycle runs the more ions you detect. One way to address this relative abundance and other is robust center log ratio transformation CLR use log scaling and finding rations between numbers

# Remove zero variance again and then run model one more time

# Check for class imbalance
table(Y_treatment)

# Use a subset of the data due to overlap between the treatments in a PCA
set.seed(123)  # For reproducibility
balanced_index <- createDataPartition(Y_treatment, p = 0.5, list = FALSE, times = 1)
X_balanced <- X_filtered[balanced_index, ]
Y_treatment_balanced <- Y_treatment[balanced_index]

# Retry OPLS-DA with balanced data
oplsda_model_treatment_balanced <- opls(X_balanced, Y_treatment_balanced, predI = 1, orthoI = NA, permI = 1000)

################################################

# Basic summary statistics
summary(combined_data)

# Check for class imbalance again
table(Y_treatment_balanced)

# Visualize distributions of some key features
# Pick a few features to visualize, e.g., the first few columns
feature_names <- colnames(X_filtered)[1:5]

# Plot histograms for the selected features by treatment
for (feature in feature_names) {
  ggplot(data = combined_data, aes_string(x = feature, fill = 'treatment')) + 
    geom_histogram(binwidth = 0.5, alpha = 0.6, position = 'identity') + 
    theme_minimal() + 
    labs(title = paste("Distribution of", feature, "by Treatment")) + 
    theme(legend.position = "top")
}
