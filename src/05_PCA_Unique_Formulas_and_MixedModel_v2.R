# ==== Loading libraries =========
rm(list=ls(all=T))

library(stringr); library(devtools);  library("plyr")
library("readr");  library(readxl);library(crayon); library(vegan)
# Load in necessary libraries first
library(reshape2)
library(ggpubr) # For to combine plots
library(dplyr) # For reorganization
library(stringr) # For string manipulation
library(tidyr)
library(vegan)
library(car)
library(lmtest) # For Breusch-Pagan test
library(mgcv)
library(lme4)
library(ggrepel)

# ==== Defining paths and working directories ======
github = 'C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/ECA_DOM_Thermodynamics/'
data_path = paste0(github,'Data/')
figure_path = paste0(github,'Figures/')

# ====== Read in data ======
data = read.csv(paste0(data_path,'Medians_of Median_molecular_properties_per_site_and_treatment_unique_formulas.csv'))

row.names(data) = paste0(data$site,'_',data$Treatment)

sample_data = read_csv(paste0(github,'EC_Data_Package/Sample_Data/EC_Sediment_Sample_Data_Summary.csv'),comment = '#', na = c('N/A', -9999)) %>%
  slice(-(1:11)) %>%  # Remove the first 11 rows
  mutate_at(vars(-Sample_Name, -Field_Name, -IGSN, -Material), as.numeric) %>%  # Convert specified columns to numeric
  mutate(
    Treatment = case_when(
      str_detect(Sample_Name, "-W") ~ "Wet",    # If Sample_Name contains "-W", assign "Wet"
      str_detect(Sample_Name, "-D") ~ "Dry",    # If Sample_Name contains "-D", assign "Dry"
      TRUE ~ NA_character_                       # Otherwise, assign NA
    )
  )

sample_data$site = gsub('-W|-D','',sample_data$Sample_Name)


# ==== Set up data =====
response_data = data %>%
  dplyr::filter(site != "EC_023") %>% 
  dplyr::filter(!(site %in% c("EC_012","EC_011","EC_053","EC_057","EC_052"))) %>%
#Remove EC_011, EC_012, EC_023, EC_052, EC_053, and EC_057 for having too much water added (11 and 12), no mg/kg calculation (23), and being duplicated NEON sites (52, 53, and 57)
  dplyr::select(-site,-Treatment)

explanatory_data = data %>%
  dplyr::select(site,Treatment)

explanatory_data = merge(explanatory_data,sample_data, by = c('site','Treatment'))

explanatory_data = explanatory_data %>%
  dplyr::select(Median_62948_Final_Gravimetric_Moisture_g_per_g, site, Treatment)

explanatory_data$Treatment <- as.factor(explanatory_data$Treatment)
explanatory_data$site <- as.factor(explanatory_data$site)

explanatory_data = explanatory_data %>%
  dplyr::filter(site != "EC_023") %>% 
  dplyr::filter(!(site %in% c("EC_012","EC_011","EC_053","EC_057","EC_052")))

# ==== Data transformation ====
cube_root <- function(x) sign(x) * (abs(x))^(1/3)

cube_predictors = explanatory_data %>%
  mutate(across(where(is.numeric) & !contains("Fine_Sand"), cube_root)) %>% # cube root transform everything except fine sand
  rename_with(~paste0("cube_", .), .cols = where(is.numeric) & !contains("Fine_Sand")) # rename everything except fine sand
  
cube_response = response_data %>% 
  mutate(across(where(is.numeric), cube_root)) %>% # cube root transform 
  rename_with(where(is.numeric), .fn = ~ paste0("cube_", .x))

final_cube_predictors_scaled <- cube_predictors %>%
  mutate_at(vars(-Treatment, -site), ~ scale(.))

# ====== PCA on Median DOM properties =====
# Removing co-linear properties
# Correlation matrix
library(corrplot)
cor_matrix <- cor(cube_response)
corrplot(cor_matrix, method = "color", type = "upper", 
         addCoef.col = "black", tl.col = "black", tl.srt = 45)


cube_response = cube_response %>%
  dplyr::select(-'cube_Median_Formulas',-'cube_Median_NOSC', -cube_Median_delGcoxPerCompmol)

# Run the PCA analysis
pca_result <- prcomp(cube_response, scale. = TRUE)
summary(pca_result)

# Extract important components
pc_scores <- as.data.frame(pca_result$x[,1:3])  # First 3 PCs

# Fit mixed models to principal components
pc1_model_moisture <- lmer(PC1 ~ cube_Median_62948_Final_Gravimetric_Moisture_g_per_g + (1|site), 
                  data = cbind(final_cube_predictors_scaled, pc_scores))

summary(pc1_model_moisture)

pc2_model_moisture <- lmer(PC2 ~ cube_Median_62948_Final_Gravimetric_Moisture_g_per_g + (1|site), 
                  data = cbind(final_cube_predictors_scaled, pc_scores))


summary(pc2_model_moisture)

# With Treatment

pc1_model_treatment <- lmer(PC1 ~ Treatment + (1|site), 
                  data = cbind(final_cube_predictors_scaled, pc_scores))

summary(pc1_model_treatment)

pc2_model_treatment <- lmer(PC2 ~ Treatment + (1|site), 
                  data = cbind(final_cube_predictors_scaled, pc_scores))


summary(pc2_model_treatment)

# Compare models using AIC
AIC(pc1_model_treatment, pc1_model_moisture)
AIC(pc2_model_treatment, pc2_model_moisture)

# PERMANOVA
perm_result <- adonis2(response_data ~ Median_62948_Final_Gravimetric_Moisture_g_per_g + site, 
                       data = explanatory_data,  # use original data
                       method = "euclidean",
                       by = "margin")


perm_result

# === PCA Plot ====
loadings_df <- data.frame(
  Variable = rownames(pca_result$rotation),
  PC1 = pca_result$rotation[,1],
  PC2 = pca_result$rotation[,2]
)

# Calculate variance explained percentages
var_exp <- summary(pca_result)$importance[2,] * 100
pc1_var <- round(var_exp[1], 2)
pc2_var <- round(var_exp[2], 2)


library(vegan)

perm_result <- adonis2(cube_response ~ Treatment + site,   data = explanatory_data, 
                       method = "euclidean",
                       by = "margin")

# Extract p-value for Treatment
treatment_p <- perm_result["Treatment", "Pr(>F)"]
treatment_R2 <- perm_result["Treatment", "R2"] * 100  # Convert to percentage

pca_plot <- ggplot() +
  geom_point(data = cbind(final_cube_predictors_scaled, pc_scores),
             aes(x = PC1, y = PC2, color = Treatment),
             size = 4, alpha = 0.8) +
  geom_segment(data = loadings_df,
               aes(x = 0, y = 0, xend = PC1 * 3, yend = PC2 * 3),
               arrow = arrow(length = unit(0.3, "cm")),
               color = "black") +
  geom_text_repel(data = loadings_df,
                  aes(x = PC1 * 3, y = PC2 * 3, label = Variable),
                  size = 4,
                  box.padding = 0.3,
                  point.padding = 0.5,
                  segment.color = "grey50") +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("PERMANOVA:\nSite R² = %.1f%%, p = 0.001\nTreatment R² = %.1f%%, p = %.3f", 
                           perm_result["site", "R2"]*100,
                           perm_result["Treatment", "R2"]*100,
                           perm_result["Treatment", "Pr(>F)"]),
           hjust = 0, vjust = 1,
           fontface = "italic") +
  scale_color_manual(values = c("Dry" = "darkorange", "Wet" = "lightblue")) +
  labs(title = " ",
       x = paste0("PC1 (", pc1_var, "% Variance Explained)"),
       y = paste0("PC2 (", pc2_var, "% Variance Explained)"),
       color = "Treatment") +
  theme_bw() +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.caption = element_text(hjust = 0, face = "italic"))


# Save plots
ggsave(paste0(figure_path,"Figure_2_PCA_plots.pdf"), 
       plot = pca_plot, width = 8, height = 8)
ggsave(paste0(figure_path,"Figure_2_PCA_plots.png"), 
       plot = pca_plot, width = 8, height = 8)
# ===== Save regressions ====
capture_all_analyses <- function(treatment_models, moisture_models, permanova, file_path) {
  sink(file_path)
  
  cat("==========================================\n")
  cat("Treatment Effect Models\n")
  cat("==========================================\n")
  cat("\nPC1 with Treatment\n")
  print(summary(treatment_models$pc1))
  cat("\nPC2 with Treatment\n")
  print(summary(treatment_models$pc2))
  
  cat("\n==========================================\n")
  cat("Continuous Moisture Models\n")
  cat("==========================================\n")
  cat("\nPC1 with Moisture\n")
  print(summary(moisture_models$pc1))
  cat("\nPC2 with Moisture\n")
  print(summary(moisture_models$pc2))
  
  cat("\n==========================================\n")
  cat("PERMANOVA Results\n")
  cat("==========================================\n")
  print(permanova)
  
  sink()
}

# Use function
capture_all_analyses(
  treatment_models = list(pc1 = pc1_model_treatment, pc2 = pc2_model_treatment),
  moisture_models = list(pc1 = pc1_model_moisture, pc2 = pc2_model_moisture),
  permanova = perm_result,
  "Data/PCA_mixed_models_summary_v2.txt")
