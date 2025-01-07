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


# ==== Defining paths and working directories ======
github = 'C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/ECA_DOM_Thermodynamics/'
data_path = paste0(github,'Data/')
figure_path = paste0(github,'Figures/')

# Create directories if they don't exist
if (!dir.exists("05_GAM_Results")) {
  dir.create("05_GAM_Results")
}
if (!dir.exists("05_GAM_Results/GAM_Results")) {
  dir.create("05_GAM_Results/GAM_Results")
}
if (!dir.exists("05_GAM_Results/Diagnostic_Plots")) {
  dir.create("05_GAM_Results/Diagnostic_Plots")
}
if (!dir.exists("05_GAM_Results/Comparison_Results")) {
  dir.create("05_GAM_Results/Comparison_Results")
}

# ====== Read in data ======
data = read.csv(paste0(data_path,'Medians_of_Weighted_averages_for_molecular_properties_per_site_and_treatment_unique_formulas.csv'))
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


grainsize_data = read_csv(paste0(github,'v4_CM_SSS_Data_Package/Sample_Data/v3_CM_SSS_Sediment_Grain_Size.csv'),comment = '#', na = c('N/A', -9999)) %>%
  slice(-(1:11))%>%
  mutate_at(vars(-Sample_Name,-Field_Name,-IGSN,-Material,-Methods_Deviation), as.numeric)%>%
  dplyr::select(-Field_Name,-IGSN,-Material,-Methods_Deviation)

grainsize_data$site = grainsize_data$Sample_Name
grainsize_data$site = gsub('CM','EC',grainsize_data$site)
grainsize_data$site = gsub('_GRN','',grainsize_data$site)

grainsize_data = grainsize_data %>%
  dplyr::select(site,Percent_Fine_Sand)
# ==== Set up data =====
response_data = data %>%
  dplyr::filter(site != "EC_023") %>%  # Remove rows where site is EC_023 because it has NA in some of the explanatory data
  dplyr::select(-site,-Treatment)

explanatory_data = data %>%
  dplyr::select(site,Treatment)

explanatory_data = merge(explanatory_data,sample_data, by = c('site','Treatment'))

explanatory_data = explanatory_data %>%
  dplyr::select(-Field_Name,-Sample_Name,-IGSN, -Material,-"Median_Respiration_Rate_mg_DO_per_L_per_H", -"Median_62948_Initial_Gravimetric_Moisture_g_per_g",- "Median_Fe_mg_per_L",-"Median_ATP_nanomoles_per_L",-"Median_Extractable_NPOC_mg_per_L",-"Median_Extractable_TN_mg_per_L",-"Median_Missing_Reps")

explanatory_data = merge(explanatory_data,grainsize_data, by = 'site')

explanatory_data$Treatment <- as.factor(explanatory_data$Treatment)
explanatory_data$site <- as.factor(explanatory_data$site)

# ===== Check for co-linearity =====
# To help with model overfitting for the number of samples we have, assess co-linearity and pick between 5-8 explanatory variables to put into the model. 
# Check for multicolinearity
vif_values <- vif(lm(Median_SpC_microsiemens_per_cm ~ Median_pH + Median_Temperature_degC + 
                       Median_Respiration_Rate_mg_DO_per_kg_per_H + Median_62948_Final_Gravimetric_Moisture_g_per_g + 
                       Median_Fe_mg_per_kg + Median_ATP_picomoles_per_g + Median_Extractable_NPOC_mg_per_kg + 
                       Median_Extractable_TN_mg_per_kg + Median_01395_C_percent_per_mg + 
                       Median_01397_N_percent_per_mg + Percent_Fine_Sand, data = explanatory_data))
print(vif_values)


# Need to keep variables with vif< 5 and also keep variables with relevance to the hypothesis

library(corrplot)
corr_matrix <- cor(select(explanatory_data, Median_SpC_microsiemens_per_cm, Median_pH, 
                          Median_Temperature_degC, Median_Respiration_Rate_mg_DO_per_kg_per_H, 
                          Median_62948_Final_Gravimetric_Moisture_g_per_g, Median_Fe_mg_per_kg, 
                          Median_ATP_picomoles_per_g, Median_Extractable_NPOC_mg_per_kg, 
                          Median_Extractable_TN_mg_per_kg, Median_01395_C_percent_per_mg, 
                          Median_01397_N_percent_per_mg, Percent_Fine_Sand), use = "complete.obs")
corrplot(corr_matrix, method = "color", addCoef.col = "black", number.cex = 0.7)

# Identify and remove high VIF variables
predictors_reduced <- explanatory_data %>%
  select(-site, -Treatment,
         -Median_Respiration_Rate_mg_DO_per_kg_per_H,
         -Median_Temperature_degC,
         -Median_pH,
         -Median_SpC_microsiemens_per_cm)

lm_reduced <- lm(Median_62948_Final_Gravimetric_Moisture_g_per_g  ~ ., data = predictors_reduced)
vif_reduced <- vif(lm_reduced)
print(vif_reduced)

final_predictors <- explanatory_data %>%
  select(site,Treatment,colnames(predictors_reduced))%>%
  drop_na() # site 23 has NA in some of geochem data

final_predictors_scaled <- final_predictors %>%
  mutate_at(vars(-Treatment, -site), ~ scale(.))


# ==== Run RDA =====
# Partial RDA controlling for site

rda_partial <- rda(response_data ~ Treatment + Median_62948_Final_Gravimetric_Moisture_g_per_g + Median_Fe_mg_per_kg + Median_ATP_picomoles_per_g + Median_Extractable_NPOC_mg_per_kg + Median_Extractable_TN_mg_per_kg + Median_01395_C_percent_per_mg + Median_01397_N_percent_per_mg + Percent_Fine_Sand + Condition(site), data = final_predictors_scaled)

summary(rda_partial)
anova(rda_partial, permutations = 999)
anova(rda_partial, by = "axis", permutations = 999)  # Axis significance
anova(rda_partial, by = "terms", permutations = 999)  # Variable significance

# Plot RDA
plot(rda_partial, scaling = 2, main = "Partial RDA Biplot")

# Check for influential samples
cooks.distance <- cooks.distance(rda_partial)
plot(cooks.distance, type = "h", main = "Cooks Distance for Partial RDA")
abline(h = 4/(length(cooks.distance)), col = "red")  # Common threshold

# ==== RDA Plots ===
site_scores <- as.data.frame(scores(rda_partial, display = "sites", scaling = 2))
sites_df <- as.data.frame(site_scores)
sites_df$Treatment <- final_predictors_scaled$Treatment  # Add treatment info

biplot_scores <- as.data.frame(scores(rda_partial, display = "bp", scaling = 2))
signif_vars <- c("Median_62948_Final_Gravimetric_Moisture_g_per_g",
                 "Median_Extractable_NPOC_mg_per_kg")  #

# Filter species_scores for significant variables
species_scores_signif <- biplot_scores[rownames(biplot_scores) %in% signif_vars, ]

# Create a data frame for significant variables
species_df_signif <- as.data.frame(species_scores_signif)
species_df_signif$Variable <- rownames(species_scores_signif)

species_df_signif <- species_df_signif %>%
  mutate(Variable = case_when(
    Variable == "Median_Extractable_NPOC_mg_per_kg" ~ "Extractable NPOC",
    Variable == "Median_62948_Final_Gravimetric_Moisture_g_per_g" ~ "Gravimetric Moisture",
    TRUE ~ Variable  # Keep other variables as they are
  ))

library(vegan)
library(ggplot2)
library(ggrepel)

# Plot Partial RDA using ggplot2
partial_rda_plot <- ggplot() +
  # Plot the sites (samples)
  geom_point(data = sites_df, 
             aes(x = RDA1, y = RDA2, color = Treatment), 
             size = 3, alpha = 0.8) +
  # Plot arrows for significant variables
  geom_segment(data = species_df_signif, 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "black") +
  # Label significant variables
  geom_text_repel(data = species_df_signif, 
                  aes(x = RDA1, y = RDA2, label = Variable),
                  size = 4, 
                  box.padding = 0.3, 
                  point.padding = 0.5,
                  segment.color = "grey50") +
  # Customize colors for treatments
  scale_color_manual(values = c("Dry" = "darkorange", "Wet" = "lightblue")) +
  # Add labels and theme adjustments
  labs(title = " ",
       x = "RDA1 (33.33% Variance Explained)",
       y = "RDA2 (6.95% Variance Explained)",
       color = "Treatment") +
  theme_bw() +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, face = "bold"))

# Save as PDF
ggsave(paste0(figure_path,"Figure3_Partial-RDA_plots.pdf"), plot = partial_rda_plot, width = 8, height = 8)

# Save as PNG
ggsave(paste0(figure_path,"Figure3_Partial-RDA_plots.png"), plot = partial_rda_plot, width = 8, height = 8)

# ======== RDA plot with arrows scaled =====
library(ggplot2)
library(ggrepel)

scaling_factor <- 5  # Change this as necessary to clearly visualize the arrows

partial_rda_plot <- ggplot() +
  # Plot the sites (samples)
  geom_point(data = sites_df, 
             aes(x = RDA1, y = RDA2, color = Treatment), 
             size = 3, alpha = 0.8) +
  # Plot arrows for significant variables with scaling
  geom_segment(data = species_df_signif, 
               aes(x = 0, y = 0, xend = scaling_factor * RDA1, yend = scaling_factor * RDA2),
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "black") +
  # Label significant variables with scaling
  geom_text_repel(data = species_df_signif, 
                  aes(x = scaling_factor * RDA1, y = scaling_factor * RDA2, label = Variable),
                  size = 4, 
                  box.padding = 0.3, 
                  point.padding = 0.5,
                  segment.color = "grey50") +
  # Customize colors for treatments
  scale_color_manual(values = c("Dry" = "darkorange", "Wet" = "lightblue")) +
  # Add labels and theme adjustments
  labs(title = "Partial RDA Plot",
       x = "RDA1 (33.33% Variance Explained)",
       y = "RDA2 (6.95% Variance Explained)",
       color = "Treatment") +
  theme_bw() +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, face = "bold"))

# Print the plot
print(partial_rda_plot)

# ====== Clean up and load data again ====
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
# ==== Defining paths and working directories 
github = 'C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/ECA_DOM_Thermodynamics/'
data_path = paste0(github,'Data/')
figure_path = paste0(github,'Figures/')
# ====== Read in data
data = read.csv(paste0(data_path,'Medians_of_Weighted_averages_for_molecular_properties_per_site_and_treatment_unique_formulas.csv'))
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

grainsize_data = read_csv(paste0(github,'v4_CM_SSS_Data_Package/Sample_Data/v3_CM_SSS_Sediment_Grain_Size.csv'),comment = '#', na = c('N/A', -9999)) %>%
  slice(-(1:11))%>%
  mutate_at(vars(-Sample_Name,-Field_Name,-IGSN,-Material,-Methods_Deviation), as.numeric)%>%
  dplyr::select(-Field_Name,-IGSN,-Material,-Methods_Deviation)

grainsize_data$site = grainsize_data$Sample_Name
grainsize_data$site = gsub('CM','EC',grainsize_data$site)
grainsize_data$site = gsub('_GRN','',grainsize_data$site)

grainsize_data = grainsize_data %>%
  dplyr::select(site,Percent_Fine_Sand)
# ==== Set up data 
data = data %>%
  dplyr::filter(site != "EC_023") # Remove rows where site is EC_023 because it has NA in some of the explanatory data

data = merge(data,sample_data, by = c('site','Treatment'))

data = data %>%
  dplyr::select(-Field_Name,-Sample_Name,-IGSN, -Material, -"Median_62948_Initial_Gravimetric_Moisture_g_per_g",- "Median_Fe_mg_per_L",-"Median_ATP_nanomoles_per_L",-"Median_Extractable_NPOC_mg_per_L",-"Median_Extractable_TN_mg_per_L",-"Median_Missing_Reps")

data = merge(data,grainsize_data, by = 'site')

data$Treatment <- as.factor(data$Treatment)
data$site <- as.factor(data$site)
combined_data = data
# ===== Checking data skewness in untransformed data before GAM ====
library(dplyr)
library(ggplot2)
library(e1071)      # For skewness calculation
library(reshape2)   # For data reshaping (if needed)


# Identify numeric columns
numeric_vars <- data %>%
  select(where(is.numeric)) %>%
  names()

# View the numeric variables
print(numeric_vars)
# Calculate skewness for each numeric variable
skewness_values <- sapply(data[numeric_vars], skewness, na.rm = TRUE)

# Create a data frame for easy viewing
skewness_df <- data.frame(
  Variable = names(skewness_values),
  Skewness = skewness_values
)

# Print skewness values
print(skewness_df)

# Highlight variables with high skewness
high_skew <- skewness_df %>%
  filter(Skewness > 1 | Skewness < -1)

print("Variables with High Skewness:")
print(high_skew)

library(gridExtra)
high_skew_vars = high_skew$Variable
# Function to plot histograms for variables with high skewness
plot_high_skew <- function(data, high_skew_vars) {
  plots <- list()
  for (var in high_skew_vars$Variable) {
    p <- ggplot(data, aes_string(x = var)) +
      geom_histogram(aes(y = ..density..), bins = 30, fill = "salmon", color = "black", alpha = 0.7) +
      geom_density(color = "darkblue", size = 1) +
      labs(title = paste("Histogram and Density of", var),
           x = var,
           y = "Density") +
      theme_minimal()
    plots[[var]] <- p
  }
  
  # Arrange plots in a grid
  do.call(grid.arrange, c(plots, ncol = 2))
}

# Plot histograms for highly skewed variables
plot_high_skew(data, high_skew)

# ====== Cube root transform data and check skewness ====
cube_root <- function(x) sign(x) * (abs(x))^(1/3)

cube_data = data %>% 
  mutate(across(where(is.numeric), cube_root)) %>% # cube root transform data
  rename_with(where(is.numeric), .fn = ~ paste0("cube_", .x))

# Identify numeric columns
numeric_vars <- cube_data %>%
  select(where(is.numeric)) %>%
  names()
# Calculate skewness for each numeric variable
skewness_values <- sapply(cube_data[numeric_vars], skewness, na.rm = TRUE)

# Create a data frame for easy viewing
skewness_df <- data.frame(
  Variable = names(skewness_values),
  Skewness = skewness_values
)

# Print skewness values
print(skewness_df)

# Note that there are mild improvements on the skewness of the data for some of the DOM properties but NOSC becomes significantly worse

# Highlight variables with high skewness
high_skew <- skewness_df %>%
  filter(Skewness > 1 | Skewness < -1)

print("Variables with High Skewness:")
print(high_skew)

library(gridExtra)
high_skew_vars = high_skew$Variable
# Function to plot histograms for variables with high skewness
plot_high_skew <- function(data, high_skew_vars) {
  plots <- list()
  for (var in high_skew_vars$Variable) {
    p <- ggplot(data, aes_string(x = var)) +
      geom_histogram(aes(y = ..density..), bins = 30, fill = "salmon", color = "black", alpha = 0.7) +
      geom_density(color = "darkblue", size = 1) +
      labs(title = paste("Histogram and Density of", var),
           x = var,
           y = "Density") +
      theme_minimal()
    plots[[var]] <- p
  }
  
  # Arrange plots in a grid
  do.call(grid.arrange, c(plots, ncol = 2))
}

# Plot histograms for highly skewed variables
plot_high_skew(cube_data, high_skew)

#  ===== Check model for transformed or untransformed data ====
# Fit models using transformed and untransformed predictors
# ==== AI_mod ====
model_transformed <- gam(cube_Median_Weighted_Avg_AI_mod ~ s(cube_Median_62948_Final_Gravimetric_Moisture_g_per_g) + 
                           s(cube_Median_Extractable_NPOC_mg_per_kg) + 
                           ti(cube_Median_62948_Final_Gravimetric_Moisture_g_per_g, cube_Median_Extractable_NPOC_mg_per_kg) + 
                           s(site, bs = "re"), data = cube_data)

model_original <- gam(Median_Weighted_Avg_AI_mod ~ s(Median_62948_Final_Gravimetric_Moisture_g_per_g) + 
                        s(Median_Extractable_NPOC_mg_per_kg) + 
                        ti(Median_62948_Final_Gravimetric_Moisture_g_per_g, Median_Extractable_NPOC_mg_per_kg) + 
                        s(site, bs = "re"), data = combined_data)

# Compare AIC
model_comparison_file <- "05_GAM_Results/Comparison_Results/AIC_comparison_AI_mod.csv"
model_comparison <- AIC(model_transformed, model_original)
write.csv(model_comparison, model_comparison_file)

summary_transformed_file <- "05_GAM_Results/GAM_Results/summary_transformed_AI_mod.txt"
summary_original_file <- "05_GAM_Results/GAM_Results/summary_original_AI_mod.txt"
diag_transformed_file <- "05_GAM_Results/GAM_Results/diagnostics_transformed_AI_mod.txt"
diag_original_file <- "05_GAM_Results/GAM_Results/diagnostics_original_AI_mod.txt"

sink(summary_transformed_file)
print(summary(model_transformed))
sink()

sink(summary_original_file)
print(summary(model_original))
sink()

sink(diag_transformed_file)
gam.check(model_transformed)
sink()

sink(diag_original_file)
gam.check(model_original)
sink()

plot_file_original <- "05_GAM_Results/Diagnostic_Plots/diagnostic_plot_original_AI_mod.png"
plot_file_transformed <- "05_GAM_Results/Diagnostic_Plots/diagnostic_plot_transformed_AI_mod.png"

png(plot_file_original, width = 2000, height = 1500, res=300)
plot(model_original, pages = 1)
dev.off()

png(plot_file_transformed, width = 2000, height = 1500, res=300)
plot(model_transformed, pages = 1)
dev.off()


# Function to create residual plots
plot_residuals <- function(model, model_name) {
  residuals <- resid(model)
  fitted_values <- fitted(model)
  
  residual_plot <- ggplot(data = NULL, aes(x = fitted_values, y = residuals)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(title = paste("Residuals vs Fitted -", model_name),
         x = "Fitted Values",
         y = "Residuals") +
    theme_minimal()
  
  plot_file <- paste0("05_GAM_Results/Diagnostic_Plots/residual_plot_", model_name, ".png")
  ggsave(plot_file, plot = residual_plot, width = 10, height = 8, units = "in")
}

# Apply the plotting function to each model
plot_residuals(model_original, "original_AI_mod")
plot_residuals(model_transformed, "transformed_AI_mod")


# Function to perform Breusch-Pagan test and export results
homoscedasticity_tests <- function(model, model_name) {
  bp_test <- bptest(model)
  results_file <- paste0("05_GAM_Results/GAM_Results/BP_test_", model_name, ".csv")
  results <- data.frame(
    Model = model_name,
    BP_p_value = bp_test$p.value
  )
  write.csv(results, results_file, row.names = FALSE)
}

# Apply the homoscedasticity test function to each model
homoscedasticity_tests(model_original, "original_AI_mod")
homoscedasticity_tests(model_transformed, "transformed_AI_mod")


# Plotting smooth functions and saving plots
plot_smooth_terms <- function(model, model_name) {
  plot_file <- paste0("05_GAM_Results/Diagnostic_Plots/smooth_plot_", model_name, ".png")
  png(plot_file, width = 2000, height = 1500, res = 300)
  plot(model, pages = 1, all.terms = TRUE)
  dev.off()
}

plot_smooth_terms(model_original, "original_AI_mod")
plot_smooth_terms(model_transformed, "transformed_AI_mod")

# ==== NOSC ====
model_transformed <- gam(cube_Median_Weighted_Avg_NOSC ~ s(cube_Median_62948_Final_Gravimetric_Moisture_g_per_g) + 
                           s(cube_Median_Extractable_NPOC_mg_per_kg) + 
                           ti(cube_Median_62948_Final_Gravimetric_Moisture_g_per_g, cube_Median_Extractable_NPOC_mg_per_kg) + 
                           s(site, bs = "re"), data = cube_data)

model_original <- gam(Median_Weighted_Avg_NOSC ~ s(Median_62948_Final_Gravimetric_Moisture_g_per_g) + 
                        s(Median_Extractable_NPOC_mg_per_kg) + 
                        ti(Median_62948_Final_Gravimetric_Moisture_g_per_g, Median_Extractable_NPOC_mg_per_kg) + 
                        s(site, bs = "re"), data = combined_data)

# Compare AIC
model_comparison_file <- "05_GAM_Results/Comparison_Results/AIC_comparison_NOSC.csv"
model_comparison <- AIC(model_transformed, model_original)
write.csv(model_comparison, model_comparison_file)

summary_transformed_file <- "05_GAM_Results/GAM_Results/summary_transformed_NOSC.txt"
summary_original_file <- "05_GAM_Results/GAM_Results/summary_original_NOSC.txt"
diag_transformed_file <- "05_GAM_Results/GAM_Results/diagnostics_transformed_NOSC.txt"
diag_original_file <- "05_GAM_Results/GAM_Results/diagnostics_original_NOSC.txt"

sink(summary_transformed_file)
print(summary(model_transformed))
sink()

sink(summary_original_file)
print(summary(model_original))
sink()

sink(diag_transformed_file)
gam.check(model_transformed)
sink()

sink(diag_original_file)
gam.check(model_original)
sink()

plot_file_original <- "05_GAM_Results/Diagnostic_Plots/diagnostic_plot_original_NOSC.png"
plot_file_transformed <- "05_GAM_Results/Diagnostic_Plots/diagnostic_plot_transformed_NOSC.png"

png(plot_file_original, width = 2000, height = 1500, res=300)
plot(model_original, pages = 1)
dev.off()

png(plot_file_transformed, width = 2000, height = 1500, res=300)
plot(model_transformed, pages = 1)
dev.off()

# Apply the plotting function to each model
plot_residuals(model_original, "original_NOSC")
plot_residuals(model_transformed, "transformed_NOSC")

# Apply the homoscedasticity test function to each model
homoscedasticity_tests(model_original, "original_NOSC")
homoscedasticity_tests(model_transformed, "transformed_NOSC")

# Plotting smooth functions and saving plots
plot_smooth_terms <- function(model, model_name) {
  plot_file <- paste0("05_GAM_Results/Diagnostic_Plots/smooth_plot_", model_name, ".png")
  png(plot_file, width = 2000, height = 1500, res = 300)
  plot(model, pages = 1, all.terms = TRUE)
  dev.off()
}

plot_smooth_terms(model_original, "original_NOSC")
plot_smooth_terms(model_transformed, "transformed_NOSC")


# ==== DBE ====
model_transformed <- gam(cube_Median_Weighted_Avg_DBE ~ s(cube_Median_62948_Final_Gravimetric_Moisture_g_per_g) + 
                           s(cube_Median_Extractable_NPOC_mg_per_kg) + 
                           ti(cube_Median_62948_Final_Gravimetric_Moisture_g_per_g, cube_Median_Extractable_NPOC_mg_per_kg) + 
                           s(site, bs = "re"), data = cube_data)

model_original <- gam(Median_Weighted_Avg_DBE ~ s(Median_62948_Final_Gravimetric_Moisture_g_per_g) + 
                        s(Median_Extractable_NPOC_mg_per_kg) + 
                        ti(Median_62948_Final_Gravimetric_Moisture_g_per_g, Median_Extractable_NPOC_mg_per_kg) + 
                        s(site, bs = "re"), data = combined_data)

# Compare AIC
model_comparison_file <- "05_GAM_Results/Comparison_Results/AIC_comparison_DBE.csv"
model_comparison <- AIC(model_transformed, model_original)
write.csv(model_comparison, model_comparison_file)

summary_transformed_file <- "05_GAM_Results/GAM_Results/summary_transformed_DBE.txt"
summary_original_file <- "05_GAM_Results/GAM_Results/summary_original_DBE.txt"
diag_transformed_file <- "05_GAM_Results/GAM_Results/diagnostics_transformed_DBE.txt"
diag_original_file <- "05_GAM_Results/GAM_Results/diagnostics_original_DBE.txt"

sink(summary_transformed_file)
print(summary(model_transformed))
sink()

sink(summary_original_file)
print(summary(model_original))
sink()

sink(diag_transformed_file)
gam.check(model_transformed)
sink()

sink(diag_original_file)
gam.check(model_original)
sink()

plot_file_original <- "05_GAM_Results/Diagnostic_Plots/diagnostic_plot_original_DBE.png"
plot_file_transformed <- "05_GAM_Results/Diagnostic_Plots/diagnostic_plot_transformed_DBE.png"

png(plot_file_original, width = 2000, height = 1500, res=300)
plot(model_original, pages = 1)
dev.off()

png(plot_file_transformed, width = 2000, height = 1500, res=300)
plot(model_transformed, pages = 1)
dev.off()

# Apply the plotting function to each model
plot_residuals(model_original, "original_DBE")
plot_residuals(model_transformed, "transformed_DBE")

# Apply the homoscedasticity test function to each model
homoscedasticity_tests(model_original, "original_DBE")
homoscedasticity_tests(model_transformed, "transformed_DBE")

# Plotting smooth functions and saving plots
plot_smooth_terms <- function(model, model_name) {
  plot_file <- paste0("05_GAM_Results/Diagnostic_Plots/smooth_plot_", model_name, ".png")
  png(plot_file, width = 2000, height = 1500, res = 300)
  plot(model, pages = 1, all.terms = TRUE)
  dev.off()
}

plot_smooth_terms(model_original, "original_DBE")
plot_smooth_terms(model_transformed, "transformed_DBE")

# ==== delGcoxPerCmol =====
model_transformed <- gam(cube_Median_Weighted_Avg_delGcoxPerCmol ~ s(cube_Median_62948_Final_Gravimetric_Moisture_g_per_g) + 
                           s(cube_Median_Extractable_NPOC_mg_per_kg) + 
                           ti(cube_Median_62948_Final_Gravimetric_Moisture_g_per_g, cube_Median_Extractable_NPOC_mg_per_kg) + 
                           s(site, bs = "re"), data = cube_data)

model_original <- gam(Median_Weighted_Avg_delGcoxPerCmol ~ s(Median_62948_Final_Gravimetric_Moisture_g_per_g) + 
                        s(Median_Extractable_NPOC_mg_per_kg) + 
                        ti(Median_62948_Final_Gravimetric_Moisture_g_per_g, Median_Extractable_NPOC_mg_per_kg) + 
                        s(site, bs = "re"), data = combined_data)

# Compare AIC
model_comparison_file <- "05_GAM_Results/Comparison_Results/AIC_comparison_delGcoxPerCmol.csv"
model_comparison <- AIC(model_transformed, model_original)
write.csv(model_comparison, model_comparison_file)

summary_transformed_file <- "05_GAM_Results/GAM_Results/summary_transformed_delGcoxPerCmol.txt"
summary_original_file <- "05_GAM_Results/GAM_Results/summary_original_delGcoxPerCmol.txt"
diag_transformed_file <- "05_GAM_Results/GAM_Results/diagnostics_transformed_delGcoxPerCmol.txt"
diag_original_file <- "05_GAM_Results/GAM_Results/diagnostics_original_delGcoxPerCmol.txt"

sink(summary_transformed_file)
print(summary(model_transformed))
sink()

sink(summary_original_file)
print(summary(model_original))
sink()

sink(diag_transformed_file)
gam.check(model_transformed)
sink()

sink(diag_original_file)
gam.check(model_original)
sink()

plot_file_original <- "05_GAM_Results/Diagnostic_Plots/diagnostic_plot_original_delGcoxPerCmol.png"
plot_file_transformed <- "05_GAM_Results/Diagnostic_Plots/diagnostic_plot_transformed_delGcoxPerCmol.png"

png(plot_file_original, width = 2000, height = 1500, res=300)
plot(model_original, pages = 1)
dev.off()

png(plot_file_transformed, width = 2000, height = 1500, res=300)
plot(model_transformed, pages = 1)
dev.off()

# Apply the plotting function to each model
plot_residuals(model_original, "original_delGcoxPerCmol")
plot_residuals(model_transformed, "transformed_delGcoxPerCmol")

# Apply the homoscedasticity test function to each model
homoscedasticity_tests(model_original, "original_delGcoxPerCmol")
homoscedasticity_tests(model_transformed, "transformed_delGcoxPerCmol")

# Plotting smooth functions and saving plots
plot_smooth_terms <- function(model, model_name) {
  plot_file <- paste0("05_GAM_Results/Diagnostic_Plots/smooth_plot_", model_name, ".png")
  png(plot_file, width = 2000, height = 1500, res = 300)
  plot(model, pages = 1, all.terms = TRUE)
  dev.off()
}

plot_smooth_terms(model_original, "original_delGcoxPerCmol")
plot_smooth_terms(model_transformed, "transformed_delGcoxPerCmol")

# ==== delGcoxPerCompmol ====

model_transformed <- gam(cube_Median_Weighted_Avg_delGcoxPerCompmol ~ s(cube_Median_62948_Final_Gravimetric_Moisture_g_per_g) + 
                           s(cube_Median_Extractable_NPOC_mg_per_kg) + 
                           ti(cube_Median_62948_Final_Gravimetric_Moisture_g_per_g, cube_Median_Extractable_NPOC_mg_per_kg) + 
                           s(site, bs = "re"), data = cube_data)

model_original <- gam(Median_Weighted_Avg_delGcoxPerCompmol ~ s(Median_62948_Final_Gravimetric_Moisture_g_per_g) + 
                        s(Median_Extractable_NPOC_mg_per_kg) + 
                        ti(Median_62948_Final_Gravimetric_Moisture_g_per_g, Median_Extractable_NPOC_mg_per_kg) + 
                        s(site, bs = "re"), data = combined_data)

# Compare AIC
model_comparison_file <- "05_GAM_Results/Comparison_Results/AIC_comparison_delGcoxPerCompmol.csv"
model_comparison <- AIC(model_transformed, model_original)
write.csv(model_comparison, model_comparison_file)

summary_transformed_file <- "05_GAM_Results/GAM_Results/summary_transformed_delGcoxPerCompmol.txt"
summary_original_file <- "05_GAM_Results/GAM_Results/summary_original_delGcoxPerCompmol.txt"
diag_transformed_file <- "05_GAM_Results/GAM_Results/diagnostics_transformed_delGcoxPerCompmol.txt"
diag_original_file <- "05_GAM_Results/GAM_Results/diagnostics_original_delGcoxPerCompmol.txt"

sink(summary_transformed_file)
print(summary(model_transformed))
sink()

sink(summary_original_file)
print(summary(model_original))
sink()

sink(diag_transformed_file)
gam.check(model_transformed)
sink()

sink(diag_original_file)
gam.check(model_original)
sink()

plot_file_original <- "05_GAM_Results/Diagnostic_Plots/diagnostic_plot_original_delGcoxPerCompmol.png"
plot_file_transformed <- "05_GAM_Results/Diagnostic_Plots/diagnostic_plot_transformed_delGcoxPerCompmol.png"

png(plot_file_original, width = 2000, height = 1500, res=300)
plot(model_original, pages = 1)
dev.off()

png(plot_file_transformed, width = 2000, height = 1500, res=300)
plot(model_transformed, pages = 1)
dev.off()

# Apply the plotting function to each model
plot_residuals(model_original, "original_delGcoxPerCompmol")
plot_residuals(model_transformed, "transformed_delGcoxPerCompmol")

# Apply the homoscedasticity test function to each model
homoscedasticity_tests(model_original, "original_delGcoxPerCompmol")
homoscedasticity_tests(model_transformed, "transformed_delGcoxPerCompmol")

# Plotting smooth functions and saving plots
plot_smooth_terms <- function(model, model_name) {
  plot_file <- paste0("05_GAM_Results/Diagnostic_Plots/smooth_plot_", model_name, ".png")
  png(plot_file, width = 2000, height = 1500, res = 300)
  plot(model, pages = 1, all.terms = TRUE)
  dev.off()
}

plot_smooth_terms(model_original, "original_delGcoxPerCompmol")
plot_smooth_terms(model_transformed, "transformed_delGcoxPerCompmol")

# ==== Lambda ====
model_transformed <- gam(cube_Median_Weighted_Avg_Lambda ~ s(cube_Median_62948_Final_Gravimetric_Moisture_g_per_g) + 
                           s(cube_Median_Extractable_NPOC_mg_per_kg) + 
                           ti(cube_Median_62948_Final_Gravimetric_Moisture_g_per_g, cube_Median_Extractable_NPOC_mg_per_kg) + 
                           s(site, bs = "re"), data = cube_data)

model_original <- gam(Median_Weighted_Avg_Lambda ~ s(Median_62948_Final_Gravimetric_Moisture_g_per_g) + 
                        s(Median_Extractable_NPOC_mg_per_kg) + 
                        ti(Median_62948_Final_Gravimetric_Moisture_g_per_g, Median_Extractable_NPOC_mg_per_kg) + 
                        s(site, bs = "re"), data = combined_data)

# Compare AIC
model_comparison_file <- "05_GAM_Results/Comparison_Results/AIC_comparison_Lambda.csv"
model_comparison <- AIC(model_transformed, model_original)
write.csv(model_comparison, model_comparison_file)

summary_transformed_file <- "05_GAM_Results/GAM_Results/summary_transformed_Lambda.txt"
summary_original_file <- "05_GAM_Results/GAM_Results/summary_original_Lambda.txt"
diag_transformed_file <- "05_GAM_Results/GAM_Results/diagnostics_transformed_Lambda.txt"
diag_original_file <- "05_GAM_Results/GAM_Results/diagnostics_original_Lambda.txt"

sink(summary_transformed_file)
print(summary(model_transformed))
sink()

sink(summary_original_file)
print(summary(model_original))
sink()

sink(diag_transformed_file)
gam.check(model_transformed)
sink()

sink(diag_original_file)
gam.check(model_original)
sink()

plot_file_original <- "05_GAM_Results/Diagnostic_Plots/diagnostic_plot_original_Lambda.png"
plot_file_transformed <- "05_GAM_Results/Diagnostic_Plots/diagnostic_plot_transformed_Lambda.png"

png(plot_file_original, width = 2000, height = 1500, res=300)
plot(model_original, pages = 1)
dev.off()

png(plot_file_transformed, width = 2000, height = 1500, res=300)
plot(model_transformed, pages = 1)
dev.off()

# Apply the plotting function to each model
plot_residuals(model_original, "original_Lambda")
plot_residuals(model_transformed, "transformed_Lambda")

# Apply the homoscedasticity test function to each model
homoscedasticity_tests(model_original, "original_Lambda")
homoscedasticity_tests(model_transformed, "transformed_Lambda")

# Plotting smooth functions and saving plots
plot_smooth_terms <- function(model, model_name) {
  plot_file <- paste0("05_GAM_Results/Diagnostic_Plots/smooth_plot_", model_name, ".png")
  png(plot_file, width = 2000, height = 1500, res = 300)
  plot(model, pages = 1, all.terms = TRUE)
  dev.off()
}

plot_smooth_terms(model_original, "original_Lambda")
plot_smooth_terms(model_transformed, "transformed_Lambda")


# ==== Testing GAM models default vs REML ====
# GAM models with default GCV method
ggam_models <- list(
  ggam_ai_mod = gam(Median_Weighted_Avg_AI_mod ~ s(Median_62948_Final_Gravimetric_Moisture_g_per_g) + 
                      s(Median_Extractable_NPOC_mg_per_kg) + 
                      ti(Median_62948_Final_Gravimetric_Moisture_g_per_g, Median_Extractable_NPOC_mg_per_kg) + 
                      s(site, bs = "re"), data = combined_data),
  ggam_nosc = gam(Median_Weighted_Avg_NOSC ~ s(Median_62948_Final_Gravimetric_Moisture_g_per_g) + 
                    s(Median_Extractable_NPOC_mg_per_kg) + 
                    ti(Median_62948_Final_Gravimetric_Moisture_g_per_g, Median_Extractable_NPOC_mg_per_kg) + 
                    s(site, bs = "re"), data = combined_data),
  ggam_dbe = gam(Median_Weighted_Avg_DBE ~ s(Median_62948_Final_Gravimetric_Moisture_g_per_g) + 
                   s(Median_Extractable_NPOC_mg_per_kg) + 
                   ti(Median_62948_Final_Gravimetric_Moisture_g_per_g, Median_Extractable_NPOC_mg_per_kg) + 
                   s(site, bs = "re"), data = combined_data),
  ggam_delGcoxPerCmol = gam(Median_Weighted_Avg_delGcoxPerCmol ~ s(Median_62948_Final_Gravimetric_Moisture_g_per_g) + 
                              s(Median_Extractable_NPOC_mg_per_kg) + 
                              ti(Median_62948_Final_Gravimetric_Moisture_g_per_g, Median_Extractable_NPOC_mg_per_kg) + 
                              s(site, bs = "re"), data = combined_data),
  ggam_delGcoxPerCompmol = gam(Median_Weighted_Avg_delGcoxPerCompmol ~ s(Median_62948_Final_Gravimetric_Moisture_g_per_g) + 
                                 s(Median_Extractable_NPOC_mg_per_kg) + 
                                 ti(Median_62948_Final_Gravimetric_Moisture_g_per_g, Median_Extractable_NPOC_mg_per_kg) + 
                                 s(site, bs = "re"), data = combined_data),
  ggam_lambda = gam(Median_Weighted_Avg_Lambda ~ s(Median_62948_Final_Gravimetric_Moisture_g_per_g) + 
                      s(Median_Extractable_NPOC_mg_per_kg) + 
                      ti(Median_62948_Final_Gravimetric_Moisture_g_per_g, Median_Extractable_NPOC_mg_per_kg) + 
                      s(site, bs = "re"), data = combined_data)
)

# GAM models with REML method
gam_models <- list(
  gam_ai_mod = gam(Median_Weighted_Avg_AI_mod ~ s(Median_62948_Final_Gravimetric_Moisture_g_per_g) + 
                     s(Median_Extractable_NPOC_mg_per_kg) + 
                     ti(Median_62948_Final_Gravimetric_Moisture_g_per_g, Median_Extractable_NPOC_mg_per_kg) + 
                     s(site, bs = "re"), data = combined_data, method = "REML"),
  gam_nosc = gam(Median_Weighted_Avg_NOSC ~ s(Median_62948_Final_Gravimetric_Moisture_g_per_g) + 
                   s(Median_Extractable_NPOC_mg_per_kg) + 
                   ti(Median_62948_Final_Gravimetric_Moisture_g_per_g, Median_Extractable_NPOC_mg_per_kg) + 
                   s(site, bs = "re"), data = combined_data, method = "REML"),
  gam_dbe = gam(Median_Weighted_Avg_DBE ~ s(Median_62948_Final_Gravimetric_Moisture_g_per_g) + 
                  s(Median_Extractable_NPOC_mg_per_kg) + 
                  ti(Median_62948_Final_Gravimetric_Moisture_g_per_g, Median_Extractable_NPOC_mg_per_kg) + 
                  s(site, bs = "re"), data = combined_data, method = "REML"),
  gam_delGcoxPerCmol = gam(Median_Weighted_Avg_delGcoxPerCmol ~ s(Median_62948_Final_Gravimetric_Moisture_g_per_g) + 
                             s(Median_Extractable_NPOC_mg_per_kg) + 
                             ti(Median_62948_Final_Gravimetric_Moisture_g_per_g, Median_Extractable_NPOC_mg_per_kg) + 
                             s(site, bs = "re"), data = combined_data, method = "REML"),
  gam_delGcoxPerCompmol = gam(Median_Weighted_Avg_delGcoxPerCompmol ~ s(Median_62948_Final_Gravimetric_Moisture_g_per_g) + 
                                s(Median_Extractable_NPOC_mg_per_kg) + 
                                ti(Median_62948_Final_Gravimetric_Moisture_g_per_g, Median_Extractable_NPOC_mg_per_kg) + 
                                s(site, bs = "re"), data = combined_data, method = "REML"),
  gam_lambda = gam(Median_Weighted_Avg_Lambda ~ s(Median_62948_Final_Gravimetric_Moisture_g_per_g) + 
                     s(Median_Extractable_NPOC_mg_per_kg) + 
                     ti(Median_62948_Final_Gravimetric_Moisture_g_per_g, Median_Extractable_NPOC_mg_per_kg) + 
                     s(site, bs = "re"), data = combined_data, method = "REML")
)

# Extracting model names
model_names <- c("AI_mod", "NOSC", "DBE", "delGcoxPerCmol", "delGcoxPerCompmol", "Lambda")

# Save summaries, diagnostics, and plots for REML models
for (i in seq_along(gam_models)) {
  summary_file <- paste0("05_GAM_Results/GAM_Results/summary_", model_names[i], "_REML.txt")
  diag_file <- paste0("05_GAM_Results/GAM_Results/diagnostics_", model_names[i], "_REML.txt")
  plot_file <- paste0("05_GAM_Results/Diagnostic_Plots/diagnostic_plot_", model_names[i], "_REML.png")
  
  sink(summary_file)
  print(summary(gam_models[[i]]))
  sink()
  
  sink(diag_file)
  gam.check(gam_models[[i]])
  sink()
  
  png(plot_file, width = 2000, height = 1500, res = 300)
  plot(gam_models[[i]], pages = 1)
  dev.off()
}

# Function to perform Breusch-Pagan test and export results
homoscedasticity_test_and_save <- function(model, model_name) {
  bp_test <- bptest(model)
  results <- data.frame(
    Model = model_name,
    BP_p_value = bp_test$p.value
  )
  results_file <- paste0("05_GAM_Results/GAM_Results/BP_test_", model_name, ".csv")
  write.csv(results, results_file, row.names = FALSE)
}

# Apply the test and save results for REML models
for (i in seq_along(gam_models)) {
  homoscedasticity_test_and_save(gam_models[[i]], model_names[i])
}

# Plotting smooth terms for REML models
for (i in seq_along(gam_models)) {
  plot_file <- paste0("05_GAM_Results/Diagnostic_Plots/smooth_plot_", model_names[i], "_REML.png")
  png(plot_file, width = 2000, height = 1500, res = 300)
  plot(gam_models[[i]], pages = 1, all.terms = TRUE)
  dev.off()
}

# AIC comparison for default vs REML method
comparison_results <- data.frame(
  Model = model_names,
  AIC_default = sapply(ggam_models, AIC),
  AIC_REML = sapply(gam_models, AIC)
)
write.csv(comparison_results, "05_GAM_Results/Comparison_Results/AIC_comparison_default_vs_REML.csv")

# Conduct and save diagnostic checks for both default and REML models
for (i in seq_along(ggam_models)) {
  gam.check(ggam_models[[i]])
  gam.check(gam_models[[i]])
}

# File paths
summary_files <- paste0("05_GAM_Results/GAM_Results/summary_", model_names, "_REML.txt")
diag_files <- paste0("05_GAM_Results/GAM_Results/diagnostics_", model_names, "_REML.txt")
plot_files <- paste0("05_GAM_Results/Diagnostic_Plots/diagnostic_plot_", model_names, "_REML.png")

# Save summaries and diagnostics
for (i in seq_along(gam_models)) {
  summary_sink <- summary_files[i]
  diag_sink <- diag_files[i]
  plot_sink <- plot_files[i]
  
  sink(summary_sink)
  print(summary(gam_models[[i]]))
  sink()
  
  sink(diag_sink)
  gam.check(gam_models[[i]])
  sink()
  
  png(plot_sink, width = 2000, height = 1500, res = 300)
  plot(gam_models[[i]], pages = 1)
  dev.off()
}

# Comparison of AIC for default vs REML
comparison_results <- data.frame(
  Model = model_names,
  AIC_default = sapply(ggam_models, AIC),
  AIC_REML = sapply(gam_models, AIC)
)
write.csv(comparison_results, "05_GAM_Results/Comparison_Results/AIC_comparison_default_vs_REML.csv")

gam.check(ggam_models$ggam_ai_mod)
gam.check(gam_models$gam_ai_mod)
gam.check(ggam_models$ggam_delGcoxPerCmol)
gam.check(gam_models$gam_delGcoxPerCmol)
gam.check(gam_models$ggam_models$ggam_delGcoxPerCompmol)
gam.check(gam_models$gam_delGcoxPerCompmol)
gam.check(gam_models$ggam_models$ggam_nosc)
gam.check(gam_models$gam_nosc)
gam.check(ggam_models$ggam_dbe)
gam.check(gam_models$gam_dbe)
gam.check(ggam_models$ggam_lambda)
gam.check(gam_models$gam_lambda)

