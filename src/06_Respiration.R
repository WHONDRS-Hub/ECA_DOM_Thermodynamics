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

# Create directories if they don't exist
if (!dir.exists("06_GAM_Results")) {
  dir.create("06_GAM_Results")
}
if (!dir.exists("06_GAM_Results/GAM_Results")) {
  dir.create("06_GAM_Results/GAM_Results")
}
if (!dir.exists("06_GAM_Results/Diagnostic_Plots")) {
  dir.create("06_GAM_Results/Diagnostic_Plots")
}
if (!dir.exists("06_GAM_Results/Comparison_Results")) {
  dir.create("06_GAM_Results/Comparison_Results")
}

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
# ====== Cube root transform data ====
cube_root <- function(x) sign(x) * (abs(x))^(1/3)

cube_data = data %>% 
  mutate(across(where(is.numeric), cube_root)) %>% # cube root transform data
  rename_with(where(is.numeric), .fn = ~ paste0("cube_", .x))

#  ===== Check model for transformed or untransformed data ====
# Fit models using transformed and untransformed predictors

# Define the refined list of predictors
cube_predictors <- c("cube_Median_Weighted_Avg_AI_mod", "cube_Median_Weighted_Avg_NOSC", "cube_Median_Weighted_Avg_DBE",
                "cube_Median_Weighted_Avg_delGcoxPerCmol", "cube_Median_Weighted_Avg_delGcoxPerCompmol",
                "cube_Median_Weighted_Avg_Lambda", "cube_Median_Extractable_NPOC_mg_per_kg")

predictors <- c("Median_Weighted_Avg_AI_mod", "Median_Weighted_Avg_NOSC", "Median_Weighted_Avg_DBE",
                     "Median_Weighted_Avg_delGcoxPerCmol", "Median_Weighted_Avg_delGcoxPerCompmol",
                     "Median_Weighted_Avg_Lambda", "Median_Extractable_NPOC_mg_per_kg")

# Fit null model
cube_null_model <- gam(cube_Median_Respiration_Rate_mg_DO_per_kg_per_H ~ 1 + s(site, bs = "re"), data = cube_data)

null_model <- gam(Median_Respiration_Rate_mg_DO_per_kg_per_H ~ 1 + s(site, bs = "re"), data = data)

# Initialize list of models with the null model
gam_models <- list(null_model)
names(gam_models) <- "Null"

cube_gam_models <- list(cube_null_model)
names(cube_gam_models) <- "Null_cube"

# Fit models with each predictor and interaction with moisture
for (predictor in predictors) {
  formula <- as.formula(paste("Median_Respiration_Rate_mg_DO_per_kg_per_H ~",
                              paste("te(Median_62948_Final_Gravimetric_Moisture_g_per_g, ", predictor, ")"),
                              "+ s(site, bs = 're')"))
  model <- gam(formula, data = data, method = "REML")
  gam_models[[predictor]] <- model
}

for (predictor in cube_predictors) {
  formula <- as.formula(paste("cube_Median_Respiration_Rate_mg_DO_per_kg_per_H ~",
                              paste("te(cube_Median_62948_Final_Gravimetric_Moisture_g_per_g, ", predictor, ")"),
                              "+ s(site, bs = 're')"))
  model <- gam(formula, data = cube_data, method = "REML")
  cube_gam_models[[predictor]] <- model
}

# Compare models using AIC
model_comparison <- sapply(gam_models, AIC)
model_comparison_cube <- sapply(cube_gam_models, AIC)

# Print model comparison results
print(model_comparison)
print(model_comparison_cube)

# Compare AIC
model_comparison_file <- "06_GAM_Results/Comparison_Results/AIC_comparison_transformed_vs_unstrasnsformed.csv"
model_comparison_data <- as.data.frame(c(model_comparison,model_comparison_cube))
colnames(model_comparison_data) = 'AIC'
write.csv(model_comparison_data, model_comparison_file)

# ==== Export Summary Files =====
for (model_name in names(gam_models)) {
  summary_file <- paste0("06_GAM_Results/GAM_Results/summary_", model_name, ".txt")
  sink(summary_file)
  print(summary(gam_models[[model_name]]))
  sink()
}

for (model_name in names(cube_gam_models)) {
  summary_file <- paste0("06_GAM_Results/GAM_Results/summary_", model_name, ".txt")
  sink(summary_file)
  print(summary(cube_gam_models[[model_name]]))
  sink()
}

# ===== Generate Residual plots and diagnostic plots ====
generate_plots <- function(model, model_name, prefix) {
  # Residual plot
  png(paste0("06_GAM_Results/Diagnostic_Plots/", model_name, "_", prefix, "_residuals.png"), width = 3000, height = 2250, res = 300)
  plot(model, residuals = TRUE, pch = 1)
  dev.off()
  
  # Smooth plot
  png(paste0("06_GAM_Results/Diagnostic_Plots/", model_name, "_", prefix, "_smooth.png"), width = 3000, height = 2250, res = 300)
  plot(model, pages = 1, all.terms = T)
  dev.off()
  
  # QQ plot for residuals
  png(paste0("06_GAM_Results/Diagnostic_Plots/", model_name, "_", prefix, "_qq.png"), width = 3000, height = 2250, res = 300)
  qq.gam(model, rep = 50, level = 1, cex = 5.5)
  dev.off()
  
  # Homoscedasticity
  residuals <- resid(model)
  fitted <- fitted(model)
  png(paste0("06_GAM_Results/Diagnostic_Plots/", model_name, "_", prefix, "_homoscedasticity.png"), width = 3000, height = 2250, res = 300)
  plot(fitted, residuals)
  abline(h = 0, col = "red")
  dev.off()
}

# Generate plots for untransformed models
for (model_name in names(gam_models)) {
  generate_plots(gam_models[[model_name]], model_name, "untransformed")
}

# Generate plots for cube root transformed models
for (model_name in names(cube_gam_models)) {
  generate_plots(cube_gam_models[[model_name]], model_name, "transformed")
}

# Function to perform Breusch-Pagan test and export results
homoscedasticity_tests <- function(model, model_name, prefix) {
  bp_test <- bptest(model)
  results_file <- paste0("06_GAM_Results/GAM_Results/", model_name, "_", prefix,"_BP_test.csv")
  results <- data.frame(
    Model = model_name,
    BP_p_value = bp_test$p.value
  )
  write.csv(results, results_file, row.names = FALSE)
}

# Apply the homoscedasticity test function to each model

for (model_name in names(cube_gam_models)) {
  homoscedasticity_tests(cube_gam_models[[model_name]], model_name, "transformed")
}

for (model_name in names(gam_models)) {
  homoscedasticity_tests(gam_models[[model_name]], model_name, "untransformed")
}

#  ===== Check default vs REML ====
# Transformed data worked better so pnly tested with cube root data

# Define the refined list of predictors
cube_predictors <- c("cube_Median_Weighted_Avg_AI_mod", "cube_Median_Weighted_Avg_NOSC", "cube_Median_Weighted_Avg_DBE",
                     "cube_Median_Weighted_Avg_delGcoxPerCmol", "cube_Median_Weighted_Avg_delGcoxPerCompmol",
                     "cube_Median_Weighted_Avg_Lambda", "cube_Median_Extractable_NPOC_mg_per_kg")

# Fit null model
cube_null_model <- gam(cube_Median_Respiration_Rate_mg_DO_per_kg_per_H ~ 1 + s(site, bs = "re"), data = cube_data)

# Initialize list of models with the null model
gam_models_REML <- list(cube_null_model)
names(gam_models_REML) <- "Null"

cube_gam_models <- list(cube_null_model)
names(cube_gam_models) <- "Null_cube"

# Fit models with each predictor and interaction with moisture
for (predictor in cube_predictors) {
  formula <- as.formula(paste("cube_Median_Respiration_Rate_mg_DO_per_kg_per_H ~",
                              paste("te(cube_Median_62948_Final_Gravimetric_Moisture_g_per_g, ", predictor, ")"),
                              "+ s(site, bs = 're')"))
  model <- gam(formula, data = cube_data)
  gam_models_REML[[predictor]] <- model
}

for (predictor in cube_predictors) {
  formula <- as.formula(paste("cube_Median_Respiration_Rate_mg_DO_per_kg_per_H ~",
                              paste("te(cube_Median_62948_Final_Gravimetric_Moisture_g_per_g, ", predictor, ")"),
                              "+ s(site, bs = 're')"))
  model <- gam(formula, data = cube_data, method = "REML")
  cube_gam_models[[predictor]] <- model
}

# Compare models using AIC
model_comparison <- sapply(gam_models_REML, AIC)
model_comparison_cube <- sapply(cube_gam_models, AIC)

# Print model comparison results
print(model_comparison)
print(model_comparison_cube)

# Compare AIC
model_comparison_file <- "06_GAM_Results/Comparison_Results/AIC_comparison_REML_vs_Default.csv"
# Compare models using AIC
model_comparison_REML <- sapply(gam_models_REML, AIC)
model_comparison_default <- sapply(cube_gam_models, AIC)

# Combine model comparisons into a data frame
model_comparison_data <- data.frame(
  Model = c(names(model_comparison_REML), names(model_comparison_default)),
  AIC = c(model_comparison_REML, model_comparison_default)
)

write.csv(model_comparison_data, model_comparison_file)


gam.check(gam_models_REML$Null)
gam.check(cube_gam_models$Null_cube)
gam.check(gam_models_REML$cube_Median_Weighted_Avg_AI_mod)
gam.check(cube_gam_models$cube_Median_Weighted_Avg_AI_mod)
gam.check(gam_models_REML$cube_Median_Weighted_Avg_NOSC)
gam.check(cube_gam_models$cube_Median_Weighted_Avg_NOSC)
gam.check(gam_models_REML$cube_Median_Weighted_Avg_DBE)
gam.check(cube_gam_models$cube_Median_Weighted_Avg_DBE)
gam.check(gam_models_REML$cube_Median_Weighted_Avg_delGcoxPerCmol)
gam.check(cube_gam_models$cube_Median_Weighted_Avg_delGcoxPerCmol)
gam.check(gam_models_REML$cube_Median_Weighted_Avg_delGcoxPerCompmol)
gam.check(cube_gam_models$cube_Median_Weighted_Avg_delGcoxPerCompmol)

# ==== Export Summary Files =====
for (model_name in names(gam_models)) {
  summary_file <- paste0("06_GAM_Results/GAM_Results/summary_", model_name, ".txt")
  sink(summary_file)
  print(summary(gam_models[[model_name]]))
  sink()
}

for (model_name in names(cube_gam_models)) {
  summary_file <- paste0("06_GAM_Results/GAM_Results/summary_", model_name, ".txt")
  sink(summary_file)
  print(summary(cube_gam_models[[model_name]]))
  sink()
}

# ===== Generate Residual plots and diagnostic plots ====

# Generate plots for REML models
for (model_name in names(gam_models_REML)) {
  generate_plots(gam_models_REML[[model_name]], model_name, "transformed_REML")
}


# Function to perform Breusch-Pagan test and export results
homoscedasticity_tests <- function(model, model_name, prefix) {
  bp_test <- bptest(model)
  results_file <- paste0("06_GAM_Results/GAM_Results/", model_name, "_", prefix,"_BP_test.csv")
  results <- data.frame(
    Model = model_name,
    BP_p_value = bp_test$p.value
  )
  write.csv(results, results_file, row.names = FALSE)
}

# Apply the homoscedasticity test function to each model

for (model_name in names(gam_models_REML)) {
  homoscedasticity_tests(gam_models_REML[[model_name]], model_name, "transformed_REML")
}

# ======== Final GAM ======
rm(gam_models);rm(gam_models_REML);rm(cube_gam_models); rm(null_model); rm(cube_null_model)

# Based on AIC, Convergence and simplicity it seems that the transformed model with default settings is the way to go

# Fit null model
null_model <- gam(cube_Median_Respiration_Rate_mg_DO_per_kg_per_H ~ 1 + s(site, bs = "re"), data = cube_data)

# Initialize list of models with the null model
gam_models <- list(null_model)
names(gam_models) <- "Null"

# Fit models with each predictor and interaction with moisture
for (predictor in cube_predictors) {
  formula <- as.formula(paste("cube_Median_Respiration_Rate_mg_DO_per_kg_per_H ~",
                              paste("te(cube_Median_62948_Final_Gravimetric_Moisture_g_per_g, ", predictor, ")"),
                              "+ s(site, bs = 're')"))
  model <- gam(formula, data = cube_data)
  gam_models[[predictor]] <- model
}

# Compare models using AIC
model_comparison <- sapply(gam_models, AIC)

# Print model comparison results
print(model_comparison)

# Identify the best model based on the lowest AIC
best_model_name <- names(model_comparison)[which.min(model_comparison)]
best_model <- gam_models[[best_model_name]]
cat("Best model based on the lowest AIC:", best_model_name, "\n")

# Summary of the best model
best_summary <- summary(best_model)
print(best_summary)

# Perform ANOVA to compare null model and best model
anova_results <- anova(null_model, best_model, test = "Chisq")
print(anova_results)

# Visualize the interaction effect
par(mfrow=c(1,1))
plot(best_model, pages = 1, rug = TRUE, seWithMean = TRUE)

# ===== Interaction Plots ====

library(mgcv)
vis.gam(best_model, view = c("cube_Median_62948_Final_Gravimetric_Moisture_g_per_g", 
                             "cube_Median_Weighted_Avg_delGcoxPerCmol"), 
        plot.type = "contour", color = "topo", 
        main = "Interaction between Moisture and delGcoxPerCmol on Respiration Rate")

library(mgcv)
library(viridis)

pdf(paste0(figure_path,"Figure4_Interaction_plot.pdf"), width = 3000, height = 2250, res = 300)
vis.gam(
  best_model, 
  view = c("cube_Median_62948_Final_Gravimetric_Moisture_g_per_g", "cube_Median_Weighted_Avg_delGcoxPerCmol"), 
  plot.type = "contour", 
  color = 'heat',
  xlab = expression(paste("Cubic Root Median Gravimetric Moisture ", (g/g))),
  ylab =   expression(paste("Cubic Root Median ", Delta, G[cox],~(kJ~Cmol^{-1}))),
  main = ''
)
dev.off()

pdf(paste0(figure_path,"Figure4_Interaction_plot_color2.pdf"))
# Create the vis.gam plot with adjusted font size and margins
vis.gam(
  best_model, 
  view = c("cube_Median_62948_Final_Gravimetric_Moisture_g_per_g", "cube_Median_Weighted_Avg_delGcoxPerCmol"), 
  plot.type = "contour", 
  color = 'topo',
  xlab = expression(paste("Cubic Root Median Gravimetric Moisture ", (g/g))),
  ylab =   expression(paste("Cubic Root Median ", Delta, G[cox],~(kJ~Cmol^{-1}))),
  main = ''
)
dev.off()



png(paste0(figure_path,"Figure4_Interaction_plot.png"), width = 3000, height = 2250, res = 300)
par(mar = c(5, 6, 4, 2) + 0.1)  # Adjust the left margin to be larger

# Create the vis.gam plot with adjusted font size and margins
vis.gam(
  best_model, 
  view = c("cube_Median_62948_Final_Gravimetric_Moisture_g_per_g", "cube_Median_Weighted_Avg_delGcoxPerCmol"), 
  plot.type = "contour", 
  color = "heat",
  xlab = "",  # Leave empty to manually set the label
  ylab = "",  # Leave empty to manually set the label
  main = '',
  cex.lab = 1.2,   # Adjusts the font size of the axis labels
  cex.axis = 1.0   # Adjusts the font size of the tick labels
)

# Manually add axis labels with adjusted positions
mtext(side = 1, text = expression(paste("Cubic Root Median Gravimetric Moisture ", (g/g))), line = 3, cex = 1.2)
mtext(side = 2, text = expression(paste("Cubic Root Median ", Delta, G[cox], ~(kJ~Cmol^{-1}))), line = 4, cex = 1.2)

dev.off()


png(paste0(figure_path,"Figure4_Interaction_plot_color2.png"), width = 3300, height = 2250, res = 300)
par(mar = c(5, 6, 4, 2) + 0.1)  # Adjust the left margin to be larger

# Create the vis.gam plot with adjusted font size and margins
vis.gam(
  best_model, 
  view = c("cube_Median_62948_Final_Gravimetric_Moisture_g_per_g", "cube_Median_Weighted_Avg_delGcoxPerCmol"), 
  plot.type = "contour", 
  color = "topo",
  xlab = "",  # Leave empty to manually set the label
  ylab = "",  # Leave empty to manually set the label
  main = '',
  cex.lab = 1.2,   # Adjusts the font size of the axis labels
  cex.axis = 1.0   # Adjusts the font size of the tick labels
)

# Manually add axis labels with adjusted positions
mtext(side = 1, text = expression(paste("Cubic Root Median Gravimetric Moisture ", (g/g))), line = 3, cex = 1.2)
mtext(side = 2, text = expression(paste("Cubic Root Median ", Delta, G[cox], ~(kJ~Cmol^{-1}))), line = 4, cex = 1.2)

dev.off()

# ====== Other plot options for interaction effects =====
library(ggplot2)
library(mgcv)

# Generate a grid of moisture and delGcoxPerCmol values
moisture_seq <- seq(min(cube_data$cube_Median_62948_Final_Gravimetric_Moisture_g_per_g), 
                    max(cube_data$cube_Median_62948_Final_Gravimetric_Moisture_g_per_g), 
                    length.out = 100)
delGcox_seq <- seq(min(cube_data$cube_Median_Weighted_Avg_delGcoxPerCmol), 
                   max(cube_data$cube_Median_Weighted_Avg_delGcoxPerCmol), 
                   length.out = 100)

# Create a new dataset for prediction
prediction_grid <- expand.grid(
  cube_Median_62948_Final_Gravimetric_Moisture_g_per_g = moisture_seq,
  cube_Median_Weighted_Avg_delGcoxPerCmol = delGcox_seq,
  site = levels(cube_data$site)[1]  # Use a reference site or average
)

# Predict respiration rates
prediction_grid$Predicted_Respiration <- predict(best_model, newdata = prediction_grid)

# Plot the interaction using a heatmap
ggplot(prediction_grid, aes(x = cube_Median_62948_Final_Gravimetric_Moisture_g_per_g, 
                            y = cube_Median_Weighted_Avg_delGcoxPerCmol, 
                            fill = Predicted_Respiration)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Predicted Respiration") +
  labs(title = "Interaction Effect of Moisture and delGcoxPerCmol on Respiration Rates",
       x = "Cube-Root Transformed Moisture (g/g)",
       y = "Cube-Root Transformed delGcoxPerCmol") +
  theme_minimal()
