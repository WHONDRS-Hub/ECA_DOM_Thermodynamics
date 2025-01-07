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
library(tidyverse);library(corrplot);library(ggpubr);library(ggpmisc);library(factoextra);library(stringr)
# ==== Defining paths and working directories ======
github = 'C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/ECA_DOM_Thermodynamics/'
data_path = paste0(github,'Data/')
figure_path = paste0(github,'Figures/')
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
data = data %>%
  dplyr::filter(site != "EC_023") # Remove rows where site is EC_023 because it has NA in some of the explanatory data

data = merge(data,sample_data, by = c('site','Treatment'))

data = data %>%
  dplyr::select(-Field_Name,-Sample_Name,-IGSN, -Material, -"Median_62948_Initial_Gravimetric_Moisture_g_per_g",- "Median_Fe_mg_per_L",-"Median_ATP_nanomoles_per_L",-"Median_Extractable_NPOC_mg_per_L",-"Median_Extractable_TN_mg_per_L",-"Median_Missing_Reps")

data = merge(data,grainsize_data, by = 'site')

data$Treatment <- as.factor(data$Treatment)
data$site <- as.factor(data$site)
# ====== Diagnostic plots ======
# Load necessary libraries
library(mgcv)
library(dplyr)
library(purrr)
library(broom)
library(ggplot2)
library(tidyr)
library(rsample) # For cross-validation
library(yardstick) # For performance metrics
library(forcats) # For factor handling


library(ggplot2)

# Original respiration variable
ggplot(data.frame(resp = all_data$Median_Respiration_Rate_mg_DO_per_kg_per_H), aes(x = resp)) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.5) +
  labs(title = "Histogram of Original Respiration Variable") +
  theme_minimal()

# Transformed respiration variable (cube root)
ggplot(data.frame(resp = all_data$cube_Median_Respiration_Rate_mg_DO_per_kg_per_H), aes(x = resp)) +
  geom_histogram(bins = 30, fill = "green", alpha = 0.5) +
  labs(title = "Histogram of Transformed Respiration Variable (Cube Root)") +
  theme_minimal()

# QQ plot for original variable
ggplot(data.frame(res = all_data$Median_Respiration_Rate_mg_DO_per_kg_per_H), 
       aes(sample = res)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "QQ Plot of Original Respiration Variable") +
  theme_minimal()

# QQ plot for transformed variable
ggplot(data.frame(res = all_data$cube_Median_Respiration_Rate_mg_DO_per_kg_per_H), 
       aes(sample = res)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "QQ Plot of Transformed Respiration Variable") +
  theme_minimal()

# Shapiro-Wilk test for original variable
shapiro.test(all_data$cube_Median_Respiration_Rate_mg_DO_per_kg_per_H)

# Shapiro-Wilk test for transformed variable
shapiro.test(all_data$cube_Median_Respiration_Rate_mg_DO_per_kg_per_H)

# Fit model with original variable
model_original <- gam(Median_Respiration_Rate_mg_DO_per_kg_per_H ~ te(Median_62948_Final_Gravimetric_Moisture_g_per_g, Median_Extractable_NPOC_mg_per_kg) + s(site, bs = 're'), data = all_data)

# Fit model with transformed variable
model_transformed <- gam(cube_Median_Respiration_Rate_mg_DO_per_kg_per_H ~ te(Median_62948_Final_Gravimetric_Moisture_g_per_g, Median_Extractable_NPOC_mg_per_kg) + s(site, bs = 're'), data = all_data)


# Extract residuals from both models
residuals_original <- resid(model_original)
residuals_transformed <- resid(model_transformed)

# Histogram of residuals for both models
ggplot(data.frame(res = residuals_original), aes(x = res)) +
  geom_histogram(bins = 30, fill = "red", alpha = 0.5) +
  labs(title = "Histogram of Residuals from Original Model") +
  theme_minimal()

ggplot(data.frame(res = residuals_transformed), aes(x = res)) +
  geom_histogram(bins = 30, fill = "orange", alpha = 0.5) +
  labs(title = "Histogram of Residuals from Transformed Model") +
  theme_minimal()

# QQ plot of residuals for the original model
ggplot(data.frame(res = residuals_original), aes(sample = res)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "QQ Plot of Residuals from Original Model") +
  theme_minimal()

# QQ plot of residuals for the transformed model
ggplot(data.frame(res = residuals_transformed), aes(sample = res)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "QQ Plot of Residuals from Transformed Model") +
  theme_minimal()

# Shapiro-Wilk test for residuals
shapiro.test(residuals_original)
shapiro.test(residuals_transformed)



# Define a function to fit GAM and perform diagnostics
fit_gam <- function(data, response_var, predictors, random_effect, family = "gaussian") {
  # Print current response variable being processed
  cat("Fitting GAM for:", response_var, "\n")
  
  # Check response variable
  response_data <- data[[response_var]]
  cat("Class:", class(response_data), "\n")
  cat("Range:", range(response_data, na.rm = TRUE), "\n")
  cat("Number of non-NA values:", sum(!is.na(response_data)), "\n\n")
  
  # Construct the formula
  formula <- as.formula(
    paste(
      response_var, "~",
      paste(predictors, collapse = " + "), "+",
      paste0("s(", random_effect, ", bs = 're')")
    )
  )
  
  # Fit the GAM
  tryCatch({
    gam_model <- gam(formula, data = data, family = family, method = "REML")
    
    # Summary of the model
    model_summary <- summary(gam_model)
    
    # Extract residuals and fitted values
    residuals <- resid(gam_model, type = "pearson")
    fitted_vals <- fitted(gam_model)
    
    # Create diagnostics plot to visualize residuals against fitted values
    p1 <- ggplot(data = NULL, aes(x = fitted_vals, y = residuals)) +
      geom_point(alpha = 0.5) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      labs(title = paste("Residuals vs Fitted for", response_var),
           x = "Fitted Values",
           y = "Pearson Residuals") +
      theme_minimal()
    
    # QQ Plot
    p2 <- ggplot(data = NULL, aes(sample = residuals)) +
      stat_qq(alpha = 0.5) +
      stat_qq_line(color = "red") +
      labs(title = paste("QQ Plot for", response_var),
           x = "Theoretical Quantiles",
           y = "Sample Quantiles") +
      theme_minimal()
    
    # Diagnostic plots for each predictor variable
    for (pred in predictors) {
      pred_effect <- predict(gam_model, newdata = data, type = "terms")[, pred]
      
      # Partial effects plot
      p_effect <- ggplot(data, aes_string(x = pred, y = pred_effect)) +
        geom_point() +
        geom_smooth(method = "loess") +
        labs(title = paste("Partial Effect of", pred, "on", response_var),
             y = "Partial Effect",
             x = pred) +
        theme_minimal()
      
      # Save the partial effect plot
      ggsave(filename = paste0("GAM_Plots/Partial_Effect_", pred, "_", response_var, ".png"),
             plot = p_effect, width = 7, height = 5)
    }
    
    # Return a list of outputs
    list(
      model = gam_model,
      summary = model_summary,
      residuals_vs_fitted = p1,
      qq_plot = p2
    )
  }, error = function(e) {
    message(paste("Error in fitting GAM for", response_var, ":", e$message))
    return(NULL)
  })
}

# Define response variables (including your cube respiration variable)
response_vars <- c(
  "Median_Weighted_Avg_AI_mod",
  "Median_Weighted_Avg_NOSC",
  "Median_Weighted_Avg_DBE",
  "Median_Weighted_Avg_delGcoxPerCmol",
  "Median_Weighted_Avg_delGcoxPerCompmol",
  "Median_Weighted_Avg_Lambda"
)

# Define predictor variables including cubed respiration
predictor_vars <- c(
  "Median_62948_Final_Gravimetric_Moisture_g_per_g",
  "Median_Extractable_NPOC_mg_per_kg",
  "cube_Median_Respiration_Rate_mg_DO_per_kg_per_H"
)

# Create your existing data preparation process
combined_data_clean <- all_data %>%
  select(all_of(response_vars), all_of(predictor_vars), site) %>%
  mutate(across(all_of(c(response_vars, predictor_vars)), as.numeric)) %>%
  mutate(site = factor(site)) %>%
  drop_na()

gam_results <- map(response_vars, ~ fit_gam(
  data = combined_data_clean,
  response_var = .x,
  predictors = predictor_vars,  # You can use 'predictors' instead, if defined
  random_effect = random_effect,
  family = response_families[[.x]]
))

# Iterate to plot the defined p1 and p2 plots for all responses
for (name in names(gam_results)) {
  result <- gam_results[[name]]
  if (!is.null(result)) {
    ggsave(
      filename = paste0("GAM_Plots/Residuals_vs_Fitted_", name, ".png"),
      plot = result$residuals_vs_fitted,
      width = 7, height = 5
    )
    
    ggsave(
      filename = paste0("GAM_Plots/QQ_Plot_", name, ".png"),
      plot = result$qq_plot,
      width = 7, height = 5
    )
  }
}

# Note: Ensure "GAM_Plots" directory exists before running.

# ======== GAM ======
# Transformation for normalization is cube root - have to cube root then add sign back to value to make it positive or negative
cube_root <- function(x) sign(x) * (abs(x))^(1/3)

cube_data = data %>% 
  mutate(across(where(is.numeric), cube_root)) %>% # cube root transform data
  rename_with(where(is.numeric), .fn = ~ paste0("cube_", .x))

# Define the refined list of predictors
predictors <- c("cube_Median_Weighted_Avg_AI_mod", "cube_Median_Weighted_Avg_NOSC", "cube_Median_Weighted_Avg_DBE",
                "cube_Median_Weighted_Avg_delGcoxPerCmol", "cube_Median_Weighted_Avg_delGcoxPerCompmol",
                "cube_Median_Weighted_Avg_Lambda", "cube_Median_Extractable_NPOC_mg_per_kg")

# Fit null model
null_model <- gam(cube_Median_Respiration_Rate_mg_DO_per_kg_per_H ~ 1 + s(site, bs = "re"), data = cube_data)

# Initialize list of models with the null model
gam_models <- list(null_model)
names(gam_models) <- "Null"

# Fit models with each predictor and interaction with moisture
for (predictor in predictors) {
  formula <- as.formula(paste("cube_Median_Respiration_Rate_mg_DO_per_kg_per_H ~",
                              paste("te(cube_Median_62948_Final_Gravimetric_Moisture_g_per_g, ", predictor, ")"),
                              "+ s(site, bs = 're')"))
  model <- gam(formula, data = cube_data, method = "REML")
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

# ===== Cross Validation ====
library(caret)
# Function to fit GAM model with interactions
fit_gam_model <- function(train_data) {
  formula <- cube_Median_Respiration_Rate_mg_DO_per_kg_per_H ~ te(cube_Median_62948_Final_Gravimetric_Moisture_g_per_g, 
                                                                  cube_Median_Weighted_Avg_delGcoxPerCmol) + 
    s(site, bs = "re")
  model <- gam(formula, data = train_data, method = "REML")
  return(model)
}

# Ensure that all levels are considered in the model
cube_data$site <- factor(cube_data$site)

# Create stratified k-fold cross-validation
set.seed(123)  # For reproducibility

folds <- createFolds(cube_data$site, k = 10, list = TRUE, returnTrain = TRUE)
# Perform cross-validation
cv_results <- lapply(folds, function(train_indices) {
  train_data <- cube_data[train_indices, ]
  test_data <- cube_data[-train_indices, ]
  
  # Ensure test data factor levels match those of the training data
  test_data$site <- factor(test_data$site, levels = levels(train_data$site))
  
  # Fit model on training data
  model <- fit_gam_model(train_data)
  
  # Predict on test data
  predictions <- predict(model, newdata = test_data, type="response")
  
  # Filter predictions and observations for only present factor levels
  valid_indices <- which(!is.na(predictions))
  valid_predictions <- predictions[valid_indices]
  valid_observations <- test_data$cube_Median_Respiration_Rate_mg_DO_per_kg_per_H[valid_indices]
  
  # Calculate RMSE using valid predictions
  rmse <- sqrt(mean((valid_observations - valid_predictions)^2))
  return(rmse)
})

# Mean and standard deviation of RMSE across folds
mean_rmse <- mean(unlist(cv_results))
sd_rmse <- sd(unlist(cv_results))

cat("Mean RMSE from cross-validation:", mean_rmse, "\n")
cat("SD of RMSE from cross-validation:", sd_rmse, "\n")


library(mgcv)
vis.gam(best_model, view = c("cube_Median_62948_Final_Gravimetric_Moisture_g_per_g", 
                             "cube_Median_Weighted_Avg_delGcoxPerCmol"), 
        plot.type = "contour", color = "heat", 
        main = "Interaction between Moisture and delGcoxPerCmol on Respiration Rate")

library(mgcv)
library(viridis)

pdf(paste0(figure_path,"Figure4_Interaction_plot.pdf"))
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

png(paste0(figure_path,"Figure4_Interaction_plot.png"))
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

# ===== Testing GAMM ====
library(mgcv)

gamm_respiration <- gam(cube_Median_Respiration_Rate_mg_DO_per_kg_per_H ~ 
                          te(cube_Median_62948_Final_Gravimetric_Moisture_g_per_g, 
                             cube_Median_Weighted_Avg_delGcoxPerCmol) + 
                          te(cube_Median_Weighted_Avg_AI_mod, cube_Median_Weighted_Avg_Lambda) + 
                          s(site, bs = "re"),
                        data = cube_data,
                        method = "REML")

summary(gamm_respiration)
plot(gamm_respiration, pages = 1, shade = TRUE)


# Transformation for normalization is cube root - have to cube root then add sign back to value to make it positive or negative
cube_root <- function(x) sign(x) * (abs(x))^(1/3)

cube_data = data %>% 
  mutate(across(where(is.numeric), cube_root)) %>% # cube root transform data
  rename_with(where(is.numeric), .fn = ~ paste0("cube_", .x))

# Define the refined list of predictors
predictors <- c("Median_Weighted_Avg_AI_mod", "Median_Weighted_Avg_NOSC", "Median_Weighted_Avg_DBE",
                "Median_Weighted_Avg_delGcoxPerCmol", "Median_Weighted_Avg_delGcoxPerCompmol",
                "Median_Weighted_Avg_Lambda", "Median_Extractable_NPOC_mg_per_kg")

# Fit null model
null_model <- gam(Median_Respiration_Rate_mg_DO_per_kg_per_H ~ 1 + s(site, bs = "re"), data = dat)

# Initialize list of models with the null model
gam_models <- list(null_model)
names(gam_models) <- "Null"

# Fit models with each predictor and interaction with moisture
for (predictor in predictors) {
  formula <- as.formula(paste("Median_Respiration_Rate_mg_DO_per_kg_per_H ~",
                              paste("te(Median_62948_Final_Gravimetric_Moisture_g_per_g, ", predictor, ")"),
                              "+ s(site, bs = 're')"))
  model <- gam(formula, data = dat, method = "REML")
  gam_models[[predictor]] <- model
}

# Compare models using AIC
model_comparison <- sapply(gam_models, AIC)

# Print model comparison results
print(model_comparison)

# Identify the best model ba("Best model based on the lowest AIC:", best_model_name, "\n")sed on the lowest AIC
best_model_name <- names(model_comparison)[which.min(model_comparison)]
best_model <- gam_models[[best_model_name]]
cat("Best model based on the lowest AIC:", best_model_name, "\n")


# Summary of the best model
best_summary <- summary(best_model)
print(best_summary)

# Perform ANOVA to compare null model and best model

# ====== GAM with diagnostics ======
# Load necessary libraries
library(mgcv)
library(caret)
library(ggplot2)
library(viridis)
library(Metrics)  # For additional metrics
# Ensure all libraries are installed; install if necessary
# install.packages(c("mgcv", "caret", "ggplot2", "viridis", "Metrics"))

# ===== Data Transformation =====
# Transformation for normalization: cube root with sign preservation
cube_root <- function(x) sign(x) * (abs(x))^(1/3)
cube_data <- data %>%
  mutate(across(where(is.numeric), cube_root)) %>% # Apply cube root transformation
  rename_with(where(is.numeric), .fn = ~ paste0("cube_", .x))

# Define the refined list of predictors
predictors <- c(
  "cube_Median_Weighted_Avg_AI_mod", "cube_Median_Weighted_Avg_NOSC", 
  "cube_Median_Weighted_Avg_DBE", "cube_Median_Weighted_Avg_delGcoxPerCmol", 
  "cube_Median_Weighted_Avg_delGcoxPerCompmol", "cube_Median_Weighted_Avg_Lambda", 
  "cube_Median_Extractable_NPOC_mg_per_kg"
)

predictors <- c(
  "Median_Weighted_Avg_AI_mod", "Median_Weighted_Avg_NOSC", 
  "Median_Weighted_Avg_DBE", "Median_Weighted_Avg_delGcoxPerCmol", 
  "Median_Weighted_Avg_delGcoxPerCompmol", "Median_Weighted_Avg_Lambda", 
  "Median_Extractable_NPOC_mg_per_kg"
)

all_data = merge(cube_data,combined_data, by = c('site','Treatment'))

# Load necessary libraries
library(mgcv)
library(dplyr)
library(purrr)
library(broom)
library(ggplot2)
library(tidyr)
library(rsample) # For cross-validation
library(yardstick) # For performance metrics
library(forcats) # For factor handling

# Define a function to fit GAM and perform diagnostics
fit_gam <- function(data, response_var, predictors, random_effect, family = "gaussian") {
  # Print current response variable being processed
  cat("Fitting GAM for:", response_var, "\n")
  
  # Check response variable
  response_data <- data[[response_var]]
  cat("Class:", class(response_data), "\n")
  cat("Range:", range(response_data, na.rm = TRUE), "\n")
  cat("Number of non-NA values:", sum(!is.na(response_data)), "\n\n")
  
  # Construct the formula
  formula <- as.formula(
    paste(
      response_var, "~",
      paste(predictors, collapse = " + "), "+",
      paste0("s(", random_effect, ", bs = 're')")
    )
  )
  
  # Fit the GAM
  tryCatch({
    gam_model <- gam(formula, data = data, family = family, method = "REML")
    
    # Summary of the model
    model_summary <- summary(gam_model)
    
    # Extract residuals and fitted values
    residuals <- resid(gam_model, type = "pearson")
    fitted_vals <- fitted(gam_model)
    
    # Diagnostic Plots
    # 1. Residuals vs Fitted
    p1 <- ggplot(data = NULL, aes(x = fitted_vals, y = residuals)) +
      geom_point(alpha = 0.5) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      labs(title = paste("Residuals vs Fitted for", response_var),
           x = "Fitted Values",
           y = "Pearson Residuals") +
      theme_minimal()
    
    # 2. QQ Plot
    p2 <- ggplot(data = NULL, aes(sample = residuals)) +
      stat_qq(alpha = 0.5) +
      stat_qq_line(color = "red") +
      labs(title = paste("QQ Plot for", response_var),
           x = "Theoretical Quantiles",
           y = "Sample Quantiles") +
      theme_minimal()
    
    # Return a list of outputs
    list(
      model = gam_model,
      summary = model_summary,
      residuals_vs_fitted = p1,
      qq_plot = p2
    )
  }, error = function(e) {
    message(paste("Error in fitting GAM for", response_var, ":", e$message))
    return(NULL)
  })
}

# Define a function to perform k-fold cross-validation on GAM models
cv_gam <- function(data, response_var, predictors, random_effect, family = "gaussian", k = 5) {
  # Print current response variable being processed
  cat("Performing", k, "-fold Cross-Validation for:", response_var, "\n")
  
  # Define the formula
  formula <- as.formula(
    paste(
      response_var, "~",
      paste(predictors, collapse = " + "), "+",
      paste0("s(", random_effect, ", bs = 're')")
    )
  )
  
  # Create k-fold cross-validation splits
  set.seed(123)  # For reproducibility
  cv_splits <- vfold_cv(data, v = k)
  
  # Define a helper function to fit model and predict on the test set
  fit_and_predict <- function(split) {
    # Extract training and testing data
    train_data <- analysis(split)
    test_data <- assessment(split)
    
    # Ensure 'site' as a factor with levels from the training set
    train_data <- train_data %>%
      mutate(site = fct_explicit_na(factor(site), na_level = "Unknown"))
    
    test_data <- test_data %>%
      mutate(site = ifelse(!(site %in% levels(train_data$site)), "Unknown", as.character(site)))
    
    # Fit the GAM model on training data
    gam_model <- tryCatch({
      gam(formula, data = train_data, family = family, method = "REML")
    }, error = function(e) {
      message(paste("Error fitting model:", e$message))
      return(NULL)
    })
    
    if (is.null(gam_model)) {
      return(tibble(RMSE = NA, Rsquared = NA))
    }
    
    # Predict on test data
    predictions <- tryCatch({
      predict(gam_model, newdata = test_data)
    }, warning = function(w) {
      message(paste("Warning during prediction for", response_var, ":", w$message))
      return(rep(NA, nrow(test_data)))
    }, error = function(e) {
      message(paste("Error during prediction for", response_var, ":", e$message))
      return(rep(NA, nrow(test_data)))
    })
    
    # Actual values
    actuals <- test_data[[response_var]]
    
    # Remove NA predictions and actuals
    valid_indices <- !is.na(predictions) & !is.na(actuals)
    actuals <- actuals[valid_indices]
    predictions <- predictions[valid_indices]
    
    # Check if we have enough data to calculate the metrics
    if (length(predictions) == 0) {
      return(tibble(RMSE = NA, Rsquared = NA))
    }
    
    # Calculate performance metrics
    metrics <- metric_set(rmse, rsq)(data = tibble(
      truth = actuals,
      estimate = predictions
    ))
    
    return(metrics)
  }
  
  # Apply the helper function to each split and combine results
  cv_metrics <- map_dfr(cv_splits$splits, fit_and_predict)
  
  # Calculate average metrics, excluding any NAs
  avg_metrics <- cv_metrics %>%
    summarize(
      RMSE = mean(rmse, na.rm = TRUE),
      Rsquared = mean(rsq, na.rm = TRUE)  # Ensure 'rsq' is the correct name
    )
  
  return(avg_metrics)
}

# Define response variables
response_vars <- c(
  "Median_Weighted_Avg_AI_mod",
  "Median_Weighted_Avg_NOSC",
  "Median_Weighted_Avg_DBE",
  "Median_Weighted_Avg_delGcoxPerCmol",
  "Median_Weighted_Avg_delGcoxPerCompmol",
  "Median_Weighted_Avg_Lambda"
)

# Define predictor variables
predictor_vars <- c(
  "Median_62948_Final_Gravimetric_Moisture_g_per_g",
  "Median_Extractable_NPOC_mg_per_kg"
)

# Define families
response_families <- list(
  "Median_Weighted_Avg_AI_mod" = "gaussian",
  "Median_Weighted_Avg_NOSC" = "gaussian",
  "Median_Weighted_Avg_DBE" = "gaussian",
  "Median_Weighted_Avg_delGcoxPerCmol" = "gaussian",
  "Median_Weighted_Avg_delGcoxPerCompmol" = "gaussian",
  "Median_Weighted_Avg_Lambda" = "gaussian"
)

# Define predictors with smooth terms and interactions
predictors <- c(
  "s(Median_62948_Final_Gravimetric_Moisture_g_per_g)",
  "s(Median_Extractable_NPOC_mg_per_kg)",
  "ti(Median_62948_Final_Gravimetric_Moisture_g_per_g, Median_Extractable_NPOC_mg_per_kg)"
)

# Define random effect
random_effect <- "site"

# Prepare and clean the data
combined_data_clean <- combined_data %>%
  select(all_of(response_vars), all_of(predictor_vars), site) %>%
  mutate(across(all_of(c(response_vars, predictor_vars)), as.numeric)) %>%
  mutate(site = factor(site)) %>%
  drop_na()

# Fit GAMs for all response variables with specified families
gam_results <- map(response_vars, ~ fit_gam(
  data = combined_data_clean,
  response_var = .x,
  predictors = predictors,
  random_effect = random_effect,
  family = response_families[[.x]]
))

# Name the list elements
names(gam_results) <- response_vars

# Perform cross-validation for each response variable
cv_results <- map(response_vars, ~ cv_gam(
  data = combined_data_clean,
  response_var = .x,
  predictors = predictors,
  random_effect = random_effect,
  family = response_families[[.x]],
  k = 5 # Number of folds
))

# Name the list elements
names(cv_results) <- response_vars

# Extract tidy model summaries
tidy_summaries <- gam_results %>%
  map(~ if(!is.null(.x)) tidy(.x$model) else NULL)

# Print tidy summaries
walk(names(tidy_summaries), function(name) {
  cat("\nTidy Summary for:", name, "\n")
  if (!is.null(tidy_summaries[[name]])) {
    print(tidy_summaries[[name]])
  } else {
    cat("Model fitting failed or no tidy summary available.\n")
  }
})

# Print Cross-Validation Results
walk(names(cv_results), function(name) {
  cat("\nCross-Validation Metrics for:", name, "\n")
  print(cv_results[[name]])
})

# Save diagnostic plots
if(!dir.exists("GAM_Plots")) {
  dir.create("GAM_Plots")
}

# Iterate over each model to save plots
for (name in names(gam_results)) {
  result <- gam_results[[name]]
  if (!is.null(result)) {
    ggsave(
      filename = paste0("GAM_Plots/Residuals_vs_Fitted_", name, ".png"),
      plot = result$residuals_vs_fitted,
      width = 7, height = 5
    )
    
    ggsave(
      filename = paste0("GAM_Plots/QQ_Plot_", name, ".png"),
      plot = result$qq_plot,
      width = 7, height = 5
    )
    
    # Save Partial Effects Plots
    png(filename = paste0("GAM_Plots/Partial_Effects_", name, ".png"), width = 800, height = 600)
    plot(result$model, pages = 1, scheme = 1, main = paste("Partial Effects for", name))
    dev.off()
  }
}

# Combine all tidy summaries into one dataframe with a Response column
combined_tidy <- bind_rows(
  map2(tidy_summaries, names(tidy_summaries), ~ {
    if(!is.null(.x)) {
      mutate(.x, Response = .y)
    } else {
      NULL
    }
  }),
  .id = "Model"
)

# Save the combined tidy summaries to a CSV file
write.csv(combined_tidy, "GAM_Tidy_Summaries.csv", row.names = FALSE)

# Combine cross-validation results into one dataframe with a Response column
combined_cv <- bind_rows(
  map2(cv_results, names(cv_results), ~ {
    if(!is.null(.x)) {
      mutate(.x, Response = .y)
    } else {
      NULL
    }
  }),
  .id = "Model"
)

# Save the combined cross-validation metrics to a CSV file
write.csv(combined_cv, "GAM_Cross_Validation_Metrics.csv", row.names = FALSE)

# Save all fitted GAM models to an RDS file for future use
saveRDS(gam_results, file = "GAM_Models.rds")


# ===== Bootstrapping ====
# Load necessary library
library(boot)

# Define a function to fit the GAM model and extract desired statistics
gam_function <- function(data, indices) {
  # Resample the dataset
  d <- data[indices, ] 
  
  # Fit the GAM model using resampled data
  model <- gam(cube_Median_Respiration_Rate_mg_DO_per_kg_per_H ~ 
                 te(cube_Median_62948_Final_Gravimetric_Moisture_g_per_g, 
                    cube_Median_Weighted_Avg_AI_mod) + 
                 s(site, bs = 're'), 
               data = d, method = "REML")
  
  # Return specific coefficients of interest
  return(coef(model))  # You can change this to return other statistics if desired
}

# Define the data and run bootstrapping
set.seed(123)  # For reproducibility
results <- boot(data = cube_data, statistic = gam_function, R = 1000)  # Adjust R for number of resamples

# Examine bootstrap results
print(results)
plot(results)  # This can visualize the distribution of bootstrapped coefficients

# ===== gaussian vs poisson ====

# Example of fitting models with different families
gaussian_model <- gam(cube_Median_Respiration_Rate_mg_DO_per_kg_per_H ~ 
                        te(cube_Median_62948_Final_Gravimetric_Moisture_g_per_g, 
                           cube_Median_Weighted_Avg_AI_mod) + 
                        s(site, bs = 're'), 
                      data = cube_data, family = gaussian(), method = "REML")

poisson_model <- gam(cube_Median_Respiration_Rate_mg_DO_per_kg_per_H ~ 
                       te(cube_Median_62948_Final_Gravimetric_Moisture_g_per_g, 
                          cube_Median_Weighted_Avg_AI_mod) + 
                       s(site, bs = 're'), 
                     data = cube_data, family = poisson(), method = "REML")

# Compare models
model_comparison <- AIC(gaussian_model, poisson_model)
print(model_comparison)

# ===== Transfromed vs untransfromed =====
formula <- as.formula(paste("cube_Median_Respiration_Rate_mg_DO_per_kg_per_H ~",
                            paste("te(cube_Median_62948_Final_Gravimetric_Moisture_g_per_g, ", predictor, ")"),
                            "+ s(site, bs = 're')"))
model <- gam(formula, data = cube_data, method = "REML")

gam_models[[predictor]] <- model
# Fit models using transformed and untransformed predictors
model_transformed <- gam((cube_Median_Respiration_Rate_mg_DO_per_kg_per_H) ~ 
                           te(cube_Median_62948_Final_Gravimetric_Moisture_g_per_g, 
                              cube_Median_Weighted_Avg_delGcoxPerCmol) + 
                           s(site, bs = 're'), 
                         data = all_data, method = "REML")

model_original <- gam(Median_Respiration_Rate_mg_DO_per_kg_per_H ~ 
                        te(Median_62948_Final_Gravimetric_Moisture_g_per_g, 
                           Median_Weighted_Avg_delGcoxPerCmol) + 
                        s(site, bs = 're'), 
                      data = all_data, method = "REML")

# Compare AIC
model_comparison <- AIC(model_transformed, model_original)
print(model_comparison)
