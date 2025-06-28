# ==== Loading libraries =========
rm(list=ls(all = TRUE))

# Load required packages
library(mgcv)       # For GAM models
library(ggplot2)    # For visualization
library(dplyr)      # For data manipulation
library(viridis)    # For color palettes
library(patchwork)  # For combining plots
library(readr)      # For reading CSV files
library(stringr)    # For string manipulation
library(segmented) # For segmented regression
# ====== Define paths and read in data ======
github = 'C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/ECA_DOM_Thermodynamics/'
data_path = paste0(github, 'Data/')
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


sample_data = read_csv(paste0(github, 'EC_Data_Package/Sample_Data/EC_Sediment_Sample_Data_Summary.csv'), 
                       comment = '#', na = c('N/A', -9999)) %>%
  slice(-(1:11)) %>%
  mutate(across(-Sample_Name:-Material, as.numeric)) %>%
  mutate(Treatment = case_when(
    str_detect(Sample_Name, "-W") ~ "Wet",
    str_detect(Sample_Name, "-D") ~ "Dry",
    TRUE ~ NA_character_)) %>%
  mutate(site = gsub('-W|-D', '', Sample_Name)) %>%
  dplyr::select(site, Treatment, Median_62948_Final_Gravimetric_Moisture_g_per_g, 
         Median_Respiration_Rate_mg_DO_per_kg_per_H, Median_Extractable_NPOC_mg_per_kg)

# Read the main data file
data = read.csv(paste0(data_path, 'Medians_of Median_molecular_properties_per_site_and_treatment_unique_formulas.csv'))
row.names(data) = paste0(data$site, '_', data$Treatment)

# Set up data
data <- data %>%
  filter(!(site %in% c("EC_023", "EC_011", "EC_012", "EC_052", "EC_053", "EC_057"))) %>%
  merge(sample_data, by = c('site', 'Treatment'))

# ===== Cube root transform ======

cube_root <- function(x) sign(x) * (abs(x))^(1/3)

cube_data = data %>% 
  mutate(across(where(is.numeric), cube_root)) %>% # cube root transform data
  rename_with(where(is.numeric), .fn = ~ paste0("cube_", .x))

dom_data = cube_data %>%
  dplyr::select(site, Treatment, resp_cube = cube_Median_Respiration_Rate_mg_DO_per_kg_per_H, gibbs = cube_Median_delGcoxPerCmol,
                moisture = cube_Median_62948_Final_Gravimetric_Moisture_g_per_g,
                lambda = cube_Median_Lambda)

dom_data$site = as.factor(dom_data$site)
dom_data$Treatment = as.factor(dom_data$Treatment)

dom_data$Treatment = as.factor(dom_data$Treatment)
# ======= Linear regressions with treatment ====
library(broom)  # For extracting model stats

# Split data by treatment and fit separate linear models
wet_data <- dom_data %>% filter(Treatment == "Wet")
dry_data <- dom_data %>% filter(Treatment == "Dry")

# Fit separate linear models for each treatment 
# Gibbs models
gibbs_wet_model <- lm(resp_cube ~ gibbs, data = wet_data)
gibbs_dry_model <- lm(resp_cube ~ gibbs, data = dry_data)

# Lambda models  
lambda_wet_model <- lm(resp_cube ~ lambda, data = wet_data)
lambda_dry_model <- lm(resp_cube ~ lambda, data = dry_data)
# Extract statistics from model summaries
gibbs_wet_r2 <- round(summary(gibbs_wet_model)$r.squared, 2)
gibbs_wet_p <- round(summary(gibbs_wet_model)$coefficients[2,4], 2)
gibbs_dry_r2 <- round(summary(gibbs_dry_model)$r.squared, 2)
gibbs_dry_p <- round(summary(gibbs_dry_model)$coefficients[2,4], 2)

lambda_wet_r2 <- round(summary(lambda_wet_model)$r.squared, 2)
lambda_wet_p <- round(summary(lambda_wet_model)$coefficients[2,4], 2)
lambda_dry_r2 <- round(summary(lambda_dry_model)$r.squared, 2)
lambda_dry_p <- round(summary(lambda_dry_model)$coefficients[2,4], 2)

# Function to format p-values with significance stars
format_p <- function(p_val) {
  if (p_val < 0.001) return(paste0(format(p_val, scientific = TRUE, digits = 2), "***"))
  else if (p_val < 0.01) return(paste0(round(p_val, 3), ""))
  else if (p_val < 0.05) return(paste0(round(p_val, 3), ""))
  else if (p_val < 0.1) return(paste0(round(p_val, 3), ""))
  else return(round(p_val, 3))
}

# Create dynamic labels
gibbs_wet_label <- paste0("R² = ", gibbs_wet_r2, ", p = ", format_p(gibbs_wet_p))
gibbs_dry_label <- paste0("R² = ", gibbs_dry_r2, ", p = ", format_p(gibbs_dry_p))

lambda_wet_label <- paste0("R² = ", lambda_wet_r2, ", p = ", format_p(lambda_wet_p))
lambda_dry_label <- paste0("R² = ", lambda_dry_r2, ", p = ", format_p(lambda_dry_p))

# Updated plots with dynamic labels
plot1 <- ggplot(dom_data, aes(x = gibbs, y = resp_cube, color = Treatment)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  # Dynamic colored text for each treatment
  annotate("text", x = Inf, y = -Inf, 
           label = gibbs_wet_label,
           hjust = 1.1, vjust = -0.2, size = 5, color = "#4682B4", fontface = "bold") +
  annotate("text", x = Inf, y = -Inf, 
           label = gibbs_dry_label,
           hjust = 1.1, vjust = -1.5, size = 5, color = "#8B4513", fontface = "bold") +
  scale_color_manual(values = c("Dry" = "#8B4513", "Wet" = "#4682B4")) +
  scale_fill_manual(values = c("Dry" = "#8B4513", "Wet" = "#4682B4")) +
  labs(title = "A",
       x = expression(paste("Median ", Delta, G[cox], " (kJ Cmol"^{-1},")"^(1/3))),
       y = expression(paste("Median ",~O[2],' consumption rate ', ~(mg~O[2]~kg^{-1}~h^{-1})^{1/3}))) +
  theme_bw()+
  theme(
    # Increase all font sizes and ensure axis text is black
    axis.text = element_text(size = 14, color = "black"),    # Black axis numbers
    axis.title = element_text(size = 16, color = "black"),   # Black axis titles  
    plot.title = element_text(size = 18, color = "black"),   # Black plot title
    aspect.ratio = 1
  )

plot2 <- ggplot(dom_data, aes(x = lambda, y = resp_cube, color = Treatment)) +
  geom_point(size = 3) +
  #geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  # Dynamic colored text for each treatment
  annotate("text", x = -Inf, y = -Inf,
           label = lambda_wet_label,
           hjust = -0.1, vjust = -0.1, size = 5, color = "#4682B4", fontface = "bold") +
  annotate("text", x = -Inf, y = -Inf,
           label = lambda_dry_label,
           hjust = -0.1, vjust = -1.5, size = 5, color = "#8B4513", fontface = "bold") +
  scale_color_manual(values = c("Dry" = "#8B4513", "Wet" = "#4682B4")) +
  scale_fill_manual(values = c("Dry" = "#8B4513", "Wet" = "#4682B4")) +
  labs(title = "B",
       x = expression(paste("Median ", lambda^(1/3))),
       y = expression(paste("Median ",~O[2],' consumption rate ', ~(mg~O[2]~kg^{-1}~h^{-1})^{1/3}))) +
  theme_bw()+  theme(
    # Increase all font sizes and ensure axis text is black
    axis.text = element_text(size = 14, color = "black"),    # Black axis numbers
    axis.title = element_text(size = 16, color = "black"),   # Black axis titles  
    plot.title = element_text(size = 18, color = "black"),   # Black plot title
    aspect.ratio = 1,                                        # Square aspect ratio
  )

# Combine plots
combined_plot <- plot1 + plot2 +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        legend.text = element_text(size = 14),      # Style the legend
        legend.title = element_text(size = 16))

print(combined_plot)

# Print summaries to see significance
summary(gibbs_wet_model)
summary(gibbs_dry_model)
summary(lambda_wet_model) 
summary(lambda_dry_model)

ggsave("Figures/Figure3_Thermodynamics_vs_Respiration.png", 
       combined_plot, 
       width = 12, height = 6, dpi = 300)

# Save as PDF
ggsave("Figures/Figure3_Thermodynamics_vs_Respiration.pdf", 
       combined_plot, 
       width = 12, height = 6)


# ====== GAM models =======
# Define the refined list of predictors
cube_predictors <- c("gibbs", 'lambda')

# Fit null model
cube_null_model <- gam(resp_cube ~ 1 + s(site, bs = "re"), 
                       data = dom_data, method = "REML")

# Initialize list of models with the null model
cube_gam_models <- list(cube_null_model)
names(cube_gam_models) <- "Null_cube"

# Fit models with individual predictors (no interaction with moisture)
for (predictor in cube_predictors) {
  # Model with thermodynamic property alone
  formula_solo <- as.formula(paste("resp_cube ~ s(", predictor, ")", 
                                   "+ s(site, bs = 're')"))
  model_solo <- gam(formula_solo, data = dom_data, method = "REML")
  cube_gam_models[[paste0(predictor, "_solo")]] <- model_solo
  
  # Model with moisture alone
  if (!exists("moisture_solo_model")) {
    formula_moisture <- as.formula("resp_cube ~ s(moisture) + s(site, bs = 're')")
    moisture_solo_model <- gam(formula_moisture, data = dom_data, method = "REML")
    cube_gam_models[["moisture_solo"]] <- moisture_solo_model
  }
  
  # Model with additive effects (no interaction)
  formula_additive <- as.formula(paste("resp_cube ~ s(moisture) + s(", 
                                       predictor, ") + s(site, bs = 're')"))
  model_additive <- gam(formula_additive, data = dom_data, method = "REML")
  cube_gam_models[[paste0("moisture_plus_", predictor)]] <- model_additive
  
  # Original interaction model
  formula_interaction <- as.formula(paste("resp_cube ~ te(moisture, ", 
                                          predictor, ") + s(site, bs = 're')"))
  model_interaction <- gam(formula_interaction, data = dom_data, method = "REML")
  cube_gam_models[[paste0("moisture_", predictor, "_interaction")]] <- model_interaction
}

# Create a function to extract significant terms from model summary
get_significant_terms <- function(model) {
  summary_info <- summary(model)
  if (is.null(summary_info$s.table)) {
    return("No smooth terms")
  }
  
  # Extract terms with p < 0.05
  smooth_terms <- summary_info$s.table
  sig_terms <- rownames(smooth_terms)[smooth_terms[,"p-value"] < 0.05]
  
  if (length(sig_terms) == 0) {
    return("No significant terms")
  }
  
  # Format significance levels
  sig_levels <- sapply(smooth_terms[sig_terms, "p-value"], function(p) {
    if (p < 0.001) return("***")
    else if (p < 0.01) return("**")
    else if (p < 0.05) return("*")
    else return("")
  })
  
  # Combine terms with significance
  return(paste(sig_terms, sig_levels, collapse = ", "))
}

# Create comparison dataframe with safer approach for formula extraction
model_comparison <- data.frame(
  Model = names(cube_gam_models),
  AIC = sapply(cube_gam_models, AIC),
  R_squared = sapply(cube_gam_models, function(m) summary(m)$r.sq),
  #Adj_R_squared = sapply(cube_gam_models, function(m) summary(m)$r.sq.adj),
  Dev_explained = sapply(cube_gam_models, function(m) summary(m)$dev.expl),
  stringsAsFactors = FALSE
)

# Add formula column separately to avoid issues
model_comparison$Formula <- sapply(cube_gam_models, function(m) {
  form <- deparse(formula(m))
  form <- gsub("\\s+", " ", form) # clean up white space
  return(form)
})

# Add significant terms column
model_comparison$Significant_terms <- sapply(cube_gam_models, get_significant_terms)

# Sort by AIC (lower is better)
model_comparison <- model_comparison[order(model_comparison$AIC), ]

# Print comprehensive comparison
print(model_comparison)

# Compare AIC
model_comparison_file <- "06_GAM_Results/Comparison_Results/AIC_comparison_multiple_models.csv"
write.csv(model_comparison, model_comparison_file)

# ==== Export Summary Files =====
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

# Generate plots for cube root transformed models
for (model_name in names(cube_gam_models)) {
  generate_plots(cube_gam_models[[model_name]], model_name, "transformed")
}

# Function to perform Breusch-Pagan test and export results
library(lmtest)
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


# ==== Second derivative for threshold detection ====

moisture_seq <- seq(min(dom_data$moisture), max(dom_data$moisture), length.out = 500)
reference_site <- levels(factor(dom_data$site))[1]
reference_gibbs <- median(dom_data$gibbs)

# Initial GAM model
gam_model <- gam(resp_cube ~ te(moisture, gibbs) + s(site, bs = "re"), 
                 data = dom_data, method = "REML")
summary(gam_model)

# Create a dataset for prediction along the moisture gradient
pred_data_gibbs <- data.frame(
  moisture = moisture_seq,
  gibbs = reference_gibbs,
  site = reference_site
)

# Generate predicted values
pred_data_gibbs$predicted <- predict(gam_model, newdata = pred_data_gibbs, se.fit = TRUE)
pred_data_gibbs$fit <- pred_data_gibbs$predicted$fit
pred_data_gibbs$se <- pred_data_gibbs$predicted$se.fit

# Calculate derivatives using smooth.spline (more robust than diff)
# Fit a smooth spline to the predictions
spline_fit <- smooth.spline(pred_data_gibbs$moisture, pred_data_gibbs$fit, cv = TRUE)

# First derivative - rate of change
first_deriv <- predict(spline_fit, pred_data_gibbs$moisture, deriv = 1)
pred_data_gibbs$first_derivative <- first_deriv$y

# Second derivative - change in rate of change
second_deriv <- predict(spline_fit, pred_data_gibbs$moisture, deriv = 2)
pred_data_gibbs$second_derivative <- second_deriv$y

#Find the threshold using second derivative maxima/minima
# Use absolute value to find points where curvature changes most dramatically
pred_data_gibbs$abs_second_deriv <- abs(pred_data_gibbs$second_derivative)

# Find local maxima in second derivative (points of maximum curvature)
# This is more robust than just taking the global maximum
local_max <- which(diff(sign(diff(pred_data_gibbs$abs_second_deriv))) == -2) + 1

# If no clear local maxima, use global maximum
if(length(local_max) == 0) {
  threshold_idx <- which.max(pred_data_gibbs$abs_second_deriv)
} else {
  # Get the strongest local maximum
  threshold_idx <- local_max[which.max(pred_data_gibbs$abs_second_deriv[local_max])]
}

moisture_threshold <- pred_data_gibbs$moisture[threshold_idx]

cat("Moisture threshold from second derivative method:", round(moisture_threshold, 3), "\n")

dom_data$moisture_regime <- ifelse(dom_data$moisture < moisture_threshold, 
                                   "Below Threshold", "Above Threshold")
# ==== Bootstrap for Confidence Intervals ====
n_boot <- 500
boot_thresholds <- numeric(n_boot)

cat("Running bootstrap for second derivative threshold...\n")
for(i in 1:n_boot) {
  # Bootstrap sample
  boot_idx <- sample(1:nrow(dom_data), replace = TRUE)
  boot_data <- dom_data[boot_idx, ]
  
  tryCatch({
    # Fit GAM to bootstrap sample
    boot_gam <- gam(resp_cube ~ te(moisture, gibbs) + s(site, bs="re"),
                    data = boot_data, method = "REML")
    
    # Generate predictions
    boot_preds <- predict(boot_gam, newdata = pred_data_gibbs)
    
    # Smooth the predictions
    boot_spline <- smooth.spline(pred_data_gibbs$moisture, boot_preds, cv = TRUE)
    
    # Calculate second derivative
    boot_deriv2 <- predict(boot_spline, pred_data_gibbs$moisture, deriv = 2)
    
    # Find threshold
    boot_abs_deriv2 <- abs(boot_deriv2$y)
    
    # Find local maxima in second derivative
    boot_local_max <- which(diff(sign(diff(boot_abs_deriv2))) == -2) + 1
    
    # If no clear local maxima, use global maximum
    if(length(boot_local_max) == 0) {
      boot_threshold_idx <- which.max(boot_abs_deriv2)
    } else {
      # Get the strongest local maximum
      boot_threshold_idx <- boot_local_max[which.max(boot_abs_deriv2[boot_local_max])]
    }
    
    boot_thresholds[i] <- pred_data_gibbs$moisture[boot_threshold_idx]
    
  }, error = function(e) {
    boot_thresholds[i] <- NA
  })
  
  if(i %% 50 == 0) cat(i, "bootstrap iterations completed\n")
}

# Calculate confidence intervals
valid_boots <- boot_thresholds[!is.na(boot_thresholds)]
ci_lower <- quantile(valid_boots, 0.025)
ci_upper <- quantile(valid_boots, 0.975)

cat("\nSecond Derivative Threshold:", round(moisture_threshold, 3),
    "\n95% CI:", round(ci_lower, 3), "-", round(ci_upper, 3), "\n")

# ----- Threshold Visualization for Gibbs -----

# ---- Plot A: Moisture Threshold Plot with GAM Fit ----
p1 <- ggplot() +
  # Original data points
  geom_point(data = dom_data, aes(x = moisture, y = resp_cube, color = site),
             alpha = 0.5, size = 3) +
  # GAM fit with confidence band
  geom_line(data = pred_data_gibbs, aes(x = moisture, y = fit),
            color = "blue", size = 1.2) +
  # Show confidence intervals for the fit
  geom_ribbon(data = pred_data_gibbs,
              aes(x = moisture, ymin = fit - 1.96*se, ymax = fit + 1.96*se),
              fill = "blue", alpha = 0.1) +
  # Threshold line with confidence interval
  geom_vline(xintercept = moisture_threshold, linetype = "dashed",
             color = "red", size = 1) +
  geom_rect(aes(xmin = ci_lower, xmax = ci_upper,
                ymin = -Inf, ymax = Inf),
            alpha = 0.1, fill = "red") +
  # Annotations
  annotate("text", x = moisture_threshold + 0.1,
           y = 0.1 * diff(range(dom_data$resp_cube)),
           label = paste0("Threshold: ", round(moisture_threshold, 2)),
           color = "red", hjust = 0) +
  annotate("text", x = moisture_threshold + 0.1,
           y = 0.01 * diff(range(dom_data$resp_cube)),
           label = paste0("(95% CI: ",
                          round(ci_lower, 2), "-", round(ci_upper, 2), ")"),
           color = "red", hjust = 0) +
  # Aesthetics
  scale_color_viridis_d(guide = "none") +  # Colorful sites without legend
  labs(
    title = "A",
    x = expression(Median~Gravimetric~Moisture~Content~(g/g)^{1/3}),
    y = expression(paste("Median",~O[2],' consumption rate', ~(mg~o[2]~kg^{-1}~h^{-1})^{1/3}))
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 18, color = "black"),
    aspect.ratio = 1
  )


# ---- Plot B: DOM Thermodynamics Split by Moisture Regime ----
library(ggplot2)
library(lmerTest)

# Get your model results (simpler than before)
below_mixed <- lmer(resp_cube ~ gibbs + (1|site),
                    data = dom_data[dom_data$moisture < moisture_threshold, ])
above_mixed <- lmer(resp_cube ~ gibbs + (1|site),
                    data = dom_data[dom_data$moisture >= moisture_threshold, ])

# Extract stats
below_r2 <- r.squaredGLMM(below_mixed)[1]  # marginal R²
above_r2 <- r.squaredGLMM(above_mixed)[1]
below_beta <- fixef(below_mixed)[2]
above_beta <- fixef(above_mixed)[2]
below_p <- summary(below_mixed)$coefficients[2,5]
above_p <- summary(above_mixed)$coefficients[2,5]

# Create significance stars
below_signif <- case_when(below_p < 0.001 ~ "***",
                          below_p < 0.01 ~ "**", 
                          below_p < 0.05 ~ "*",
                          TRUE ~ " ")
above_signif <- case_when(above_p < 0.001 ~ "***",
                          above_p < 0.01 ~ "**",
                          above_p < 0.05 ~ "*", 
                          TRUE ~ " ")

# Simple plot using built-in geom_smooth
p2 <- ggplot(dom_data, aes(x = gibbs, y = resp_cube, color = moisture_regime)) +
  # Data points
  geom_point(alpha = 1, size = 3) +
  # Automatic regression lines with confidence intervals
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) +
  # Statistics annotations
  annotate("text", x = min(dom_data$gibbs) + 0.002,
           y = 0.1 * diff(range(dom_data$resp_cube)),
           label = paste0("R² = ", round(below_r2, 2),
                          ", β = ", round(below_beta, 2), below_signif),
           color = "darkorange", hjust = 0, fontface = "bold") +
  annotate("text", x = min(dom_data$gibbs) + 0.002,
           y = 0.2 * diff(range(dom_data$resp_cube)),
           label = paste0("R² = ", round(above_r2, 2),
                          ", β = ", round(above_beta, 2), above_signif),
           color = "lightblue", hjust = 0, fontface = "bold") +
  # Styling
  scale_color_manual(values = c("Below Threshold" = "darkorange",
                                "Above Threshold" = "lightblue")) +
  labs(
    title = "B",
    x = expression(paste(Median~Delta, G[cox], " (kJ Cmol"^{-1},")"^(1/3))),
    y = expression(paste("Median",~O[2],' consumption rate', ~(mg~o[2]~kg^{-1}~h^{-1})^(1/3))),
    color = " "
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 18, color = "black"),
    aspect.ratio = 1
  )
# Combine plots
combined_plot <- p1 + p2 +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        legend.text = element_text(size = 14),      # Style the legend
        legend.title = element_text(size = 16))

print(combined_plot)


ggsave("Figures/Figure4_Thermodynamics_vs_Respiration.png", 
       combined_plot, 
       width = 12, height = 6, dpi = 300)

# Save as PDF
ggsave("Figures/Figure4_Thermodynamics_vs_Respiration.pdf", 
       combined_plot, 
       width = 12, height = 6)
