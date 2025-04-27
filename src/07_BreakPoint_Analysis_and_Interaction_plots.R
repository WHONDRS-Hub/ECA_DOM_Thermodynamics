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
  theme(aspect.ratio = 1)

# ---- Plot B: DOM Thermodynamics Split by Moisture Regime ----

library(lme4)

# Mixed models with site as random effect
below_mixed <- lmer(resp_cube ~ gibbs + (1|site), 
                    data = dom_data[dom_data$moisture < moisture_threshold, ])
above_mixed <- lmer(resp_cube ~ gibbs + (1|site), 
                    data = dom_data[dom_data$moisture >= moisture_threshold, ])

# Summarize models
below_summary <- summary(below_mixed)
above_summary <- summary(above_mixed)

# Extract fixed-effect coefficients (these are your adjusted betas)
below_beta <- fixef(below_mixed)[2]
above_beta <- fixef(above_mixed)[2]

# Get p-values (using lmerTest approach)
library(lmerTest)
below_mixed_t <- lmer(resp_cube ~ gibbs + (1|site), 
                      data = dom_data[dom_data$moisture < moisture_threshold, ])
above_mixed_t <- lmer(resp_cube ~ gibbs + (1|site), 
                      data = dom_data[dom_data$moisture >= moisture_threshold, ])

below_p <- summary(below_mixed_t)$coefficients[2,5]
above_p <- summary(above_mixed_t)$coefficients[2,5]

# Calculate conditional R-squared (variance explained by fixed + random effects)
library(MuMIn)
below_r2 <- r.squaredGLMM(below_mixed)[1] # Use [2] for conditional R²
above_r2 <- r.squaredGLMM(above_mixed)[1] # Use [2] for conditional R²

# Print comparison with adjusted statistics
cat("DOM-Respiration Relationship by Moisture Regime (Adjusted for Site Effects):\n")
cat("Below threshold (", round(moisture_threshold, 2), "):\n", sep="")
cat("  Adjusted Coefficient:", round(below_beta, 3),
    "\n  Marginal R-squared:", round(below_r2, 3),
    "\n  p-value:", format.pval(below_p, digits = 3), "\n\n")
cat("Above threshold (", round(moisture_threshold, 2), "):\n", sep="")
cat("  Adjusted Coefficient:", round(above_beta, 3),
    "\n  Marginal R-squared:", round(above_r2, 3),
    "\n  p-value:", format.pval(above_p, digits = 3), "\n\n")

# Create prediction frames from mixed models
# Since site is a random effect, we'll use the average site effect
# Create sequence for predictions
gibbs_range <- range(dom_data$gibbs, na.rm = TRUE)
gibbs_seq <- seq(gibbs_range[1], gibbs_range[2], length.out = 100)

# Now create new_data for predictions
new_data <- data.frame(gibbs = gibbs_seq)

# For below threshold
below_pred_fixed <- fixef(below_mixed)[1] + fixef(below_mixed)[2] * gibbs_seq
# For confidence intervals, we'll use bootMer to account for uncertainty
library(boot)
bootfun <- function(model) {
  predict(model, newdata = new_data, re.form = NA)
}
below_boot <- bootMer(below_mixed, bootfun, nsim = 500)
below_ci <- t(apply(below_boot$t, 2, quantile, c(0.025, 0.975)))

# Repeat for above threshold
above_pred_fixed <- fixef(above_mixed)[1] + fixef(above_mixed)[2] * gibbs_seq
above_boot <- bootMer(above_mixed, bootfun, nsim = 500)
above_ci <- t(apply(above_boot$t, 2, quantile, c(0.025, 0.975)))

# Create dataframes for plotting
below_frame <- data.frame(
  gibbs = gibbs_seq,
  predicted = below_pred_fixed,
  lower = below_ci[,1],
  upper = below_ci[,2],
  regime = "Below Threshold"
)

above_frame <- data.frame(
  gibbs = gibbs_seq,
  predicted = above_pred_fixed,
  lower = above_ci[,1],
  upper = above_ci[,2],
  regime = "Above Threshold"
)

pred_combined <- rbind(below_frame, above_frame)


# Create significance indicators
below_signif <- ifelse(below_p < 0.05, 
                       ifelse(below_p < 0.01, 
                              ifelse(below_p < 0.001, "***", "**"), "*"), "ns")

above_signif <- ifelse(above_p < 0.05, 
                       ifelse(above_p < 0.01, 
                              ifelse(above_p < 0.001, "***", "**"), "*"), "ns")

# Plot with adjusted statistics
p2 <- ggplot() +
  # Data points
  geom_point(data = dom_data,
             aes(x = gibbs, y = resp_cube, color = moisture_regime),
             alpha = 1, size = 3) +
  # Model fits with confidence intervals
  geom_line(data = pred_combined,
            aes(x = gibbs, y = predicted, color = regime),
            size = 1) +
  geom_ribbon(data = pred_combined,
              aes(x = gibbs, ymin = lower, ymax = upper, fill = regime),
              alpha = 0.5) +
  # Annotations with site-adjusted statistics
  annotate("text", x = min(dom_data$gibbs) + 0.002,
           y =  0.1 * diff(range(dom_data$resp_cube)),
           label = paste0("Below: R² = ", round(below_r2, 2),
                          ", β = ", round(below_beta, 2), below_signif),
           color = "darkorange", hjust = 0, fontface = "bold") +
  annotate("text", x = min(dom_data$gibbs) + 0.002,
           y =  0.2 * diff(range(dom_data$resp_cube)),
           label = paste0("Above: R² = ", round(above_r2, 2),
                          ", β = ", round(above_beta, 2), above_signif),
           color = "lightblue", hjust = 0, fontface = "bold") +
  # Theme and labels
  scale_color_manual(values = c("Below Threshold" = "darkorange",
                                "Above Threshold" = "lightblue")) +
  scale_fill_manual(values = c("Below Threshold" = "darkorange",
                               "Above Threshold" = "lightblue")) +
  labs(
    title = "B",
    x = expression(paste(Median~Delta, G[cox], " (kJ Cmol"^{-1},")"^(1/3))),
    y = expression(paste("Median",~O[2],' consumption rate', ~(mg~o[2]~kg^{-1}~h^{-1})^(1/3))),
    color = " ",
    fill = " "
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'top', aspect.ratio = 1)

# ---- Plot C: Interaction Heatmap with Contours ----
# Create prediction grid for interaction plot
moisture_grid <- seq(min(dom_data$moisture), max(dom_data$moisture), length.out = 50)
gibbs_grid <- seq(min(dom_data$gibbs), max(dom_data$gibbs), length.out = 50)
grid_data <- expand.grid(moisture = moisture_grid, gibbs = gibbs_grid, site = reference_site)

# Generate predictions from GAM model
grid_data$predicted <- predict(gam_model, newdata = grid_data)

# Calculate contour breaks for the plot
pred_range <- range(grid_data$predicted)
contour_breaks <- seq(pred_range[1], pred_range[2], length.out = 10)

# Create interaction heatmap
p3 <- ggplot(grid_data, aes(x = moisture, y = gibbs)) +
  # Create filled contours with custom colormap
  geom_contour_filled(aes(z = predicted),
                      breaks = contour_breaks,
                      alpha = 0.8) +
  # Add contour lines
  geom_contour(aes(z = predicted),
               breaks = contour_breaks,
               color = "black",
               alpha = 0.5,
               size = 0.3) +
  # Add threshold line
  geom_vline(xintercept = moisture_threshold,
             linetype = "dashed",
             color = "black",
             size = 1) +
  # Add data points
  geom_point(data = dom_data,
             aes(size = abs(resp_cube)),
             shape = 21,
             color = "black",
             fill = "white",
             alpha = 0.7) +
  # Use a color scale similar to 'topo' in vis.gam
  scale_fill_viridis_d(option = "plasma",
                       name = expression(paste("Predicted ",~O[2]," Consumption (Cube Root)"))) +
  scale_size_continuous(name = expression(paste("Observed",~O[2], " Consumption (Cube Root)")),
                        range = c(1, 5)) +
  # Properly label the axes with expressions
  labs(
    title = "A",
    x = expression(paste("Median Gravimetric Moisture ", (g/g)^(1/3))),
    y = expression(paste(Median~Delta, G[cox], " (kJ Cmol"^{-1},")"^(1/3))),
  ) +
  # Theme settings to match vis.gam aesthetics
  theme_bw(base_size = 12) +
  theme(
    aspect.ratio = 1,
    panel.grid = element_blank(),
    legend.position = "right"
  )


# ===== SI plots =====
# ---- SI Plot 1: First and Second Derivatives ----
# Plot First Derivative
si_p1a <- ggplot(pred_data_gibbs, aes(x = moisture, y = first_derivative)) +
  geom_line(size = 1, color = "darkgreen") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = moisture_threshold, color = "red", 
             linetype = "dashed", size = 1) +
  labs(
    title = "A",
    x = expression(Median~Gravimetric~Moisture~Content~(g/g)^{1/3}),
    y = "First Derivative (Rate of Change)"
  ) +
  theme_bw(base_size = 12) + 
  theme(aspect.ratio = 1)

# Plot Second Derivative with threshold
si_p1b <- ggplot(pred_data_gibbs, aes(x = moisture, y = second_derivative)) +
  geom_line(size = 1, color = "purple") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = moisture_threshold, color = "red", 
             linetype = "dashed", size = 1) +
  geom_point(aes(x = moisture[threshold_idx], y = second_derivative[threshold_idx]), 
             color = "red", size = 3) +
  labs(
    title = "B",
    x = expression(Median~Gravimetric~Moisture~Content~(g/g)^{1/3}),
    y = "Second Derivative (Curvature)"
  ) +
  theme_bw(base_size = 12) + 
  theme(aspect.ratio = 1)

# ---- SI Plot 2: Bootstrap Distribution ----
si_p2 <- ggplot(data.frame(threshold = valid_boots), aes(x = threshold)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = moisture_threshold, color = "red", size = 1) +
  geom_vline(xintercept = ci_lower, color = "darkred", linetype = "dashed") +
  geom_vline(xintercept = ci_upper, color = "darkred", linetype = "dashed") +
  annotate("text", x = moisture_threshold - 0.2, 
           y = max(hist(valid_boots, breaks = 30, plot = FALSE)$counts) * 0.9,
           label = paste0("Threshold: ", round(moisture_threshold, 2)), 
           color = "red", hjust = 0) +
  annotate("text", x =  moisture_threshold - 0.2, 
           y = max(hist(valid_boots, breaks = 30, plot = FALSE)$counts) * 0.8,
           label = paste0("95% CI: [", round(ci_lower, 2), ", ", round(ci_upper, 2), "]"), 
           color = "red", hjust = 0) +
  labs(
    title = "C",
    x = expression(Moisture~Threshold~(g/g)^{1/3}),
    y = "Frequency",
    subtitle = paste0(length(valid_boots), " valid bootstrap samples out of ", n_boot)
  ) +
  theme_bw(base_size = 12) +
  theme(aspect.ratio = 1)

# ---- SI Plot 3: Model Validation ----
# Add predictions to dom_data if they're not already there
if(!"predicted" %in% colnames(dom_data)) {
  dom_data$predicted <- predict(gam_model, newdata = dom_data)
}

si_p3 <- ggplot(dom_data, aes(x = predicted, y = resp_cube, color = moisture_regime)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  stat_smooth(method = "lm", se = FALSE, linetype = "solid") +
  scale_color_manual(values = c("Below Threshold" = "darkorange", 
                                "Above Threshold" = "lightblue"), name = "Moisture Regime") +
  labs(
    title = "D",
    x = expression(paste("Predicted",~O[2],' consumption rate', ~(mg~o[2]~kg^{-1}~h^{-1})^{1/3})),
    y = expression(paste("Observed",~O[2],' consumption rate', ~(mg~o[2]~kg^{-1}~h^{-1})^{1/3}))
  ) +
  theme_bw(base_size = 12) +
  theme(aspect.ratio = 1, legend.position = "none")

# ---- SI Plot 4: Statistical Testing ----
dom_data$moisture_below <- pmin(dom_data$moisture, moisture_threshold)
dom_data$moisture_above <- pmax(0, dom_data$moisture - moisture_threshold)

# Model with single tensor product vs. threshold model with separate tensor products
model_single <- gam(resp_cube ~ te(moisture, gibbs) + s(site, bs="re"),
                    data = dom_data, method = "REML")

model_threshold <- gam(resp_cube ~ te(moisture_below, gibbs) + te(moisture_above, gibbs) + s(site, bs="re"),
                       data = dom_data, method = "REML")

# Compare models
anova_result <- anova(model_single, model_threshold, test = "F")

# Create text for statistical test
stat_text <- paste0(
  "Statistical Test for Threshold Significance\n\n",
  "Method: F-test comparing nested GAM models with tensor product smooths\n",
  "F-statistic: ", round(anova_result$F[2], 3), "\n",
  "p-value: ", format.pval(anova_result$`p-value`[2], digits = 3), "\n",
  "Significant at α=0.05: ", ifelse(anova_result$`p-value`[2] < 0.05, "Yes", "No"), "\n\n",
  "A significant test indicates that modeling separate\n",
  "moisture-DOM thermodynamics interactions before and after\n",
  "the threshold is statistically preferred.\n\n",
  "Deviance explained (single model): ", 
  round(summary(model_single)$dev.expl * 100, 1), "%\n",
  "Deviance explained (threshold model): ", 
  round(summary(model_threshold)$dev.expl * 100, 1), "%"
)

si_p4 <- ggplot() + 
  annotate("text", x = 0.5, y = 0.5, label = stat_text, size = 5) +
  theme_void() +
  labs(title = "E") +
  theme(plot.title = element_text(hjust = 0, size = 14, face = "bold"))
# ---- SI Plot 5: Residual Diagnostics ----
# Calculate residuals
dom_data$residuals <- residuals(gam_model)

# Create residuals plot by moisture with threshold
si_p5 <- ggplot(dom_data, aes(x = moisture, y = residuals, color = moisture_regime)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = moisture_threshold, color = "red", linetype = "dashed") +
  geom_smooth(method = "loess", se = TRUE, alpha = 0.2) +
  scale_color_manual(values = c("Below Threshold" = "darkorange", 
                                "Above Threshold" = "lightblue"), name = "Moisture Regime") +
  labs(
    title = "F",
    x = expression(Median~Gravimetric~Moisture~Content~(g/g)^{1/3}),
    y = "Model Residuals"
  ) +
  theme_bw(base_size = 12) +
  theme(aspect.ratio = 1, legend.position = "none")

# ---- SI Plot 6: Methodology Diagram ----
# Create a simplified diagram explaining the second derivative method
method_steps <- data.frame(
  x = 1:4,
  y = rep(1, 4),
  label = c("1. Generate predictions\nacross moisture gradient",
            "2. Calculate first and\nsecond derivatives\nusing smooth splines",
            "3. Identify threshold at\nmaximum curvature\n(second derivative)",
            "4. Bootstrap for\nconfidence intervals")
)

arrows_df <- data.frame(
  x = 1:3,
  xend = 2:4,
  y = rep(1, 3),
  yend = rep(1, 3)
)

si_p6 <- ggplot() +
  geom_point(data = method_steps, aes(x = x, y = y), size = 15, color = "steelblue", alpha = 0.7) +
  geom_text(data = method_steps, aes(x = x, y = y, label = label), size = 3.5, color = "white") +
  geom_segment(data = arrows_df, aes(x = x + 0.2, y = y, xend = xend - 0.2, yend = yend),
               arrow = arrow(length = unit(0.3, "cm")), size = 1) +
  labs(title = "G. Second Derivative Method for Threshold Detection") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.margin = margin(20, 20, 20, 20))

# Arrange SI plots
si_p1 <- si_p1a + si_p1b  + si_p2 + plot_layout(ncol = 3)

si_row2 <- si_p3 + si_p5 + plot_layout(ncol = 2)

si_row3 <-   si_p4 

# Save the SI plots
ggsave(paste0(figure_path,"SI_Figure1_Gibbs.png"), si_p1, width = 12, height = 16, dpi = 300)
ggsave(paste0(figure_path,"SI_Figure2_Gibbs.png"), si_row2, width = 12, height = 16, dpi = 300)
ggsave(paste0(figure_path,"SI_Figure3_Gibbs.png"), si_row3, width = 12, height = 16, dpi = 300)


ggsave(paste0(figure_path,"SI_Figure1_Gibbs.pdf"), si_p1, width = 12, height = 16)
ggsave(paste0(figure_path,"SI_Figure2_Gibbs.pdf"), si_row2, width = 12, height = 16)
ggsave(paste0(figure_path,"SI_Figure3_Gibbs.pdf"), si_row3, width = 12, height = 16)

# ====== NOW FOR LAMBDA ====
# ==== Second derivative for threshold detection ====

moisture_seq <- seq(min(dom_data$moisture), max(dom_data$moisture), length.out = 500)
reference_site <- levels(factor(dom_data$site))[1]
reference_lambda <- median(dom_data$lambda)

# Initial GAM model
gam_model <- gam(resp_cube ~ te(moisture, lambda) + s(site, bs = "re"), 
                 data = dom_data, method = "REML")
summary(gam_model)

# Create a dataset for prediction along the moisture gradient
pred_data_lambda <- data.frame(
  moisture = moisture_seq,
  lambda = reference_lambda,
  site = reference_site
)

# Generate predicted values
pred_data_lambda$predicted <- predict(gam_model, newdata = pred_data_lambda, se.fit = TRUE)
pred_data_lambda$fit <- pred_data_lambda$predicted$fit
pred_data_lambda$se <- pred_data_lambda$predicted$se.fit

# Calculate derivatives using smooth.spline (more robust than diff)
# Fit a smooth spline to the predictions
spline_fit <- smooth.spline(pred_data_lambda$moisture, pred_data_lambda$fit, cv = TRUE)

# First derivative - rate of change
first_deriv <- predict(spline_fit, pred_data_lambda$moisture, deriv = 1)
pred_data_lambda$first_derivative <- first_deriv$y

# Second derivative - change in rate of change
second_deriv <- predict(spline_fit, pred_data_lambda$moisture, deriv = 2)
pred_data_lambda$second_derivative <- second_deriv$y

#Find the threshold using second derivative maxima/minima
# Use absolute value to find points where curvature changes most dramatically
pred_data_lambda$abs_second_deriv <- abs(pred_data_lambda$second_derivative)

# Find local maxima in second derivative (points of maximum curvature)
# This is more robust than just taking the global maximum
local_max <- which(diff(sign(diff(pred_data_lambda$abs_second_deriv))) == -2) + 1

# If no clear local maxima, use global maximum
if(length(local_max) == 0) {
  threshold_idx <- which.max(pred_data_lambda$abs_second_deriv)
} else {
  # Get the strongest local maximum
  threshold_idx <- local_max[which.max(pred_data_lambda$abs_second_deriv[local_max])]
}

moisture_threshold <- pred_data_lambda$moisture[threshold_idx]

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
    boot_gam <- gam(resp_cube ~ te(moisture, lambda) + s(site, bs="re"),
                    data = boot_data, method = "REML")
    
    # Generate predictions
    boot_preds <- predict(boot_gam, newdata = pred_data_lambda)
    
    # Smooth the predictions
    boot_spline <- smooth.spline(pred_data_lambda$moisture, boot_preds, cv = TRUE)
    
    # Calculate second derivative
    boot_deriv2 <- predict(boot_spline, pred_data_lambda$moisture, deriv = 2)
    
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
    
    boot_thresholds[i] <- pred_data_lambda$moisture[boot_threshold_idx]
    
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

# ----- Threshold Visualization for lambda -----

# ---- Plot A: Moisture Threshold Plot with GAM Fit ----
p1_lambda <- ggplot() +
  # Original data points
  geom_point(data = dom_data, aes(x = moisture, y = resp_cube, color = site),
             alpha = 0.5, size = 3) +
  # GAM fit with confidence band
  geom_line(data = pred_data_lambda, aes(x = moisture, y = fit),
            color = "blue", size = 1.2) +
  # Show confidence intervals for the fit
  geom_ribbon(data = pred_data_lambda,
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
    title = "C",
    x = expression(Median~Gravimetric~Moisture~Content~(g/g)^{1/3}),
    y = expression(paste("Median",~O[2],' consumption rate', ~(mg~o[2]~kg^{-1}~h^{-1})^{1/3}))
  ) +
  theme_bw(base_size = 12) +
  theme(aspect.ratio = 1)

# ---- Plot B: DOM Thermodynamics Split by Moisture Regime ----

library(lme4)

# Mixed models with site as random effect
below_mixed <- lmer(resp_cube ~ lambda + (1|site), 
                    data = dom_data[dom_data$moisture < moisture_threshold, ])
above_mixed <- lmer(resp_cube ~ lambda + (1|site), 
                    data = dom_data[dom_data$moisture >= moisture_threshold, ])

# Summarize models
below_summary <- summary(below_mixed)
above_summary <- summary(above_mixed)

# Extract fixed-effect coefficients (these are your adjusted betas)
below_beta <- fixef(below_mixed)[2]
above_beta <- fixef(above_mixed)[2]

# Get p-values (using lmerTest approach)
library(lmerTest)
below_mixed_t <- lmer(resp_cube ~ lambda + (1|site), 
                      data = dom_data[dom_data$moisture < moisture_threshold, ])
above_mixed_t <- lmer(resp_cube ~ lambda + (1|site), 
                      data = dom_data[dom_data$moisture >= moisture_threshold, ])

below_p <- summary(below_mixed_t)$coefficients[2,5]
above_p <- summary(above_mixed_t)$coefficients[2,5]

# Calculate conditional R-squared (variance explained by fixed + random effects)
library(MuMIn)
below_r2 <- r.squaredGLMM(below_mixed)[1] # Use [2] for conditional R²
above_r2 <- r.squaredGLMM(above_mixed)[1] # Use [2] for conditional R²

# Print comparison with adjusted statistics
cat("DOM-Respiration Relationship by Moisture Regime (Adjusted for Site Effects):\n")
cat("Below threshold (", round(moisture_threshold, 2), "):\n", sep="")
cat("  Adjusted Coefficient:", round(below_beta, 3),
    "\n  Marginal R-squared:", round(below_r2, 3),
    "\n  p-value:", format.pval(below_p, digits = 3), "\n\n")
cat("Above threshold (", round(moisture_threshold, 2), "):\n", sep="")
cat("  Adjusted Coefficient:", round(above_beta, 3),
    "\n  Marginal R-squared:", round(above_r2, 3),
    "\n  p-value:", format.pval(above_p, digits = 3), "\n\n")

# Create prediction frames from mixed models
# Since site is a random effect, we'll use the average site effect
lambda_range <- range(dom_data$lambda, na.rm = TRUE)
lambda_seq <- seq(lambda_range[1], lambda_range[2], length.out = 100)

new_data <- data.frame(lambda = lambda_seq)

# For below threshold
below_pred_fixed <- fixef(below_mixed)[1] + fixef(below_mixed)[2] * lambda_seq
# For confidence intervals, we'll use bootMer to account for uncertainty
library(boot)
bootfun <- function(model) {
  predict(model, newdata = new_data, re.form = NA)
}
below_boot <- bootMer(below_mixed, bootfun, nsim = 500)
below_ci <- t(apply(below_boot$t, 2, quantile, c(0.025, 0.975)))

# Repeat for above threshold
above_pred_fixed <- fixef(above_mixed)[1] + fixef(above_mixed)[2] * lambda_seq
above_boot <- bootMer(above_mixed, bootfun, nsim = 500)
above_ci <- t(apply(above_boot$t, 2, quantile, c(0.025, 0.975)))

# Create dataframes for plotting
below_frame <- data.frame(
  lambda = lambda_seq,
  predicted = below_pred_fixed,
  lower = below_ci[,1],
  upper = below_ci[,2],
  regime = "Below Threshold"
)

above_frame <- data.frame(
  lambda = lambda_seq,
  predicted = above_pred_fixed,
  lower = above_ci[,1],
  upper = above_ci[,2],
  regime = "Above Threshold"
)

pred_combined <- rbind(below_frame, above_frame)


# Create significance indicators
below_signif <- ifelse(below_p < 0.05, 
                       ifelse(below_p < 0.01, 
                              ifelse(below_p < 0.001, "***", "**"), "*"), "ns")

above_signif <- ifelse(above_p < 0.05, 
                       ifelse(above_p < 0.01, 
                              ifelse(above_p < 0.001, "***", "**"), "*"), "ns")

# Plot with adjusted statistics
p2_lambda <- ggplot() +
  # Data points
  geom_point(data = dom_data,
             aes(x = lambda, y = resp_cube, color = moisture_regime),
             alpha = 1, size = 3) +
  # Model fits with confidence intervals
  geom_line(data = pred_combined,
            aes(x = lambda, y = predicted, color = regime),
            size = 1) +
  geom_ribbon(data = pred_combined,
              aes(x = lambda, ymin = lower, ymax = upper, fill = regime),
              alpha = 0.5) +
  # Annotations with site-adjusted statistics
  annotate("text", x = min(dom_data$lambda) + 0.002,
           y =  0.1 * diff(range(dom_data$resp_cube)),
           label = paste0("Below: R² = ", round(below_r2, 2),
                          ", β = ", round(below_beta, 2), below_signif),
           color = "darkorange", hjust = 0, fontface = "bold") +
  annotate("text", x = min(dom_data$lambda) + 0.002,
           y =  0.2 * diff(range(dom_data$resp_cube)),
           label = paste0("Above: R² = ", round(above_r2, 2),
                          ", β = ", round(above_beta, 2), above_signif),
           color = "lightblue", hjust = 0, fontface = "bold") +
  # Theme and labels
  scale_color_manual(values = c("Below Threshold" = "darkorange",
                                "Above Threshold" = "lightblue")) +
  scale_fill_manual(values = c("Below Threshold" = "darkorange",
                               "Above Threshold" = "lightblue")) +
  labs(
    title = "B",
    x = expression(paste(Median~lambda^(1/3))),
    y = expression(paste("Median",~O[2],' consumption rate', ~(mg~o[2]~kg^{-1}~h^{-1})^(1/3))),
    color = " ",
    fill = " "
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'top', aspect.ratio = 1)

# ---- Plot C: Interaction Heatmap with Contours ----
# Create prediction grid for interaction plot
moisture_grid <- seq(min(dom_data$moisture), max(dom_data$moisture), length.out = 50)
lambda_grid <- seq(min(dom_data$lambda), max(dom_data$lambda), length.out = 50)
grid_data <- expand.grid(moisture = moisture_grid, lambda = lambda_grid, site = reference_site)

# Generate predictions from GAM model
grid_data$predicted <- predict(gam_model, newdata = grid_data)

# Calculate contour breaks for the plot
pred_range <- range(grid_data$predicted)
contour_breaks <- seq(pred_range[1], pred_range[2], length.out = 10)

# Create interaction heatmap
p3_lambda <- ggplot(grid_data, aes(x = moisture, y = lambda)) +
  # Create filled contours with custom colormap
  geom_contour_filled(aes(z = predicted),
                      breaks = contour_breaks,
                      alpha = 0.8) +
  # Add contour lines
  geom_contour(aes(z = predicted),
               breaks = contour_breaks,
               color = "black",
               alpha = 0.5,
               size = 0.3) +
  # Add threshold line
  geom_vline(xintercept = moisture_threshold,
             linetype = "dashed",
             color = "black",
             size = 1) +
  # Add data points
  geom_point(data = dom_data,
             aes(size = abs(resp_cube)),
             shape = 21,
             color = "black",
             fill = "white",
             alpha = 0.7) +
  # Use a color scale similar to 'topo' in vis.gam
  scale_fill_viridis_d(option = "plasma",
                       name = expression(paste("Predicted ",~O[2]," Consumption (Cube Root)"))) +
  scale_size_continuous(name = expression(paste("Observed",~O[2], " Consumption (Cube Root)")),
                        range = c(1, 5)) +
  # Properly label the axes with expressions
  labs(
    title = "B",
    x = expression(paste("Median Gravimetric Moisture ", (g/g)^(1/3))),
    y = expression(paste(Median~lambda^(1/3))),
  ) +
  # Theme settings to match vis.gam aesthetics
  theme_bw(base_size = 12) +
  theme(
    aspect.ratio = 1,
    panel.grid = element_blank(),
    legend.position = "right"
  )


# ===== SI plots =====
# ---- SI Plot 1: First and Second Derivatives ----
# Plot First Derivative
si_p1a_lambda <- ggplot(pred_data_lambda, aes(x = moisture, y = first_derivative)) +
  geom_line(size = 1, color = "darkgreen") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = moisture_threshold, color = "red", 
             linetype = "dashed", size = 1) +
  labs(
    title = "A",
    x = expression(Median~Gravimetric~Moisture~Content~(g/g)^{1/3}),
    y = "First Derivative (Rate of Change)"
  ) +
  theme_bw(base_size = 12) + 
  theme(aspect.ratio = 1)

# Plot Second Derivative with threshold
si_p1b_lambda <- ggplot(pred_data_lambda, aes(x = moisture, y = second_derivative)) +
  geom_line(size = 1, color = "purple") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = moisture_threshold, color = "red", 
             linetype = "dashed", size = 1) +
  geom_point(aes(x = moisture[threshold_idx], y = second_derivative[threshold_idx]), 
             color = "red", size = 3) +
  labs(
    title = "B",
    x = expression(Median~Gravimetric~Moisture~Content~(g/g)^{1/3}),
    y = "Second Derivative (Curvature)"
  ) +
  theme_bw(base_size = 12) + 
  theme(aspect.ratio = 1)

# ---- SI Plot 2: Bootstrap Distribution ----
si_p2_lambda <- ggplot(data.frame(threshold = valid_boots), aes(x = threshold)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = moisture_threshold, color = "red", size = 1) +
  geom_vline(xintercept = ci_lower, color = "darkred", linetype = "dashed") +
  geom_vline(xintercept = ci_upper, color = "darkred", linetype = "dashed") +
  annotate("text", x = moisture_threshold - 0.2, 
           y = max(hist(valid_boots, breaks = 30, plot = FALSE)$counts) * 0.9,
           label = paste0("Threshold: ", round(moisture_threshold, 2)), 
           color = "red", hjust = 0) +
  annotate("text", x =  moisture_threshold - 0.2, 
           y = max(hist(valid_boots, breaks = 30, plot = FALSE)$counts) * 0.8,
           label = paste0("95% CI: [", round(ci_lower, 2), ", ", round(ci_upper, 2), "]"), 
           color = "red", hjust = 0) +
  labs(
    title = "C",
    x = expression(Moisture~Threshold~(g/g)^{1/3}),
    y = "Frequency",
    subtitle = paste0(length(valid_boots), " valid bootstrap samples out of ", n_boot)
  ) +
  theme_bw(base_size = 12) +
  theme(aspect.ratio = 1)

# ---- SI Plot 3: Model Validation ----
# Add predictions to dom_data if they're not already there
if(!"predicted" %in% colnames(dom_data)) {
  dom_data$predicted <- predict(gam_model, newdata = dom_data)
}

si_p3_lambda <- ggplot(dom_data, aes(x = predicted, y = resp_cube, color = moisture_regime)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  stat_smooth(method = "lm", se = FALSE, linetype = "solid") +
  scale_color_manual(values = c("Below Threshold" = "darkorange", 
                                "Above Threshold" = "lightblue"), name = "Moisture Regime") +
  labs(
    title = "D",
    x = expression(paste("Predicted",~O[2],' consumption rate', ~(mg~o[2]~kg^{-1}~h^{-1})^{1/3})),
    y = expression(paste("Observed",~O[2],' consumption rate', ~(mg~o[2]~kg^{-1}~h^{-1})^{1/3}))
  ) +
  theme_bw(base_size = 12) +
  theme(aspect.ratio = 1, legend.position = "none")

# ---- SI Plot 4: Statistical Testing ----
# Create text summary for statistical test comparing GAM models
dom_data$moisture_below <- pmin(dom_data$moisture, moisture_threshold)
dom_data$moisture_above <- pmax(0, dom_data$moisture - moisture_threshold)

# Model with single tensor product vs. threshold model with separate tensor products
model_single <- gam(resp_cube ~ te(moisture, lambda) + s(site, bs="re"),
                    data = dom_data, method = "REML")

model_threshold <- gam(resp_cube ~ te(moisture_below, lambda) + te(moisture_above, gibbs) + s(site, bs="re"),
                       data = dom_data, method = "REML")

# Compare models
anova_result <- anova(model_single, model_threshold, test = "F")

# Create text for statistical test
stat_text <- paste0(
  "Statistical Test for Threshold Significance\n\n",
  "Method: F-test comparing nested GAM models with tensor product smooths\n",
  "F-statistic: ", round(anova_result$F[2], 3), "\n",
  "p-value: ", format.pval(anova_result$`p-value`[2], digits = 3), "\n",
  "Significant at α=0.05: ", ifelse(anova_result$`p-value`[2] < 0.05, "Yes", "No"), "\n\n",
  "A significant test indicates that modeling separate\n",
  "moisture-DOM thermodynamics interactions before and after\n",
  "the threshold is statistically preferred.\n\n",
  "Deviance explained (single model): ", 
  round(summary(model_single)$dev.expl * 100, 1), "%\n",
  "Deviance explained (threshold model): ", 
  round(summary(model_threshold)$dev.expl * 100, 1), "%"
)

si_p4_lambda <- ggplot() + 
  annotate("text", x = 0.5, y = 0.5, label = stat_text, size = 5) +
  theme_void() +
  labs(title = "E") +
  theme(plot.title = element_text(hjust = 0, size = 14, face = "bold"))

# ---- SI Plot 5: Residual Diagnostics ----
# Calculate residuals
dom_data$residuals <- residuals(gam_model)

# Create residuals plot by moisture with threshold
si_p5_lambda <- ggplot(dom_data, aes(x = moisture, y = residuals, color = moisture_regime)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = moisture_threshold, color = "red", linetype = "dashed") +
  geom_smooth(method = "loess", se = TRUE, alpha = 0.2) +
  scale_color_manual(values = c("Below Threshold" = "darkorange", 
                                "Above Threshold" = "lightblue"), name = "Moisture Regime") +
  labs(
    title = "F",
    x = expression(Median~Gravimetric~Moisture~Content~(g/g)^{1/3}),
    y = "Model Residuals"
  ) +
  theme_bw(base_size = 12) +
  theme(aspect.ratio = 1, legend.position = "none")

# ---- SI Plot 6: Methodology Diagram ----
# Create a simplified diagram explaining the second derivative method
method_steps <- data.frame(
  x = 1:4,
  y = rep(1, 4),
  label = c("1. Generate predictions\nacross moisture gradient",
            "2. Calculate first and\nsecond derivatives\nusing smooth splines",
            "3. Identify threshold at\nmaximum curvature\n(second derivative)",
            "4. Bootstrap for\nconfidence intervals")
)

arrows_df <- data.frame(
  x = 1:3,
  xend = 2:4,
  y = rep(1, 3),
  yend = rep(1, 3)
)

si_p6_lambda <- ggplot() +
  geom_point(data = method_steps, aes(x = x, y = y), size = 15, color = "steelblue", alpha = 0.7) +
  geom_text(data = method_steps, aes(x = x, y = y, label = label), size = 3.5, color = "white") +
  geom_segment(data = arrows_df, aes(x = x + 0.2, y = y, xend = xend - 0.2, yend = yend),
               arrow = arrow(length = unit(0.3, "cm")), size = 1) +
  labs(title = "G. Second Derivative Method for Threshold Detection") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.margin = margin(20, 20, 20, 20))

# Arrange SI plots
si_p1_lambda <- si_p1a_lambda + si_p1b_lambda  + si_p2_lambda + plot_layout(ncol = 3)

si_row2_lambda <- si_p3_lambda + si_p5_lambda + plot_layout(ncol = 2)

si_row3_lambda <-   si_p4_lambda 

# Save the SI plots
ggsave(paste0(figure_path,"SI_Figure1_lambda.png"), si_p1_lambda, width = 12, height = 16, dpi = 300)
ggsave(paste0(figure_path,"SI_Figure2_lambda.png"), si_row2_lambda, width = 12, height = 16, dpi = 300)
ggsave(paste0(figure_path,"SI_Figure3_lambda.png"), si_row3_lambda, width = 12, height = 16, dpi = 300)


ggsave(paste0(figure_path,"SI_Figure1_lambda.pdf"), si_p1_lambda, width = 12, height = 16)
ggsave(paste0(figure_path,"SI_Figure2_lambda.pdf"), si_row2_lambda, width = 12, height = 16)
ggsave(paste0(figure_path,"SI_Figure3_lambda.pdf"), si_row3_lambda, width = 12, height = 16)

# ----- Combine all six plots into one figure -----
final_plot <- (p1 | p2 ) / (p1_lambda | p2_lambda ) + plot_layout(nrow = 2)
print(final_plot)

# Save the combined plot
ggsave(paste0(figure_path,"Figure2_combined_gibbs_lambda_threshold_plots.png"), final_plot, width = 18, height = 10, dpi = 300)

ggsave(paste0(figure_path,"Figure2_combined_gibbs_lambda_threshold_plots.pdf"), final_plot, width = 18, height = 10)

# Interaction 
interaction_plot <- (p3 | p3_lambda )
print(interaction_plot)

# Save the combined plot
ggsave(paste0(figure_path,"Figure3_Interaction_plots.png"), interaction_plot, width = 18, height = 10, dpi = 300)

ggsave(paste0(figure_path,"Figure3_Interaction_plots.pdf"), interaction_plot, width = 18, height = 10)
