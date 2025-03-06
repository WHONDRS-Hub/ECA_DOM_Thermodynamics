# ==== Loading libraries =========
rm(list=ls(all=T))

# Comprehensive DOM-Moisture Threshold Analysis
# ------------------------------------------------------
# This code identifies and visualizes critical moisture thresholds
# where DOM thermodynamic properties influence respiration differently

# Load required packages
library(mgcv)       # For GAM models
library(ggplot2)    # For visualization
library(dplyr)      # For data manipulation
library(viridis)    # For color palettes
library(patchwork)  # For combining plots
library(gridExtra)  # For arranging multiple plots

# ==== Defining paths and working directories 
github = 'C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/ECA_DOM_Thermodynamics/'
data_path = paste0(github,'Data/')
figure_path = paste0(github,'Figures/')

# ====== Read in data
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

sample_data = sample_data %>%
  dplyr::select(site,Treatment,Median_62948_Final_Gravimetric_Moisture_g_per_g,Median_Respiration_Rate_mg_DO_per_kg_per_H,Median_Extractable_NPOC_mg_per_kg)

# ==== Set up data 
data = data %>%
  dplyr::filter(site != "EC_023") %>%
  dplyr::filter(!(site %in% c("EC_011", "EC_012", "EC_052", "EC_053", "EC_057")))

data = merge(data,sample_data, by = c('site','Treatment'))

data$Treatment <- as.factor(data$Treatment)
data$site <- as.factor(data$site)
combined_data = data

# ====== Cube root transform data ====
cube_root <- function(x) sign(x) * (abs(x))^(1/3)

cube_data = data %>% 
  mutate(across(where(is.numeric), cube_root)) %>% # cube root transform data
  rename_with(where(is.numeric), .fn = ~ paste0("cube_", .x))

dom_data = cube_data %>%
  dplyr::select(site, Treatment, resp_cube = cube_Median_Respiration_Rate_mg_DO_per_kg_per_H, gibbs = cube_Median_delGcoxPerCmol,
                moisture = cube_Median_62948_Final_Gravimetric_Moisture_g_per_g,
                lambda = cube_Median_Lambda)

# ==== Exlore linear regressions =====
p1 <- ggplot(dom_data, aes(x = moisture, y = resp_cube, color = gibbs)) +
  geom_point(alpha = 0.7) +
  scale_color_viridis_c(limits = c(-5.26, -5.22), name = "Gibbs") +
  labs(title = "Moisture vs. Respiration (colored by Gibbs)") +
  theme_bw()

p2 <- ggplot(dom_data, aes(x = moisture, y = resp_cube, color = lambda)) +
  geom_point(alpha = 0.7) +
  scale_color_viridis_c(limits = c(0.32, 0.35), name = "Lambda") +
  labs(title = "Moisture vs. Respiration (colored by Lambda)") +
  theme_bw()

combined_plot <- p1 / p2

# Save as PNG
ggsave("exploratory_linear_regressions.png", combined_plot, width = 10, height = 8, dpi = 300)

# Save as PDF
ggsave("exploratory_linear_regressions.pdf", combined_plot, width = 10, height = 8)

# ===== Moisture thresholds for Gibbs ======
# Fit the improved GAM model with tensor interaction and site random effect
# This properly accounts for the complex interaction and site-specific variation
gam_model <- gam(resp_cube ~ te(moisture, gibbs) + s(site, bs="re"), 
                 data = dom_data, method = "REML")

# Summarize model
summary_gam <- summary(gam_model)
cat("Model summary:\n")
print(summary_gam)
summary(gam_model)$dev.expl * 100  # Gets percentage of deviance explained
sqrt(mean((dom_data$resp_cube - fitted(gam_model))^2))  # RMSE


# ==== Moisture Threshold Identification ====

# Create a sequence along the moisture range for threshold detection
moisture_seq <- seq(min(dom_data$moisture), max(dom_data$moisture), length.out = 100)
reference_site <- levels(factor(dom_data$site))[1]
reference_gibbs <- median(dom_data$gibbs)

# Create prediction frame for the effect of moisture at median Gibbs
pred_frame <- data.frame(
  moisture = moisture_seq,
  gibbs = reference_gibbs,
  site = reference_site
)

# Generate predictions
pred_frame$predicted <- predict(gam_model, newdata = pred_frame)
# Calculate first derivative to find rate of change
pred_frame$derivative <- c(NA, diff(pred_frame$predicted) / diff(pred_frame$moisture))

# Calculate second derivative to find where slope changes most rapidly
# Initialize with NAs
pred_frame$second_deriv <- NA  

# Get positions where first derivative is available (not NA)
valid_deriv <- which(!is.na(pred_frame$derivative))

# We need at least 2 points with valid derivatives to calculate second derivative
if(length(valid_deriv) >= 2) {
  # Calculate second derivative only where we have valid first derivatives
  second_deriv_values <- diff(pred_frame$derivative[valid_deriv]) / 
    diff(pred_frame$moisture[valid_deriv])
  
  # Place results in correct positions (skip first valid derivative position)
  pred_frame$second_deriv[valid_deriv[-1]] <- second_deriv_values
}

# Find the threshold (maximum absolute second derivative)
valid_rows <- which(!is.na(pred_frame$second_deriv))
if(length(valid_rows) > 0) {
  threshold_idx <- valid_rows[which.max(abs(pred_frame$second_deriv[valid_rows]))]
  moisture_threshold <- pred_frame$moisture[threshold_idx]
} else {
  # Fallback if no valid second derivatives are found
  warning("No valid second derivatives found. Using alternative threshold method.")
  # Find where first derivative crosses zero or changes most
  valid_first_deriv <- which(!is.na(pred_frame$derivative))
  if(length(valid_first_deriv) > 0) {
    # Find where derivative changes sign or has maximum change
    deriv_changes <- diff(sign(pred_frame$derivative[valid_first_deriv]))
    if(any(deriv_changes != 0)) {
      # Find where derivative changes sign
      sign_change_idx <- valid_first_deriv[which(deriv_changes != 0)[1] + 1]
      moisture_threshold <- pred_frame$moisture[sign_change_idx]
    } else {
      # Otherwise use maximum absolute derivative
      max_deriv_idx <- valid_first_deriv[which.max(abs(pred_frame$derivative[valid_first_deriv]))]
      moisture_threshold <- pred_frame$moisture[max_deriv_idx]
    }
  } else {
    # Last resort - use median moisture
    moisture_threshold <- median(dom_data$moisture)
    warning("Using median moisture as threshold due to analysis limitations.")
  }
}

# Bootstrap to get confidence interval
n_boot <- 100
boot_thresholds <- numeric(n_boot)

cat("Bootstrap thresholds (this may take a moment)...\n")
for(i in 1:n_boot) {
  # Bootstrap sample
  boot_idx <- sample(1:nrow(dom_data), replace = TRUE)
  boot_data <- dom_data[boot_idx, ]
  
  # Fit GAM to bootstrap sample
  tryCatch({
    boot_model <- gam(resp_cube ~ te(moisture, gibbs) + s(site, bs="re"), 
                      data = boot_data, method = "REML")
    
    # Generate predictions
    boot_preds <- predict(boot_model, newdata = pred_frame)
    
    # Calculate first derivative
    boot_deriv <- c(NA, diff(boot_preds) / diff(pred_frame$moisture))
    
    # Initialize second derivative
    boot_deriv2 <- rep(NA, length(boot_deriv))
    
    # Calculate second derivative properly
    valid_deriv <- which(!is.na(boot_deriv))
    if(length(valid_deriv) >= 2) {
      second_deriv_values <- diff(boot_deriv[valid_deriv]) / 
        diff(pred_frame$moisture[valid_deriv])
      boot_deriv2[valid_deriv[-1]] <- second_deriv_values
    }
    
    # Find threshold
    valid_boot <- which(!is.na(boot_deriv2))
    if(length(valid_boot) > 0) {
      max_idx <- which.max(abs(boot_deriv2[valid_boot]))
      boot_thresholds[i] <- pred_frame$moisture[valid_boot[max_idx]]
    } else if(length(valid_deriv) > 0) {
      # Fallback to first derivative if no valid second derivatives
      max_idx <- which.max(abs(boot_deriv[valid_deriv]))
      boot_thresholds[i] <- pred_frame$moisture[valid_deriv[max_idx]]
    } else {
      boot_thresholds[i] <- NA
    }
  }, error = function(e) {
    boot_thresholds[i] <- NA
  })
  
  # Progress indicator
  if(i %% 10 == 0) cat(".")
  if(i %% 50 == 0) cat(" ", i, "\n")
}
cat("\n")

# Calculate confidence interval
valid_boots <- boot_thresholds[!is.na(boot_thresholds)]
ci_lower <- quantile(valid_boots, 0.025)
ci_upper <- quantile(valid_boots, 0.975)

cat("Moisture threshold:", round(moisture_threshold, 3), 
    "\n95% CI:", round(ci_lower, 3), "-", round(ci_upper, 3), "\n\n")

# ------------------------------------------------------
# 3. Visualizing Main Threshold Effect
# ------------------------------------------------------

# Main threshold visualization
p1 <- ggplot() +
  # Original data points
  geom_point(data = dom_data, aes(x = moisture, y = resp_cube, color = site), 
             alpha = 0.4, size = 1.5) +
  # GAM fit
  geom_line(data = pred_frame, aes(x = moisture, y = predicted),
            color = "blue", size = 1.2) +
  # Threshold line with confidence interval
  geom_vline(xintercept = moisture_threshold, linetype = "dashed", 
             color = "red", size = 1) +
  geom_rect(aes(xmin = ci_lower, xmax = ci_upper, 
                ymin = -Inf, ymax = Inf),
            alpha = 0.1, fill = "red") +
  # Annotations
  annotate("text", x = moisture_threshold + 0.1, 
           y = 0,
           label = paste0("Threshold: ", round(moisture_threshold, 2)), 
           color = "red", hjust = 0) +
  annotate("text", x = moisture_threshold + 0.1, 
           y = -1,
           label = paste0("(95% CI: ", 
                          round(ci_lower, 2), "-", round(ci_upper, 2), ")"), 
           color = "red", hjust = 0) +
  # Aesthetics
  scale_color_viridis_d(guide = "none") +  # Colorful sites without legend
  labs(
    title = "A",
    x = "Cube Root Gravimetric Moisture Content (g/g)",
    y = expression(paste("Cube Root Median Respiration Rate", ~(mg~DO~kg^{-1}~h^{-1})))
  ) +
  theme_bw(base_size = 12)+ theme(aspect.ratio = 1)

# Rate of change visualization
p2 <- ggplot() +
  geom_line(data = pred_frame[!is.na(pred_frame$derivative), ], 
            aes(x = moisture, y = derivative),
            size = 1, color = "darkblue") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
  geom_vline(xintercept = moisture_threshold, linetype = "dashed", 
             color = "red", size = 1) +
  geom_rect(aes(xmin = ci_lower, xmax = ci_upper, 
                ymin = -Inf, ymax = Inf),
            alpha = 0.1, fill = "red") +
  labs(
    title = "A",
    x = "Cube Root Gravimetric Moisture Content (g/g)",
    y = "Slope (dRespiration/dMoisture)"
  ) +
  theme_bw(base_size = 12)

# ------------------------------------------------------
# 4. Split Analysis by Threshold
# ------------------------------------------------------

# Create categorical variable for before/after threshold
dom_data$moisture_regime <- ifelse(dom_data$moisture < moisture_threshold, 
                                   "Below Threshold", "Above Threshold")

# Fit separate models for each regime
below_model <- lm(resp_cube ~ gibbs, 
                  data = dom_data[dom_data$moisture < moisture_threshold, ])
above_model <- lm(resp_cube ~ gibbs, 
                  data = dom_data[dom_data$moisture >= moisture_threshold, ])

# Summarize models
below_summary <- summary(below_model)
above_summary <- summary(above_model)

# Print comparison
cat("DOM-Respiration Relationship by Moisture Regime:\n")
cat("Below threshold (", round(moisture_threshold, 2), "):\n", sep="")
cat("  Coefficient:", round(coef(below_model)[2], 3), 
    "\n  R-squared:", round(below_summary$r.squared, 3), 
    "\n  p-value:", format.pval(below_summary$coefficients[2, 4], digits = 3), "\n\n")

cat("Above threshold (", round(moisture_threshold, 2), "):\n", sep="")
cat("  Coefficient:", round(coef(above_model)[2], 3), 
    "\n  R-squared:", round(above_summary$r.squared, 3), 
    "\n  p-value:", format.pval(above_summary$coefficients[2, 4], digits = 3), "\n\n")

# Create prediction frames for visualization
gibbs_range <- range(dom_data$gibbs)
gibbs_seq <- seq(gibbs_range[1], gibbs_range[2], length.out = 100)

below_preds <- predict(below_model, 
                       newdata = data.frame(gibbs = gibbs_seq), 
                       interval = "confidence")
above_preds <- predict(above_model, 
                       newdata = data.frame(gibbs = gibbs_seq), 
                       interval = "confidence")

below_frame <- data.frame(
  gibbs = gibbs_seq,
  predicted = below_preds[, "fit"],
  lower = below_preds[, "lwr"],
  upper = below_preds[, "upr"],
  regime = "Below Threshold"
)

above_frame <- data.frame(
  gibbs = gibbs_seq,
  predicted = above_preds[, "fit"],
  lower = above_preds[, "lwr"],
  upper = above_preds[, "upr"],
  regime = "Above Threshold"
)

pred_combined <- rbind(below_frame, above_frame)

# DOM thermodynamics split visualization
p3 <- ggplot() +
  # Data points
  geom_point(data = dom_data, 
             aes(x = gibbs, y = resp_cube, color = moisture_regime),
             alpha = 0.6, size = 2) +
  # Model fits with confidence intervals
  geom_line(data = pred_combined, 
            aes(x = gibbs, y = predicted, color = regime),
            size = 1) +
  geom_ribbon(data = pred_combined, 
              aes(x = gibbs, ymin = lower, ymax = upper, fill = regime),
              alpha = 0.2) +
  # Annotations
  annotate("text", x = min(dom_data$gibbs) + 0.002, 
           y = 1,
           label = paste0("Below: R² = ", round(below_summary$r.squared, 2),
                          ", β = ", round(coef(below_model)[2], 2)),
           color = "darkorange", hjust = 0) +
  annotate("text", x = min(dom_data$gibbs) + 0.002, 
           y = 2,
           label = paste0("Above: R² = ", round(above_summary$r.squared, 2),
                          ", β = ", round(coef(above_model)[2], 2)),
           color = "darkblue", hjust = 0) +
  # Theme and labels
  scale_color_manual(values = c("Below Threshold" = "darkorange", 
                                "Above Threshold" = "darkblue")) +
  scale_fill_manual(values = c("Below Threshold" = "darkorange", 
                               "Above Threshold" = "darkblue")) +
  labs(
    title = "B",
    x = expression(paste("Cube Root Median ", Delta, G[cox], ~(kJ~Cmol^{-1}))),
    y = expression(paste("Cube Root Median Respiration Rate ", ~(mg~DO~kg^{-1}~h^{-1}))),
    color = " ",
    fill = " "
  ) +
  theme_bw(base_size = 12)+  theme(legend.position = 'top')+ theme(aspect.ratio = 1)

# ------------------------------------------------------
# 5. 2D Interaction Visualization
# ------------------------------------------------------

# Create a detailed prediction grid
moisture_grid <- seq(min(dom_data$moisture), max(dom_data$moisture), length.out = 40)
gibbs_grid <- seq(min(dom_data$gibbs), max(dom_data$gibbs), length.out = 40)
grid_data <- expand.grid(moisture = moisture_grid, gibbs = gibbs_grid)

# Add reference site for prediction
grid_data$site <- reference_site

# Generate predictions
grid_data$predicted <- predict(gam_model, newdata = grid_data)

# Heat map visualization of the interaction
p4 <- ggplot(grid_data, aes(x = moisture, y = gibbs, fill = predicted)) +
  geom_tile() +
  # Threshold line
  geom_vline(xintercept = moisture_threshold, linetype = "dashed", 
             color = "black", size = 1) +
  # Aesthetics
  scale_fill_viridis_c(name = "Predicted\nCube Root Respiration",
                       option = "plasma") +
  labs(
    title = "C",
    x = "Cube Root Gravimetric Moisture Content (g/g)",
    y = expression(paste("Cube Root Median ", Delta, G[cox], ~(kJ~Cmol^{-1})))
  ) +
  theme_bw(base_size = 12)+ theme(legend.position = 'top') +
theme(aspect.ratio = 1)

# Contour plot alternative
p5 <- ggplot(grid_data, aes(x = moisture, y = gibbs, z = predicted)) +
  geom_contour_filled() +
  geom_vline(xintercept = moisture_threshold, linetype = "dashed", 
             color = "black", size = 1) +
  labs(
    title = "C",
    x = "Cube Root Gravimetric Moisture Content (g/g)",
    y = expression(paste("Cube Root Median ", Delta, G[cox], ~(kJ~Cmol^{-1}))),
    fill = "Predicted\n Cube Root Respiration"
  ) +
  theme_bw(base_size = 12)+ theme(aspect.ratio = 1)

# ------------------------------------------------------
# 6. Calculate Marginal Effect of DOM Properties
# ------------------------------------------------------

# Function to calculate marginal effect of Gibbs at different moisture levels
get_gibbs_effect <- function(moisture_val) {
  # Create data for high and low Gibbs at this moisture
  low_gibbs <- min(dom_data$gibbs)
  high_gibbs <- max(dom_data$gibbs)
  
  test_data <- data.frame(
    moisture = rep(moisture_val, 2),
    gibbs = c(low_gibbs, high_gibbs),
    site = rep(reference_site, 2)
  )
  
  # Predict respiration
  preds <- predict(gam_model, newdata = test_data)
  
  # Calculate marginal effect
  effect <- (preds[2] - preds[1]) / (high_gibbs - low_gibbs)
  return(effect)
}

# Calculate marginal effect across moisture range
effect_data <- data.frame(
  moisture = moisture_seq,
  gibbs_effect = sapply(moisture_seq, get_gibbs_effect)
)

# Visualize changing marginal effect
p6 <- ggplot(effect_data, aes(x = moisture, y = gibbs_effect)) +
  geom_line(size = 1.2, color = "purple") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
  geom_vline(xintercept = moisture_threshold, linetype = "dashed", 
             color = "red", size = 1) +
  labs(
    title = "D",
    x = "Cube Root Median Gravimetric Moisture Content (g/g)",
    y = "Effect of Gibbs on Respiration"
  ) +
  theme_bw(base_size = 12)+ theme(aspect.ratio = 1)

# ------------------------------------------------------
# 7. Formal Statistical Test of Interaction
# ------------------------------------------------------

# Test if DOM-respiration relationship differs between regimes using parametric model
interaction_test <- lm(resp_cube ~ gibbs * moisture_regime + site, data = dom_data)
interaction_summary <- summary(interaction_test)

cat("Statistical test for difference in DOM-respiration relationship\n")
cat("----------------------------------------------------------------\n")
cat("Interaction model results:\n")
print(interaction_summary)

# Extract interaction term significance
interaction_pval <- interaction_summary$coefficients["gibbs:moisture_regimeBelow Threshold", 4]
cat("\nInteraction p-value:", format.pval(interaction_pval, digits = 3), "\n")

# ------------------------------------------------------
# 8. Combine and Save Plots
# ------------------------------------------------------

# Combine main plots
main_plot <- (p1 | p3) / (p4 | p6)  # Using patchwork for elegant layout

# Add title to combined plot
main_plot <- main_plot + 
  plot_annotation(
    title = " ",
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5))
  )

# Display combined plot
print(main_plot)

# Save individual plots
ggsave("moisture_threshold_main.png", p1, width = 8, height = 6, dpi = 300)
ggsave("dom_response_by_regime.png", p3, width = 8, height = 6, dpi = 300)
ggsave("moisture_dom_interaction.png", p4, width = 8, height = 6, dpi = 300)
ggsave("marginal_effect_plot.png", p6, width = 8, height = 6, dpi = 300)

# Save combined plot
ggsave("combined_threshold_plots.png", main_plot, width = 12, height = 10, dpi = 300)

# Save individual plots as PDF
ggsave("moisture_threshold_main.pdf", p1, width = 8, height = 6)
ggsave("dom_response_by_regime.pdf", p3, width = 8, height = 6)
ggsave("moisture_dom_interaction.pdf", p4, width = 8, height = 6)
ggsave("marginal_effect_plot.pdf", p6, width = 8, height = 6)

# Save combined plot as PDF
ggsave("combined_threshold_plots.pdf", main_plot, width = 12, height = 10)
cat("\nAnalysis complete! All plots saved.\n")

# ------------------------------------------------------
# ===== Additional for Lambda =======
# ------------------------------------------------------
# Fit the improved GAM model with tensor interaction and site random effect
# This properly accounts for the complex interaction and site-specific variation
gam_model <- gam(resp_cube ~ te(moisture, lambda) + s(site, bs="re"), 
                 data = dom_data, method = "REML")

# Summarize model
summary_gam <- summary(gam_model)
cat("Model summary:\n")
print(summary_gam)
summary(gam_model)$dev.expl * 100  # Gets percentage of deviance explained
sqrt(mean((dom_data$resp_cube - fitted(gam_model))^2))  # RMSE

# ==== Moisture Threshold Identification ====

# Create a sequence along the moisture range for threshold detection
moisture_seq <- seq(min(dom_data$moisture), max(dom_data$moisture), length.out = 100)
reference_site <- levels(factor(dom_data$site))[1]
reference_lambda <- median(dom_data$lambda)

# Create prediction frame for the effect of moisture at median lambda
pred_frame <- data.frame(
  moisture = moisture_seq,
  lambda = reference_lambda,
  site = reference_site
)

# Generate predictions
pred_frame$predicted <- predict(gam_model, newdata = pred_frame)
# Calculate first derivative to find rate of change
pred_frame$derivative <- c(NA, diff(pred_frame$predicted) / diff(pred_frame$moisture))

# Calculate second derivative to find where slope changes most rapidly
# Initialize with NAs
pred_frame$second_deriv <- NA  

# Get positions where first derivative is available (not NA)
valid_deriv <- which(!is.na(pred_frame$derivative))

# We need at least 2 points with valid derivatives to calculate second derivative
if(length(valid_deriv) >= 2) {
  # Calculate second derivative only where we have valid first derivatives
  second_deriv_values <- diff(pred_frame$derivative[valid_deriv]) / 
    diff(pred_frame$moisture[valid_deriv])
  
  # Place results in correct positions (skip first valid derivative position)
  pred_frame$second_deriv[valid_deriv[-1]] <- second_deriv_values
}

# Find the threshold (maximum absolute second derivative)
valid_rows <- which(!is.na(pred_frame$second_deriv))
if(length(valid_rows) > 0) {
  threshold_idx <- valid_rows[which.max(abs(pred_frame$second_deriv[valid_rows]))]
  moisture_threshold <- pred_frame$moisture[threshold_idx]
} else {
  # Fallback if no valid second derivatives are found
  warning("No valid second derivatives found. Using alternative threshold method.")
  # Find where first derivative crosses zero or changes most
  valid_first_deriv <- which(!is.na(pred_frame$derivative))
  if(length(valid_first_deriv) > 0) {
    # Find where derivative changes sign or has maximum change
    deriv_changes <- diff(sign(pred_frame$derivative[valid_first_deriv]))
    if(any(deriv_changes != 0)) {
      # Find where derivative changes sign
      sign_change_idx <- valid_first_deriv[which(deriv_changes != 0)[1] + 1]
      moisture_threshold <- pred_frame$moisture[sign_change_idx]
    } else {
      # Otherwise use maximum absolute derivative
      max_deriv_idx <- valid_first_deriv[which.max(abs(pred_frame$derivative[valid_first_deriv]))]
      moisture_threshold <- pred_frame$moisture[max_deriv_idx]
    }
  } else {
    # Last resort - use median moisture
    moisture_threshold <- median(dom_data$moisture)
    warning("Using median moisture as threshold due to analysis limitations.")
  }
}

# Bootstrap to get confidence interval
n_boot <- 100
boot_thresholds <- numeric(n_boot)

cat("Bootstrap thresholds (this may take a moment)...\n")
for(i in 1:n_boot) {
  # Bootstrap sample
  boot_idx <- sample(1:nrow(dom_data), replace = TRUE)
  boot_data <- dom_data[boot_idx, ]
  
  # Fit GAM to bootstrap sample
  tryCatch({
    boot_model <- gam(resp_cube ~ te(moisture, lambda) + s(site, bs="re"), 
                      data = boot_data, method = "REML")
    
    # Generate predictions
    boot_preds <- predict(boot_model, newdata = pred_frame)
    
    # Calculate first derivative
    boot_deriv <- c(NA, diff(boot_preds) / diff(pred_frame$moisture))
    
    # Initialize second derivative
    boot_deriv2 <- rep(NA, length(boot_deriv))
    
    # Calculate second derivative properly
    valid_deriv <- which(!is.na(boot_deriv))
    if(length(valid_deriv) >= 2) {
      second_deriv_values <- diff(boot_deriv[valid_deriv]) / 
        diff(pred_frame$moisture[valid_deriv])
      boot_deriv2[valid_deriv[-1]] <- second_deriv_values
    }
    
    # Find threshold
    valid_boot <- which(!is.na(boot_deriv2))
    if(length(valid_boot) > 0) {
      max_idx <- which.max(abs(boot_deriv2[valid_boot]))
      boot_thresholds[i] <- pred_frame$moisture[valid_boot[max_idx]]
    } else if(length(valid_deriv) > 0) {
      # Fallback to first derivative if no valid second derivatives
      max_idx <- which.max(abs(boot_deriv[valid_deriv]))
      boot_thresholds[i] <- pred_frame$moisture[valid_deriv[max_idx]]
    } else {
      boot_thresholds[i] <- NA
    }
  }, error = function(e) {
    boot_thresholds[i] <- NA
  })
  
  # Progress indicator
  if(i %% 10 == 0) cat(".")
  if(i %% 50 == 0) cat(" ", i, "\n")
}
cat("\n")

# Calculate confidence interval
valid_boots <- boot_thresholds[!is.na(boot_thresholds)]
ci_lower <- quantile(valid_boots, 0.025)
ci_upper <- quantile(valid_boots, 0.975)

cat("Moisture threshold:", round(moisture_threshold, 3), 
    "\n95% CI:", round(ci_lower, 3), "-", round(ci_upper, 3), "\n\n")

# ------------------------------------------------------
# 3. Visualizing Main Threshold Effect
# ------------------------------------------------------

# Main threshold visualization
p1 <- ggplot() +
  # Original data points
  geom_point(data = dom_data, aes(x = moisture, y = resp_cube, color = site), 
             alpha = 0.4, size = 1.5) +
  # GAM fit
  geom_line(data = pred_frame, aes(x = moisture, y = predicted),
            color = "blue", size = 1.2) +
  # Threshold line with confidence interval
  geom_vline(xintercept = moisture_threshold, linetype = "dashed", 
             color = "red", size = 1) +
  geom_rect(aes(xmin = ci_lower, xmax = ci_upper, 
                ymin = -Inf, ymax = Inf),
            alpha = 0.1, fill = "red") +
  # Annotations
  annotate("text", x = moisture_threshold + 0.1,
           y = 0,
           label = paste0("Threshold = ", round(moisture_threshold, 2)),color = "red", hjust = 0) +
  annotate("text", x = moisture_threshold + 0.1,
           y = -1,
           label = paste0(" (95% CI: ", 
                          round(ci_lower, 2), "-", round(ci_upper, 2), ")"),color = "red", hjust = 0) +
  # Aesthetics
  scale_color_viridis_d(guide = "none") +  # Colorful sites without legend
  labs(
    title = "A",
    x = "Cube Root Gravimetric Moisture Content (g/g)",
    y = expression(paste("Cube Root Median Respiration Rate ", ~(mg~DO~kg^{-1}~h^{-1}))),
  ) +
  theme_bw(base_size = 12)+  theme(aspect.ratio = 1)

# Rate of change visualization
p2 <- ggplot() +
  geom_line(data = pred_frame[!is.na(pred_frame$derivative), ], 
            aes(x = moisture, y = derivative),
            size = 1, color = "darkblue") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
  geom_vline(xintercept = moisture_threshold, linetype = "dashed", 
             color = "red", size = 1) +
  geom_rect(aes(xmin = ci_lower, xmax = ci_upper, 
                ymin = -Inf, ymax = Inf),
            alpha = 0.1, fill = "red") +
  labs(
    title = "Rate of Change in Respiration with Moisture",
    x = "Cube Root Gravimetric Moisture Content (g/g)",
    y = "Slope (dRespiration/dMoisture)"
  ) +
  theme_bw(base_size = 12)

# ------------------------------------------------------
# 4. Split Analysis by Threshold
# ------------------------------------------------------

# Create categorical variable for before/after threshold
dom_data$moisture_regime <- ifelse(dom_data$moisture < moisture_threshold, 
                                   "Below Threshold", "Above Threshold")

# Fit separate models for each regime
below_model <- lm(resp_cube ~ lambda, 
                  data = dom_data[dom_data$moisture < moisture_threshold, ])
above_model <- lm(resp_cube ~ lambda, 
                  data = dom_data[dom_data$moisture >= moisture_threshold, ])

# Summarize models
below_summary <- summary(below_model)
above_summary <- summary(above_model)

# Print comparison
cat("DOM-Respiration Relationship by Moisture Regime:\n")
cat("Below threshold (", round(moisture_threshold, 2), "):\n", sep="")
cat("  Coefficient:", round(coef(below_model)[2], 3), 
    "\n  R-squared:", round(below_summary$r.squared, 3), 
    "\n  p-value:", format.pval(below_summary$coefficients[2, 4], digits = 3), "\n\n")

cat("Above threshold (", round(moisture_threshold, 2), "):\n", sep="")
cat("  Coefficient:", round(coef(above_model)[2], 3), 
    "\n  R-squared:", round(above_summary$r.squared, 3), 
    "\n  p-value:", format.pval(above_summary$coefficients[2, 4], digits = 3), "\n\n")

# Create prediction frames for visualization
lambda_range <- range(dom_data$lambda)
lambda_seq <- seq(lambda_range[1], lambda_range[2], length.out = 100)

below_preds <- predict(below_model, 
                       newdata = data.frame(lambda = lambda_seq), 
                       interval = "confidence")
above_preds <- predict(above_model, 
                       newdata = data.frame(lambda = lambda_seq), 
                       interval = "confidence")

below_frame <- data.frame(
  lambda = lambda_seq,
  predicted = below_preds[, "fit"],
  lower = below_preds[, "lwr"],
  upper = below_preds[, "upr"],
  regime = "Below Threshold"
)

above_frame <- data.frame(
  lambda = lambda_seq,
  predicted = above_preds[, "fit"],
  lower = above_preds[, "lwr"],
  upper = above_preds[, "upr"],
  regime = "Above Threshold"
)

pred_combined <- rbind(below_frame, above_frame)

# DOM thermodynamics split visualization
p3 <- ggplot() +
  # Data points
  geom_point(data = dom_data, 
             aes(x = lambda, y = resp_cube, color = moisture_regime),
             alpha = 0.6, size = 2) +
  # Model fits with confidence intervals
  geom_line(data = pred_combined, 
            aes(x = lambda, y = predicted, color = regime),
            size = 1) +
  geom_ribbon(data = pred_combined, 
              aes(x = lambda, ymin = lower, ymax = upper, fill = regime),
              alpha = 0.2) +
  # Annotations
  annotate("text", x = min(dom_data$lambda) + 0.002, 
           y = max(dom_data$resp_cube) + 1,
           label = paste0("Below: R² = ", round(below_summary$r.squared, 2),
                          ", β = ", round(coef(below_model)[2], 2)),
           color = "darkorange", hjust = 0) +
  annotate("text", x = min(dom_data$lambda) + 0.002, 
           y = max(dom_data$resp_cube) + 2,
           label = paste0("Above: R² = ", round(above_summary$r.squared, 2),
                          ", β = ", round(coef(above_model)[2], 2)),
           color = "darkblue", hjust = 0) +
  # Theme and labels
  scale_color_manual(values = c("Below Threshold" = "darkorange", 
                                "Above Threshold" = "darkblue")) +
  scale_fill_manual(values = c("Below Threshold" = "darkorange", 
                               "Above Threshold" = "darkblue")) +
  labs(
    title = "B",
    x = "Cube Root Lambda",
    y = expression(paste("Cube Root Median Respiration Rate ", ~(mg~DO~kg^{-1}~h^{-1}))),
    color = " ",
    fill = " "
  ) + 
  theme_bw(base_size = 12)+
  theme(legend.position = 'top')+  theme(aspect.ratio = 1)

# ------------------------------------------------------
# 5. 2D Interaction Visualization
# ------------------------------------------------------

# Create a detailed prediction grid
moisture_grid <- seq(min(dom_data$moisture), max(dom_data$moisture), length.out = 40)
lambda_grid <- seq(min(dom_data$lambda), max(dom_data$lambda), length.out = 40)
grid_data <- expand.grid(moisture = moisture_grid, lambda = lambda_grid)

# Add reference site for prediction
grid_data$site <- reference_site

# Generate predictions
grid_data$predicted <- predict(gam_model, newdata = grid_data)

# Heat map visualization of the interaction
p4 <- ggplot(grid_data, aes(x = moisture, y = lambda, fill = predicted)) +
  geom_tile() +
  # Threshold line
  geom_vline(xintercept = moisture_threshold, linetype = "dashed", 
             color = "black", size = 1) +
  # Aesthetics
  scale_fill_viridis_c(name = "Predicted\nCube Root Respiration",
                       option = "plasma") +
  labs(
    title = "C",
    x = "Cube Root Gravimetric Moisture Content",
    y = "Cube Root Lambda"
  ) +
  theme_bw(base_size = 12)+
  theme(legend.position = 'top')+  theme(aspect.ratio = 1)


# Contour plot alternative
p5 <- ggplot(grid_data, aes(x = moisture, y = lambda, z = predicted)) +
  geom_contour_filled() +
  geom_vline(xintercept = moisture_threshold, linetype = "dashed", 
             color = "black", size = 1) +
  labs(
    title = "C",
    x = "Cube Root Gravimetric Moisture Content",
    y = "Cube Root Lambda",
    fill = "Predicted\nCube Root Respiration"
  ) +
  theme_bw(base_size = 12) +  theme(aspect.ratio = 1)

# ------------------------------------------------------
# 6. Calculate Marginal Effect of DOM Properties
# ------------------------------------------------------

# Function to calculate marginal effect of lambda at different moisture levels
get_lambda_effect <- function(moisture_val) {
  # Create data for high and low lambda at this moisture
  low_lambda <- min(dom_data$lambda)
  high_lambda <- max(dom_data$lambda)
  
  test_data <- data.frame(
    moisture = rep(moisture_val, 2),
    lambda = c(low_lambda, high_lambda),
    site = rep(reference_site, 2)
  )
  
  # Predict respiration
  preds <- predict(gam_model, newdata = test_data)
  
  # Calculate marginal effect
  effect <- (preds[2] - preds[1]) / (high_lambda - low_lambda)
  return(effect)
}

# Calculate marginal effect across moisture range
effect_data <- data.frame(
  moisture = moisture_seq,
  lambda_effect = sapply(moisture_seq, get_lambda_effect)
)

# Visualize changing marginal effect
p6 <- ggplot(effect_data, aes(x = moisture, y = lambda_effect)) +
  geom_line(size = 1.2, color = "purple") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
  geom_vline(xintercept = moisture_threshold, linetype = "dashed", 
             color = "red", size = 1) +
  labs(
    title = "D",
    x = "Cube Root Gravimetric Moisture Content (g/g)",
    y = "Effect of Lambda on Respiration"
  ) +
  theme_bw(base_size = 12)+  theme(aspect.ratio = 1)

# ------------------------------------------------------
# 7. Formal Statistical Test of Interaction
# ------------------------------------------------------

# Test if DOM-respiration relationship differs between regimes using parametric model
interaction_test <- lm(resp_cube ~ lambda * moisture_regime + site, data = dom_data)
interaction_summary <- summary(interaction_test)

cat("Statistical test for difference in DOM-respiration relationship\n")
cat("----------------------------------------------------------------\n")
cat("Interaction model results:\n")
print(interaction_summary)

# Extract interaction term significance
interaction_pval <- interaction_summary$coefficients["lambda:moisture_regimeBelow Threshold", 4]
cat("\nInteraction p-value:", format.pval(interaction_pval, digits = 3), "\n")

# ------------------------------------------------------
# 8. Combine and Save Plots
# ------------------------------------------------------

# Combine main plots
main_plot <- (p1 | p3) / (p4 | p6)  # Using patchwork for elegant layout

# Add title to combined plot
main_plot <- main_plot + 
  plot_annotation(
    theme = theme(plot.title = element_text(size = 8, hjust = 0.5))
  )

# Display combined plot
print(main_plot)

# Save individual plots
ggsave("moisture_threshold_main_lambda.png", p1, width = 8, height = 6, dpi = 300)
ggsave("dom_response_by_regime_lambda.png", p3, width = 8, height = 6, dpi = 300)
ggsave("moisture_dom_interaction_lambda.png", p4, width = 8, height = 6, dpi = 300)
ggsave("marginal_effect_plot_lambda.png", p6, width = 8, height = 6, dpi = 300)

# Save combined plot
ggsave("combined_threshold_plots_lambda.png", main_plot, width = 12, height = 10, dpi = 300)

# Save individual plots as PDF
ggsave("moisture_threshold_main_lambda.pdf", p1, width = 8, height = 6)
ggsave("dom_response_by_regime_lambda.pdf", p3, width = 8, height = 6)
ggsave("moisture_dom_interaction_lambda.pdf", p4, width = 8, height = 6)
ggsave("marginal_effect_plot_lambda.pdf", p6, width = 8, height = 6)

# Save combined plot as PDF
ggsave("combined_threshold_plots_lambda.pdf", main_plot, width = 12, height = 10)
cat("\nAnalysis complete! All plots saved.\n")
