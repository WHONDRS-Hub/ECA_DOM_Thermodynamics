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
# ======= Linear regressions with treatment ====
library(broom)  # For extracting model stats

# Split data by treatment and fit separate linear models
wet_data <- dom_data %>% filter(Treatment == "Wet")
dry_data <- dom_data %>% filter(Treatment == "Dry")

# Fit separate linear models for each treatment (no random effects needed)
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
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  # Dynamic colored text for each treatment
  annotate("text", x = Inf, y = -Inf, 
           label = gibbs_wet_label,
           hjust = 1.1, vjust = -0.2, size = 4, color = "#4682B4", fontface = "bold") +
  annotate("text", x = Inf, y = -Inf, 
           label = gibbs_dry_label,
           hjust = 1.1, vjust = -1.5, size = 4, color = "#8B4513", fontface = "bold") +
  scale_color_manual(values = c("Dry" = "#8B4513", "Wet" = "#4682B4")) +
  scale_fill_manual(values = c("Dry" = "#8B4513", "Wet" = "#4682B4")) +
  labs(title = "A",
       x = expression(paste("Median ", Delta, G[cox], " (kJ Cmol"^{-1},")"^(1/3))),
       y = expression(paste("Median ",~O[2],' consumption rate ', ~(mg~O[2]~kg^{-1}~h^{-1})^{1/3}))) +
  theme_bw()

plot2 <- ggplot(dom_data, aes(x = lambda, y = resp_cube, color = Treatment)) +
  geom_point(size = 3, alpha = 0.7) +
  #geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  # Dynamic colored text for each treatment
  annotate("text", x = -Inf, y = -Inf,
           label = lambda_wet_label,
           hjust = -0.1, vjust = -0.1, size = 4, color = "#4682B4", fontface = "bold") +
  annotate("text", x = -Inf, y = -Inf,
           label = lambda_dry_label,
           hjust = -0.1, vjust = -1.5, size = 4, color = "#8B4513", fontface = "bold") +
  scale_color_manual(values = c("Dry" = "#8B4513", "Wet" = "#4682B4")) +
  scale_fill_manual(values = c("Dry" = "#8B4513", "Wet" = "#4682B4")) +
  labs(title = "B",
       x = expression(paste("Median ", lambda^(1/3))),
       y = expression(paste("Median ",~O[2],' consumption rate ', ~(mg~O[2]~kg^{-1}~h^{-1})^{1/3}))) +
  theme_bw()

# Combine plots
combined_plot <- plot1 + plot2 + 
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

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

# =====
# Simple reality check
cor(dom_data$moisture, dom_data$gibbs)  # Are predictors correlated?
plot(dom_data$moisture, dom_data$gibbs)  # Visual check

# How much variance does each predictor explain alone?
simple_moisture <- lmer(resp_cube ~ moisture + (1|site), data = dom_data)
simple_gibbs <- lmer(resp_cube ~ gibbs + (1|site), data = dom_data)

r.squaredGLMM(simple_moisture)
r.squaredGLMM(simple_gibbs)
##############################
# make sure treatment is a factor
dom_data$Treatment <- factor(dom_data$Treatment, levels = c("Dry","Wet"))
dom_data$gibbs  <- scale(dom_data$gibbs,   center=TRUE, scale=FALSE)

# main‐effects only
mod_cat       <- lmer(resp_cube ~ gibbs + Treatment + (1 | site),
                      data = dom_data)

# if you want to allow the gibbs‐effect to differ by treatment:
mod_cat_int   <- lmer(resp_cube ~ gibbs * Treatment + (1 | site),
                      data = dom_data)

summary(mod_cat)
summary(mod_cat_int)
anova(mod_cat_int)
AIC(mod_cat)
AIC(mod_cat_int)

# centre the moisture variable to aid interpretation
dom_data$moist_c <- scale(dom_data$moisture, center = TRUE, scale = FALSE)

# main‐effects only
mod_cont      <- lmer(resp_cube ~ gibbs + moisture + (1 | site),
                      data = dom_data)

# or with a gibbs × moisture interaction
mod_cont_int  <- lmer(resp_cube ~ gibbs * moisture + (1 | site),
                      data = dom_data)

summary(mod_cont)
summary(mod_cont_int)     
anova(mod_cont_int)
AIC(mod_cont)
AIC(mod_cont_int)
#####
# 1) load packages
library(mgcv)    # for gam()
library(gratia)  # for derivatives() & nicer plots (optional)
library(ggplot2) # for custom plotting

# 2) fit the tensor‐product GAM, REML + automatic smoothing‐penalty selection
mod_gam <- gam(
  resp_cube ~ 
    te(moisture, gibbs,    # tensor product smooth in 2D
       bs = c("tp","tp"),   # thin‐plate bases on each margin
       k  = c(10, 10)) +    # max edf per margin
    s(site, bs="re"),      # random intercept for site
  data   = dom_data,
  method = "REML",
  select = TRUE           # allows mgcv to shrink terms if not needed
)

# 3) check for basic over‐fit diagnostics
gam.check(mod_gam)
#   • look at k‐index for each smooth: should be >= 1
#   • look at QQ‐plot and residuals vs. linear predictor

# 4) visualize the fitted surface
# 4a) contour plot
vis.gam(mod_gam,
        view     = c("moisture","gibbs"),
        plot.type= "contour",
        color    = "topo",
        main     = "Contour of te(moisture, gibbs)")

# 4b) perspective plot
vis.gam(mod_gam,
        view      = c("moisture","gibbs"),
        plot.type = "persp",
        theta     = 30, phi = 30,
        main      = "3D surface of resp_cube")

# 5) estimate ∂resp/∂moisture at a fixed gibbs (say its median),
#    to see where the slope changes most rapidly
gibbs_med <- median(dom_data$gibbs)
moist_seq <- seq(min(dom_data$moisture),
                 max(dom_data$moisture), length=200)
newdat <- data.frame(moisture = moist_seq,
                     gibbs    = gibbs_med,
                     site     = NA)  # NA or any existing level; random effects average to zero

# use gratia::derivatives() to get the partial derivative w.r.t. moisture
deriv_df <- derivatives(mod_gam,
                        newdata = newdat,
                        term    = "moisture",
                        type    = "central")

# plot the derivatives
ggplot(deriv_df, aes(x = moisture, est)) +
  geom_line() +
  geom_hline(yintercept = 0, lty = 2) +
  labs(x = "moisture",
       y = "∂resp_cube/∂moisture",
       title = "Slope of resp_cube vs. moisture at median gibbs")

# you can inspect where the derivative peaks or crosses zero

#––– 0) prerequisite
library(lme4)    # for lmer()
library(nlme)    # for lme()
library(segmented)
library(nlme)
library(segmented)

#––– 1) center both predictors once and for all
dom_data$gibbs_c  <- scale(dom_data$gibbs,   center=TRUE, scale=FALSE)
dom_data$moist_c  <- scale(dom_data$moisture,center=TRUE, scale=FALSE)

# store the centering constants for back‐transformation
mean_gibbs  <- attr(dom_data$gibbs_c,  "scaled:center")
mean_moist  <- attr(dom_data$moist_c,  "scaled:center")

#––– 2) re‐fit the continuous × continuous interaction LMM with centered vars
mod_ci <- lmer(resp_cube ~ gibbs_c * moist_c + (1 | site),
               data = dom_data,
               REML = TRUE)
summary(mod_ci)
AIC(mod_ci)
# this AIC you already know was lowest of your four linear models

# 1) center moisture (and gibbs if you want) up‐front
dom_data$moist_c <- scale(dom_data$moisture, center = TRUE, scale = FALSE)

# 2) fit the base LMM (with the random intercept for site)
fm0_c <- lme(resp_cube ~ moist_c + gibbs, 
             random = ~1 | site, 
             data   = dom_data, 
             method = "REML")

# 3) choose a starting value for the breakpoint on the centered scale
psi_start <- 0   # i.e. start at mean(moisture)

# 4) call segmented() *without* a random= argument

seg_c <- segmented(
  fm0_c,
  seg.Z   = ~ moist_c,           # break on moisture
  psi     = list(moist_c = psi_start),
  control = seg.control(display=TRUE)
)



summary(seg_c)

#––– 4) extract & back‐transform the breakpoint
psi_c      <- coef(seg_c)$psi["moist_c","Est."]
se_psi_c   <- coef(seg_c)$psi["moist_c","St.Err"]
# back to original moisture units:
psi_orig   <- psi_c + mean_moist
se_psi_orig<- se_psi_c           # SE stays the same on the original scale
cat("Threshold in original moisture =", round(psi_orig,2),
    "±", round(se_psi_orig,2), "\n")

#––– 5) compare AICs if you like
AIC(seg_c)    # might work directly
# or
logLik_seg <- logLik(seg_c)
AIC_seg     <- -2*as.numeric(logLik_seg) + 2*length(logLik_seg@.Data)
cat("AIC(mod_ci)=", AIC(mod_ci),
    "  AIC(segmented)=", AIC_seg, "\n")

#––– 6) plot the segmented fit
plot(seg_c, 
     xlab="Centered moisture",
     ylab="resp_cube",
     main=sprintf("Piecewise fit, break at %.2f (orig=%.2f)", psi_c, psi_orig))
points(dom_data$moist_c, dom_data$resp_cube, pch=16, cex=0.6)
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

# ===== CROSS-VALIDATION FOR YOUR GAM MODELS =====

# Function to perform k-fold cross-validation for GAM with site random effects
cv_gam <- function(model_formula, data, k = 5, n_reps = 10) {
  set.seed(123)
  
  cv_results <- matrix(NA, nrow = n_reps, ncol = k)
  
  for(rep in 1:n_reps) {
    # Create folds - stratified by site to ensure representation
    sites <- unique(data$site)
    folds <- sample(rep(1:k, length.out = nrow(data)))
    
    for(i in 1:k) {
      tryCatch({
        # Split data
        train_data <- data[folds != i, ]
        test_data <- data[folds == i, ]
        
        # Fit model on training data
        train_model <- gam(model_formula, data = train_data, method = "REML")
        
        # Predict on test data
        predictions <- predict(train_model, newdata = test_data)
        
        # Calculate R-squared
        ss_res <- sum((test_data$resp_cube - predictions)^2, na.rm = TRUE)
        ss_tot <- sum((test_data$resp_cube - mean(test_data$resp_cube))^2, na.rm = TRUE)
        cv_results[rep, i] <- 1 - (ss_res / ss_tot)
        
      }, error = function(e) {
        cv_results[rep, i] <- NA
      })
    }
  }
  
  # Return summary statistics
  valid_results <- cv_results[!is.na(cv_results)]
  return(list(
    mean_r2 = mean(valid_results),
    sd_r2 = sd(valid_results),
    median_r2 = median(valid_results),
    all_results = cv_results
  ))
}

# ===== CROSS-VALIDATE YOUR KEY MODELS =====

formula_interaction <- resp_cube ~ te(moisture, gibbs) + s(site, bs = "re")
cv_interaction <- cv_gam(formula_interaction, dom_data, k = 5, n_reps = 10)

# Additive model for comparison
formula_additive <- resp_cube ~ s(moisture) + s(gibbs) + s(site, bs = "re")
cv_additive <- cv_gam(formula_additive, dom_data, k = 5, n_reps = 10)

# Moisture only model
formula_moisture <- resp_cube ~ s(moisture) + s(site, bs = "re")
cv_moisture <- cv_gam(formula_moisture, dom_data, k = 5, n_reps = 10)

# ===== LEAVE-ONE-SITE-OUT CROSS-VALIDATION =====
# This is important because of your site random effect

loso_cv <- function(model_formula, data) {
  sites <- unique(data$site)
  loso_results <- numeric(length(sites))
  
  for(i in 1:length(sites)) {
    tryCatch({
      # Leave one site out
      train_data <- data[data$site != sites[i], ]
      test_data <- data[data$site == sites[i], ]
      
      # Fit model
      train_model <- gam(model_formula, data = train_data, method = "REML")
      
      # Predict
      predictions <- predict(train_model, newdata = test_data)
      
      # Calculate R-squared
      ss_res <- sum((test_data$resp_cube - predictions)^2, na.rm = TRUE)
      ss_tot <- sum((test_data$resp_cube - mean(test_data$resp_cube))^2, na.rm = TRUE)
      loso_results[i] <- 1 - (ss_res / ss_tot)
      
    }, error = function(e) {
      loso_results[i] <- NA
    })
  }
  
  return(list(
    site_r2 = setNames(loso_results, sites),
    mean_r2 = mean(loso_results, na.rm = TRUE),
    sd_r2 = sd(loso_results, na.rm = TRUE)
  ))
}

# Run leave-one-site-out for interaction model
loso_interaction <- loso_cv(formula_interaction, dom_data)

# ===== COMPILE RESULTS =====

validation_results <- data.frame(
  Model = c("Interaction (te)", "Additive", "Moisture only"),
  Original_R2 = c(
    summary(cube_gam_models[["moisture_gibbs_interaction"]])$r.sq,
    summary(cube_gam_models[["moisture_plus_gibbs"]])$r.sq,
    summary(cube_gam_models[["moisture_solo"]])$r.sq
  ),
  Original_DevExpl = c(
    summary(cube_gam_models[["moisture_gibbs_interaction"]])$dev.expl,
    summary(cube_gam_models[["moisture_plus_gibbs"]])$dev.expl,
    summary(cube_gam_models[["moisture_solo"]])$dev.expl
  ),
  CV_R2_mean = c(cv_interaction$mean_r2, cv_additive$mean_r2, cv_moisture$mean_r2),
  CV_R2_sd = c(cv_interaction$sd_r2, cv_additive$sd_r2, cv_moisture$sd_r2),
  stringsAsFactors = FALSE
)

# Add LOSO results for interaction model
validation_results$LOSO_R2[1] <- loso_interaction$mean_r2
validation_results$LOSO_R2_sd[1] <- loso_interaction$sd_r2

print("=== CROSS-VALIDATION RESULTS ===")
print(validation_results)

# Export results
write.csv(validation_results, "06_GAM_Results/Comparison_Results/cross_validation_results.csv")

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


# === Compare linear mixed effects models with GAM ====
library(mgcv)  # for GAM
library(lme4)  # for lmer
library(ggplot2)
library(gridExtra)

# Define response variable
response_var <- "cube_Median_Respiration_Rate_mg_DO_per_kg_per_H"

# CORRECTED: Create formulas properly for GAM
# You need to use as.formula() or write the formula directly

# LINEAR MIXED MODEL (what we've been doing)
lmer_model <- lmer(
  cube_Median_Respiration_Rate_mg_DO_per_kg_per_H ~ Treatment + cube_Median_ATP_picomoles_per_g + cube_Median_62948_Final_Gravimetric_Moisture_g_per_g + (1|site),
  data = cube_df
)

# GAM VERSION - Linear terms (CORRECTED)
cube_df$Treatment <- as.factor(cube_df$Treatment)
cube_df$site <- as.factor(cube_df$site)

gam_linear <- gam(
  cube_Median_Respiration_Rate_mg_DO_per_kg_per_H ~ Treatment + cube_Median_62948_Final_Gravimetric_Moisture_g_per_g + s(site, bs='re'),
  data = cube_df,
  method = "REML"
)

# GAM VERSION - With smooth terms (non-linear) (CORRECTED)
gam_smooth <- gam(
  cube_Median_Respiration_Rate_mg_DO_per_kg_per_H ~ Treatment + s(cube_Median_ATP_picomoles_per_g) + s(cube_Median_62948_Final_Gravimetric_Moisture_g_per_g) + s(site, bs='re'),
  data = cube_df,
  method = "REML"
)

# Compare models
comparison_table <- data.frame(
  Model = c("LMER (Linear Mixed)", "GAM (Linear)", "GAM (Smooth)"),
  AIC = c(AIC(lmer_model), AIC(gam_linear), AIC(gam_smooth)),
  BIC = c(BIC(lmer_model), BIC(gam_linear), BIC(gam_smooth)),
  R2_marginal = c(
    r.squaredGLMM(lmer_model)[1],
    summary(gam_linear)$r.sq,
    summary(gam_smooth)$r.sq
  ),
  Deviance_Explained = c(
    NA,
    summary(gam_linear)$dev.expl,
    summary(gam_smooth)$dev.expl
  )
)

print("=== MODEL COMPARISON ===")
print(comparison_table)

# Check if smooth terms are needed
cat("\n=== GAM SMOOTH SUMMARY ===\n")
print(summary(gam_smooth))
cat("\nLook for 'edf' values > 1, which suggest non-linearity is needed\n")

# Model summaries
cat("\n=== LMER SUMMARY ===\n")
print(summary(lmer_model))

cat("\n=== GAM LINEAR SUMMARY ===\n")
print(summary(gam_linear))

# Check residual patterns
par(mfrow=c(2,3))

# LMER diagnostics
plot(lmer_model, main="LMER: Residuals vs Fitted")
qqnorm(residuals(lmer_model), main="LMER: Q-Q Plot")
qqline(residuals(lmer_model))

# GAM diagnostics
gam.check(gam_linear, main="GAM Linear")
gam.check(gam_smooth, main="GAM Smooth")

# Plot smooth functions if they exist
if(length(gam_smooth$smooth) > 0) {
  par(mfrow=c(2,2))
  plot(gam_smooth, pages=1, main="GAM Smooth Functions")
}

# Residual comparison plots
cube_df$lmer_resid <- residuals(lmer_model)
cube_df$gam_linear_resid <- residuals(gam_linear)
cube_df$gam_smooth_resid <- residuals(gam_smooth)

cube_df$lmer_fitted <- fitted(lmer_model)
cube_df$gam_linear_fitted <- fitted(gam_linear)
cube_df$gam_smooth_fitted <- fitted(gam_smooth)

# Create comparison plots
p1 <- ggplot(cube_df, aes(x=lmer_fitted, y=lmer_resid)) + 
  geom_point() + geom_smooth(se=FALSE) + 
  labs(title="LMER Residuals", x="Fitted Values", y="Residuals") +
  geom_hline(yintercept=0, linetype="dashed", color="red")

p2 <- ggplot(cube_df, aes(x=gam_linear_fitted, y=gam_linear_resid)) + 
  geom_point() + geom_smooth(se=FALSE) + 
  labs(title="GAM Linear Residuals", x="Fitted Values", y="Residuals") +
  geom_hline(yintercept=0, linetype="dashed", color="red")

p3 <- ggplot(cube_df, aes(x=gam_smooth_fitted, y=gam_smooth_resid)) + 
  geom_point() + geom_smooth(se=FALSE) + 
  labs(title="GAM Smooth Residuals", x="Fitted Values", y="Residuals") +
  geom_hline(yintercept=0, linetype="dashed", color="red")

# Fitted vs Observed
p4 <- ggplot(cube_df, aes(x=lmer_fitted, y=cube_df[[response_var]])) + 
  geom_point() + geom_abline(slope=1, intercept=0, color="red") +
  labs(title="LMER: Fitted vs Observed", x="Fitted", y="Observed")

p5 <- ggplot(cube_df, aes(x=gam_smooth_fitted, y=cube_df[[response_var]])) + 
  geom_point() + geom_abline(slope=1, intercept=0, color="red") +
  labs(title="GAM Smooth: Fitted vs Observed", x="Fitted", y="Observed")

# Display plots
grid.arrange(p1, p2, p3, ncol=3)
grid.arrange(p4, p5, ncol=2)

# Print interpretation
cat("\n=== INTERPRETATION GUIDE ===\n")
cat("1. Compare AIC/BIC: Lower is better\n")
cat("2. Compare R²: Higher is better\n") 
cat("3. GAM edf values:\n")
cat("   - edf ≈ 1: Linear relationship (LMER is fine)\n")
cat("   - edf > 1: Non-linear relationship (GAM might be better)\n")
cat("4. Residual plots should show random scatter around 0\n")
cat("5. If GAM doesn't improve much, stick with LMER for simplicity\n")