# ============================================================
# BEST-OF-BOTH WORLDS: Bootstrapped Random Forest (robust small N)
# ------------------------------------------------------------
# What this script does:
#   1) Builds effect sizes (Wet/Dry ratios) exactly like your pipeline
#   2) Uses your site_annotation_num predictors (ATP, SSA, C, N, sand, Fe, moisture)
#   3) Models BOTH outcomes as log-ratios:
#        y_gibbs  = log(effect_delGcoxPerCmol)
#        y_lambda = log(effect_Lambda)
#   4) Bootstraps times = B:
#        - TRAIN = bootstrap sample (with replacement)
#        - TEST  = out-of-bootstrap sites (not sampled in TRAIN)  [external-ish validation]
#        - Also computes OOB metrics from ranger on TRAIN          [stable internal check]
#        - Median imputation for predictors using TRAIN medians (no leakage)
#   5) Returns:
#        - OOB RMSE / OOB R2 distributions
#        - TEST RMSE / TEST R2 distributions (when n_test >= min_test_n)
#        - Bootstrapped permutation importance + 95% CI
#   6) Saves CSVs + PDF plots per outcome
#
# Why “best-of-both”:
#   - OOB metrics are stable even when N is tiny
#   - Out-of-bootstrap TEST metrics approximate generalization but can be noisy;
#     we gate them with min_test_n so we don’t overinterpret tiny test sets.
# ============================================================

rm(list=ls(all=T))

suppressPackageStartupMessages({
  library(tidyverse)
  library(ranger)
  library(ggplot2)
})

set.seed(123)

# ----------------------------
# User settings
# ----------------------------
figure_path <- "Figures/"
dir.create(figure_path, showWarnings = FALSE, recursive = TRUE)

B <- 1000             # bootstrap replicates (e.g., 500–2000)
TREES <- 5000         # trees per forest (more = more stable importance)
MIN_NODE <- 5         # a bit conservative for small N
MIN_TEST_N <- 5       # only compute "test" metrics if test set has >= this many sites
SEED_BASE <- 123

# ----------------------------
# 1) Load / prep your data
# ----------------------------
data <- utils::read.csv("Data/Medians_of_Median_molecular_properties_per_site_and_treatment_unique_formulas.csv")
row.names(data) <- paste0(data$site, "_", data$Treatment)

# You can keep these if needed elsewhere, but predictors come from sample_data below
field_metadata <- utils::read.csv("EC_Data_Package/EC_Field_Metadata.csv") %>%
  dplyr::select(site = Parent_ID, State)

ecoregions <- utils::read.csv("Data/EC_Ecoregions.csv") %>%
  dplyr::select(site = Parent_ID, Ecoregion = Name) %>%
  dplyr::mutate(Ecoregion = stringr::str_trim(stringr::str_extract(Ecoregion, "(?<=- ).*")))

field_metadata <- base::merge(field_metadata, ecoregions, by = "site")

# Sample data providing the site_annotation_num predictors
sample_data <- readr::read_csv(
  "v4_CM_SSS_Data_Package_1/Sample_Data/v3_CM_SSS_Sediment_Sample_Data_Summary.csv",
  comment = "#",
  na = c("N/A", "-9999", -9999)
) %>%
  dplyr::slice(-(1:11)) %>%
  dplyr::mutate(dplyr::across(-c(Sample_Name, Field_Name, IGSN, Material), as.numeric)) %>%
  dplyr::select(
    site = Sample_Name,
    Mean_ATP_picomoles_per_g,
    Mean_Specific_Surface_Area_m2_per_g,
    C_percent_per_mg = `01395_C_percent_per_mg`,
    N_percent_per_mg = `01397_N_percent_per_mg`,
    Percent_Tot_Sand,
    Mean_Fe_mg_per_kg,
    Mean_Gravimetric_Moisture_g_per_g = Mean_62948_Gravimetric_Moisture_g_per_g
  )

# Standardize site ids like your code
sample_data$site <- gsub("CM", "EC", sample_data$site)
sample_data$site <- gsub("_Sediment", "", sample_data$site)

# DOM subset
dom_data <- data %>%
  dplyr::filter(site != "EC_023") %>%
  dplyr::filter(!(site %in% c("EC_012", "EC_011", "EC_053", "EC_057", "EC_052"))) %>%
  dplyr::select(site, Treatment, Median_delGcoxPerCmol, Median_Lambda)

sites <- dom_data %>% dplyr::distinct(site)

# merge metadata just to mirror your pipeline (not used in predictors)
site_metadata <- base::merge(sample_data, field_metadata, by = "site")
site_metadata <- base::merge(site_metadata, sites, by = "site")

# ----------------------------
# 2) Compute treatment effect sizes (Wet/Dry ratios)
# ----------------------------
treatment_effects <- dom_data %>%
  dplyr::filter(Treatment %in% c("Wet", "Dry")) %>%
  dplyr::select(site, Treatment, dplyr::starts_with("Median_")) %>%
  tidyr::pivot_wider(
    id_cols = site,
    names_from = Treatment,
    values_from = dplyr::starts_with("Median_")
  ) %>%
  dplyr::mutate(
    effect_delGcoxPerCmol = (Median_delGcoxPerCmol_Wet / Median_delGcoxPerCmol_Dry),
    effect_Lambda         = (Median_Lambda_Wet / Median_Lambda_Dry)
  ) %>%
  dplyr::select(site, effect_delGcoxPerCmol, effect_Lambda)

# ----------------------------
# 3) Predictors (your requested site_annotation_num)
# ----------------------------
site_annotation_num <- site_metadata %>%
  dplyr::select(
    site,
    Mean_ATP_picomoles_per_g,
    Mean_Specific_Surface_Area_m2_per_g,
    C_percent_per_mg,
    N_percent_per_mg,
    Percent_Tot_Sand,
    Mean_Fe_mg_per_kg,
    Mean_Gravimetric_Moisture_g_per_g
  )

predictors <- c(
  "Mean_ATP_picomoles_per_g",
  "Mean_Specific_Surface_Area_m2_per_g",
  "C_percent_per_mg",
  "N_percent_per_mg",
  "Percent_Tot_Sand",
  "Mean_Fe_mg_per_kg",
  "Mean_Gravimetric_Moisture_g_per_g"
)

# ----------------------------
# 4) Modeling df + outcomes (log ratios)
# ----------------------------
model_df <- treatment_effects %>%
  dplyr::left_join(site_annotation_num, by = "site") %>%
  dplyr::mutate(
    effect_delGcoxPerCmol = as.numeric(effect_delGcoxPerCmol),
    effect_Lambda         = as.numeric(effect_Lambda),
    dplyr::across(dplyr::all_of(predictors), as.numeric),
    y_gibbs  = ifelse(effect_delGcoxPerCmol > 0, log(effect_delGcoxPerCmol), NA_real_),
    y_lambda = ifelse(effect_Lambda > 0,         log(effect_Lambda),         NA_real_)
  )

cat("\n===== DATA CHECK =====\n")
cat("Rows (sites):", nrow(model_df), "\n")
cat("Complete cases (gibbs): ",  sum(complete.cases(model_df[, c("y_gibbs",  predictors)])), "\n")
cat("Complete cases (lambda): ", sum(complete.cases(model_df[, c("y_lambda", predictors)])), "\n\n")

# ----------------------------
# 5) Metrics helpers
# ----------------------------
rmse <- function(y_true, y_pred) {
  sqrt(mean((y_true - y_pred)^2, na.rm = TRUE))
}

r2_sse <- function(y_true, y_pred) {
  ok <- is.finite(y_true) & is.finite(y_pred)
  y_true <- y_true[ok]; y_pred <- y_pred[ok]
  if (length(y_true) < 2) return(NA_real_)
  if (!is.finite(sd(y_true)) || sd(y_true) == 0) return(NA_real_)
  sse <- sum((y_true - y_pred)^2)
  sst <- sum((y_true - mean(y_true))^2)
  1 - (sse / sst)
}

# ----------------------------
# 6) Bootstrapped RF (OOB + out-of-bootstrap TEST)
# ----------------------------
rf_bootstrap_best <- function(df, outcome, predictors,
                              times = 1000, seed = 1,
                              trees = 5000, mtry = NULL,
                              min_node_size = 5,
                              min_test_n = 5) {
  
  set.seed(seed)
  
  dat <- df %>%
    dplyr::select(site, dplyr::all_of(c(outcome, predictors))) %>%
    dplyr::mutate(dplyr::across(-site, as.numeric)) %>%
    dplyr::filter(!is.na(.data[[outcome]]))
  
  # drop predictors that are constant in the full data (helps stability)
  sds <- sapply(dat[predictors], stats::sd, na.rm = TRUE)
  predictors_use <- predictors[is.finite(sds) & sds > 0]
  if (length(predictors_use) < 2) stop("Too few non-constant predictors.")
  
  n <- nrow(dat)
  p <- length(predictors_use)
  if (is.null(mtry)) mtry <- max(1, floor(sqrt(p)))
  
  metrics_list <- vector("list", times)
  imp_list <- vector("list", times)
  
  for (b in seq_len(times)) {
    
    train_idx <- sample.int(n, size = n, replace = TRUE)
    test_idx  <- setdiff(seq_len(n), unique(train_idx))
    
    train <- dat[train_idx, , drop = FALSE]
    test  <- if (length(test_idx) > 0) dat[test_idx, , drop = FALSE] else NULL
    
    # Median impute predictors using TRAIN medians (no leakage)
    med <- sapply(train[, predictors_use, drop = FALSE], function(x) stats::median(x, na.rm = TRUE))
    for (v in predictors_use) {
      train[[v]][is.na(train[[v]])] <- med[[v]]
      if (!is.null(test)) test[[v]][is.na(test[[v]])] <- med[[v]]
    }
    
    # If training outcome has 0 variance, skip
    if (!is.finite(stats::sd(train[[outcome]])) || stats::sd(train[[outcome]]) == 0) {
      metrics_list[[b]] <- tibble::tibble(
        boot_id = b,
        n_train = nrow(train),
        n_test  = ifelse(is.null(test), 0, nrow(test)),
        oob_rmse = NA_real_, oob_r2 = NA_real_,
        test_rmse = NA_real_, test_r2 = NA_real_
      )
      imp_list[[b]] <- tibble::tibble(
        boot_id = b,
        term = predictors_use,
        importance = NA_real_
      )
      next
    }
    
    fit <- ranger::ranger(
      dependent.variable.name = outcome,
      data = train %>% dplyr::select(-site),
      num.trees = trees,
      mtry = mtry,
      min.node.size = min_node_size,
      importance = "permutation",
      oob.error = TRUE,
      seed = seed + b
    )
    
    # OOB metrics from ranger (stable)
    oob_rmse <- sqrt(as.numeric(fit$prediction.error))  # OOB MSE -> RMSE
    oob_r2   <- as.numeric(fit$r.squared)
    
    # Importance
    imp_vec <- unlist(fit$variable.importance, use.names = TRUE)
    imp_list[[b]] <- tibble::tibble(
      boot_id = b,
      term = names(imp_vec),
      importance = as.numeric(imp_vec)
    )
    
    # TEST metrics (out-of-bootstrap), only if test set is large enough
    test_rmse_val <- NA_real_
    test_r2_val   <- NA_real_
    
    if (!is.null(test) && nrow(test) >= min_test_n) {
      pred <- predict(fit, data = test[, predictors_use, drop = FALSE])$predictions
      y_true <- test[[outcome]]
      test_rmse_val <- rmse(y_true, pred)
      test_r2_val   <- r2_sse(y_true, pred)
    }
    
    metrics_list[[b]] <- tibble::tibble(
      boot_id = b,
      n_train = nrow(train),
      n_test  = ifelse(is.null(test), 0, nrow(test)),
      oob_rmse = oob_rmse,
      oob_r2   = oob_r2,
      test_rmse = test_rmse_val,
      test_r2   = test_r2_val
    )
  }
  
  metrics <- dplyr::bind_rows(metrics_list)
  imps <- dplyr::bind_rows(imp_list)
  
  imp_summary <- imps %>%
    dplyr::group_by(term) %>%
    dplyr::summarise(
      median_importance = stats::median(importance, na.rm = TRUE),
      mean_importance   = mean(importance, na.rm = TRUE),
      lo95 = stats::quantile(importance, 0.025, na.rm = TRUE),
      hi95 = stats::quantile(importance, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(median_importance))
  
  # Summaries (OOB always; TEST only on finite values)
  perf_summary <- tibble::tibble(
    outcome = outcome,
    n = nrow(dat),
    p_used = length(predictors_use),
    
    oob_r2_med = stats::median(metrics$oob_r2, na.rm = TRUE),
    oob_r2_lo  = stats::quantile(metrics$oob_r2, 0.025, na.rm = TRUE),
    oob_r2_hi  = stats::quantile(metrics$oob_r2, 0.975, na.rm = TRUE),
    
    oob_rmse_med = stats::median(metrics$oob_rmse, na.rm = TRUE),
    oob_rmse_lo  = stats::quantile(metrics$oob_rmse, 0.025, na.rm = TRUE),
    oob_rmse_hi  = stats::quantile(metrics$oob_rmse, 0.975, na.rm = TRUE),
    
    test_r2_med = stats::median(metrics$test_r2, na.rm = TRUE),
    test_r2_lo  = stats::quantile(metrics$test_r2, 0.025, na.rm = TRUE),
    test_r2_hi  = stats::quantile(metrics$test_r2, 0.975, na.rm = TRUE),
    
    test_rmse_med = stats::median(metrics$test_rmse, na.rm = TRUE),
    test_rmse_lo  = stats::quantile(metrics$test_rmse, 0.025, na.rm = TRUE),
    test_rmse_hi  = stats::quantile(metrics$test_rmse, 0.975, na.rm = TRUE),
    
    test_n_finite = sum(is.finite(metrics$test_r2) & is.finite(metrics$test_rmse))
  )
  
  list(
    outcome = outcome,
    predictors_used = predictors_use,
    metrics = metrics,
    importance_summary = imp_summary,
    perf_summary = perf_summary,
    settings = list(times=times, trees=trees, mtry=mtry, min_node_size=min_node_size, min_test_n=min_test_n)
  )
}

# ----------------------------
# 7) Run models (Gibbs + Lambda)
# ----------------------------
res_gibbs  <- rf_bootstrap_best(model_df, "y_gibbs",  predictors,
                                times = B, seed = SEED_BASE + 101,
                                trees = TREES, min_node_size = MIN_NODE,
                                min_test_n = MIN_TEST_N)

res_lambda <- rf_bootstrap_best(model_df, "y_lambda", predictors,
                                times = B, seed = SEED_BASE + 202,
                                trees = TREES, min_node_size = MIN_NODE,
                                min_test_n = MIN_TEST_N)

cat("\n===== PERFORMANCE SUMMARY (OOB + TEST when test >= MIN_TEST_N) =====\n")
print(res_gibbs$perf_summary)
print(res_lambda$perf_summary)

# ----------------------------
# 8) Plots
# ----------------------------
plot_importance <- function(res, top_n = 20) {
  res$importance_summary %>%
    dplyr::slice_head(n = top_n) %>%
    ggplot2::ggplot(ggplot2::aes(x = reorder(term, median_importance), y = median_importance)) +
    ggplot2::geom_point() +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lo95, ymax = hi95), width = 0.0) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = NULL,
      y = "Permutation importance (bootstrapped median, 95% CI)",
      title = paste("Bootstrapped RF importance:", res$outcome)
    )
}

plot_metric_box <- function(res, cols = c("oob_r2","oob_rmse","test_r2","test_rmse")) {
  res$metrics %>%
    dplyr::select(dplyr::all_of(cols)) %>%
    tidyr::pivot_longer(dplyr::everything(), names_to = "metric", values_to = "value") %>%
    dplyr::filter(is.finite(value)) %>%
    ggplot2::ggplot(ggplot2::aes(x = metric, y = value)) +
    ggplot2::geom_boxplot(outlier.alpha = 0.25) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = paste("Bootstrap performance distributions:", res$outcome),
      x = NULL, y = NULL
    )
}

p_imp_g <- plot_importance(res_gibbs, top_n = 20)
p_imp_l <- plot_importance(res_lambda, top_n = 20)

p_perf_g <- plot_metric_box(res_gibbs)
p_perf_l <- plot_metric_box(res_lambda)

ggplot2::ggsave(paste0(figure_path, "RF_best_importance_y_gibbs.pdf"),  p_imp_g, width = 7, height = 5)
ggplot2::ggsave(paste0(figure_path, "RF_best_importance_y_lambda.pdf"), p_imp_l, width = 7, height = 5)
ggplot2::ggsave(paste0(figure_path, "RF_best_perf_y_gibbs.pdf"),        p_perf_g, width = 7, height = 4)
ggplot2::ggsave(paste0(figure_path, "RF_best_perf_y_lambda.pdf"),       p_perf_l, width = 7, height = 4)

# ----------------------------
# 9) Save outputs
# ----------------------------
dir.create("rf_bootstrap_outputs_best", showWarnings = FALSE)

readr::write_csv(model_df, "rf_bootstrap_outputs_best/model_df_for_sharing.csv")

readr::write_csv(res_gibbs$importance_summary,  "rf_bootstrap_outputs_best/importance_summary_y_gibbs.csv")
readr::write_csv(res_lambda$importance_summary, "rf_bootstrap_outputs_best/importance_summary_y_lambda.csv")

readr::write_csv(res_gibbs$metrics,  "rf_bootstrap_outputs_best/metrics_y_gibbs.csv")
readr::write_csv(res_lambda$metrics, "rf_bootstrap_outputs_best/metrics_y_lambda.csv")

readr::write_csv(res_gibbs$perf_summary,  "rf_bootstrap_outputs_best/perf_summary_y_gibbs.csv")
readr::write_csv(res_lambda$perf_summary, "rf_bootstrap_outputs_best/perf_summary_y_lambda.csv")

saveRDS(list(gibbs=res_gibbs, lambda=res_lambda),
        file = paste0(figure_path, "RF_bootstrap_best_all_results.rds"))

cat("\nSaved outputs to:\n")
cat("  Figures/ (PDFs + RDS)\n")
cat("  rf_bootstrap_outputs_best/ (CSVs)\n")
cat("DONE.\n")
