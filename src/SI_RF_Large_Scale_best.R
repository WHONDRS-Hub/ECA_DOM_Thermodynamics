# ============================================================
# BEST-OF-BOTH WORLDS: Bootstrapped Random Forest (Large-scale predictors)
# ------------------------------------------------------------
# Predictors = ALL numeric "large-scale" variables from EC_Field_Geospatial_NPP_ET.csv
# Targets (ratios Wet/Dry), modeled as log-ratios:
#   y_gibbs  = log(effect_delGcoxPerCmol)
#   y_lambda = log(effect_Lambda)
#
# Bootstrap procedure (robust small N):
#   - TRAIN = bootstrap sample (with replacement)
#   - TEST  = out-of-bootstrap rows (not sampled)  [more conservative generalization]
#   - Also compute ranger OOB metrics on TRAIN     [stable internal check]
#   - Median imputation of predictors using TRAIN medians (no leakage)
#   - Drop zero-variance predictors
#   - Only compute TEST metrics if n_test >= MIN_TEST_N
#
# Outputs:
#   Figures/
#     RF_LS_best_perf_summary_y_gibbs.csv
#     RF_LS_best_perf_summary_y_lambda.csv
#     RF_LS_best_metrics_y_gibbs.csv
#     RF_LS_best_metrics_y_lambda.csv
#     RF_LS_best_importance_summary_y_gibbs.csv
#     RF_LS_best_importance_summary_y_lambda.csv
#     RF_LS_best_importance_y_gibbs.pdf
#     RF_LS_best_importance_y_lambda.pdf
#     RF_LS_best_perf_y_gibbs.pdf
#     RF_LS_best_perf_y_lambda.pdf
#     RF_LS_best_model_df_for_sharing.csv
#     RF_LS_best_all_results.rds
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
TREES <- 5000         # trees per forest
MIN_NODE <- 5         # conservative for small N
MIN_TEST_N <- 5       # compute TEST metrics only if test set has >= this many sites
SEED_BASE <- 123

# ----------------------------
# 1) Load / prep your data
# ----------------------------
data <- utils::read.csv("Data/Medians_of_Median_molecular_properties_per_site_and_treatment_unique_formulas.csv")
row.names(data) <- paste0(data$site, "_", data$Treatment)

large_data <- utils::read.csv("Data/EC_Field_Geospatial_NPP_PET_VGC.csv")

# match your earlier naming swap (site/code columns)
names(large_data)[2] <- "code"
names(large_data)[1] <- "site"

# DOM subset (same filters as your pipeline)
dom_data <- data %>%
  dplyr::filter(site != "EC_023") %>%
  dplyr::filter(!(site %in% c("EC_012","EC_011","EC_053","EC_057","EC_052"))) %>%
  dplyr::select(site, Treatment, Median_delGcoxPerCmol, Median_Lambda)

sites <- dom_data %>% dplyr::distinct(site)

# Explanatory data aligned to DOM sites
explanatory_data <- base::merge(sites, large_data, by = "site") %>%
  dplyr::distinct()

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
# 3) Build modeling dataframe (large-scale predictors only)
# ----------------------------
# Predictors = all numeric cols in explanatory_data excluding code
predictors_all <- explanatory_data %>%
  dplyr::select(dplyr::where(is.numeric)) %>%
  dplyr::select(-dplyr::any_of(c("code"))) %>%
  names()

# Model df: join predictors + effects; use log ratios as outcomes
model_df <- treatment_effects %>%
  dplyr::left_join(explanatory_data, by = "site") %>%
  dplyr::mutate(
    effect_delGcoxPerCmol = as.numeric(effect_delGcoxPerCmol),
    effect_Lambda         = as.numeric(effect_Lambda),
    dplyr::across(dplyr::all_of(predictors_all), as.numeric),
    y_gibbs  = ifelse(effect_delGcoxPerCmol > 0, log(effect_delGcoxPerCmol), NA_real_),
    y_lambda = ifelse(effect_Lambda > 0,         log(effect_Lambda),         NA_real_)
  )

# remove non-predictor numeric columns that might sneak in
# (explanatory_data has code removed; site is non-numeric so ok)
predictors <- predictors_all

cat("\n===== DATA CHECK (Large-scale predictors) =====\n")
cat("Rows (sites):", nrow(model_df), "\n")
cat("Num predictors (before filtering):", length(predictors), "\n")
cat("Complete cases (gibbs): ",  sum(complete.cases(model_df[, c("y_gibbs",  predictors)])), "\n")
cat("Complete cases (lambda): ", sum(complete.cases(model_df[, c("y_lambda", predictors)])), "\n\n")

# ----------------------------
# 4) Metrics helpers
# ----------------------------
rmse <- function(y_true, y_pred) sqrt(mean((y_true - y_pred)^2, na.rm = TRUE))

r2_sse <- function(y_true, y_pred) {
  ok <- is.finite(y_true) & is.finite(y_pred)
  y_true <- y_true[ok]; y_pred <- y_pred[ok]
  if (length(y_true) < 2) return(NA_real_)
  if (!is.finite(stats::sd(y_true)) || stats::sd(y_true) == 0) return(NA_real_)
  sse <- sum((y_true - y_pred)^2)
  sst <- sum((y_true - mean(y_true))^2)
  1 - (sse / sst)
}

# ----------------------------
# 5) Bootstrapped RF (OOB + out-of-bootstrap TEST)
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
  
  # Drop predictors constant in the full dataset
  sds <- sapply(dat[predictors], stats::sd, na.rm = TRUE)
  predictors_use <- predictors[is.finite(sds) & sds > 0]
  if (length(predictors_use) < 2) stop(paste0("Too few non-constant predictors for ", outcome))
  
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
    
    # Skip if training outcome has 0 variance
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
    
    # OOB metrics (stable)
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
# 6) Run models (Gibbs + Lambda)
# ----------------------------
res_gibbs  <- rf_bootstrap_best(model_df, "y_gibbs",  predictors,
                                times = B, seed = SEED_BASE + 101,
                                trees = TREES, min_node_size = MIN_NODE,
                                min_test_n = MIN_TEST_N)

res_lambda <- rf_bootstrap_best(model_df, "y_lambda", predictors,
                                times = B, seed = SEED_BASE + 202,
                                trees = TREES, min_node_size = MIN_NODE,
                                min_test_n = MIN_TEST_N)

cat("\n===== PERFORMANCE SUMMARY (Large-scale predictors) =====\n")
print(res_gibbs$perf_summary)
print(res_lambda$perf_summary)

# ----------------------------
# 7) Plots
# ----------------------------
plot_importance <- function(res, top_n = 20) {
  res$importance_summary %>%
    dplyr::slice_head(n = min(top_n, nrow(res$importance_summary))) %>%
    ggplot2::ggplot(ggplot2::aes(x = reorder(term, median_importance), y = median_importance)) +
    ggplot2::geom_point() +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lo95, ymax = hi95), width = 0.0) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = NULL,
      y = "Permutation importance (bootstrapped median, 95% CI)",
      title = paste("Bootstrapped RF importance:", res$outcome, "(large-scale predictors)")
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
      title = paste("Bootstrap performance distributions:", res$outcome, "(large-scale predictors)"),
      x = NULL, y = NULL
    )
}

p_imp_g <- plot_importance(res_gibbs, top_n = 20)
p_imp_l <- plot_importance(res_lambda, top_n = 20)

p_perf_g <- plot_metric_box(res_gibbs)
p_perf_l <- plot_metric_box(res_lambda)

ggplot2::ggsave(paste0(figure_path, "RF_LS_best_importance_y_gibbs.pdf"),  p_imp_g, width = 8, height = 6)
ggplot2::ggsave(paste0(figure_path, "RF_LS_best_importance_y_lambda.pdf"), p_imp_l, width = 8, height = 6)
ggplot2::ggsave(paste0(figure_path, "RF_LS_best_perf_y_gibbs.pdf"),        p_perf_g, width = 8, height = 5)
ggplot2::ggsave(paste0(figure_path, "RF_LS_best_perf_y_lambda.pdf"),       p_perf_l, width = 8, height = 5)

# ----------------------------
# 8) Save outputs
# ----------------------------
utils::write.csv(model_df, paste0(figure_path, "RF_LS_best_model_df_for_sharing.csv"), row.names = FALSE)

utils::write.csv(res_gibbs$perf_summary,  paste0(figure_path, "RF_LS_best_perf_summary_y_gibbs.csv"), row.names = FALSE)
utils::write.csv(res_lambda$perf_summary, paste0(figure_path, "RF_LS_best_perf_summary_y_lambda.csv"), row.names = FALSE)

utils::write.csv(res_gibbs$metrics,  paste0(figure_path, "RF_LS_best_metrics_y_gibbs.csv"), row.names = FALSE)
utils::write.csv(res_lambda$metrics, paste0(figure_path, "RF_LS_best_metrics_y_lambda.csv"), row.names = FALSE)

utils::write.csv(res_gibbs$importance_summary,  paste0(figure_path, "RF_LS_best_importance_summary_y_gibbs.csv"), row.names = FALSE)
utils::write.csv(res_lambda$importance_summary, paste0(figure_path, "RF_LS_best_importance_summary_y_lambda.csv"), row.names = FALSE)

saveRDS(list(gibbs=res_gibbs, lambda=res_lambda),
        file = paste0(figure_path, "RF_LS_best_all_results.rds"))

cat("\nSaved outputs to Figures/:\n")
cat("  RF_LS_best_* (CSVs, PDFs, RDS)\n")
cat("DONE.\n")
