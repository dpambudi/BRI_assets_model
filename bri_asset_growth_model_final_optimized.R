# ============================================================================
# UNIFIED ANALYSIS: BRI ASSET GROWTH & DIAGNOSTIC FRAMEWORK
# ============================================================================
# Diagnostic & Comparative Evaluation of Bounded-Growth (Sigmoid) vs 
# Unbounded Trend Models (Linear, Quadratic, Exponential) with Full Residual
# Diagnostics, 11-Step Incremental OOS Validation, and OOS Parameter Tracking.
# ============================================================================

# --- PART 0: ENVIRONMENT SETUP & LIBRARIES ---
rm(list = ls())
graphics.off()

if (basename(Sys.getenv("RSTUDIO_USER_IDENTITY")) != "") {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

required_packages <- c("ggplot2", "minpack.lm", "stats", "patchwork")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# ============================================================================
# PART 1: DATA INGESTION & SETUP
# ============================================================================
if (!file.exists("DataBRI.csv")) {
  stop("Critical Error: 'DataBRI.csv' not found in working directory.")
}

df <- read.csv("DataBRI.csv", sep = ";")
df <- na.omit(df[, c("t", "Year", "Q", "Asset")])
df$Label <- paste0(df$Year, "-", df$Q)

# NLS Optimization Settings
control_nls <- nls.lm.control(maxiter = 1024, maxfev = 50000)

# Master Fitting Function for All 6 Models
fit_all_models <- function(data) {
  max_asset <- max(data$Asset)
  min_asset <- min(data$Asset)
  
  models <- list(
    Linear    = lm(Asset ~ t, data = data),
    Quadratic = lm(Asset ~ t + I(t^2), data = data),
    Exp       = nlsLM(Asset ~ A0 * exp(r * t), data = data,
                      start = list(A0 = min_asset, r = 0.04), control = control_nls),
    Log       = nlsLM(Asset ~ K / (1 + C * exp(-r * t)), data = data,
                      start = list(K = max_asset * 1.5, C = 5, r = 0.04), control = control_nls),
    Gom       = nlsLM(Asset ~ K * exp(-B * exp(-k_g * t)), data = data,
                      start = list(K = max_asset * 1.5, B = 5, k_g = 0.04), control = control_nls),
    Ric       = nlsLM(Asset ~ K / (1 + C * exp(-r * t))^(1/v), data = data,
                      start = list(K = max_asset * 1.5, C = 5, r = 0.04, v = 1), control = control_nls)
  )
  return(models)
}

# Statistical Error Metrics Helper
metric_fun <- function(actual, predicted) {
  c(MAE = mean(abs(actual - predicted)),
    RMSE = sqrt(mean((actual - predicted)^2)),
    MAPE = mean(abs((actual - predicted) / actual)) * 100)
}

# Helper to Extract Optimization Iterations and Convergence Status
get_conv_info <- function(model) {
  if (inherits(model, "lm")) {
    return(list(iters = 1, status = "Linear OLS (Closed-form)"))
  } else if (inherits(model, "nls")) {
    iters <- model$convInfo$finIter
    is_conv <- model$convInfo$isConv
    status_msg <- ifelse(is_conv, "Converged", paste("Not Converged:", model$convInfo$stopMessage))
    return(list(iters = iters, status = status_msg))
  } else {
    return(list(iters = NA, status = "Unknown"))
  }
}

# Parameter Extraction Helper (Includes Iteration & Convergence Tracking)
extract_params <- function(model, name) {
  coef_sum <- summary(model)$coefficients
  ci <- tryCatch(confint(model), error = function(e) matrix(NA, nrow = nrow(coef_sum), ncol = 2))
  if (name == "Richards") ci <- matrix(NA, nrow = nrow(coef_sum), ncol = 2) # Override Richards CI if unstable
  
  conv_info <- get_conv_info(model)
  
  data.frame(
    Model = name,
    Parameter = rownames(coef_sum),
    Estimate = coef_sum[, "Estimate"],
    Std_Error = coef_sum[, "Std. Error"],
    t_value = coef_sum[, "t value"],
    p_value = coef_sum[, "Pr(>|t|)"],
    CI_Lower = ci[, 1],
    CI_Upper = ci[, 2],
    Iterations = conv_info$iters,
    Convergence_Status = conv_info$status
  )
}

model_names <- c("Linear", "Quadratic", "Exponential", "Verhulst", "Gompertz", "Richards")

# ============================================================================
# PART 2: INCREMENTAL (1-SAMPLE STEP) OOS VALIDATION & OOS PARAMETER LOGGING
# ============================================================================
# Training N ranges from 38 to 48 (11 iterations), Testing M ranges from 15 down to 5
val_results <- data.frame()
oos_param_logs <- data.frame()

for (tr_size in 38:48) {
  ts_start <- tr_size + 1
  te_end   <- 53
  ts_length <- te_end - tr_size
  
  df_train <- subset(df, t <= tr_size)
  df_test  <- subset(df, t >= ts_start & t <= te_end)
  
  models_oos <- fit_all_models(df_train)
  
  preds <- list(
    predict(models_oos$Linear, newdata = df_test),
    predict(models_oos$Quadratic, newdata = df_test),
    predict(models_oos$Exp, newdata = df_test),
    predict(models_oos$Log, newdata = df_test),
    predict(models_oos$Gom, newdata = df_test),
    predict(models_oos$Ric, newdata = df_test)
  )
  
  # 1. Log OOS Forecast Accuracy Metrics
  for (j in seq_along(model_names)) {
    m <- metric_fun(df_test$Asset, preds[[j]])
    val_results <- rbind(val_results, data.frame(
      Training_N = tr_size,
      Testing_M = ts_length,
      Split = paste0(tr_size, "/", ts_length),
      Model = model_names[j], 
      MAE = m["MAE"], 
      RMSE = m["RMSE"], 
      MAPE = m["MAPE"]
    ))
    
    # 2. Log Model Parameter Estimates for Every OOS Iteration
    p_df <- extract_params(models_oos[[j]], model_names[j])
    p_df$Training_N <- tr_size
    p_df$Testing_M  <- ts_length
    p_df$Split      <- paste0(tr_size, "/", ts_length)
    oos_param_logs  <- rbind(oos_param_logs, p_df)
  }
}

# Write Full Incremental OOS Metrics & OOS Parameter Estimates
write.csv(val_results, "Export_MultiWindow_Validation.csv", row.names = FALSE)
write.csv(oos_param_logs, "Export_OOS_Parameter_Estimates_All_Iter.csv", row.names = FALSE)

# Compute Aggregated OOS Validation Summary (Mean & Median Across 11 Windows)
val_summary <- aggregate(cbind(MAE, RMSE, MAPE) ~ Model, data = val_results, mean)
val_summary <- val_summary[order(val_summary$MAPE), ]
write.csv(val_summary, "Export_MultiWindow_Validation_Summary.csv", row.names = FALSE)


# ============================================================================
# PART 3: RECURSIVE PARAMETER STABILITY TEST (TRACKING K & t* FOR ALL SIGMOIDS)
# ============================================================================
rec_results <- data.frame()

for (n_sample in 38:53) {
  df_sub <- subset(df, t <= n_sample)
  m_sub <- fit_all_models(df_sub)
  
  # 1. Verhulst Parameters & Inflection Point (t*)
  K_log <- coef(m_sub$Log)["K"]
  C_log <- coef(m_sub$Log)["C"]
  r_log <- coef(m_sub$Log)["r"]
  t_star_log <- log(C_log) / r_log
  
  # 2. Gompertz Parameters & Inflection Point (t*)
  K_gom <- coef(m_sub$Gom)["K"]
  B_gom <- coef(m_sub$Gom)["B"]
  kg_gom <- coef(m_sub$Gom)["k_g"]
  t_star_gom <- log(B_gom) / kg_gom
  
  # 3. Richards Parameters & Inflection Point (t*)
  K_ric <- coef(m_sub$Ric)["K"]
  C_ric <- coef(m_sub$Ric)["C"]
  r_ric <- coef(m_sub$Ric)["r"]
  v_ric <- coef(m_sub$Ric)["v"]
  
  # Analytical Formula for Richards Inflection Point: ln(C / v) / r
  t_star_ric <- tryCatch({
    val <- log(C_ric / v_ric) / r_ric
    if (is.nan(val) || is.infinite(val) || val > 150 || val < 0) NA else val
  }, error = function(e) NA)
  
  rec_results <- rbind(rec_results, data.frame(
    Sample_Size = n_sample,
    Cutoff_Period = df_sub$Label[nrow(df_sub)],
    Verhulst_K = K_log,
    Verhulst_t_star = t_star_log,
    Gompertz_K = K_gom,
    Gompertz_t_star = t_star_gom,
    Richards_K = K_ric,
    Richards_t_star = t_star_ric
  ))
}

write.csv(rec_results, "Export_Recursive_Parameter_Stability.csv", row.names = FALSE)


# ============================================================================
# PART 4: FULL SAMPLE CALIBRATION (N = 53) & SELECTION METRICS
# ============================================================================
full_models <- fit_all_models(df)

# Parameter Estimates Export for Full Sample
param_table_full <- do.call(rbind, mapply(extract_params, full_models, model_names, SIMPLIFY = FALSE))
write.csv(param_table_full, "Export_Model_Parameter_Estimates.csv", row.names = FALSE)

# Information Criteria Assessment
rss_vals <- sapply(full_models, function(m) sum(residuals(m)^2))
aic_vals <- sapply(full_models, AIC)
bic_vals <- sapply(full_models, BIC)
k_params <- c(2, 3, 2, 3, 3, 4)

# Track Iterations & Status for Selection Table Summary
iters_vals <- sapply(full_models, function(m) get_conv_info(m)$iters)
conv_status_vals <- sapply(full_models, function(m) get_conv_info(m)$status)

metrics_ic <- data.frame(
  Model = model_names, 
  Parameters = k_params,
  Iterations = iters_vals,
  Convergence = conv_status_vals,
  RSS = rss_vals, 
  AIC = aic_vals, 
  BIC = bic_vals
)
metrics_ic$Delta_AIC <- metrics_ic$AIC - min(metrics_ic$AIC)
metrics_ic$Delta_BIC <- metrics_ic$BIC - min(metrics_ic$BIC)
metrics_ic$Akaike_Weight <- exp(-0.5 * metrics_ic$Delta_AIC) / sum(exp(-0.5 * metrics_ic$Delta_AIC))
metrics_ic <- metrics_ic[order(metrics_ic$AIC), ]

write.csv(metrics_ic, "Export_AIC_BIC_Metrics.csv", row.names = FALSE)


# ============================================================================
# PART 5: RESIDUAL DIAGNOSTICS & SEASONALITY TESTS (ALL MODELS)
# ============================================================================
resids_list <- lapply(full_models, residuals)
names(resids_list) <- model_names

diag_summary <- data.frame()
for (i in seq_along(model_names)) {
  m_name <- model_names[i]
  r <- resids_list[[i]]
  
  shapiro_p  <- shapiro.test(r)$p.value
  ljung_lag4 <- Box.test(r, lag = 4, type = "Ljung-Box")$p.value
  ljung_lag8 <- Box.test(r, lag = 8, type = "Ljung-Box")$p.value
  
  # ACF lag-4
  acf_obj <- acf(r, lag.max = 8, plot = FALSE)
  acf_lag4 <- as.numeric(acf_obj$acf[5])
  
  diag_summary <- rbind(diag_summary, data.frame(
    Model = m_name,
    Shapiro_Normality_p = shapiro_p,
    LJungBox_Lag4_p = ljung_lag4,
    ACF_lag4 = acf_lag4,
    LJungBox_Lag8_p = ljung_lag8
  ))
}

write.csv(diag_summary, "Export_Residual_Diagnostics.csv", row.names = FALSE)



# ============================================================================
# PART 6: PLOTTING FIGURE 1 - RECURSIVE PARAMETER TRAJECTORY
# ============================================================================
fig_traj_K <- ggplot(rec_results, aes(x = Sample_Size)) +
  geom_line(aes(y = Verhulst_K / 1e9, color = "Verhulst Log"), linewidth = 1.1) +
  geom_line(aes(y = Gompertz_K / 1e9, color = "Gompertz"), linewidth = 1.1, linetype = "dashed") +
  geom_line(aes(y = Richards_K / 1e9, color = "Richards"), linewidth = 1.1, linetype = "dotdash") +
  scale_color_manual(values = c("Verhulst Log" = "#ff7f00", "Gompertz" = "#4daf4a", "Richards" = "#a65628")) +
  labs(title = "Figure 1A. Trajectory of Estimated Carrying Capacity (K)",
       x = "Sample Size (n)", y = "Carrying Capacity K (IDR Quadrillion)", color = "Model") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "right", plot.title = element_text(face = "bold"))

fig_traj_t <- ggplot(rec_results, aes(x = Sample_Size)) +
  geom_line(aes(y = Verhulst_t_star, color = "Verhulst Log"), linewidth = 1.1) +
  geom_line(aes(y = Gompertz_t_star, color = "Gompertz"), linewidth = 1.1, linetype = "dashed") +
  geom_line(aes(y = Richards_t_star, color = "Richards"), linewidth = 1.1, linetype = "dotdash", na.rm = TRUE) +
  scale_color_manual(values = c("Verhulst Log" = "#ff7f00", "Gompertz" = "#4daf4a", "Richards" = "#a65628")) +
  labs(title = "Figure 1B. Trajectory of Estimated Inflection Quarter (t*)",
       x = "Sample Size (n)", y = "Inflection Point (t*)", color = "Model") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "right", plot.title = element_text(face = "bold"))

fig_stability <- fig_traj_K / fig_traj_t + plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave("Figure_1_Recursive_Parameter_Stability.png", plot = fig_stability, width = 8.5, height = 7, dpi = 350)


# ============================================================================
# PART 7: PLOTTING FIGURE 2 - MULTI-PANEL TRAJECTORIES (REPRESENTATIVE OOS & FULL SAMPLE)
# ============================================================================
generate_panel_plot <- function(train_n, panel_title, is_full_sample = FALSE) {
  panel_data <- data.frame(t = 1:53, Label = df$Label)
  
  if (is_full_sample) {
    m <- fit_all_models(df)
    panel_data$Linear_Fit    <- predict(m$Linear, newdata = panel_data)
    panel_data$Quad_Fit      <- predict(m$Quadratic, newdata = panel_data)
    panel_data$Exp_Fit       <- predict(m$Exp, newdata = panel_data)
    panel_data$Log_Fit       <- predict(m$Log, newdata = panel_data)
    panel_data$Gom_Fit       <- predict(m$Gom, newdata = panel_data)
    panel_data$Ric_Fit       <- predict(m$Ric, newdata = panel_data)
    panel_data$Data_Group    <- "Full Sample"
    panel_data$Asset         <- df$Asset
  } else {
    df_train <- subset(df, t <= train_n)
    m <- fit_all_models(df_train)
    
    panel_data$Linear_Fit    <- predict(m$Linear, newdata = panel_data)
    panel_data$Quad_Fit      <- predict(m$Quadratic, newdata = panel_data)
    panel_data$Exp_Fit       <- predict(m$Exp, newdata = panel_data)
    panel_data$Log_Fit       <- predict(m$Log, newdata = panel_data)
    panel_data$Gom_Fit       <- predict(m$Gom, newdata = panel_data)
    panel_data$Ric_Fit       <- predict(m$Ric, newdata = panel_data)
    
    panel_data <- merge(panel_data, df[, c("t", "Asset")], by = "t", all.x = TRUE)
    panel_data$Data_Group    <- ifelse(panel_data$t <= train_n, "Training", "Testing")
  }
  
  model_colors <- c(
    "Linear"          = "#e41a1c",
    "Quadratic"       = "#377eb8",
    "Exponential"     = "#984ea3",
    "Verhulst Log"    = "#ff7f00",
    "Gompertz"        = "#4daf4a",
    "Richards"        = "#a65628"
  )
  
  p <- ggplot(panel_data, aes(x = t))
  
  if (!is_full_sample) {
    p <- p + geom_vline(xintercept = train_n + 0.5, linetype = "dashed", color = "grey40", linewidth = 0.7)
  }
  
  p <- p + 
    geom_line(aes(y = Linear_Fit / 1e9, color = "Linear"), linewidth = 0.8, linetype = "dotted") +
    geom_line(aes(y = Quad_Fit / 1e9, color = "Quadratic"), linewidth = 1.1) +
    geom_line(aes(y = Exp_Fit / 1e9, color = "Exponential"), linewidth = 0.8, linetype = "twodash") +
    geom_line(aes(y = Log_Fit / 1e9, color = "Verhulst Log"), linewidth = 1.0, linetype = "longdash") +
    geom_line(aes(y = Gom_Fit / 1e9, color = "Gompertz"), linewidth = 1.0) +
    geom_line(aes(y = Ric_Fit / 1e9, color = "Richards"), linewidth = 0.8, linetype = "dotdash") +
    geom_point(aes(y = Asset / 1e9, fill = Data_Group, shape = Data_Group), size = 2, color = "black") +
    scale_x_continuous(breaks = seq(1, 53, by = 8), labels = panel_data$Label[seq(1, 53, by = 8)]) +
    scale_y_continuous(labels = function(x) sprintf("%.2f", x / 1e9)) +
    scale_color_manual(values = model_colors) +
    scale_fill_manual(values = c("Training" = "grey80", "Testing" = "#e41a1c", "Full Sample" = "black")) +
    scale_shape_manual(values = c("Training" = 21, "Testing" = 23, "Full Sample" = 21)) +
    labs(title = panel_title, x = "Timeline (Quarters)", y = "Assets (IDR Quadrillion)", 
         color = "Models", fill = "Data Partition", shape = "Data Partition") +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.border = element_rect(color = "grey80", fill = NA),
      legend.position = "bottom"
    )
  
  return(p)
}

p1 <- generate_panel_plot(38, "(A) OOS Validation: 38/15 Split (n = 38)")
p2 <- generate_panel_plot(43, "(B) OOS Validation: 43/10 Split (n = 43)")
p3 <- generate_panel_plot(48, "(C) OOS Validation: 48/5 Split (n = 48)")
p4 <- generate_panel_plot(53, "(D) Full Sample Calibration (N = 53)", is_full_sample = TRUE)

combined_plot <- (p1 + p2) / (p3 + p4) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

ggsave("Figure_2_Trajectories_OOS_and_FullSample.png", plot = combined_plot, 
       width = 12, height = 9, dpi = 350)


# ============================================================================
# PART 8: PLOTTING FIGURE 3 - NORMAL Q-Q PLOTS (ALL 6 MODELS)
# ============================================================================
df_qq <- data.frame()

for (name in model_names) {
  r <- resids_list[[name]]
  df_qq <- rbind(df_qq, data.frame(
    Residuals = r,
    Model = name
  ))
}

df_qq$Model <- factor(df_qq$Model, levels = model_names)

fig_qq <- ggplot(df_qq, aes(sample = Residuals / 1e6)) +
  stat_qq(color = "#2b8cbe", size = 1.8) +
  stat_qq_line(color = "#e41a1c", linewidth = 0.8, linetype = "dashed") +
  facet_wrap(~ Model, scales = "free_y", ncol = 3) +
  labs(title = "Figure 3. Normal Q-Q Plots of Residuals Across All Models",
       x = "Theoretical Quantiles", y = "Sample Quantiles (IDR Trillion)") +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

ggsave("Figure_3_Normal_QQ_Plots_All_Models.png", plot = fig_qq, width = 10, height = 6.5, dpi = 350)


# ============================================================================
# PART 9: PLOTTING FIGURE 4 - RESIDUAL ACF WITH EXPLICIT QUARTERLY LAGS (ALL MODELS)
# ============================================================================
get_acf_df <- function(resid, name, max_lag = 8) {
  acf_obj <- acf(resid, lag.max = max_lag, plot = FALSE)
  data.frame(
    Lag_Num = acf_obj$lag[-1],
    Lag_Label = paste0("Q", acf_obj$lag[-1]),
    ACF = acf_obj$acf[-1],
    Model = name
  )
}

df_acf_all <- data.frame()
for (name in model_names) {
  df_acf_all <- rbind(df_acf_all, get_acf_df(resids_list[[name]], name))
}

df_acf_all$Model <- factor(df_acf_all$Model, levels = model_names)

ci_bound <- 1.96 / sqrt(nrow(df))

fig_acf <- ggplot(df_acf_all, aes(x = Lag_Num, y = ACF)) +
  geom_vline(xintercept = c(4, 8), linetype = "dotted", color = "grey40", linewidth = 0.6) +
  geom_hline(yintercept = c(-ci_bound, ci_bound), linetype = "solid", color = "#e41a1c", alpha = 0.45, linewidth = 0.8) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_bar(stat = "identity", fill = "darkslategrey", width = 0.3) +
  scale_x_continuous(breaks = 1:8, labels = paste0("Q", 1:8)) +
  facet_wrap(~ Model, ncol = 3) +
  labs(title = "Figure 4. Residual Autocorrelation Functions (ACF) with Quarterly Lags",
       subtitle = "Red solid lines denote 95% confidence limits (±0.269, transparent). Vertical dotted lines highlight annual lag cycles (Q4 & Q8).",
       x = "Autocorrelation Lag Horizon (Quarters)", y = "ACF Coefficient") +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

ggsave("Figure_4_Residual_ACF_All_Models.png", plot = fig_acf, width = 10, height = 6.5, dpi = 350)

print("=== ALL COMPUTATIONS & FIGURE GENERATIONS COMPLETED SUCCESSFULLY ===")

