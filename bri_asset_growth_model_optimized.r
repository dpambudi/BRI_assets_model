# ============================================================================
# BRI Asset - Optimized Version
# ============================================================================
# Core Goal: Full Statistical Computation, Adjudication Matrices, Kinetics,
#            CSV Table Exports, and 7 Independent High-Res Figure Exports.
# ============================================================================

# --- PART 0: ENVIRONMENT CLEANING & WORKING DIRECTORY INITIALIZATION ---
rm(list = ls())
graphics.off()
if (basename(Sys.getenv("RSTUDIO_USER_IDENTITY")) != "") {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

# Load Required Libraries Safely
required_packages <- c("ggplot2", "minpack.lm", "stats", "patchwork")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# ============================================================================
# PART 1: DATA INGESTION & SUBSETTING (2013Q1 - 2026Q1)
# ============================================================================
if (!file.exists("DataBRI.csv")) {
  stop("Critical Error: 'DataBRI.csv' not found in the working directory.")
}

df <- read.csv("DataBRI.csv", sep = ";")
df <- na.omit(df[, c("t", "Year", "Q", "Asset")])

# ============================================================================
# HELPER FUNCTIONS FOR MODELING & EVALUATION
# ============================================================================
control_nls <- nls.lm.control(maxiter = 1024, maxfev = 50000)

fit_models <- function(data) {
  max_asset <- max(data$Asset)
  min_asset <- min(data$Asset)
  
  models <- list(
    Exp = nlsLM(Asset ~ A0 * exp(r * t), data = data,
                start = list(A0 = min_asset, r = 0.04), control = control_nls),
    Log = nlsLM(Asset ~ K / (1 + C * exp(-r * t)), data = data,
                start = list(K = max_asset * 1.5, C = 5, r = 0.04), control = control_nls),
    Gom = nlsLM(Asset ~ K * exp(-B * exp(-k_g * t)), data = data,
                start = list(K = max_asset * 1.5, B = 5, k_g = 0.04), control = control_nls),
    Ric = nlsLM(Asset ~ K / (1 + C * exp(-r * t))^(1/v), data = data,
                start = list(K = max_asset * 1.5, C = 5, r = 0.04, v = 1), control = control_nls)
  )
  return(models)
}

metric_fun <- function(actual, predicted) {
  c(MAE = mean(abs(actual - predicted)),
    RMSE = sqrt(mean((actual - predicted)^2)),
    MAPE = mean(abs((actual - predicted) / actual)) * 100)
}

extract_params <- function(model, name) {
  coef_sum <- summary(model)$coefficients
  ci <- tryCatch(confint(model), error = function(e) matrix(NA, nrow = nrow(coef_sum), ncol = 2))
  
  if(name == "Richards") ci <- matrix(NA, nrow = nrow(coef_sum), ncol = 2) # Manual override for Richards per original
  
  data.frame(
    Model = name,
    Iterations = model$convInfo$finIter,
    Status = ifelse(model$convInfo$isConv, "Converged", "Not Converged"),
    Parameter = rownames(coef_sum),
    Estimate = coef_sum[, "Estimate"],
    Std_Error = coef_sum[, "Std. Error"],
    t_value = coef_sum[, "t value"],
    p_value = coef_sum[, "Pr(>|t|)"],
    CI_Lower = ci[, 1],
    CI_Upper = ci[, 2]
  )
}

# ============================================================================
# PART 2: VALIDATION LOOP
# ============================================================================
splits <- list(list(tr = 38, ts = 39, te = 53), list(tr = 43, ts = 44, te = 53), list(tr = 48, ts = 49, te = 53))

val_results <- data.frame()
saved_models <- list()
panels <- list()

model_names <- c("Exponential", "Verhulst", "Gompertz", "Richards")

for(i in seq_along(splits)) {
  s <- splits[[i]]
  df_train <- subset(df, t <= s$tr)
  df_test  <- subset(df, t >= s$ts & t <= s$te)
  
  models <- fit_models(df_train)
  saved_models[[i]] <- c(models, list(train_end = s$tr))
  
  # Predictions
  preds <- list(
    predict(models$Exp, newdata = df_test),
    predict(models$Log, newdata = df_test),
    predict(models$Gom, newdata = df_test),
    predict(models$Ric, newdata = df_test)
  )
  
  # Metrics
  for (j in 1:4) {
    m <- metric_fun(df_test$Asset, preds[[j]])
    val_results <- rbind(val_results, data.frame(
      Split = paste0(s$tr, "/", s$te - s$tr),
      Model = model_names[j], MAE = m["MAE"], RMSE = m["RMSE"], MAPE = m["MAPE"]
    ))
  }
  
  # Panel Data Creation
  panel <- data.frame(t = 1:max(df$t))
  panel$Exp_Fit_Train <- predict(models$Exp, newdata = panel)
  panel$Log_Fit_Train <- predict(models$Log, newdata = panel)
  panel$Gom_Fit_Train <- predict(models$Gom, newdata = panel)
  panel$Ric_Fit_Train <- predict(models$Ric, newdata = panel)
  
  panel <- merge(panel, df[, c("t", "Year", "Q", "Asset")], by = "t", all.x = TRUE)
  panel$Label <- paste0(panel$Year, "-", panel$Q)
  panel$Data_Group <- ifelse(panel$t <= s$tr, "Training", "Testing")
  panels[[i]] <- panel
}

write.csv(val_results, "Export_MultiWindow_Validation.csv", row.names = FALSE)
write.csv(aggregate(. ~ Model, data = val_results[, -1], mean), "Export_MultiWindow_Validation_Summary.csv", row.names = FALSE)

# 48-window Parameter Export
tr48_params <- do.call(rbind, mapply(extract_params, saved_models[[3]][1:4], model_names, SIMPLIFY = FALSE))
write.csv(tr48_params, "Export_Training48_Parameter_Estimates.csv", row.names = FALSE)

# ============================================================================
# FUNCTION: VALIDATION FIGURES
# ============================================================================
create_val_plot <- function(panel, train_end, fig_title, file_name) {
  fig <- ggplot(panel, aes(x = t)) +
    geom_vline(xintercept = train_end + 0.5, linetype = "dashed", color = "purple", linewidth = 0.8) +
    geom_line(aes(y = Exp_Fit_Train, color = "Exponential"), linewidth = 0.9, linetype = "twodash") +
    geom_line(aes(y = Log_Fit_Train, color = "Verhulst"), linewidth = 1.3) +
    geom_line(aes(y = Gom_Fit_Train, color = "Gompertz"), linewidth = 1.0, linetype = "longdash") +
    geom_line(aes(y = Ric_Fit_Train, color = "Richards"), linewidth = 0.9, linetype = "dotted") +
    geom_point(aes(y = Asset, fill = Data_Group, shape = Data_Group), size = 2.4, color = "black") +
    annotate("text", x = train_end / 2, y = max(panel$Asset, na.rm = TRUE) * 0.85, 
             label = paste0("Training Phase\n(n = ", train_end, ")"), color = "grey30", fontface = "italic", size = 3.5) +
    annotate("text", x = (train_end + max(panel$t)) / 2, y = max(panel$Asset, na.rm = TRUE) * 1.05, 
             label = paste0("Testing Phase\n(m = ", max(panel$t) - train_end, ")"), color = "purple", fontface = "bold", size = 3.5) +
    scale_x_continuous(breaks = seq(1, 53, by = 4), labels = panel$Label[seq(1, 53, by = 4)]) +
    scale_y_continuous(labels = function(x) sprintf("%.2f", x / 1e9)) +
    scale_color_manual(values = c("Exponential" = "#e34a33", "Verhulst" = "#ff7f00", "Gompertz" = "#4daf4a", "Richards" = "#984ea3")) +
    scale_fill_manual(values = c("Training" = "grey80", "Testing" = "#e41a1c")) +
    scale_shape_manual(values = c("Training" = 21, "Testing" = 23)) +
    labs(title = fig_title, x = "Timeline Parameters (Quarters)", y = "Total Assets (Quadrillion IDR)", color = "Models", fill = "Data Split", shape = "Data Split") +
    theme_minimal(base_size = 11) + theme(plot.title = element_text(face = "bold"), axis.text.x = element_text(angle = 45, hjust = 1), panel.border = element_rect(color = "grey85", fill = NA))
  
  ggsave(file_name, plot = fig, width = 8.5, height = 5.5, dpi = 350)
}

letters_fig <- c("A", "B", "C")
for (i in 1:3) {
  create_val_plot(panels[[i]], splits[[i]]$tr, 
                  paste0("Figure 1", letters_fig[i], ". Out-of-sample validation (", splits[[i]]$tr, "-quarter training window)."), 
                  paste0("Figure_1", letters_fig[i], ".png"))
}

# ============================================================================
# PART 3: FULL DATA CALIBRATION & ESTIMATION (t = 1 to 53)
# ============================================================================
full_models <- fit_models(df)

# Parameter Extraction & Export
param_table_full <- do.call(rbind, mapply(extract_params, full_models, model_names, SIMPLIFY = FALSE))
write.csv(param_table_full, "Export_Model_Parameter_Estimates.csv", row.names = FALSE)

# AIC/BIC & RSS
rss_vals <- sapply(full_models, function(m) sum(residuals(m)^2))
aic_vals <- sapply(full_models, AIC)
bic_vals <- sapply(full_models, BIC)

metrics_ic <- data.frame(
  Model = model_names, Parameters = c(2, 3, 3, 4),
  RSS = rss_vals, AIC = aic_vals, BIC = bic_vals
)
metrics_ic$Delta_AIC <- metrics_ic$AIC - min(metrics_ic$AIC)
metrics_ic$Akaike_Weight <- exp(-0.5 * metrics_ic$Delta_AIC) / sum(exp(-0.5 * metrics_ic$Delta_AIC))
metrics_ic <- metrics_ic[order(metrics_ic$AIC), ]
write.csv(metrics_ic, "Export_AIC_BIC_Metrics.csv", row.names = FALSE)

# Diagnostics
resids <- lapply(full_models, residuals)
shapiro_tests <- sapply(resids, function(r) shapiro.test(r)$p.value)
ljung_tests <- sapply(resids, function(r) Box.test(r, lag = 8, type = "Ljung-Box")$p.value)

# ============================================================================
# PART 4: SYSTEM KINETICS & THRESHOLDS HORIZON
# ============================================================================
K_val <- coef(full_models$Log)["K"]; C_val <- coef(full_models$Log)["C"]; r_val <- coef(full_models$Log)["r"]
K_gom <- coef(full_models$Gom)["K"]; B_gom <- coef(full_models$Gom)["B"]; kg_gom <- coef(full_models$Gom)["k_g"]
se_K <- summary(full_models$Log)$coefficients["K", "Std. Error"]

t_kin <- seq(1, 72, by = 0.1)
A_log <- K_val / (1 + C_val * exp(-r_val * t_kin))
A_gom <- K_gom * exp(-B_gom * exp(-kg_gom * t_kin))

df_kin <- data.frame(
  t = t_kin,
  Velocity_Verhulst = r_val * A_log * (1 - (A_log / K_val)),
  Velocity_Gompertz = kg_gom * A_gom * B_gom * exp(-kg_gom * t_kin),
  Acc_Verhulst = (r_val^2) * A_log * (1 - (A_log / K_val)) * (1 - (2 * A_log / K_val)),
  Acc_Gompertz = (kg_gom^2) * A_gom * B_gom * exp(-kg_gom * t_kin) * (B_gom * exp(-kg_gom * t_kin) - 1)
)

t_inf_log <- log(C_val) / r_val; yr_inf_log <- 2013 + (t_inf_log - 1) / 4
t_inf_gom <- log(B_gom) / kg_gom; yr_inf_gom <- 2013 + (t_inf_gom - 1) / 4

# Vectorized Saturation Thresholds
thresholds <- c(0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 0.99)
t_star_log <- (log(C_val) - log((1 - thresholds) / thresholds)) / r_val
t_star_gom <- (log(B_gom) - log(-log(thresholds))) / kg_gom

generate_threshold_df <- function(t_star, thresh_vals, model_name) {
  quarter_index <- ceiling(t_star) - 1
  data.frame(
    Threshold = thresh_vals * 100, t_star = t_star,
    Year = 2013 + (quarter_index %/% 4), Q = (quarter_index %% 4) + 1,
    Status = ifelse(t_star <= 53, "Passed", "Projected Horizon"), Model = model_name
  )
}

thresh_log_df <- generate_threshold_df(t_star_log, thresholds, "Verhulst")
thresh_gom_df <- generate_threshold_df(t_star_gom, thresholds, "Gompertz")

write.csv(thresh_log_df[, -6], "Export_Verhulst_Saturation_Thresholds.csv", row.names = FALSE)
write.csv(thresh_gom_df[, -6], "Export_Gompertz_Saturation_Thresholds.csv", row.names = FALSE)
write.csv(rbind(thresh_log_df, thresh_gom_df), "Export_Saturation_Thresholds_Comparison.csv", row.names = FALSE)

sat_summary <- data.frame(
  Metric = c("Carrying Capacity", "Current Saturation (%)", "Inflection Time", "Inflection Year"),
  Verhulst = c(K_val, (df$Asset[df$t==53]/K_val)*100, t_inf_log, yr_inf_log),
  Gompertz = c(K_gom, (df$Asset[df$t==53]/K_gom)*100, t_inf_gom, yr_inf_gom)
)
write.csv(sat_summary, "Export_Verhulst_Gompertz_Comparison.csv", row.names = FALSE)

# ============================================================================
# PART 5: FIGURES 2 - 7
# ============================================================================

# Prepare Plot Data (Fig 2)
df_plot <- data.frame(t = 1:72, Year = c(df$Year, rep(2026, 3), rep(2027:2030, each = 4)), 
                      Q = c(df$Q, "Q2", "Q3", "Q4", rep(paste0("Q", 1:4), 4)))
df_plot <- merge(df_plot, df[, c("t", "Asset")], by = "t", all.x = TRUE)
df_plot$Label <- paste0(df_plot$Year, "-", df_plot$Q)

df_plot$Exp_Fit <- predict(full_models$Exp, newdata = df_plot)
df_plot$Log_Fit <- predict(full_models$Log, newdata = df_plot)
df_plot$Gom_Fit <- predict(full_models$Gom, newdata = df_plot)
df_plot$Ric_Fit <- predict(full_models$Ric, newdata = df_plot)
df_plot$Log_Low <- (K_val - 1.96 * se_K) / (1 + C_val * exp(-r_val * df_plot$t))
df_plot$Log_High <- (K_val + 1.96 * se_K) / (1 + C_val * exp(-r_val * df_plot$t))

# FIGURE 2
fig2 <- ggplot(df_plot, aes(x = t)) +
  geom_ribbon(data = subset(df_plot, t > 53), aes(ymin = Log_Low, ymax = Log_High), fill = "#3182bd", alpha = 0.15) +
  geom_hline(yintercept = K_val, linetype = "dashed", color = "red", linewidth = 0.7) +
  geom_hline(yintercept = K_gom, color = "#4daf4a", linetype = "dashed", linewidth = 0.8) +
  geom_line(aes(y = Exp_Fit, color = "Exponential"), linewidth = 0.9, linetype = "twodash") +
  geom_line(aes(y = Log_Fit, color = "Verhulst Log"), linewidth = 1.3) +
  geom_line(aes(y = Gom_Fit, color = "Gompertz"), linewidth = 0.9, linetype = "longdash") +
  geom_line(aes(y = Ric_Fit, color = "Richards"), linewidth = 0.9, linetype = "dotted") +
  geom_point(data = subset(df_plot, t <= 53), aes(y = Asset), shape = 21, size = 2.2, fill = "grey40", color = "black") +
  geom_vline(xintercept = 53.5, linetype = "dotted", color = "grey40", linewidth = 0.8) +
  annotate("text", x=18, y=K_val*1.05, label=paste0("Verhulst K ≈ ", round(K_val/1e9, 3), "QIDR"), color="red", fontface="bold", size=3.5) +
  annotate("text", x=18, y=K_gom*1.03, label=paste0("Gompertz K ≈ ", round(K_gom/1e9, 3)," QIDR"), color="#4daf4a", size=3.2, hjust=0) +
  scale_x_continuous(breaks = seq(1, 72, by = 4), labels = df_plot$Label[seq(1, 72, by = 4)]) +
  scale_y_continuous(labels = function(x) sprintf("%.2f", x / 1e9)) +
  scale_color_manual(values = c("Exponential" = "#e34a33", "Verhulst Log" = "#ff7f00", "Gompertz" = "#4daf4a", "Richards" = "#984ea3")) +
  labs(title = "Figure 2. Macro structural projections horizon.", x = "Timeline Parameters (Quarters)", y = "Total Assets (QIDR)", color = "Models") +
  theme_minimal(base_size = 11) + theme(plot.title = element_text(face="bold"), axis.text.x = element_text(angle=45, hjust=1), panel.border = element_rect(color="grey85", fill=NA))
ggsave("Figure_2.png", plot = fig2, width = 8.5, height = 5.5, dpi = 350)

# FIGURE 3 & 4: Residuals & QQ
df_res <- rbind(
  data.frame(Fitted = predict(full_models$Log), Residuals = resids$Log, Model = "Verhulst"),
  data.frame(Fitted = predict(full_models$Gom), Residuals = resids$Gom, Model = "Gompertz")
)
fig3 <- ggplot(df_res, aes(x = Fitted / 1e9, y = Residuals / 1e9)) + geom_point(color = "#3182bd", alpha = 0.75, size = 2.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#e34a33", linewidth = 0.9) + facet_wrap(~Model) +
  labs(title = "Figure 3. Fitted-versus-residuals diagnostic plots.", x = "Fitted Assets (QIDR)", y = "Residuals (QIDR)") + theme_bw(base_size = 11)
ggsave("Figure_3.png", plot = fig3, width = 8, height = 4.8, dpi = 350)

fig4 <- ggplot(df_res, aes(sample = Residuals)) + stat_qq(color = "#984ea3", size = 2.5) + stat_qq_line(color = "black") + facet_wrap(~Model) +
  labs(title = "Figure 4. Normal Q-Q plots of model residuals.", x = "Theoretical Quantiles", y = "Sample Quantiles") + theme_bw()
ggsave("Figure_4.png", plot = fig4, width = 8, height = 4.8, dpi = 350)

# FIGURE 5: ACF
acf_data <- function(resid, name) { a <- acf(resid, lag.max=8, plot=FALSE); data.frame(Lag=a$lag[-1], ACF=a$acf[-1], Model=name) }
df_acf <- rbind(acf_data(resids$Log, "Verhulst"), acf_data(resids$Gom, "Gompertz"))
ci_bound <- 1.96 / sqrt(nrow(df))

fig5 <- ggplot(df_acf, aes(x = Lag, y = ACF)) + geom_bar(stat = "identity", fill = "darkslategrey", width = 0.25) +
  geom_hline(yintercept = c(-ci_bound, ci_bound), linetype = "dashed", color = "blue", linewidth = 0.7) +
  geom_hline(yintercept = 0, color = "black") + facet_wrap(~Model) + scale_x_continuous(breaks=1:8, labels=paste0("Q",1:8)) +
  labs(title = "Figure 5. Residual autocorrelation functions (ACF).", x = "Lag Horizon (Quarters)", y = "ACF Coefficient") + theme_bw()
ggsave("Figure_5.png", plot = fig5, width = 8, height = 4.5, dpi = 350)

# FIGURE 6: Velocity (Using patchwork)
emp_vel <- data.frame(t = df$t[-1], Velocity = diff(df$Asset))
fig6A <- ggplot(df_kin) + geom_line(aes(t, Velocity_Verhulst/1e6, color="Verhulst"), linewidth=1.3) +
  geom_line(aes(t, Velocity_Gompertz/1e6, color="Gompertz"), linewidth=1.3) +
  geom_vline(xintercept=t_inf_log, color="#ff7f00", linetype="dotted") + geom_vline(xintercept=t_inf_gom, color="#4daf4a", linetype="dotted") +
  scale_color_manual(values=c("Verhulst"="#ff7f00", "Gompertz"="#4daf4a")) + labs(title="(A) Model-Implied Growth", x="t", y="dA/dt (TIDR)") + theme_minimal()

fig6B <- ggplot(emp_vel, aes(t, Velocity/1e6)) + 
  geom_line(data=df_kin, aes(t, Velocity_Verhulst/1e6), color="#ff7f00", alpha=0.5, linewidth=1.1) +
  geom_line(data=df_kin, aes(t, Velocity_Gompertz/1e6), color="#4daf4a", alpha=0.5, linewidth=1.1) +
  geom_line(color="grey60") + geom_point() + geom_vline(xintercept=t_inf_log, color="#ff7f00", linetype="dotted") + geom_vline(xintercept=t_inf_gom, color="#4daf4a", linetype="dotted") +
  labs(title="(B) Empirical Asset Growth", x="t", y="ΔAssets (TIDR)") + theme_minimal()

ggsave("Figure_6.png", plot = fig6A + fig6B + plot_annotation(title="Figure 6. Empirical and model-implied growth velocity profiles."), width = 10.5, height = 4.5, dpi = 350)

# FIGURE 7: Acceleration
fig7 <- ggplot(df_kin) + geom_line(aes(t, Acc_Verhulst, color="Verhulst"), linewidth=1.3) + geom_line(aes(t, Acc_Gompertz, color="Gompertz"), linewidth=1.3) +
  geom_hline(yintercept=0, color="grey40", linetype="dashed") + geom_vline(xintercept=t_inf_log, color="#ff7f00", linetype="dotted") +
  geom_vline(xintercept=t_inf_gom, color="#4daf4a", linetype="dotted") + scale_color_manual(values=c("Verhulst"="#ff7f00", "Gompertz"="#4daf4a")) +
  labs(title="Figure 7. Growth acceleration comparison (d²A/dt²).", x="t", y="Acceleration", color="Model") + theme_minimal()
ggsave("Figure_7.png", fig7, width = 7.5, height = 4.5, dpi = 350)

print("--- PROCESS COMPLETED SUCCESSFULLY: ALL TABLES EXPORTED AND PLOTS SAVED ---")

