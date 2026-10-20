# ============================================================================
# MASTER_STRUCTURAL_CALIBRATION_TOTAL_SYSTEM.R
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
required_packages <- c("ggplot2", "minpack.lm", "stats")
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

df_training <- subset(df, t <= 48)   # 2013Q1 - 2024Q4 (n = 48)
df_testing  <- subset(df, t >= 49 & t <= 53) # 2025Q1 - 2026Q1 (m = 5)

# ============================================================================
# PART 2: OUT-OF-SAMPLE TRAINING PHASE (t = 1 to 48)
# ============================================================================
# 1. Exponential
model_exp_tr <- lm(log(Asset) ~ t, data = df_training)
# 2. Verhulst Logistic
model_log_tr <- nlsLM(Asset ~ K / (1 + C * exp(-r * t)), data = df_training,
                      start = list(K = max(df_training$Asset) * 1.5, C = 5, r = 0.08))
# 3. Gompertz
model_gom_tr <- nlsLM(Asset ~ K * exp(-B * exp(-k_g * t)), data = df_training,
                      start = list(K = max(df_training$Asset) * 1.5, B = 2, k_g = 0.04))
# 4. Richards
model_ric_tr <- nlsLM(Asset ~ K / (1 + C * exp(-r * t))^(1/v), data = df_training,
                      start = list(K = max(df_training$Asset) * 1.5, C = 5, r = 0.08, v = 1.0),
                      control = nls.lm.control(maxiter = 1000))
# 5. Bass Diffusion
model_bas_tr <- nlsLM(Asset ~ m_b * ((1 - exp(-(p + q) * t)) / (1 + (q / p) * exp(-(p + q) * t))), 
                      data = df_training,
                      start = list(m_b = max(df_training$Asset) * 1.5, p = 0.005, q = 0.05),
                      control = nls.lm.control(maxiter = 1000))
# Verhulst Baseline Extraction
K_val_tr <- as.numeric(coef(model_log_tr)["K"])
C_val_tr <- as.numeric(coef(model_log_tr)["C"])
r_val_tr <- as.numeric(coef(model_log_tr)["r"])
se_K_tr  <- summary(model_log_tr)$coefficients["K", "Std. Error"]

K_val_tr
C_val_tr
r_val_tr
se_K_tr

summary(model_exp_tr) 
confint(model_exp_tr)
summary(model_log_tr) 
confint(model_log_tr)
summary(model_gom_tr) 
confint(model_gom_tr)
summary(model_ric_tr) 
summary(model_bas_tr) 



# Out-of-Sample Predictions (t = 49 to 53)
act_test <- df_testing$Asset
pred_exp_test <- exp(coef(model_exp_tr)[1] + coef(model_exp_tr)[2] * df_testing$t)
pred_log_test <- coef(model_log_tr)["K"] / (1 + coef(model_log_tr)["C"] * exp(-coef(model_log_tr)["r"] * df_testing$t))
pred_gom_test <- coef(model_gom_tr)["K"] * exp(-coef(model_gom_tr)["B"] * exp(-coef(model_gom_tr)["k_g"] * df_testing$t))
pred_ric_test <- coef(model_ric_tr)["K"] / (1 + coef(model_ric_tr)["C"] * exp(-coef(model_ric_tr)["r"] * df_testing$t))^(1/coef(model_ric_tr)["v"])
pred_bas_test <- coef(model_bas_tr)["m_b"] * ((1 - exp(-(coef(model_bas_tr)["p"] + coef(model_bas_tr)["q"]) * df_testing$t)) / (1 + (coef(model_bas_tr)["q"] / coef(model_bas_tr)["p"]) * exp(-(coef(model_bas_tr)["p"] + coef(model_bas_tr)["q"]) * df_testing$t)))

# Metrics Evaluation
mape_exp <- mean(abs((act_test - pred_exp_test) / act_test)) * 100
mape_log <- mean(abs((act_test - pred_log_test) / act_test)) * 100
mape_gom <- mean(abs((act_test - pred_gom_test) / act_test)) * 100
mape_ric <- mean(abs((act_test - pred_ric_test) / act_test)) * 100
mape_bas <- mean(abs((act_test - pred_bas_test) / act_test)) * 100

rmse_exp <- sqrt(mean((act_test - pred_exp_test)^2))
rmse_log <- sqrt(mean((act_test - pred_log_test)^2))
rmse_gom <- sqrt(mean((act_test - pred_gom_test)^2))
rmse_ric <- sqrt(mean((act_test - pred_ric_test)^2))
rmse_bas <- sqrt(mean((act_test - pred_bas_test)^2))

# Setup Figure 1 Panel Data
df_panelA <- data.frame(t = 1:53)
df_panelA$Exp_Fit_Train <- exp(coef(model_exp_tr)[1] + coef(model_exp_tr)[2] * df_panelA$t)
df_panelA$Log_Fit_Train <- coef(model_log_tr)["K"] / (1 + coef(model_log_tr)["C"] * exp(-coef(model_log_tr)["r"] * df_panelA$t))
df_panelA$Gom_Fit_Train <- coef(model_gom_tr)["K"] * exp(-coef(model_gom_tr)["B"] * exp(-coef(model_gom_tr)["k_g"] * df_panelA$t))
df_panelA$Ric_Fit_Train <- coef(model_ric_tr)["K"] / (1 + coef(model_ric_tr)["C"] * exp(-coef(model_ric_tr)["r"] * df_panelA$t))^(1/coef(model_ric_tr)["v"])
df_panelA$Bas_Fit_Train <- coef(model_bas_tr)["m_b"] * ((1 - exp(-(coef(model_bas_tr)["p"] + coef(model_bas_tr)["q"]) * df_panelA$t)) / (1 + (coef(model_bas_tr)["q"] / coef(model_bas_tr)["p"]) * exp(-(coef(model_bas_tr)["p"] + coef(model_bas_tr)["q"]) * df_panelA$t)))

df_panelA <- merge(df_panelA, df[, c("t", "Year", "Q", "Asset")], by = "t", all.x = TRUE)
df_panelA$Label <- paste0(df_panelA$Year, "-", df_panelA$Q)
df_panelA$Data_Group <- ifelse(df_panelA$t <= 48, "Training Data (2013-2024)", "Testing Data (2025-2026Q1)")

# ============================================================================
# PART 3: FULL DATA CALIBRATION & ESTIMATION (t = 1 to 53)
# ============================================================================
model_exp_full <- lm(log(Asset) ~ t, data = df)
model_log_full <- nlsLM(Asset ~ K / (1 + C * exp(-r * t)), data = df,
                        start = list(K = max(df$Asset) * 1.3, C = 5, r = 0.06))
model_gom_full <- nlsLM(Asset ~ K * exp(-B * exp(-k_g * t)), data = df,
                        start = list(K = max(df$Asset) * 1.3, B = 2, k_g = 0.03))
model_ric_full <- nlsLM(Asset ~ K / (1 + C * exp(-r * t))^(1/v), data = df,
                        start = list(K = max(df$Asset) * 1.3, C = 5, r = 0.06, v = 1.0),
                        control = nls.lm.control(maxiter = 1000))
model_bas_full <- nlsLM(Asset ~ m_b * ((1 - exp(-(p + q) * t)) / (1 + (q / p) * exp(-(p + q) * t))), 
                        data = df,
                        start = list(m_b = max(df$Asset) * 1.3, p = 0.005, q = 0.04),
                        control = nls.lm.control(maxiter = 1000))

# Verhulst Baseline Extraction
K_val <- as.numeric(coef(model_log_full)["K"])
C_val <- as.numeric(coef(model_log_full)["C"])
r_val <- as.numeric(coef(model_log_full)["r"])
se_K  <- summary(model_log_full)$coefficients["K", "Std. Error"]

K_val
C_val
r_val
se_K

summary(model_exp_full) 
confint(model_exp_full)
summary(model_log_full) 
confint(model_log_full)
summary(model_gom_full) 
confint(model_gom_full)
summary(model_ric_full) 
confint(model_ric_full)
summary(model_bas_full)
confint(model_bas_full)

rss_exp <- sum(residuals(model_exp_full)^2)
rss_log <- sum(residuals(model_log_full)^2)
rss_gom <- sum(residuals(model_gom_full)^2)
rss_ric <- sum(residuals(model_ric_full)^2)
rss_bas <- sum(residuals(model_bas_full)^2)

# AIC/BIC Adjudication Matrix
metrics_ic <- data.frame(
  Model = c("Exponential", "Verhulst Logistic", "Gompertz", "Richards", "Bass Diffusion"),
  Parameters = c(2, 3, 3, 4, 3),
  RSS = c(
    rss_exp,
    rss_log,
    rss_gom,
    rss_ric,
    rss_bas
  ),
  AIC = c(AIC(model_exp_full), AIC(model_log_full), AIC(model_gom_full), AIC(model_ric_full), AIC(model_bas_full)),
  BIC = c(BIC(model_exp_full), BIC(model_log_full), BIC(model_gom_full), BIC(model_ric_full), BIC(model_bas_full))
)


# Residual Diagnostics
res_exp <- residuals(model_exp_full)
res_log <- residuals(model_log_full)
res_gom <- residuals(model_gom_full)
res_ric <- residuals(model_ric_full)
res_bas <- residuals(model_bas_full)

sw_exp <- shapiro.test(res_exp)
sw_log <- shapiro.test(res_log)
sw_gom <- shapiro.test(res_gom)
sw_ric <- shapiro.test(res_ric)
sw_bas <- shapiro.test(res_bas)
lb_exp_8 <- Box.test(res_exp, lag = 8, type = "Ljung-Box")
lb_log_8 <- Box.test(res_log, lag = 8, type = "Ljung-Box")
lb_gom_8 <- Box.test(res_gom, lag = 8, type = "Ljung-Box")
lb_ric_8 <- Box.test(res_ric, lag = 8, type = "Ljung-Box")
lb_bas_8 <- Box.test(res_bas, lag = 8, type = "Ljung-Box")

sw_exp$p.value
sw_log$p.value
sw_gom$p.value
sw_ric$p.value
sw_bas$p.value

lb_exp_8$p.value
lb_log_8$p.value
lb_gom_8$p.value
lb_ric_8$p.value
lb_bas_8$p.value

# Setup Figure 2 Horizon Data (2026Q2 - 2030Q4)
t_future <- 54:72
Year_future <- c(rep(2026, 3), rep(2027:2030, each = 4))
Q_future <- c("Q2", "Q3", "Q4", rep(paste0("Q", 1:4), 4))

df_plot <- data.frame(t = 1:72, Year = c(df$Year, Year_future), Q = c(df$Q, Q_future))
df_plot$Exp_Fit <- exp(coef(model_exp_full)[1] + coef(model_exp_full)[2] * df_plot$t)
df_plot$Log_Fit <- K_val / (1 + C_val * exp(-r_val * df_plot$t))
df_plot$Gom_Fit <- coef(model_gom_full)["K"] * exp(-coef(model_gom_full)["B"] * exp(-coef(model_gom_full)["k_g"] * df_plot$t))
df_plot$Ric_Fit <- coef(model_ric_full)["K"] / (1 + coef(model_ric_full)["C"] * exp(-coef(model_ric_full)["r"] * df_plot$t))^(1/coef(model_ric_full)["v"])
df_plot$Bas_Fit <- coef(model_bas_full)["m_b"] * ((1 - exp(-(coef(model_bas_full)["p"] + coef(model_bas_full)["q"]) * df_plot$t)) / (1 + (coef(model_bas_full)["q"] / coef(model_bas_full)["p"]) * exp(-(coef(model_bas_full)["p"] + coef(model_bas_full)["q"]) * df_plot$t)))
df_plot$Log_Low  <- (K_val - 1.96 * se_K) / (1 + C_val * exp(-r_val * df_plot$t))
df_plot$Log_High <- (K_val + 1.96 * se_K) / (1 + C_val * exp(-r_val * df_plot$t))
df_plot <- merge(df_plot, df[, c("t", "Asset")], by = "t", all.x = TRUE)
df_plot$Label <- paste0(df_plot$Year, "-", df_plot$Q)

# ============================================================================
# PART 4: SYSTEM KINETICS & THRESHOLDS HORIZON
# ============================================================================
t_kin <- seq(1, 72, by = 0.1)
A_kin <- K_val / (1 + C_val * exp(-r_val * t_kin))
vel_kin <- r_val * A_kin * (1 - (A_kin / K_val))
acc_kin <- (r_val^2) * A_kin * (1 - (A_kin / K_val)) * (1 - (2 * A_kin / K_val))
df_kin <- data.frame(t = t_kin, Velocity = vel_kin, Acceleration = acc_kin)

t_inflection <- log(C_val) / r_val
yr_inflection <- 2013 + (t_inflection - 1) / 4
t_for_fraction <- function(f) { (log(C_val) - log((1 - f) / f)) / r_val }
thresholds <- c(0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 0.99)
threshold_results <- data.frame(Threshold = numeric(), t_star = numeric(), Year = numeric(), Q = numeric(), Status = character())

for (f in thresholds) {
  t_star <- t_for_fraction(f)
  yr_cal <- 2013 + (t_star - 1) / 4
  q_star <- ceiling(((t_star - 1) %% 4) + 1)
  status <- ifelse(t_star <= 53, "Passed", "Projected Horizon")
  threshold_results <- rbind(threshold_results, data.frame(Threshold = f*100, t_star = t_star, Year = floor(yr_cal), Q = q_star, Status = status))
}

# CSV Exports
write.csv(metrics_ic, "Export_AIC_BIC_Metrics.csv", row.names = FALSE)
write.csv(threshold_results, "Export_Saturation_Thresholds.csv", row.names = FALSE)

# Console Output Report
cat("\n============================================================\n")
cat("          METRIC AND DIAGNOSTIC STATISTICAL REPORT          \n")
cat("============================================================\n")
cat("Out-of-Sample Accuracy Matrices (t=49 to t=53):\n")
cat(sprintf("  - Exponential Model MAPE : %.2f%% | RMSE: %.2e\n", mape_exp, rmse_exp))
cat(sprintf("  - Verhulst Logistic MAPE : %.2f%% | RMSE: %.2e\n", mape_log, rmse_log))
cat(sprintf("  - Gompertz Model MAPE    : %.2f%% | RMSE: %.2e\n", mape_gom, rmse_gom))
cat(sprintf("  - Richards Model MAPE    : %.2f%% | RMSE: %.2e\n", mape_ric, rmse_ric))
cat(sprintf("  - Bass Diffusion MAPE    : %.2f%% | RMSE: %.2e\n\n", mape_bas, rmse_bas))

cat("Information Criteria (Parsimony Adjudication Matrix):\n")
print(metrics_ic)
cat("\n")

cat(sprintf("Verhulst Full-Sample Parameter Estimates:\n"))
cat(sprintf("  - Carrying Capacity (K)  : IDR %.2f Trillion\n", K_val))
cat(sprintf("  - Growth Parameter (r)   : %.4f\n", r_val))
cat(sprintf("  - Integration Const (C)  : %.4f\n", C_val))
cat(sprintf("  - Saturation Ratio (2026Q1): %.2f%%\n\n", (df$Asset[df$t == 53] / K_val) * 100))

cat(sprintf("Residual Diagnostic Probabilities:\n"))
cat(sprintf("  - Exponential Shapiro-Wilk (Normality) p-val : %.5f\n", sw_exp$p.value))
cat(sprintf("  - Verhulst Logistic Shapiro-Wilk Normality p : %.5f\n", sw_log$p.value))
cat(sprintf("  - Exponential Ljung-Box Serial Corr (lag=8) p: %.5f\n", lb_exp_8$p.value))
cat(sprintf("  - Verhulst Logistic Ljung-Box Serial Corr p  : %.5f\n", lb_log_8$p.value))
cat("============================================================\n")
cat("--- TIME-TO-SATURATION ANALYSIS HORIZON ---\n")
print(threshold_results)
cat("============================================================\n\n")

# ============================================================================
# PART 5: GENERATING AND SAVING THE 7 INDEPENDENT FIGURES
# ============================================================================

# --- FIGURE 1: OUT-OF-SAMPLE PREDICTIVE VERIFICATION ---
brks_1 <- seq(1, 53, by = 4)
lbls_1 <- df_panelA$Label[brks_1]

fig1 <- ggplot(df_panelA, aes(x = t)) +
  geom_vline(xintercept = 48.5, linetype = "dashed", color = "purple", linewidth = 0.8) +
  geom_line(aes(y = Exp_Fit_Train, color = "Exponential Trend"), linewidth = 0.9, linetype = "twodash") +
  geom_line(aes(y = Log_Fit_Train, color = "Verhulst Logistic"), linewidth = 1.3) +
  geom_line(aes(y = Gom_Fit_Train, color = "Gompertz Curve"), linewidth = 0.9, linetype = "longdash") +
  geom_line(aes(y = Ric_Fit_Train, color = "Richards Function"), linewidth = 0.9, linetype = "dotted") +
  geom_line(aes(y = Bas_Fit_Train, color = "Bass Diffusion"), linewidth = 0.9, linetype = "dotdash") +
  geom_point(aes(y = Asset, fill = Data_Group, shape = Data_Group), size = 2.5, color = "black") +
  annotate("text", x = 24, y = max(df_panelA$Asset, na.rm=TRUE)*0.85, label = "Training Phase (n = 48)", color = "grey30", fontface = "italic", size = 3.5) +
  annotate("text", x = 51, y = max(df_panelA$Asset, na.rm=TRUE)*1.1, label = "Testing Phase\n(m = 5)", color = "purple", fontface = "bold", size = 3.5) +
  scale_x_continuous(breaks = brks_1, labels = lbls_1) +
  scale_y_continuous(labels = function(x) paste0("IDR ", round(x), " T")) +
  scale_color_manual(values = c("Exponential Trend" = "#e34a33", "Verhulst Logistic" = "#3182bd", "Gompertz Curve" = "#4daf4a", "Richards Function" = "#984ea3", "Bass Diffusion" = "#ff7f00")) +
  scale_fill_manual(values = c("Training Data (2013-2024)" = "grey80", "Testing Data (2025-2026Q1)" = "#e41a1c")) +
  scale_shape_manual(values = c("Training Data (2013-2024)" = 21, "Testing Data (2025-2026Q1)" = 23)) +
  labs(title = "Figure 1. Out-of-sample predictive verification (multi-model vs. hold-out data).", x = "Timeline Parameters (Quarters)", y = "Total Banking Assets", color = "Models", fill = "Data Split", shape = "Data Split") +
  theme_minimal(base_size = 11) + theme(plot.title = element_text(face = "bold"), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.border = element_rect(color = "grey85", fill = NA))
ggsave("Figure_1.png", plot = fig1, width = 8.5, height = 5.5, dpi = 350)

# --- FIGURE 2: MACRO STRUCTURAL PROJECTIONS HORIZON ---
brks_2 <- seq(1, 72, by = 4)
lbls_2 <- df_plot$Label[brks_2]

fig2 <- ggplot(df_plot, aes(x = t)) +
  geom_ribbon(data = subset(df_plot, t > 53), aes(ymin = Log_Low, ymax = Log_High), fill = "#3182bd", alpha = 0.15) +
  geom_hline(yintercept = K_val, linetype = "dashed", color = "red", linewidth = 0.7) +
  geom_line(aes(y = Exp_Fit, color = "Exponential Trend"), linewidth = 0.9, linetype = "twodash") +
  geom_line(aes(y = Log_Fit, color = "Verhulst Logistic"), linewidth = 1.3) +
  geom_line(aes(y = Gom_Fit, color = "Gompertz Curve"), linewidth = 0.9, linetype = "longdash") +
  geom_line(aes(y = Ric_Fit, color = "Richards Function"), linewidth = 0.9, linetype = "dotted") +
  geom_line(aes(y = Bas_Fit, color = "Bass Diffusion"), linewidth = 0.9, linetype = "dotdash") +
  geom_point(data = subset(df_plot, t <= 53), aes(y = Asset), shape = 21, size = 2.2, fill = "grey40", color = "black") +
  geom_vline(xintercept = 53.5, linetype = "dotted", color = "grey40", linewidth = 0.8) +
  annotate("text", x = 26, y = K_val * 1.06, label = paste("Verhulst Saturation Horizon (K) ≈ IDR", round(K_val, 2), "Trillion"), color = "red", fontface = "bold", size = 3.5) +
  annotate("text", x = 63, y = K_val * 0.2, label = "Projection Space\n(2026Q2-2030Q4)", color = "grey30", fontface = "italic", size = 3.5) +
  scale_x_continuous(breaks = brks_2, labels = lbls_2) +
  scale_y_continuous(labels = function(x) paste0("IDR ", round(x), " T")) +
  scale_color_manual(values = c("Exponential Trend" = "#e34a33", "Verhulst Logistic" = "#3182bd", "Gompertz Curve" = "#4daf4a", "Richards Function" = "#984ea3", "Bass Diffusion" = "#ff7f00")) +
  labs(title = "Figure 2. Macro structural projections horizon (full dataset structural calibration).", x = "Timeline Parameters (Quarters)", y = "Total Banking Assets", color = "Models") +
  theme_minimal(base_size = 11) + theme(plot.title = element_text(face = "bold"), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.border = element_rect(color = "grey85", fill = NA))
ggsave("Figure_2.png", plot = fig2, width = 8.5, height = 5.5, dpi = 350)

# --- FIGURE 3: FITTED VS RESIDUALS ---
df_res <- data.frame(Fitted = predict(model_log_full), Residuals = res_log)
fig3 <- ggplot(df_res, aes(x = Fitted, y = Residuals)) +
  geom_point(color = "#3182bd", alpha = 0.75, size = 2.5) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "#e34a33", linewidth = 0.9) +
  labs(title = "Figure 3. Fitted versus residuals diagnostic plot.", x = "Fitted Assets (IDR Trillion)", y = "Residuals") + 
  theme_bw(base_size = 11) + theme(plot.title = element_text(face="bold"))
ggsave("Figure_3.png", plot = fig3, width = 7, height = 5, dpi = 350)

# --- FIGURE 4: Q-Q PLOT ---
fig4 <- ggplot(df_res, aes(sample = Residuals)) +
  stat_qq(color = "#984ea3", alpha = 0.75, size = 2.5) + stat_qq_line(color = "black", linewidth = 0.7) +
  labs(title = "Figure 4. Normal quantile-quantile (Q-Q) plot.", x = "Theoretical Quantiles", y = "Sample Quantiles") + 
  theme_bw(base_size = 11) + theme(plot.title = element_text(face="bold"))
ggsave("Figure_4.png", plot = fig4, width = 7, height = 5, dpi = 350)

# --- FIGURE 5: ACF BOUNDS ---
acf_calc <- acf(res_log, lag.max = 8, plot = FALSE)
df_acf <- data.frame(Lag = acf_calc$lag[-1], ACF = acf_calc$acf[-1])
ci_bound <- 1.96 / sqrt(nrow(df))

fig5 <- ggplot(df_acf, aes(x = Lag, y = ACF)) +
  geom_bar(stat = "identity", fill = "darkslategrey", width = 0.25) +
  geom_hline(yintercept = c(-ci_bound, ci_bound), linetype = "dashed", color = "blue", linewidth = 0.7) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  scale_x_continuous(breaks = 1:8) +
  labs(title = "Figure 5. Autocorrelation function (ACF) error bounds.", x = "Lag Horizon (Quarters)", y = "ACF Coefficient") + 
  theme_bw(base_size = 11) + theme(plot.title = element_text(face="bold"))
ggsave("Figure_5.png", plot = fig5, width = 7, height = 4.5, dpi = 350)

# --- FIGURE 6: VELOCITY PROFILE ---
fig6 <- ggplot(df_kin, aes(x = t, y = Velocity)) +
  geom_line(color = "#ff7f00", linewidth = 1.3) +
  geom_vline(xintercept = t_inflection, linetype = "dotted", color = "red", linewidth = 0.9) +
  labs(title = "Figure 6. Growth velocity profile (dA/dt).", x = "Timeline Parameters (t)", y = "Velocity (dAsset/dt per Quarter)") + 
  theme_minimal(base_size = 11) + theme(panel.border = element_rect(color="grey90", fill=NA), plot.title = element_text(face="bold"))
ggsave("Figure_6.png", plot = fig6, width = 7.5, height = 4.5, dpi = 350)

# --- FIGURE 7: ACCELERATION PROFILE ---
fig7 <- ggplot(df_kin, aes(x = t, y = Acceleration)) +
  geom_line(color = "deeppink4", linewidth = 1.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = t_inflection, linetype = "dotted", color = "red", linewidth = 0.9) +
  labs(title = "Figure 7. Growth acceleration profile (d²A/dt²).", x = "Timeline Parameters (t)", y = "Acceleration (d²Asset/dt²)") + 
  theme_minimal(base_size = 11) + theme(panel.border = element_rect(color="grey90", fill=NA), plot.title = element_text(face="bold"))
ggsave("Figure_7.png", plot = fig7, width = 7.5, height = 4.5, dpi = 350)

print("--- PROCESS COMPLETED SUCCESSFULLY: ALL TABLES EXPORTED AND PLOTS SAVED ---")

