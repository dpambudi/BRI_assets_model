# ============================================================================
# BRI Asset, this code is a prototype, not optimized
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

# ============================================================
# SPLITS
# ============================================================

splits <- list(
  list(train_end = 38, test_start = 39, test_end = 53),
  list(train_end = 43, test_start = 44, test_end = 53),
  list(train_end = 48, test_start = 49, test_end = 53)
)

validation_results <- data.frame()

# store fitted models
saved_models <- list()

# VALIDATION LOOP
metric_fun <- function(actual, predicted){
  
  mae <- mean(abs(actual - predicted))
  
  rmse <- sqrt(
    mean((actual - predicted)^2)
  )
  
  mape <- mean(
    abs((actual - predicted) / actual)
  ) * 100
  
  c(
    MAE = mae,
    RMSE = rmse,
    MAPE = mape
  )
}

for(i in seq_along(splits)) {

  train_end  <- splits[[i]]$train_end
  test_start <- splits[[i]]$test_start
  test_end   <- splits[[i]]$test_end
  
  df_training <- subset(df, t <= train_end)
  df_testing  <- subset(df, t >= test_start & t <= test_end)
  
  # ----------------------------------------------------------
  # FIT MODELS
  # ----------------------------------------------------------
  
  model_exp <- nlsLM(
    Asset ~ A0 * exp(r * t),
    data = df_training,
    start = list(
      A0 = min(df_training$Asset),
      r  = 0.04
    ),
    control = nls.lm.control(
      maxiter = 1024,
      maxfev = 50000
    )
  )
  
  model_log <- nlsLM(
    Asset ~ K / (1 + C * exp(-r * t)),
    data = df_training,
    start = list(
      K = max(df_training$Asset) * 1.5,
      C = 5,
      r = 0.04
    ),
    control = nls.lm.control(
      maxiter = 1024,
      maxfev = 50000
    )
  )
  
  model_gom <- nlsLM(
    Asset ~ K * exp(-B * exp(-k_g * t)),
    data = df_training,
    start = list(
      K   = max(df_training$Asset) * 1.5,
      B   = 5,
      k_g = 0.04
    ),
    control = nls.lm.control(
      maxiter = 1024,
      maxfev = 50000
    )
  )
  
  model_ric <- nlsLM(
    Asset ~ K / (1 + C * exp(-r * t))^(1/v),
    data = df_training,
    start = list(
      K = max(df_training$Asset) * 1.5,
      C = 5,
      r = 0.04,
      v = 1
    ),
    control = nls.lm.control(
      maxiter = 1024,
      maxfev  = 50000
    )
  )
  
  # ----------------------------------------------------------
  # TEST DATA
  # ----------------------------------------------------------
  
  act_test <- df_testing$Asset
  
  pred_exp <- coef(model_exp)["A0"] *
    exp(coef(model_exp)["r"] * df_testing$t)
  
  pred_log <- coef(model_log)["K"] /
    (1 + coef(model_log)["C"] *
       exp(-coef(model_log)["r"] * df_testing$t))
  
  pred_gom <- coef(model_gom)["K"] *
    exp(-coef(model_gom)["B"] *
          exp(-coef(model_gom)["k_g"] * df_testing$t))
  
  pred_ric <- coef(model_ric)["K"] /
    (1 + coef(model_ric)["C"] *
       exp(-coef(model_ric)["r"] * df_testing$t)) ^
    (1 / coef(model_ric)["v"])
  
  # ----------------------------------------------------------
  # METRICS
  # ----------------------------------------------------------
  
  m_exp <- metric_fun(act_test, pred_exp)
  m_log <- metric_fun(act_test, pred_log)
  m_gom <- metric_fun(act_test, pred_gom)
  m_ric <- metric_fun(act_test, pred_ric)
  
  validation_results <- rbind(
    validation_results,
    data.frame(
      Split = paste0(
        train_end,
        "/",
        test_end - train_end
      ),
      Model = c(
        "Exponential",
        "Verhulst",
        "Gompertz",
        "Richards"
      ),
      MAE = c(
        m_exp["MAE"],
        m_log["MAE"],
        m_gom["MAE"],
        m_ric["MAE"]
      ),
      RMSE = c(
        m_exp["RMSE"],
        m_log["RMSE"],
        m_gom["RMSE"],
        m_ric["RMSE"]
      ),
      MAPE = c(
        m_exp["MAPE"],
        m_log["MAPE"],
        m_gom["MAPE"],
        m_ric["MAPE"]
      )
    )
  )
  
  # SAVE MODELS
  saved_models[[i]] <- list(
    exp = model_exp,
    log = model_log,
    gom = model_gom,
    ric = model_ric,
    train_end = train_end
  )
}

write.csv(
  validation_results,
  "Export_MultiWindow_Validation.csv",
  row.names = FALSE
)

summary_validation <-
  aggregate(
    cbind(MAE,RMSE,MAPE) ~ Model,
    data = validation_results,
    FUN = mean
  )

write.csv(
  summary_validation,
  "Export_MultiWindow_Validation_Summary.csv",
  row.names = FALSE
)

print(validation_results)
print(summary_validation)



# ============================================================
# TRAINING-SAMPLE PARAMETER ESTIMATES
# Final validation window (48/5 split)
# ============================================================

model_exp_tr_final <- saved_models[[3]]$exp
model_log_tr_final <- saved_models[[3]]$log
model_gom_tr_final <- saved_models[[3]]$gom
model_ric_tr_final <- saved_models[[3]]$ric

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

get_iter <- function(model){
  model$convInfo$finIter
}

get_status <- function(model){
  if(model$convInfo$isConv){
    "Converged"
  } else {
    "Not Converged"
  }
}

get_stop <- function(model){
  model$convInfo$stopMessage
}

# ------------------------------------------------------------
# Exponential
# ------------------------------------------------------------

coef_exp <- summary(model_exp_tr_final)$coefficients
ci_exp <- confint(model_exp_tr_final)

tbl_exp <- data.frame(
  Model = "Exponential",
  Iterations = get_iter(model_exp_tr_final),
  Status = get_status(model_exp_tr_final),
  Parameter = rownames(coef_exp),
  Estimate = coef_exp[, "Estimate"],
  Std_Error = coef_exp[, "Std. Error"],
  t_value = coef_exp[, "t value"],
  p_value = coef_exp[, "Pr(>|t|)"],
  CI_Lower = ci_exp[,1],
  CI_Upper = ci_exp[,2]
)

# ------------------------------------------------------------
# Verhulst
# ------------------------------------------------------------

coef_log <- summary(model_log_tr_final)$coefficients
ci_log <- confint(model_log_tr_final)

tbl_log <- data.frame(
  Model = "Verhulst Logistic",
  Iterations = get_iter(model_log_tr_final),
  Status = get_status(model_log_tr_final),
  Parameter = rownames(coef_log),
  Estimate = coef_log[, "Estimate"],
  Std_Error = coef_log[, "Std. Error"],
  t_value = coef_log[, "t value"],
  p_value = coef_log[, "Pr(>|t|)"],
  CI_Lower = ci_log[,1],
  CI_Upper = ci_log[,2]
)

# ------------------------------------------------------------
# Gompertz
# ------------------------------------------------------------

coef_gom <- summary(model_gom_tr_final)$coefficients
ci_gom <- confint(model_gom_tr_final)

tbl_gom <- data.frame(
  Model = "Gompertz",
  Iterations = get_iter(model_gom_tr_final),
  Status = get_status(model_gom_tr_final),
  Parameter = rownames(coef_gom),
  Estimate = coef_gom[, "Estimate"],
  Std_Error = coef_gom[, "Std. Error"],
  t_value = coef_gom[, "t value"],
  p_value = coef_gom[, "Pr(>|t|)"],
  CI_Lower = ci_gom[,1],
  CI_Upper = ci_gom[,2]
)

# ------------------------------------------------------------
# Richards
# ------------------------------------------------------------

coef_ric <- summary(model_ric_tr_final)$coefficients
ci_ric <- confint(model_ric_tr_final)

tbl_ric <- data.frame(
  Model = "Richards",
  Iterations = get_iter(model_ric_tr_final),
  Status = get_status(model_ric_tr_final),
  Parameter = rownames(coef_ric),
  Estimate = coef_ric[, "Estimate"],
  Std_Error = coef_ric[, "Std. Error"],
  t_value = coef_ric[, "t value"],
  p_value = coef_ric[, "Pr(>|t|)"],
  CI_Lower = NA,
  CI_Upper = NA
)

# ------------------------------------------------------------
# Merge
# ------------------------------------------------------------

parameter_table_tr48 <- rbind(
  tbl_exp,
  tbl_log,
  tbl_gom,
  tbl_ric
)

print(parameter_table_tr48)

# ------------------------------------------------------------
# Export
# ------------------------------------------------------------

write.csv(
  parameter_table_tr48,
  "Export_Training48_Parameter_Estimates.csv",
  row.names = FALSE
)




# ============================================================
# PANEL DATA CREATOR

create_panel_data <- function(
    model_exp,
    model_log,
    model_gom,
    model_ric,
    train_end,
    df){
  
  panel <- data.frame(
    t = 1:max(df$t)
  )
  
  panel$Exp_Fit_Train <-
    coef(model_exp)["A0"] *
    exp(
      coef(model_exp)["r"] * panel$t
    )
  
  panel$Log_Fit_Train <-
    coef(model_log)["K"] /
    (
      1 +
        coef(model_log)["C"] *
        exp(-coef(model_log)["r"] * panel$t)
    )
  
  panel$Gom_Fit_Train <-
    coef(model_gom)["K"] *
    exp(
      -coef(model_gom)["B"] *
        exp(-coef(model_gom)["k_g"] * panel$t)
    )
  
  panel$Ric_Fit_Train <-
    coef(model_ric)["K"] /
    (
      1 +
        coef(model_ric)["C"] *
        exp(-coef(model_ric)["r"] * panel$t)
    )^(1/coef(model_ric)["v"])
  
  panel <- merge(
    panel,
    df[,c("t","Year","Q","Asset")],
    by = "t",
    all.x = TRUE
  )
  
  panel$Label <-
    paste0(
      panel$Year,
      "-",
      panel$Q
    )
  
  panel$Data_Group <-
    ifelse(
      panel$t <= train_end,
      "Training",
      "Testing"
    )
  
  return(panel)
}

# ============================================================
# PANELS
# ============================================================

df_panelA <- create_panel_data(
  saved_models[[1]]$exp,
  saved_models[[1]]$log,
  saved_models[[1]]$gom,
  saved_models[[1]]$ric,
  train_end = 38,
  df = df
)

df_panelB <- create_panel_data(
  saved_models[[2]]$exp,
  saved_models[[2]]$log,
  saved_models[[2]]$gom,
  saved_models[[2]]$ric,
  train_end = 43,
  df = df
)

df_panelC <- create_panel_data(
  saved_models[[3]]$exp,
  saved_models[[3]]$log,
  saved_models[[3]]$gom,
  saved_models[[3]]$ric,
  train_end = 48,
  df = df
)

# ============================================================
# FUNCTION : CREATE VALIDATION FIGURE
# ============================================================

create_validation_plot <- function(
    df_panel,
    train_end,
    fig_title,
    file_name){
  
  #test_start <- train_end + 1
  
  train_mid <- train_end / 2
  test_mid  <- (train_end + max(df_panel$t)) / 2
  
  fig <- ggplot(df_panel, aes(x = t)) +
    
    geom_vline(
      xintercept = train_end + 0.5,
      linetype = "dashed",
      color = "purple",
      linewidth = 0.8
    ) +
    
    geom_line(
      aes(y = Exp_Fit_Train,
          color = "Exponential"),
      linewidth = 0.9,
      linetype = "twodash"
    ) +
    
    geom_line(
      aes(y = Log_Fit_Train,
          color = "Verhulst"),
      linewidth = 1.3
    ) +
    
    geom_line(
      aes(y = Gom_Fit_Train,
          color = "Gompertz"),
      linewidth = 1.0,
      linetype = "longdash"
    ) +
    
    geom_line(
      aes(y = Ric_Fit_Train,
          color = "Richards"),
      linewidth = 0.9,
      linetype = "dotted"
    ) +
    
    geom_point(
      aes(
        y = Asset,
        fill = Data_Group,
        shape = Data_Group
      ),
      size = 2.4,
      color = "black"
    ) +
    
    annotate(
      "text",
      x = train_mid,
      y = max(df_panel$Asset, na.rm = TRUE) * 0.85,
      label = paste0(
        "Training Phase\n(n = ",
        train_end,
        ")"
      ),
      color = "grey30",
      fontface = "italic",
      size = 3.5
    ) +
    
    annotate(
      "text",
      x = test_mid,
      y = max(df_panel$Asset, na.rm = TRUE) * 1.05,
      label = paste0(
        "Testing Phase\n(m = ",
        max(df_panel$t) - train_end,
        ")"
      ),
      color = "purple",
      fontface = "bold",
      size = 3.5
    ) +
    
    scale_x_continuous(
      breaks = seq(1, 53, by = 4),
      labels = df_panel$Label[seq(1, 53, by = 4)]
    ) +
    
    scale_y_continuous(
      labels = function(x)
        sprintf("%.2f", x / 1e9)
    ) +
    
    scale_color_manual(
      values = c(
        "Exponential" = "#e34a33",
        "Verhulst"    = "#ff7f00",
        "Gompertz"    = "#4daf4a",
        "Richards"    = "#984ea3"
      )
    ) +
    
    scale_fill_manual(
      values = c(
        "Training" = "grey80",
        "Testing"  = "#e41a1c"
      )
    ) +
    
    scale_shape_manual(
      values = c(
        "Training" = 21,
        "Testing"  = 23
      )
    ) +
    
    labs(
      title = fig_title,
      x = "Timeline Parameters (Quarters)",
      y = "Total Assets (Quadrillion IDR)",
      color = "Models",
      fill = "Data Split",
      shape = "Data Split"
    ) +
    
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(
        angle = 45,
        vjust = 1,
        hjust = 1
      ),
      panel.border = element_rect(
        color = "grey85",
        fill = NA
      )
    )
  
  ggsave(
    file_name,
    plot = fig,
    width = 8.5,
    height = 5.5,
    dpi = 350
  )
  
  return(fig)
}

fig1A <- create_validation_plot(
  df_panelA,
  train_end = 38,
  fig_title =
    "Figure 1A. Out-of-sample validation (38-quarter training window).",
  file_name =
    "Figure_1A.png"
)

fig1B <- create_validation_plot(
  df_panelB,
  train_end = 43,
  fig_title =
    "Figure 1B. Out-of-sample validation (43-quarter training window).",
  file_name =
    "Figure_1B.png"
)

fig1C <- create_validation_plot(
  df_panelC,
  train_end = 48,
  fig_title =
    "Figure 1C. Out-of-sample validation (48-quarter training window).",
  file_name =
    "Figure_1C.png"
)

# ============================================================================
# PART 3: FULL DATA CALIBRATION & ESTIMATION (t = 1 to 53)
# ============================================================================

control_nls <- nls.lm.control(
  maxiter = 1024,
  maxfev  = 50000
)

model_exp_full <- nlsLM(
  Asset ~ A0 * exp(r * t),
  data = df,
  start = list(
    A0 = min(df$Asset),
    r = 0.04
  ),
  control = control_nls
)
model_log_full <- nlsLM(Asset ~ K / (1 + C * exp(-r * t)), data = df,
                        start = list(K = max(df$Asset) * 1.5, C = 5, r = 0.04),
                        control = control_nls)
model_gom_full <- nlsLM(Asset ~ K * exp(-B * exp(-k_g * t)), data = df,
                        start = list(K = max(df$Asset) * 1.5, B = 5, k_g = 0.04),
                        control = control_nls)
model_ric_full <- nlsLM(Asset ~ K / (1 + C * exp(-r * t))^(1/v), data = df,
                        start = list(K = max(df$Asset) * 1.5, C = 5, r = 0.04, v = 1.0),
                        control = control_nls)

# Verhulst Baseline Extraction
K_val <- as.numeric(coef(model_log_full)["K"])
C_val <- as.numeric(coef(model_log_full)["C"])
r_val <- as.numeric(coef(model_log_full)["r"])
se_K  <- summary(model_log_full)$coefficients["K", "Std. Error"]
se_C <- summary(model_log_full)$coefficients["C", "Std. Error"]
se_r <- summary(model_log_full)$coefficients["r", "Std. Error"]

K_val
C_val
r_val
se_K
se_C 
se_r 

# Gompertz Parameter Extraction
K_gom <- as.numeric(coef(model_gom_full)["K"])
B_gom <- as.numeric(coef(model_gom_full)["B"])
kg_gom <- as.numeric(coef(model_gom_full)["k_g"])

se_K_gom <- summary(model_gom_full)$coefficients["K", "Std. Error"]
se_B_gom <- summary(model_gom_full)$coefficients["B", "Std. Error"]
se_kg_gom <- summary(model_gom_full)$coefficients["k_g", "Std. Error"]

K_gom
B_gom
kg_gom

#===============
summary(model_exp_full) 
confint(model_exp_full)
summary(model_log_full) 
confint(model_log_full)
summary(model_gom_full) 
confint(model_gom_full)
summary(model_ric_full) 
confint(model_ric_full)


# ============================================================
# PARAMETER ESTIMATION TABLE
# ============================================================



# Exponential
coef_exp <- summary(model_exp_full)$coefficients
ci_exp   <- confint(model_exp_full)

tbl_exp <- data.frame(
  Model = "Exponential",
  Iterations = get_iter(model_exp_full),
  Status = get_status(model_exp_full),
  Parameter = rownames(coef_exp),
  Estimate = coef_exp[, "Estimate"],
  Std_Error = coef_exp[, "Std. Error"],
  t_value = coef_exp[, "t value"],
  p_value = coef_exp[, "Pr(>|t|)"],
  CI_Lower = ci_exp[,1],
  CI_Upper = ci_exp[,2]
)

# ------------------------------------------------------------
# Verhulst
# ------------------------------------------------------------

coef_log <- summary(model_log_full)$coefficients
ci_log   <- confint(model_log_full)

tbl_log <- data.frame(
  Model = "Verhulst Logistic",
  Iterations = get_iter(model_log_full),
  Status = get_status(model_log_full),
  Parameter = rownames(coef_log),
  Estimate = coef_log[, "Estimate"],
  Std_Error = coef_log[, "Std. Error"],
  t_value = coef_log[, "t value"],
  p_value = coef_log[, "Pr(>|t|)"],
  CI_Lower = ci_log[,1],
  CI_Upper = ci_log[,2]
)

# ------------------------------------------------------------
# Gompertz
# ------------------------------------------------------------

coef_gom <- summary(model_gom_full)$coefficients
ci_gom   <- confint(model_gom_full)

tbl_gom <- data.frame(
  Model = "Gompertz",
  Iterations = get_iter(model_gom_full),
  Status = get_status(model_gom_full),
  Parameter = rownames(coef_gom),
  Estimate = coef_gom[, "Estimate"],
  Std_Error = coef_gom[, "Std. Error"],
  t_value = coef_gom[, "t value"],
  p_value = coef_gom[, "Pr(>|t|)"],
  CI_Lower = ci_gom[,1],
  CI_Upper = ci_gom[,2]
)

# ------------------------------------------------------------
# Richards
# ------------------------------------------------------------

coef_ric <- summary(model_ric_full)$coefficients

tbl_ric <- data.frame(
  Model = "Richards",
  Iterations = get_iter(model_ric_full),
  Status = get_status(model_ric_full),
  Parameter = rownames(coef_ric),
  Estimate = coef_ric[, "Estimate"],
  Std_Error = coef_ric[, "Std. Error"],
  t_value = coef_ric[, "t value"],
  p_value = coef_ric[, "Pr(>|t|)"],
  CI_Lower = NA,
  CI_Upper = NA
)

# ------------------------------------------------------------
# MERGE ALL
# ------------------------------------------------------------

parameter_table <- rbind(
  tbl_exp,
  tbl_log,
  tbl_gom,
  tbl_ric
)

print(parameter_table)

# ------------------------------------------------------------
# EXPORT CSV
# ------------------------------------------------------------

write.csv(
  parameter_table,
  "Export_Model_Parameter_Estimates.csv",
  row.names = FALSE
)

convergence_table <- data.frame(
  Model = c(
    "Exponential",
    "Verhulst",
    "Gompertz",
    "Richards"
  ),
  Iterations = c(
    get_iter(model_exp_full),
    get_iter(model_log_full),
    get_iter(model_gom_full),
    get_iter(model_ric_full)
  ),
  Status = c(
    get_status(model_exp_full),
    get_status(model_log_full),
    get_status(model_gom_full),
    get_status(model_ric_full)
  ),
  Stop_Message = c(
    get_stop(model_exp_full),
    get_stop(model_log_full),
    get_stop(model_gom_full),
    get_stop(model_ric_full)
  )
)

print(convergence_table)

rss_exp <- sum(residuals(model_exp_full)^2)
rss_log <- sum(residuals(model_log_full)^2)
rss_gom <- sum(residuals(model_gom_full)^2)
rss_ric <- sum(residuals(model_ric_full)^2)


# AIC/BIC Adjudication Matrix
metrics_ic <- data.frame(
  Model = c("Exponential", "Verhulst Logistic", "Gompertz", "Richards"),
  Parameters = c(2, 3, 3, 4),
  RSS = c(
    rss_exp,
    rss_log,
    rss_gom,
    rss_ric
  ),
  AIC = c(AIC(model_exp_full), AIC(model_log_full), AIC(model_gom_full), AIC(model_ric_full)),
  BIC = c(BIC(model_exp_full), BIC(model_log_full), BIC(model_gom_full), BIC(model_ric_full))
)

# DELTA AIC
best_aic <- min(metrics_ic$AIC)
metrics_ic$Delta_AIC <- metrics_ic$AIC - best_aic

# AKAIKE WEIGHTS
aic_weight_raw <- exp(-0.5 * metrics_ic$Delta_AIC)
metrics_ic$Akaike_Weight <- aic_weight_raw / sum(aic_weight_raw)

# SORT BY AIC (OPTIONAL)
metrics_ic <- metrics_ic[order(metrics_ic$AIC), ]

rownames(metrics_ic) <- NULL

# DISPLAY
print(metrics_ic)

# Residual Diagnostics
res_exp <- residuals(model_exp_full)
res_log <- residuals(model_log_full)
res_gom <- residuals(model_gom_full)
res_ric <- residuals(model_ric_full)


sw_exp <- shapiro.test(res_exp)
sw_log <- shapiro.test(res_log)
sw_gom <- shapiro.test(res_gom)
sw_ric <- shapiro.test(res_ric)

lb_exp_8 <- Box.test(res_exp, lag = 8, type = "Ljung-Box")
lb_log_8 <- Box.test(res_log, lag = 8, type = "Ljung-Box")
lb_gom_8 <- Box.test(res_gom, lag = 8, type = "Ljung-Box")
lb_ric_8 <- Box.test(res_ric, lag = 8, type = "Ljung-Box")


sw_exp$p.value
sw_log$p.value
sw_gom$p.value
sw_ric$p.value

lb_exp_8$p.value
lb_log_8$p.value
lb_gom_8$p.value
lb_ric_8$p.value


# Setup Figure 2 Horizon Data (2026Q2 - 2030Q4)
t_future <- 54:72
Year_future <- c(rep(2026, 3), rep(2027:2030, each = 4))
Q_future <- c("Q2", "Q3", "Q4", rep(paste0("Q", 1:4), 4))

df_plot <- data.frame(t = 1:72, Year = c(df$Year, Year_future), Q = c(df$Q, Q_future))
df_plot$Exp_Fit <-
  coef(model_exp_full)["A0"] *
  exp(
    coef(model_exp_full)["r"] *
      df_plot$t
  )
df_plot$Log_Fit <- K_val / (1 + C_val * exp(-r_val * df_plot$t))
df_plot$Gom_Fit <- coef(model_gom_full)["K"] * exp(-coef(model_gom_full)["B"] * exp(-coef(model_gom_full)["k_g"] * df_plot$t))
df_plot$Ric_Fit <- coef(model_ric_full)["K"] / (1 + coef(model_ric_full)["C"] * exp(-coef(model_ric_full)["r"] * df_plot$t))^(1/coef(model_ric_full)["v"])
#df_plot$Bas_Fit <- coef(model_bas_full)["m_b"] * ((1 - exp(-(coef(model_bas_full)["p"] + coef(model_bas_full)["q"]) * df_plot$t)) / (1 + (coef(model_bas_full)["q"] / coef(model_bas_full)["p"]) * exp(-(coef(model_bas_full)["p"] + coef(model_bas_full)["q"]) * df_plot$t)))
df_plot$Log_Low  <- (K_val - 1.96 * se_K) / (1 + C_val * exp(-r_val * df_plot$t))
df_plot$Log_High <- (K_val + 1.96 * se_K) / (1 + C_val * exp(-r_val * df_plot$t))
df_plot <- merge(df_plot, df[, c("t", "Asset")], by = "t", all.x = TRUE)
df_plot$Label <- paste0(df_plot$Year, "-", df_plot$Q)

# ============================================================================
# PART 4: SYSTEM KINETICS & THRESHOLDS HORIZON
# ============================================================================
t_kin <- seq(1, 72, by = 0.1)
#A_kin <- K_val / (1 + C_val * exp(-r_val * t_kin))
#vel_kin <- r_val * A_kin * (1 - (A_kin / K_val))
#acc_kin <- (r_val^2) * A_kin * (1 - (A_kin / K_val)) * (1 - (2 * A_kin / K_val))

# VERHULST KINETICS
A_log <- K_val / (1 + C_val * exp(-r_val * t_kin))
vel_log <- r_val * A_log * (1 - (A_log / K_val))
acc_log <- (r_val^2) * A_log * (1 - (A_log / K_val)) * (1 - (2 * A_log / K_val))

# GOMPERTZ KINETICS
A_gom <- K_gom * exp(-B_gom * exp(-kg_gom * t_kin))
vel_gom <- kg_gom * A_gom * B_gom * exp(-kg_gom * t_kin)
acc_gom <- (kg_gom^2) * A_gom * B_gom * exp(-kg_gom * t_kin) * (B_gom * exp(-kg_gom * t_kin) - 1)

# INFLECTION POINTS
t_inflection_log <- log(C_val) / r_val
t_inflection_gom <- log(B_gom) / kg_gom
yr_inflection_log <- 2013 + (t_inflection_log - 1) / 4
yr_inflection_gom <- 2013 + (t_inflection_gom - 1) / 4

#df_kin <- data.frame(t = t_kin, Velocity = vel_kin, Acceleration = acc_kin)
#t_inflection <- log(C_val) / r_val
#yr_inflection <- 2013 + (t_inflection - 1) / 4

# DATAFRAME
df_kin <- data.frame(
  t = t_kin,
  Velocity_Verhulst = vel_log,
  Velocity_Gompertz = vel_gom,
  Acc_Verhulst = acc_log,
  Acc_Gompertz = acc_gom
)

# ============================================================================
# SATURATION ANALYSIS
# ============================================================================

# VERHULST SATURATION
t_for_fraction_log <- function(f) {
  (log(C_val) - log((1 - f) / f)) / r_val
}

# GOMPERTZ SATURATION
t_for_fraction_gom <- function(f) {
  (log(B_gom) - log(-log(f))) / kg_gom
}

thresholds <- c(0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 0.99)

threshold_results_log <-
  data.frame(
    Threshold = numeric(),
    t_star = numeric(),
    Year = numeric(),
    Q = numeric(),
    Status = character()
  )

threshold_results_gom <-
  data.frame(
    Threshold = numeric(),
    t_star = numeric(),
    Year = numeric(),
    Q = numeric(),
    Status = character()
  )

# VERHULST
for (f in thresholds) {
  
  t_star <- t_for_fraction_log(f)
  
  yr_cal <- 2013 + (t_star - 1) / 4
  
  q_star <- ceiling(((t_star - 1) %% 4) + 1)
  
  status <- ifelse(
    t_star <= 53,
    "Passed",
    "Projected Horizon"
  )
  
  threshold_results_log <-
    rbind(
      threshold_results_log,
      data.frame(
        Threshold = f * 100,
        t_star = t_star,
        Year = floor(yr_cal),
        Q = q_star,
        Status = status
      )
    )
}

threshold_results_gom <- data.frame(
  Threshold = numeric(),
  t_star = numeric(),
  Year = numeric(),
  Q = integer(),
  Status = character(),
  stringsAsFactors = FALSE
)

# GOMPERTZ
for (f in thresholds) {
  
  t_star <- t_for_fraction_gom(f)
  quarter_index <- ceiling(t_star- 1)
  yr_cal <- 2013 + quarter_index %/% 4
  q_star<- (quarter_index %% 4) + 1
  
  status <- ifelse(
    t_star <= 53,
    "Passed",
    "Projected Horizon"
  )
  
  
  
  threshold_results_gom <-
    rbind(
      threshold_results_gom,
      data.frame(
        Threshold = f * 100,
        t_star = t_star,
        Year = floor(yr_cal),
        Q = q_star,
        Status = status
      )
    )
}

# CSV Exports
write.csv(metrics_ic, "Export_AIC_BIC_Metrics.csv", row.names = FALSE)
write.csv(threshold_results_log, "Export_Verhulst_Saturation_Thresholds.csv", row.names = FALSE)
write.csv(threshold_results_gom, "Export_Gompertz_Saturation_Thresholds.csv", row.names = FALSE)

threshold_results_log$Model <- "Verhulst"
threshold_results_gom$Model <- "Gompertz"

threshold_compare <- rbind(threshold_results_log, threshold_results_gom)
write.csv(threshold_compare, "Export_Saturation_Thresholds_Comparison.csv", row.names = FALSE)
saturation_summary <- data.frame(
  Metric = c(
    "Carrying Capacity",
    "Current Saturation (%)",
    "Inflection Time",
    "Inflection Year"
  ),
  
  Verhulst = c(
    K_val,
    (df$Asset[df$t==53]/K_val)*100,
    t_inflection_log,
    yr_inflection_log
  ),
  
  Gompertz = c(
    K_gom,
    (df$Asset[df$t==53]/K_gom)*100,
    t_inflection_gom,
    yr_inflection_gom
  )
)

write.csv(
  saturation_summary,
  "Export_Verhulst_Gompertz_Comparison.csv",
  row.names = FALSE
)

# Console Output Report


cat(sprintf("Verhulst Full-Sample Parameter Estimates:\n"))
cat(sprintf("  - Carrying Capacity (K)  : IDR %.2f\n", K_val))
cat(sprintf("  - Growth Parameter (r)   : %.4f\n", r_val))
cat(sprintf("  - Integration Const (C)  : %.4f\n", C_val))
cat(sprintf("  - Saturation Ratio (2026Q1): %.2f%%\n\n", (df$Asset[df$t == 53] / K_val) * 100))

cat("\n")
cat("============================================================\n")
cat("GOMPERTZ FULL-SAMPLE PARAMETER ESTIMATES\n")
cat("============================================================\n")

cat(sprintf(" - Carrying Capacity (K) : %.2f\n", K_gom))
cat(sprintf(" - Growth Parameter (k_g): %.5f\n", kg_gom))
cat(sprintf(" - Shape Parameter (B) : %.5f\n", B_gom))
cat(sprintf(" - Saturation Ratio (2026Q1): %.2f%%\n\n", (df$Asset[df$t == 53] / K_gom) * 100))

cat(sprintf("Residual Diagnostic Probabilities:\n"))
cat(sprintf("  - Exponential Shapiro-Wilk (Normality) p-val : %.5f\n", sw_exp$p.value))
cat(sprintf("  - Verhulst Logistic Shapiro-Wilk Normality p : %.5f\n", sw_log$p.value))
cat(sprintf("  - Gompertz Shapiro-Wilk (Normality) p-value : %.5f\n", sw_gom$p.value))
cat("\n")
cat(sprintf("  - Exponential Ljung-Box Serial Corr (lag=8) p: %.5f\n", lb_exp_8$p.value))
cat(sprintf("  - Verhulst Logistic Ljung-Box Serial Corr p  : %.5f\n", lb_log_8$p.value))
cat(sprintf("  - Gompertz Ljung-Box Serial Corr p-value            : %.5f\n", lb_gom_8$p.value))
cat("============================================================\n")
cat("--- TIME-TO-SATURATION ANALYSIS HORIZON ---\n")
print(threshold_results_log)
print(threshold_results_gom)
cat("============================================================\n\n")

# ============================================================================
# PART 5: GENERATING AND SAVING THE 7 INDEPENDENT FIGURES
# ============================================================================


# --- FIGURE 2: MACRO STRUCTURAL PROJECTIONS HORIZON ---
brks_2 <- seq(1, 72, by = 4)
lbls_2 <- df_plot$Label[brks_2]

fig2 <- ggplot(df_plot, aes(x = t)) +
  geom_ribbon(data = subset(df_plot, t > 53), aes(ymin = Log_Low, ymax = Log_High), fill = "#3182bd", alpha = 0.15) +
  geom_hline(yintercept = K_val, linetype = "dashed", color = "red", linewidth = 0.7) +
  geom_hline(yintercept = K_gom, color = "#4daf4a", linetype = "dashed",linewidth = 0.8) +
  geom_line(aes(y = Exp_Fit, color = "Exponential"), linewidth = 0.9, linetype = "twodash") +
  geom_line(aes(y = Log_Fit, color = "Verhulst Log"), linewidth = 1.3) +
  geom_line(aes(y = Gom_Fit, color = "Gompertz"), linewidth = 0.9, linetype = "longdash") +
  geom_line(aes(y = Ric_Fit, color = "Richards"), linewidth = 0.9, linetype = "dotted") +
  #geom_line(aes(y = Bas_Fit, color = "Bass Diffusion"), linewidth = 0.9, linetype = "dotdash") +
  geom_point(data = subset(df_plot, t <= 53), aes(y = Asset), shape = 21, size = 2.2, fill = "grey40", color = "black") +
  geom_vline(xintercept = 53.5, linetype = "dotted", color = "grey40", linewidth = 0.8) +
  annotate("text", x = 18, y = K_val * 1.05, label = paste0("Verhulst K ≈ ", round(K_val/1e9, 3), "QIDR"), color = "red", fontface = "bold", size = 3.5) +
  annotate("text", x = 18, y = K_gom * 1.03, label = paste0("Gompertz K ≈ ", round(K_gom / 1e9, 3)," QIDR"), color = "#4daf4a", size = 3.2, hjust = 0) +
  annotate("text", x = 63, y = K_val * 0.2, label = "Projection Space\n(2026Q2-2030Q4)", color = "grey30", fontface = "italic", size = 3.5) +
  scale_x_continuous(breaks = brks_2, labels = lbls_2) +
  scale_y_continuous(labels = function(x) sprintf("%.2f", x / 1e9)) +
  scale_color_manual(values = c("Exponential" = "#e34a33", "Verhulst Log" = "#ff7f00", "Gompertz" = "#4daf4a", "Richards" = "#984ea3", "Bass Diffusion" = "#ff7f00")) +
  labs(title = "Figure 2. Macro structural projections horizon (full dataset structural calibration).", x = "Timeline Parameters (Quarters)", y = "Total Assets (Quadrillion IDR)", color = "Models") +
  theme_minimal(base_size = 11) + theme(plot.title = element_text(face = "bold"), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.border = element_rect(color = "grey85", fill = NA))
ggsave("Figure_2.png", plot = fig2, width = 8.5, height = 5.5, dpi = 350)

# --- FIGURE 3: FITTED VS RESIDUALS ----------------------------

res_gom <- residuals(model_gom_full)

df_res_compare <- rbind(
  data.frame(
    Fitted = predict(model_log_full),
    Residuals = res_log,
    Model = "Verhulst"
  ),
  data.frame(
    Fitted = predict(model_gom_full),
    Residuals = res_gom,
    Model = "Gompertz"
  )
)

fig3 <- ggplot(
  df_res_compare,
  aes(
    x = Fitted / 1e9,
    y = Residuals / 1e9
  )
) +
  geom_point(
    color = "#3182bd",
    alpha = 0.75,
    size = 2.5
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "#e34a33",
    linewidth = 0.9
  ) +
  facet_wrap(~Model, nrow = 1) +
  labs(
    title = "Figure 3. Fitted-versus-residuals diagnostic plots.",
    x = "Fitted Assets (Quadrillion IDR)",
    y = "Residuals (Quadrillion IDR)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold")
  )

ggsave(
  "Figure_3.png",
  plot = fig3,
  width = 8,
  height = 4.8,
  dpi = 350
)

# --- FIGURE 4: Q-Q PLOTS -------------------------------------

df_qq <- rbind(
  data.frame(
    Residuals = res_log,
    Model = "Verhulst"
  ),
  data.frame(
    Residuals = res_gom,
    Model = "Gompertz"
  )
)

fig4 <- ggplot(
  df_qq,
  aes(sample = Residuals)
) +
  stat_qq(
    color = "#984ea3",
    alpha = 0.75,
    size = 2.5
  ) +
  stat_qq_line(
    color = "black",
    linewidth = 0.7
  ) +
  facet_wrap(~Model, nrow = 1) +
  labs(
    title = "Figure 4. Normal Q-Q plots of model residuals.",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold")
  )

ggsave(
  "Figure_4.png",
  plot = fig4,
  width = 8,
  height = 4.8,
  dpi = 350
)

# --- FIGURE 5: ACF COMPARISON -------------------------------

acf_log <- acf(
  res_log,
  lag.max = 8,
  plot = FALSE
)

acf_gom <- acf(
  res_gom,
  lag.max = 8,
  plot = FALSE
)

df_acf_compare <- rbind(
  data.frame(
    Lag = as.numeric(acf_log$lag[-1]),
    ACF = as.numeric(acf_log$acf[-1]),
    Model = "Verhulst"
  ),
  data.frame(
    Lag = as.numeric(acf_gom$lag[-1]),
    ACF = as.numeric(acf_gom$acf[-1]),
    Model = "Gompertz"
  )
)

acf_log$acf[which(acf_log$lag == 4)]
acf_gom$acf[which(acf_gom$lag == 4)]
ci_bound <- 1.96 / sqrt(nrow(df))
ci_bound
fig5 <- ggplot(
  df_acf_compare,
  aes(
    x = Lag,
    y = ACF
  )
) +
  geom_bar(
    stat = "identity",
    fill = "darkslategrey",
    width = 0.25
  ) +
  geom_hline(
    yintercept = c(-ci_bound, ci_bound),
    linetype = "dashed",
    color = "blue",
    linewidth = 0.7
  ) +
  geom_hline(
    yintercept = 0,
    color = "black",
    linewidth = 0.5
  ) +
  geom_vline(
    xintercept = c(4,8),
    color = "red",
    linetype = "dotted",
    linewidth = 0.8
  ) +
  facet_wrap(~Model, nrow = 1) +
  scale_x_continuous(
    breaks = 1:8,
    labels = paste0("Q", 1:8)
  ) +
  labs(
    title = "Figure 5. Residual autocorrelation functions (ACF).",
    x = "Lag Horizon (Quarters)",
    y = "ACF Coefficient"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold")
  )

ggsave(
  "Figure_5.png",
  plot = fig5,
  width = 8,
  height = 4.5,
  dpi = 350
)

# --- FIGURE 6: VELOCITY PROFILE ---
#fig6 <- ggplot(df_kin, aes(x = t, y = Velocity)) +
#  geom_line(color = "#ff7f00", linewidth = 1.3) +
#  geom_vline(xintercept = t_inflection, linetype = "dotted", color = "red", linewidth = 0.9) +
#  labs(title = "Figure 6. Growth velocity profile (dA/dt).", x = "Timeline Parameters (t)", y = "Velocity (dAsset/dt per Quarter)") + 
#  theme_minimal(base_size = 11) + theme(panel.border = element_rect(color="grey90", fill=NA), plot.title = element_text(face="bold"))
#ggsave("Figure_6.png", plot = fig6, width = 7.5, height = 4.5, dpi = 350)

# --- FIGURE 6: VELOCITY PROFILE Verhulst+Gompertz---
library(patchwork)

# ------------------------------------------------------------
# Empirical velocity
# ------------------------------------------------------------

emp_vel <- data.frame(
  t = df$t[-1],
  Velocity = diff(df$Asset)
)

fig6B <- ggplot(
  emp_vel,
  aes(x = t, y = Velocity / 1e6)
) +
  
  # model-implied velocity (background)
  geom_line(
    data = df_kin,
    aes(
      x = t,
      y = Velocity_Verhulst / 1e6
    ),
    color = "#ff7f00",
    linewidth = 1.1,
    alpha = 0.50
  ) +
  
  geom_line(
    data = df_kin,
    aes(
      x = t,
      y = Velocity_Gompertz / 1e6
    ),
    color = "#4daf4a",
    linewidth = 1.1,
    alpha = 0.50
  ) +
  
  # empirical growth line
  geom_line(
    color = "grey60",
    linewidth = 0.6
  ) +
  
  geom_point(
    color = "black",
    size = 1.8
  ) +
  
  # Verhulst peak velocity
  geom_vline(
    xintercept = t_inflection_log,
    color = "#ff7f00",
    linetype = "dotted",
    linewidth = 0.8
  ) +
  
  # Gompertz peak velocity
  geom_vline(
    xintercept = t_inflection_gom,
    color = "#4daf4a",
    linetype = "dotted",
    linewidth = 0.8
  ) +
  
  #annotate(
  #  "text",
  #  x = t_inflection_log,
  #  y = max(emp_vel$Velocity/1e6, na.rm = TRUE) * 0.95,
  #  label = "Verhulst Peak",
  #  color = "#ff7f00",
  #  angle = 90,
  #  vjust = -0.3,
  #  size = 3
  #) +
  
  #annotate(
  #  "text",
  #  x = t_inflection_gom,
  #  y = max(emp_vel$Velocity/1e6, na.rm = TRUE) * 0.95,
  #  label = "Gompertz Peak",
  #  color = "#4daf4a",
  #  angle = 90,
  #  vjust = -0.3,
  #  size = 3
  #) +
  
  labs(
    title = "(B) Empirical Asset Growth",
    x = "Timeline Parameters (t)",
    y = "ΔAssets (Trillion IDR per Quarter)"
  ) +
  
  theme_minimal(base_size = 11)

# ------------------------------------------------------------
# Model-implied velocity
# ------------------------------------------------------------

fig6A <- ggplot() +
  
  geom_line(
    data = df_kin,
    aes(
      x = t,
      y = Velocity_Verhulst / 1e6,
      color = "Verhulst"
    ),
    linewidth = 1.3
  ) +
  
  geom_line(
    data = df_kin,
    aes(
      x = t,
      y = Velocity_Gompertz / 1e6,
      color = "Gompertz"
    ),
    linewidth = 1.3
  ) +
  
  geom_vline(
    xintercept = t_inflection_log,
    color = "#ff7f00",
    linetype = "dotted"
  ) +
  
  geom_vline(
    xintercept = t_inflection_gom,
    color = "#4daf4a",
    linetype = "dotted"
  ) +
  
  scale_color_manual(
    values = c(
      "Verhulst" = "#ff7f00",
      "Gompertz" = "#4daf4a"
    )
  ) +
  
  labs(
    title = "(A) Model-Implied Growth Velocity",
    x = "Timeline Parameters (t)",
    y = "dA/dt (Trillion IDR per Quarter)",
    color = "Model"
  ) +
  
  theme_minimal(base_size = 11)

# ------------------------------------------------------------
# Combine
# ------------------------------------------------------------

fig6 <- fig6A + fig6B+
  plot_annotation(
    title = "Figure 6. Empirical and model-implied growth velocity profiles."
  )

ggsave(
  "Figure_6.png",
  plot = fig6,
  width = 10.5,
  height = 4.5,
  dpi = 350
)

# --- FIGURE 7: ACCELERATION PROFILE ---
#fig7 <- ggplot(df_kin, aes(x = t, y = Acceleration)) +
#  geom_line(color = "deeppink4", linewidth = 1.3) +
#  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
#  geom_vline(xintercept = t_inflection, linetype = "dotted", color = "red", linewidth = 0.9) +
#  labs(title = "Figure 7. Growth acceleration profile (d²A/dt²).", x = "Timeline Parameters (t)", y = "Acceleration (d²Asset/dt²)") + 
#  theme_minimal(base_size = 11) + theme(panel.border = element_rect(color="grey90", fill=NA), plot.title = element_text(face="bold"))
#ggsave("Figure_7.png", plot = fig7, width = 7.5, height = 4.5, dpi = 350)

# --- FIGURE 7: ACCELERATION PROFILE Verhulst+Gompertz---
fig7 <- ggplot() +
  
  geom_line(
    data = df_kin,
    aes(
      x = t,
      y = Acc_Verhulst,
      color = "Verhulst"
    ),
    linewidth = 1.3
  ) +
  
  geom_line(
    data = df_kin,
    aes(
      x = t,
      y = Acc_Gompertz,
      color = "Gompertz"
    ),
    linewidth = 1.3,
    #linetype = "dashed"
  ) +
  
  geom_hline(
    yintercept = 0,
    color = "grey40",
    linetype = "dashed"
  ) +
  
  geom_vline(
    xintercept = t_inflection_log,
    color = "#ff7f00",
    linetype = "dotted"
  ) +
  
  geom_vline(
    xintercept = t_inflection_gom,
    color = "#4daf4a",
    linetype = "dotted"
  ) +
  
  scale_color_manual(
    values = c(
      "Verhulst" = "#ff7f00",
      "Gompertz" = "#4daf4a"
    )
  ) +
  
  labs(
    title = "Figure 7. Growth acceleration comparison (d²A/dt²).",
    x = "Timeline Parameters (t)",
    y = "Acceleration",
    color = "Model"
  ) +
  
  theme_minimal(base_size = 11)

ggsave(
  "Figure_7.png",
  fig7,
  width = 7.5,
  height = 4.5,
  dpi = 350
)

print("--- PROCESS COMPLETED SUCCESSFULLY: ALL TABLES EXPORTED AND PLOTS SAVED ---")

