# ============================================================
# BCCA-BASED BIAS CORRECTION AND VALIDATION
# Raw GCM vs GEQM, GQM, LS, PT
# Reference: XGBoost-corrected CHIRPS
# Period: 1994–2014
# ============================================================

# ============================================================
# 1. INSTALL AND LOAD PACKAGES
# ============================================================

packages <- c("terra", "dplyr", "ggplot2", "qmap", "lubridate", "tidyr")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# ============================================================
# 2. INPUT AND OUTPUT FOLDERS
# ============================================================

# XGBoost-corrected CHIRPS folder
chirps_folder <- "D:/Desktop/HWRM_Thesis/From Gh/TM proposal/data/CHIRPS_Monthly_Abay/XGBOOST_corrected/resampled_0_25deg"

# Raw GCM folder
gcm_folder <- "D:/Desktop/HWRM_Thesis/From Gh/NEX-GDDP-CMIP6"

# Output folder
output_folder <- "D:/Desktop/HWRM_Thesis/From Gh/BCCA_Bias_Correction_Result"

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# Selected GCM name
# Change this if you want another model
selected_gcm <- "MIROC6"

start_year <- 1994
end_year <- 2014

# ============================================================
# 3. FIND TIF FILES
# ============================================================

chirps_files <- list.files(
  chirps_folder,
  pattern = "\\.(tif|tiff)$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

gcm_files_all <- list.files(
  gcm_folder,
  pattern = "\\.(tif|tiff)$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

gcm_files <- gcm_files_all[
  grepl(selected_gcm, gcm_files_all, ignore.case = TRUE)
]

cat("CHIRPS files found:", length(chirps_files), "\n")
cat("Selected GCM files found:", length(gcm_files), "\n")

if (length(chirps_files) == 0) {
  stop("ERROR: No CHIRPS tif files found.")
}

if (length(gcm_files) == 0) {
  stop("ERROR: No GCM tif files found for selected model.")
}

# ============================================================
# 4. EXTRACT YEAR-MONTH FROM FILE NAME
# ============================================================

extract_year_month <- function(file_path) {
  
  fname <- basename(file_path)
  
  match1 <- regmatches(
    fname,
    regexpr("(199[4-9]|200[0-9]|201[0-4])[_\\.-](0[1-9]|1[0-2])", fname)
  )
  
  if (length(match1) > 0 && match1 != "") {
    return(gsub("[_\\.-]", "-", match1))
  }
  
  match2 <- regmatches(
    fname,
    regexpr("(199[4-9]|200[0-9]|201[0-4])(0[1-9]|1[0-2])", fname)
  )
  
  if (length(match2) > 0 && match2 != "") {
    yr <- substr(match2, 1, 4)
    mn <- substr(match2, 5, 6)
    return(paste0(yr, "-", mn))
  }
  
  return(NA)
}

# ============================================================
# 5. CREATE FILE TABLES
# ============================================================

chirps_table <- data.frame(
  date_id = sapply(chirps_files, extract_year_month),
  chirps_file = chirps_files
)

gcm_table <- data.frame(
  date_id = sapply(gcm_files, extract_year_month),
  gcm_file = gcm_files
)

chirps_table <- chirps_table[!is.na(chirps_table$date_id), ]
gcm_table <- gcm_table[!is.na(gcm_table$date_id), ]

chirps_table <- chirps_table[
  chirps_table$date_id >= paste0(start_year, "-01") &
    chirps_table$date_id <= paste0(end_year, "-12"),
]

gcm_table <- gcm_table[
  gcm_table$date_id >= paste0(start_year, "-01") &
    gcm_table$date_id <= paste0(end_year, "-12"),
]

matched_table <- inner_join(
  chirps_table,
  gcm_table,
  by = "date_id"
)

matched_table <- matched_table[order(matched_table$date_id), ]

cat("Matched CHIRPS-GCM months:", nrow(matched_table), "\n")

if (nrow(matched_table) < 30) {
  stop("ERROR: Too few matched months. Check file names and selected_gcm.")
}

write.csv(
  matched_table,
  file.path(output_folder, "matched_CHIRPS_GCM_files.csv"),
  row.names = FALSE
)

# ============================================================
# 6. EXTRACT BASIN-AVERAGE MONTHLY RAINFALL
# ============================================================

results <- data.frame(
  date = as.Date(character()),
  date_id = character(),
  month = integer(),
  CHIRPS_XGBoost = numeric(),
  Raw_GCM = numeric()
)

for (i in 1:nrow(matched_table)) {
  
  date_id <- matched_table$date_id[i]
  cat("Processing:", date_id, "\n")
  
  chirps_r <- rast(matched_table$chirps_file[i])
  gcm_r <- rast(matched_table$gcm_file[i])
  
  # Resample GCM to CHIRPS grid if needed
  if (!compareGeom(chirps_r, gcm_r, stopOnError = FALSE)) {
    gcm_r <- resample(gcm_r, chirps_r, method = "bilinear")
  }
  
  chirps_mean <- global(chirps_r, "mean", na.rm = TRUE)[1, 1]
  gcm_mean <- global(gcm_r, "mean", na.rm = TRUE)[1, 1]
  
  results <- rbind(
    results,
    data.frame(
      date = as.Date(paste0(date_id, "-01")),
      date_id = date_id,
      month = as.numeric(substr(date_id, 6, 7)),
      CHIRPS_XGBoost = chirps_mean,
      Raw_GCM = gcm_mean
    )
  )
}

# Remove missing values
results <- results[complete.cases(results), ]

# Rainfall cannot be negative
results$CHIRPS_XGBoost[results$CHIRPS_XGBoost < 0] <- 0
results$Raw_GCM[results$Raw_GCM < 0] <- 0

# ============================================================
# 7. BIAS CORRECTION METHODS
# ============================================================

obs <- results$CHIRPS_XGBoost
raw <- results$Raw_GCM

# ------------------------------------------------------------
# 7.1 Linear Scaling, LS
# Monthly multiplicative correction
# ------------------------------------------------------------

results$LS <- NA

for (m in 1:12) {
  
  idx <- which(results$month == m)
  
  obs_mean <- mean(results$CHIRPS_XGBoost[idx], na.rm = TRUE)
  raw_mean <- mean(results$Raw_GCM[idx], na.rm = TRUE)
  
  factor <- ifelse(raw_mean == 0, 1, obs_mean / raw_mean)
  
  results$LS[idx] <- results$Raw_GCM[idx] * factor
}

# ------------------------------------------------------------
# 7.2 Power Transformation, PT
# y_corrected = a * raw^b
# Parameters fitted by minimizing RMSE against CHIRPS
# ------------------------------------------------------------

raw_safe <- pmax(raw, 0.001)

pt_objective <- function(par) {
  a <- par[1]
  b <- par[2]
  pred <- a * (raw_safe ^ b)
  sqrt(mean((obs - pred)^2, na.rm = TRUE))
}

pt_fit <- optim(
  par = c(1, 1),
  fn = pt_objective,
  method = "L-BFGS-B",
  lower = c(0.0001, 0.1),
  upper = c(100, 5)
)

a_pt <- pt_fit$par[1]
b_pt <- pt_fit$par[2]

results$PT <- a_pt * (raw_safe ^ b_pt)

# ------------------------------------------------------------
# 7.3 Gamma Quantile Mapping, GQM
# ------------------------------------------------------------

obs_gamma <- pmax(obs, 0.001)
raw_gamma <- pmax(raw, 0.001)

gqm_fit <- fitQmapDIST(
  obs = obs_gamma,
  mod = raw_gamma,
  distr = "gamma",
  qstep = 0.001
)

results$GQM <- doQmapDIST(
  fobj = gqm_fit,
  x = raw_gamma
)

results$GQM[results$GQM < 0] <- 0

# ------------------------------------------------------------
# 7.4 Generalized / Empirical Quantile Mapping, GEQM
# ------------------------------------------------------------

geqm_fit <- fitQmapQUANT(
  obs = obs,
  mod = raw,
  qstep = 0.01,
  nboot = 1,
  wet.day = TRUE,
  type = "linear"
)

results$GEQM <- doQmapQUANT(
  x = raw,
  fobj = geqm_fit
)

results$GEQM[results$GEQM < 0] <- 0

# ============================================================
# 8. PERFORMANCE METRICS
# ============================================================

NSE <- function(obs, sim) {
  1 - sum((obs - sim)^2, na.rm = TRUE) /
    sum((obs - mean(obs, na.rm = TRUE))^2, na.rm = TRUE)
}

RMSE <- function(obs, sim) {
  sqrt(mean((obs - sim)^2, na.rm = TRUE))
}

MAE <- function(obs, sim) {
  mean(abs(obs - sim), na.rm = TRUE)
}

PBIAS <- function(obs, sim) {
  100 * sum(sim - obs, na.rm = TRUE) / sum(obs, na.rm = TRUE)
}

R_value <- function(obs, sim) {
  cor(obs, sim, use = "complete.obs")
}

metrics <- data.frame(
  Method = c("Raw_GCM", "GEQM", "GQM", "LS", "PT"),
  NSE = c(
    NSE(obs, results$Raw_GCM),
    NSE(obs, results$GEQM),
    NSE(obs, results$GQM),
    NSE(obs, results$LS),
    NSE(obs, results$PT)
  ),
  RMSE = c(
    RMSE(obs, results$Raw_GCM),
    RMSE(obs, results$GEQM),
    RMSE(obs, results$GQM),
    RMSE(obs, results$LS),
    RMSE(obs, results$PT)
  ),
  MAE = c(
    MAE(obs, results$Raw_GCM),
    MAE(obs, results$GEQM),
    MAE(obs, results$GQM),
    MAE(obs, results$LS),
    MAE(obs, results$PT)
  ),
  PBIAS = c(
    PBIAS(obs, results$Raw_GCM),
    PBIAS(obs, results$GEQM),
    PBIAS(obs, results$GQM),
    PBIAS(obs, results$LS),
    PBIAS(obs, results$PT)
  ),
  r = c(
    R_value(obs, results$Raw_GCM),
    R_value(obs, results$GEQM),
    R_value(obs, results$GQM),
    R_value(obs, results$LS),
    R_value(obs, results$PT)
  )
)

metrics <- metrics %>%
  arrange(RMSE)

cat("\nPerformance metrics:\n")
print(metrics)

# ============================================================
# 9. SAVE TABLES
# ============================================================

write.csv(
  results,
  file.path(output_folder, "bias_correction_timeseries.csv"),
  row.names = FALSE
)

write.csv(
  metrics,
  file.path(output_folder, "bias_correction_performance_metrics.csv"),
  row.names = FALSE
)

# ============================================================
# 10. PREPARE DATA FOR PDF PLOT
# ============================================================

pdf_data <- results %>%
  select(
    CHIRPS_XGBoost,
    Raw_GCM,
    GEQM,
    GQM,
    LS,
    PT
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Dataset",
    values_to = "Rainfall"
  )

pdf_data$Rainfall[pdf_data$Rainfall < 0] <- 0

pdf_data$Dataset <- factor(
  pdf_data$Dataset,
  levels = c("CHIRPS_XGBoost", "Raw_GCM", "GEQM", "GQM", "LS", "PT")
)

# ============================================================
# 11. PDF / DENSITY COMPARISON PLOT
# ============================================================

p_pdf <- ggplot(
  pdf_data,
  aes(x = Rainfall, color = Dataset, fill = Dataset)
) +
  geom_density(alpha = 0.15, linewidth = 1.1) +
  labs(
    title = "PDF Comparison of Raw and Bias-Corrected GCM Rainfall",
    subtitle = paste0(
      selected_gcm,
      " vs XGBoost-corrected CHIRPS, ",
      start_year,
      "–",
      end_year
    ),
    x = "Monthly Rainfall (mm)",
    y = "Density",
    color = "Dataset",
    fill = "Dataset"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

print(p_pdf)

ggsave(
  filename = file.path(output_folder, "PDF_density_comparison.png"),
  plot = p_pdf,
  width = 10,
  height = 6,
  dpi = 300
)

# ============================================================
# 12. TIME SERIES COMPARISON PLOT
# ============================================================

ts_data <- results %>%
  select(
    date,
    CHIRPS_XGBoost,
    Raw_GCM,
    GEQM,
    GQM,
    LS,
    PT
  ) %>%
  pivot_longer(
    cols = -date,
    names_to = "Dataset",
    values_to = "Rainfall"
  )

ts_data$Dataset <- factor(
  ts_data$Dataset,
  levels = c("CHIRPS_XGBoost", "Raw_GCM", "GEQM", "GQM", "LS", "PT")
)

p_ts <- ggplot(
  ts_data,
  aes(x = date, y = Rainfall, color = Dataset)
) +
  geom_line(linewidth = 0.7) +
  labs(
    title = "Time Series Comparison of Bias Correction Methods",
    subtitle = paste0(
      selected_gcm,
      " rainfall corrected against XGBoost-corrected CHIRPS"
    ),
    x = "Date",
    y = "Monthly Rainfall (mm)",
    color = "Dataset"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

print(p_ts)

ggsave(
  filename = file.path(output_folder, "Time_series_comparison.png"),
  plot = p_ts,
  width = 12,
  height = 6,
  dpi = 300
)

# ============================================================
# 13. SELECT BEST METHOD
# ============================================================

best_method <- metrics %>%
  arrange(RMSE, desc(NSE), desc(r)) %>%
  slice(1)

summary_text <- paste0(
  "The PDF and performance statistics were used to compare the raw ",
  selected_gcm,
  " rainfall with four bias correction methods, namely GEQM, GQM, LS, and PT, ",
  "against the XGBoost-corrected CHIRPS reference data. The best-performing method was ",
  best_method$Method,
  ", with NSE = ",
  round(best_method$NSE, 3),
  ", RMSE = ",
  round(best_method$RMSE, 2),
  " mm, MAE = ",
  round(best_method$MAE, 2),
  " mm, PBIAS = ",
  round(best_method$PBIAS, 2),
  "%, and r = ",
  round(best_method$r, 3),
  "."
)

writeLines(
  summary_text,
  file.path(output_folder, "bias_correction_summary_text.txt")
)

cat("\nSummary:\n")
cat(summary_text, "\n")

cat("\nDONE SUCCESSFULLY!\n")
cat("Bias correction results saved in:\n")
cat(output_folder, "\n")