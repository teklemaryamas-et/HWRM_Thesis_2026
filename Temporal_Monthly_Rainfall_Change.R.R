# ============================================================
# TEMPORAL CHANGE IN MONTHLY RAINFALL
# Study area: Upper Blue Nile Basin, Ethiopia
# Historical baseline: XGBoost-corrected CHIRPS, 1994–2014
# Future rainfall: Final BCCA-PT downscaled MIROC6
# Future periods: 2030–2050 and 2051–2080
# Scenarios: SSP2-4.5 and SSP5-8.5
# ============================================================

# ============================================================
# 1. INSTALL AND LOAD PACKAGES
# ============================================================

packages <- c("terra", "dplyr", "ggplot2", "tidyr", "lubridate")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# ============================================================
# 2. INPUT FOLDERS
# ============================================================

# Final historical reference rainfall: XGBoost-corrected CHIRPS
chirps_folder <- "D:/Desktop/HWRM_Thesis/From Gh/XGBOOST_corrected"

# Final future rainfall: BCCA-PT downscaled MIROC6
future_downscaled_folder <- "D:/Desktop/HWRM_Thesis/From Gh/BCCA_PT_Downscaled_MIROC6_FINAL/03_Future_Downscaled"

# Output folder
output_folder <- "D:/Desktop/HWRM_Thesis/From Gh/Temporal_Monthly_Rainfall_Change_FINAL"

dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_folder, "Tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_folder, "Figures"), recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 3. BASIC SETTINGS
# ============================================================

selected_gcm <- "MIROC6"

ssp_scenarios <- c("ssp245", "ssp585")

historical_start <- 1994
historical_end <- 2014

future_periods <- data.frame(
  period = c("Near_future_2030_2050", "Far_future_2051_2080"),
  start_year = c(2030, 2051),
  end_year = c(2050, 2080)
)

month_names <- c(
  "Jan", "Feb", "Mar", "Apr", "May", "Jun",
  "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
)

# ============================================================
# 4. FUNCTION TO EXTRACT YEAR-MONTH FROM FILE NAME
# ============================================================

extract_year_month <- function(file_path) {
  
  fname <- basename(file_path)
  
  match1 <- regmatches(
    fname,
    regexpr("(19[0-9]{2}|20[0-9]{2})[_\\.-](0[1-9]|1[0-2])", fname)
  )
  
  if (length(match1) > 0 && match1 != "") {
    return(gsub("[_\\.-]", "-", match1))
  }
  
  match2 <- regmatches(
    fname,
    regexpr("(19[0-9]{2}|20[0-9]{2})(0[1-9]|1[0-2])", fname)
  )
  
  if (length(match2) > 0 && match2 != "") {
    yr <- substr(match2, 1, 4)
    mn <- substr(match2, 5, 6)
    return(paste0(yr, "-", mn))
  }
  
  return(NA)
}

# ============================================================
# 5. FUNCTION TO READ AREA-MEAN RAINFALL
# ============================================================

extract_area_mean <- function(files, dataset_name) {
  
  out <- data.frame()
  
  for (i in seq_along(files)) {
    
    date_id <- extract_year_month(files[i])
    
    if (is.na(date_id)) next
    
    r <- terra::rast(files[i])
    
    mean_value <- terra::global(r, fun = "mean", na.rm = TRUE)[1, 1]
    
    out <- rbind(
      out,
      data.frame(
        date_id = date_id,
        date = as.Date(paste0(date_id, "-01")),
        year = as.numeric(substr(date_id, 1, 4)),
        month = as.numeric(substr(date_id, 6, 7)),
        rainfall = mean_value,
        dataset = dataset_name,
        file = files[i]
      )
    )
    
    if (i %% 50 == 0) {
      cat(dataset_name, "processed:", i, "of", length(files), "\n")
    }
  }
  
  out <- out[order(out$date), ]
  
  return(out)
}

# ============================================================
# 6. READ HISTORICAL XGBOOST-CORRECTED CHIRPS
# ============================================================

chirps_files <- list.files(
  path = chirps_folder,
  pattern = "\\.(tif|tiff)$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

if (length(chirps_files) == 0) {
  stop("ERROR: No CHIRPS tif files found. Check chirps_folder path.")
}

chirps_ts <- extract_area_mean(
  files = chirps_files,
  dataset_name = "XGBoost_corrected_CHIRPS"
)

chirps_ts <- chirps_ts[
  chirps_ts$year >= historical_start &
    chirps_ts$year <= historical_end,
]

cat("Historical CHIRPS months:", nrow(chirps_ts), "\n")

if (nrow(chirps_ts) != 252) {
  warning("Expected 252 CHIRPS months for 1994–2014. Check missing files.")
}

write.csv(
  chirps_ts,
  file.path(output_folder, "Tables", "historical_CHIRPS_XGBoost_timeseries_1994_2014.csv"),
  row.names = FALSE
)

# ============================================================
# 7. HISTORICAL MONTHLY BASELINE
# ============================================================

baseline_monthly <- chirps_ts %>%
  group_by(month) %>%
  summarise(
    historical_mean_rainfall = mean(rainfall, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Month = factor(month_names[month], levels = month_names)
  )

write.csv(
  baseline_monthly,
  file.path(output_folder, "Tables", "historical_monthly_baseline_1994_2014.csv"),
  row.names = FALSE
)

cat("\nHistorical monthly baseline:\n")
print(baseline_monthly)

# ============================================================
# 8. READ FINAL BCCA-PT DOWNSCALED FUTURE MIROC6 FILES
# ============================================================

future_files_all <- list.files(
  path = future_downscaled_folder,
  pattern = "\\.(tif|tiff)$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

# Remove backup or accidental non-final files
future_files_all <- future_files_all[
  !grepl("BACKUP", future_files_all, ignore.case = TRUE) &
    !grepl("BAD", future_files_all, ignore.case = TRUE) &
    !grepl("Tables", future_files_all, ignore.case = TRUE)
]

if (length(future_files_all) == 0) {
  stop("ERROR: No final BCCA-PT future tif files found. Check future_downscaled_folder path.")
}

future_results <- data.frame()

for (ssp in ssp_scenarios) {
  
  ssp_files <- future_files_all[
    grepl(selected_gcm, future_files_all, ignore.case = TRUE) &
      grepl(ssp, future_files_all, ignore.case = TRUE) &
      grepl("BCCA_PT", future_files_all, ignore.case = TRUE)
  ]
  
  if (length(ssp_files) == 0) {
    warning(paste("No BCCA-PT future files found for", selected_gcm, ssp))
    next
  }
  
  cat("\nReading final BCCA-PT future files for:", selected_gcm, ssp, "\n")
  cat("Files found:", length(ssp_files), "\n")
  
  ssp_ts <- extract_area_mean(
    files = ssp_files,
    dataset_name = paste0(selected_gcm, "_", ssp, "_BCCA_PT")
  )
  
  ssp_ts$SSP <- ssp
  
  future_results <- rbind(future_results, ssp_ts)
}

if (nrow(future_results) == 0) {
  stop("ERROR: No future BCCA-PT time series created.")
}

cat("\nFuture BCCA-PT records:", nrow(future_results), "\n")
cat("SSP245 records:", sum(future_results$SSP == "ssp245"), "\n")
cat("SSP585 records:", sum(future_results$SSP == "ssp585"), "\n")

if (nrow(future_results) != 1224) {
  warning("Expected 1224 future monthly records: 612 for ssp245 and 612 for ssp585.")
}

# ============================================================
# 9. ASSIGN FUTURE PERIODS
# ============================================================

future_results$period <- NA

for (i in 1:nrow(future_periods)) {
  
  period_name <- future_periods$period[i]
  sy <- future_periods$start_year[i]
  ey <- future_periods$end_year[i]
  
  future_results$period[
    future_results$year >= sy &
      future_results$year <= ey
  ] <- period_name
}

future_results <- future_results[!is.na(future_results$period), ]

write.csv(
  future_results,
  file.path(output_folder, "Tables", "future_BCCA_PT_MIROC6_monthly_timeseries.csv"),
  row.names = FALSE
)

# ============================================================
# 10. SAFE FUNCTION FOR PERCENTAGE CHANGE
# ============================================================

safe_percentage_change <- function(future, historical) {
  ifelse(
    is.na(historical) | historical <= 0 | is.na(future),
    NA_real_,
    ((future - historical) / historical) * 100
  )
}

# ============================================================
# 11. CALCULATE MONTHLY FUTURE MEAN AND PERCENTAGE CHANGE
# ============================================================

future_monthly_change <- future_results %>%
  group_by(SSP, period, month) %>%
  summarise(
    future_mean_rainfall = mean(rainfall, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    baseline_monthly[, c("month", "historical_mean_rainfall")],
    by = "month"
  ) %>%
  mutate(
    percentage_change = safe_percentage_change(
      future_mean_rainfall,
      historical_mean_rainfall
    ),
    Month = factor(month_names[month], levels = month_names)
  ) %>%
  filter(is.finite(percentage_change))

write.csv(
  future_monthly_change,
  file.path(output_folder, "Tables", "monthly_percentage_change_future_vs_historical.csv"),
  row.names = FALSE
)

cat("\nMonthly percentage change table:\n")
print(future_monthly_change)

# ============================================================
# 12. PLOT MONTHLY PERCENTAGE CHANGE BARPLOT
# ============================================================

p1 <- ggplot(
  future_monthly_change,
  aes(x = Month, y = percentage_change, fill = SSP)
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "black"
  ) +
  facet_wrap(~ period, ncol = 1) +
  geom_hline(yintercept = 0, linewidth = 0.7) +
  labs(
    title = "Projected Monthly Rainfall Change over the Upper Blue Nile Basin",
    subtitle = "Final BCCA-PT downscaled MIROC6 rainfall relative to the 1994–2014 XGBoost-corrected CHIRPS baseline",
    x = "Month",
    y = "Percentage Change (%)",
    fill = "Scenario"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

print(p1)

ggsave(
  filename = file.path(output_folder, "Figures", "monthly_percentage_change_barplot.png"),
  plot = p1,
  width = 12,
  height = 8,
  dpi = 300
)

# ============================================================
# 13. BOXPLOT OF MONTHLY FUTURE CHANGES
# ============================================================

future_results_change <- future_results %>%
  left_join(
    baseline_monthly[, c("month", "historical_mean_rainfall")],
    by = "month"
  ) %>%
  mutate(
    percentage_change = safe_percentage_change(
      rainfall,
      historical_mean_rainfall
    ),
    Month = factor(month_names[month], levels = month_names)
  )

problem_rows <- future_results_change %>%
  filter(!is.finite(percentage_change) | is.na(percentage_change))

cat("\nProblematic rows removed from boxplot:", nrow(problem_rows), "\n")

write.csv(
  problem_rows,
  file.path(output_folder, "Tables", "removed_nonfinite_percentage_change_rows.csv"),
  row.names = FALSE
)

future_results_change_clean <- future_results_change %>%
  filter(is.finite(percentage_change))

write.csv(
  future_results_change_clean,
  file.path(output_folder, "Tables", "future_monthly_percentage_change_clean.csv"),
  row.names = FALSE
)

p2 <- ggplot(
  future_results_change_clean,
  aes(x = Month, y = percentage_change, fill = SSP)
) +
  geom_boxplot(
    outlier.size = 0.8,
    alpha = 0.75
  ) +
  facet_wrap(~ period, ncol = 1) +
  geom_hline(yintercept = 0, linewidth = 0.7) +
  labs(
    title = "Distribution of Projected Monthly Rainfall Change",
    subtitle = "Final BCCA-PT downscaled MIROC6 rainfall relative to the historical CHIRPS baseline",
    x = "Month",
    y = "Percentage Change (%)",
    fill = "Scenario"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

print(p2)

ggsave(
  filename = file.path(output_folder, "Figures", "monthly_percentage_change_boxplot_clean.png"),
  plot = p2,
  width = 12,
  height = 8,
  dpi = 300
)

# ============================================================
# 14. JJAS MONSOONAL CHANGE
# Upper Blue Nile main rainy season: June–September
# ============================================================

jjas_change <- future_results %>%
  filter(month %in% c(6, 7, 8, 9)) %>%
  group_by(SSP, period, year) %>%
  summarise(
    JJAS_future_rainfall = sum(rainfall, na.rm = TRUE),
    .groups = "drop"
  )

historical_jjas <- chirps_ts %>%
  filter(month %in% c(6, 7, 8, 9)) %>%
  group_by(year) %>%
  summarise(
    JJAS_historical_rainfall = sum(rainfall, na.rm = TRUE),
    .groups = "drop"
  )

historical_jjas_mean <- mean(historical_jjas$JJAS_historical_rainfall, na.rm = TRUE)

jjas_change_summary <- jjas_change %>%
  group_by(SSP, period) %>%
  summarise(
    future_JJAS_mean = mean(JJAS_future_rainfall, na.rm = TRUE),
    historical_JJAS_mean = historical_jjas_mean,
    JJAS_percentage_change = safe_percentage_change(
      future_JJAS_mean,
      historical_JJAS_mean
    ),
    .groups = "drop"
  )

write.csv(
  jjas_change_summary,
  file.path(output_folder, "Tables", "JJAS_monsoonal_rainfall_change.csv"),
  row.names = FALSE
)

cat("\nJJAS rainfall change summary:\n")
print(jjas_change_summary)

# ============================================================
# 15. SAVE SUMMARY TEXT
# ============================================================

summary_text <- paste0(
  "Temporal changes in monthly rainfall were assessed by comparing final BCCA-PT downscaled ",
  selected_gcm,
  " future rainfall projections under SSP2-4.5 and SSP5-8.5 with the historical ",
  "XGBoost-corrected CHIRPS baseline period of 1994–2014. Percentage changes were ",
  "calculated for each month for the near-future period 2030–2050 and far-future ",
  "period 2051–2080. The analysis also summarized changes during the main rainy ",
  "season over the Upper Blue Nile Basin, defined as June to September."
)

writeLines(
  summary_text,
  file.path(output_folder, "Tables", "temporal_monthly_rainfall_change_summary.txt")
)

# ============================================================
# 16. FINAL MESSAGE
# ============================================================

cat("\nDONE SUCCESSFULLY!\n")
cat("Temporal rainfall change results saved in:\n")
cat(output_folder, "\n\n")

cat("Main output tables:\n")
cat(file.path(output_folder, "Tables", "monthly_percentage_change_future_vs_historical.csv"), "\n")
cat(file.path(output_folder, "Tables", "future_monthly_percentage_change_clean.csv"), "\n")
cat(file.path(output_folder, "Tables", "JJAS_monsoonal_rainfall_change.csv"), "\n\n")

cat("Main output figures:\n")
cat(file.path(output_folder, "Figures", "monthly_percentage_change_barplot.png"), "\n")
cat(file.path(output_folder, "Figures", "monthly_percentage_change_boxplot_clean.png"), "\n")
