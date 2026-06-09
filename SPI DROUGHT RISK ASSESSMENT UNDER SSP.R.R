# ============================================================
# SPI DROUGHT RISK ASSESSMENT UNDER SSP SCENARIOS
# Study area: Upper Blue Nile Basin
# Historical baseline: XGBoost-corrected CHIRPS, 1994–2014
# Future rainfall: Final BCCA-PT downscaled MIROC6, 2030–2080
# Scenarios: SSP2-4.5 and SSP5-8.5
# SPI scales: SPI-1, SPI-3, SPI-6, SPI-12
# ============================================================

# ============================================================
# 1. LOAD PACKAGES
# ============================================================

packages <- c(
  "terra", "dplyr", "ggplot2", "tidyr",
  "zoo", "lubridate", "scales"
)

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# ============================================================
# 2. INPUT FOLDERS
# ============================================================

historical_folder <- "D:/Desktop/HWRM_Thesis/From Gh/XGBOOST_corrected"

future_folder <- "D:/Desktop/HWRM_Thesis/From Gh/BCCA_PT_Downscaled_MIROC6_FINAL/03_Future_Downscaled"

output_folder <- "D:/Desktop/HWRM_Thesis/From Gh/SPI_Drought_Risk_Assessment_FINAL"

dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_folder, "Tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_folder, "Figures"), recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 3. SETTINGS
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

# SPI accumulation periods
# SPI-1  = short-term meteorological drought
# SPI-3  = seasonal drought
# SPI-6  = medium-term agricultural/hydrological drought
# SPI-12 = long-term hydrological drought
spi_scales <- c(1, 3, 6, 12)

# ============================================================
# 4. FUNCTIONS
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

make_inventory <- function(folder, dataset_name) {
  
  files <- list.files(
    path = folder,
    pattern = "\\.(tif|tiff)$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  
  files <- files[
    !grepl("BACKUP", files, ignore.case = TRUE) &
      !grepl("BAD", files, ignore.case = TRUE) &
      !grepl("Tables", files, ignore.case = TRUE) &
      !grepl("GeoTIFF", files, ignore.case = TRUE)
  ]
  
  if (length(files) == 0) {
    stop(paste("ERROR: No valid tif files found in:", folder))
  }
  
  inv <- data.frame(
    date_id = sapply(files, extract_year_month),
    file = files,
    dataset = dataset_name
  )
  
  inv <- inv[!is.na(inv$date_id), ]
  inv$year <- as.numeric(substr(inv$date_id, 1, 4))
  inv$month <- as.numeric(substr(inv$date_id, 6, 7))
  inv$date <- as.Date(paste0(inv$date_id, "-01"))
  inv <- inv[order(inv$date), ]
  
  return(inv)
}

extract_basin_mean <- function(inv_table, dataset_name) {
  
  out <- data.frame()
  
  for (i in 1:nrow(inv_table)) {
    
    r <- terra::rast(inv_table$file[i])
    
    mean_value <- terra::global(
      r,
      fun = "mean",
      na.rm = TRUE
    )[1, 1]
    
    out <- rbind(
      out,
      data.frame(
        date_id = inv_table$date_id[i],
        date = inv_table$date[i],
        year = inv_table$year[i],
        month = inv_table$month[i],
        rainfall = mean_value,
        dataset = dataset_name,
        file = inv_table$file[i]
      )
    )
    
    if (i %% 100 == 0) {
      cat(dataset_name, "processed:", i, "of", nrow(inv_table), "\n")
    }
  }
  
  return(out)
}

add_accumulated_rainfall <- function(df, scale_months) {
  
  df <- df %>%
    arrange(date) %>%
    mutate(
      rainfall_accumulated = zoo::rollsum(
        rainfall,
        k = scale_months,
        fill = NA,
        align = "right"
      )
    )
  
  return(df)
}

fit_spi_gamma_parameters <- function(hist_df, scale_months) {
  
  hist_acc <- add_accumulated_rainfall(hist_df, scale_months)
  
  params <- data.frame()
  
  for (m in 1:12) {
    
    x <- hist_acc$rainfall_accumulated[
      hist_acc$month == m &
        is.finite(hist_acc$rainfall_accumulated)
    ]
    
    x <- x[x >= 0]
    
    if (length(x) < 10) {
      warning(paste("Too few historical values for SPI scale", scale_months, "month", m))
      next
    }
    
    zero_probability <- sum(x == 0) / length(x)
    x_pos <- x[x > 0]
    
    if (length(x_pos) < 10) {
      shape <- NA
      scale <- NA
    } else {
      mean_x <- mean(x_pos, na.rm = TRUE)
      var_x <- var(x_pos, na.rm = TRUE)
      
      shape <- (mean_x^2) / var_x
      scale <- var_x / mean_x
    }
    
    params <- rbind(
      params,
      data.frame(
        SPI_scale = scale_months,
        month = m,
        shape = shape,
        scale = scale,
        zero_probability = zero_probability
      )
    )
  }
  
  return(params)
}

calculate_spi <- function(df, params, scale_months) {
  
  df_acc <- add_accumulated_rainfall(df, scale_months)
  
  df_acc$SPI_scale <- scale_months
  df_acc$SPI <- NA_real_
  
  for (i in 1:nrow(df_acc)) {
    
    x <- df_acc$rainfall_accumulated[i]
    m <- df_acc$month[i]
    
    if (!is.finite(x)) next
    
    p_row <- params %>%
      filter(SPI_scale == scale_months, month == m)
    
    if (nrow(p_row) == 0) next
    
    shape <- p_row$shape[1]
    scale <- p_row$scale[1]
    zero_probability <- p_row$zero_probability[1]
    
    if (!is.finite(shape) | !is.finite(scale)) next
    
    if (x <= 0) {
      probability <- zero_probability
    } else {
      probability <- zero_probability +
        (1 - zero_probability) * pgamma(
          x,
          shape = shape,
          scale = scale
        )
    }
    
    probability <- pmin(pmax(probability, 0.0001), 0.9999)
    
    df_acc$SPI[i] <- qnorm(probability)
  }
  
  df_acc$SPI_class <- dplyr::case_when(
    df_acc$SPI <= -2.0 ~ "Extreme drought",
    df_acc$SPI <= -1.5 ~ "Severe drought",
    df_acc$SPI <= -1.0 ~ "Moderate drought",
    df_acc$SPI < 1.0 ~ "Near normal",
    df_acc$SPI < 1.5 ~ "Moderately wet",
    df_acc$SPI < 2.0 ~ "Very wet",
    df_acc$SPI >= 2.0 ~ "Extremely wet",
    TRUE ~ NA_character_
  )
  
  return(df_acc)
}

extract_drought_events <- function(df) {
  
  df <- df %>%
    arrange(date) %>%
    mutate(
      is_drought = ifelse(is.finite(SPI) & SPI <= -1.0, 1, 0),
      event_start = is_drought == 1 & lag(is_drought, default = 0) == 0,
      event_id = cumsum(event_start)
    )
  
  events <- df %>%
    filter(is_drought == 1) %>%
    group_by(SSP_label, period_label, SPI_scale, event_id) %>%
    summarise(
      start_date = min(date),
      end_date = max(date),
      duration_months = n(),
      minimum_SPI = min(SPI, na.rm = TRUE),
      severity = sum(abs(SPI), na.rm = TRUE),
      .groups = "drop"
    )
  
  return(events)
}

# ============================================================
# 5. READ HISTORICAL BASELINE
# ============================================================

hist_inv <- make_inventory(historical_folder, "XGBoost_CHIRPS")

hist_inv <- hist_inv %>%
  filter(year >= historical_start, year <= historical_end)

cat("Historical files:", nrow(hist_inv), "\n")

if (nrow(hist_inv) != 252) {
  warning("Expected 252 historical monthly files for 1994–2014.")
}

hist_ts <- extract_basin_mean(
  hist_inv,
  dataset_name = "Historical_XGBoost_CHIRPS"
)

write.csv(
  hist_ts,
  file.path(output_folder, "Tables", "historical_basin_mean_rainfall_1994_2014.csv"),
  row.names = FALSE
)

# ============================================================
# 6. READ FUTURE DOWNSCALED RAINFALL
# ============================================================

future_inv <- make_inventory(
  future_folder,
  dataset_name = "Future_BCCA_PT_MIROC6"
)

future_inv <- future_inv %>%
  filter(
    grepl(selected_gcm, file, ignore.case = TRUE),
    grepl("BCCA_PT", file, ignore.case = TRUE),
    year >= 2030,
    year <= 2080
  )

cat("Future files:", nrow(future_inv), "\n")
cat("SSP245:", sum(grepl("ssp245", future_inv$file, ignore.case = TRUE)), "\n")
cat("SSP585:", sum(grepl("ssp585", future_inv$file, ignore.case = TRUE)), "\n")

if (nrow(future_inv) != 1224) {
  warning("Expected 1224 future files.")
}

future_ts <- data.frame()

for (ssp in ssp_scenarios) {
  
  ssp_inv <- future_inv %>%
    filter(grepl(ssp, file, ignore.case = TRUE)) %>%
    arrange(date)
  
  if (nrow(ssp_inv) == 0) {
    warning(paste("No future files found for", ssp))
    next
  }
  
  ssp_ts <- extract_basin_mean(
    ssp_inv,
    dataset_name = paste0("Future_", selected_gcm, "_", ssp, "_BCCA_PT")
  )
  
  ssp_ts$SSP <- ssp
  ssp_ts$SSP_label <- ifelse(
    ssp == "ssp245",
    "SSP2-4.5",
    "SSP5-8.5"
  )
  
  future_ts <- rbind(future_ts, ssp_ts)
}

future_ts$period <- NA
future_ts$period_label <- NA

for (i in 1:nrow(future_periods)) {
  
  period_name <- future_periods$period[i]
  sy <- future_periods$start_year[i]
  ey <- future_periods$end_year[i]
  
  future_ts$period[
    future_ts$year >= sy &
      future_ts$year <= ey
  ] <- period_name
  
  future_ts$period_label[
    future_ts$year >= sy &
      future_ts$year <= ey
  ] <- ifelse(
    period_name == "Near_future_2030_2050",
    "Near future (2030–2050)",
    "Far future (2051–2080)"
  )
}

future_ts <- future_ts %>%
  filter(!is.na(period))

write.csv(
  future_ts,
  file.path(output_folder, "Tables", "future_basin_mean_rainfall_for_SPI.csv"),
  row.names = FALSE
)

# ============================================================
# 7. FIT SPI PARAMETERS FROM HISTORICAL BASELINE
# ============================================================

spi_params <- data.frame()

for (s in spi_scales) {
  
  params_s <- fit_spi_gamma_parameters(
    hist_ts,
    scale_months = s
  )
  
  spi_params <- rbind(spi_params, params_s)
}

write.csv(
  spi_params,
  file.path(output_folder, "Tables", "SPI_gamma_parameters_1994_2014_baseline.csv"),
  row.names = FALSE
)

# ============================================================
# 8. CALCULATE FUTURE SPI
# ============================================================

future_spi_all <- data.frame()

for (s in spi_scales) {
  
  for (ssp in ssp_scenarios) {
    
    future_ssp <- future_ts %>%
      filter(SSP == ssp) %>%
      arrange(date)
    
    spi_ssp <- calculate_spi(
      df = future_ssp,
      params = spi_params,
      scale_months = s
    )
    
    future_spi_all <- rbind(future_spi_all, spi_ssp)
  }
}

write.csv(
  future_spi_all,
  file.path(output_folder, "Tables", "future_SPI_timeseries.csv"),
  row.names = FALSE
)

# ============================================================
# 9. DROUGHT RISK SUMMARY
# ============================================================

drought_risk_summary <- future_spi_all %>%
  group_by(SSP_label, period, period_label, SPI_scale) %>%
  summarise(
    total_valid_months = sum(is.finite(SPI)),
    moderate_drought_months = sum(SPI <= -1.0, na.rm = TRUE),
    severe_drought_months = sum(SPI <= -1.5, na.rm = TRUE),
    extreme_drought_months = sum(SPI <= -2.0, na.rm = TRUE),
    drought_frequency_percent = 100 * moderate_drought_months / total_valid_months,
    severe_frequency_percent = 100 * severe_drought_months / total_valid_months,
    extreme_frequency_percent = 100 * extreme_drought_months / total_valid_months,
    mean_SPI = mean(SPI, na.rm = TRUE),
    minimum_SPI = min(SPI, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  drought_risk_summary,
  file.path(output_folder, "Tables", "future_SPI_drought_risk_summary.csv"),
  row.names = FALSE
)

cat("\nDrought risk summary:\n")
print(drought_risk_summary)

# ============================================================
# 10. DROUGHT CATEGORY SUMMARY
# ============================================================

drought_category_summary <- future_spi_all %>%
  filter(is.finite(SPI)) %>%
  group_by(SSP_label, period_label, SPI_scale, SPI_class) %>%
  summarise(
    months = n(),
    .groups = "drop"
  ) %>%
  group_by(SSP_label, period_label, SPI_scale) %>%
  mutate(
    percentage = 100 * months / sum(months)
  ) %>%
  ungroup()

write.csv(
  drought_category_summary,
  file.path(output_folder, "Tables", "future_SPI_drought_category_summary.csv"),
  row.names = FALSE
)

# ============================================================
# 11. DROUGHT EVENT ANALYSIS
# Event = continuous period with SPI <= -1.0
# Corrected version: keeps grouping variables
# ============================================================

extract_drought_events <- function(df) {
  
  df <- df %>%
    arrange(date) %>%
    mutate(
      is_drought = ifelse(is.finite(SPI) & SPI <= -1.0, 1, 0),
      event_start = is_drought == 1 & lag(is_drought, default = 0) == 0,
      event_id = cumsum(event_start)
    )
  
  events <- df %>%
    filter(is_drought == 1) %>%
    group_by(event_id) %>%
    summarise(
      start_date = min(date),
      end_date = max(date),
      duration_months = n(),
      minimum_SPI = min(SPI, na.rm = TRUE),
      severity = sum(abs(SPI), na.rm = TRUE),
      .groups = "drop"
    )
  
  return(events)
}

drought_events <- future_spi_all %>%
  filter(is.finite(SPI)) %>%
  group_by(SSP_label, period_label, SPI_scale) %>%
  group_modify(~ extract_drought_events(.x)) %>%
  ungroup()

write.csv(
  drought_events,
  file.path(output_folder, "Tables", "future_SPI_drought_events.csv"),
  row.names = FALSE
)

cat("\nDrought events extracted successfully.\n")
cat("Number of drought events:", nrow(drought_events), "\n")
# ============================================================
# 12. FIGURES
# ============================================================

# SPI time series
p_spi <- ggplot(
  future_spi_all,
  aes(x = date, y = SPI, color = SSP_label)
) +
  geom_line(linewidth = 0.35) +
  geom_hline(yintercept = -1.0, linetype = "dashed", color = "orange") +
  geom_hline(yintercept = -1.5, linetype = "dashed", color = "red") +
  geom_hline(yintercept = -2.0, linetype = "dashed", color = "darkred") +
  facet_grid(
    SPI_scale ~ period_label,
    scales = "free_x",
    labeller = labeller(
      SPI_scale = function(x) paste0("SPI-", x)
    )
  ) +
  labs(
    title = "Projected SPI drought conditions over the Upper Blue Nile Basin",
    subtitle = "SPI-1, SPI-3, SPI-6, and SPI-12 calculated using 1994–2014 XGBoost-corrected CHIRPS baseline parameters",
    x = "Year",
    y = "SPI",
    color = "Scenario"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )

ggsave(
  filename = file.path(output_folder, "Figures", "future_SPI_timeseries_SPI1_SPI3_SPI6_SPI12.png"),
  plot = p_spi,
  width = 14,
  height = 11,
  dpi = 400
)

# Drought frequency
p_freq <- ggplot(
  drought_risk_summary,
  aes(x = period_label, y = drought_frequency_percent, fill = SSP_label)
) +
  geom_col(
    position = position_dodge(width = 0.8),
    color = "black",
    width = 0.7
  ) +
  facet_wrap(
    ~ SPI_scale,
    ncol = 2,
    labeller = labeller(
      SPI_scale = function(x) paste0("SPI-", x)
    )
  ) +
  labs(
    title = "Projected drought frequency over the Upper Blue Nile Basin",
    subtitle = "Drought frequency is based on months with SPI ≤ -1.0",
    x = "Future period",
    y = "Drought frequency (%)",
    fill = "Scenario"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 10, hjust = 1),
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )

ggsave(
  filename = file.path(output_folder, "Figures", "future_SPI_drought_frequency_SPI1_SPI3_SPI6_SPI12.png"),
  plot = p_freq,
  width = 12,
  height = 8,
  dpi = 400
)

# Minimum SPI
p_min <- ggplot(
  drought_risk_summary,
  aes(x = period_label, y = minimum_SPI, fill = SSP_label)
) +
  geom_col(
    position = position_dodge(width = 0.8),
    color = "black",
    width = 0.7
  ) +
  geom_hline(yintercept = -1.0, linetype = "dashed", color = "orange") +
  geom_hline(yintercept = -1.5, linetype = "dashed", color = "red") +
  geom_hline(yintercept = -2.0, linetype = "dashed", color = "darkred") +
  facet_wrap(
    ~ SPI_scale,
    ncol = 2,
    labeller = labeller(
      SPI_scale = function(x) paste0("SPI-", x)
    )
  ) +
  labs(
    title = "Minimum projected SPI over the Upper Blue Nile Basin",
    subtitle = "Lower SPI values indicate stronger drought conditions",
    x = "Future period",
    y = "Minimum SPI",
    fill = "Scenario"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 10, hjust = 1),
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )

ggsave(
  filename = file.path(output_folder, "Figures", "future_minimum_SPI_SPI1_SPI3_SPI6_SPI12.png"),
  plot = p_min,
  width = 12,
  height = 8,
  dpi = 400
)

# ============================================================
# 13. SAVE METHOD SUMMARY
# ============================================================

summary_text <- paste0(
  "Drought risk under SSP2-4.5 and SSP5-8.5 was assessed using SPI-1, SPI-3, SPI-6, and SPI-12. ",
  "The historical XGBoost-corrected CHIRPS rainfall for 1994–2014 was used as the baseline ",
  "for fitting gamma distribution parameters. These baseline parameters were then applied ",
  "to the final BCCA-PT downscaled MIROC6 future rainfall for 2030–2080. Drought risk was ",
  "summarized using standard SPI thresholds: moderate drought (SPI <= -1.0), severe drought ",
  "(SPI <= -1.5), and extreme drought (SPI <= -2.0). The outputs include future SPI time series, ",
  "drought frequency, drought categories, minimum SPI, and drought-event characteristics."
)

writeLines(
  summary_text,
  file.path(output_folder, "Tables", "SPI_drought_risk_method_summary.txt")
)

# ============================================================
# 14. FINAL MESSAGE
# ============================================================

cat("\nDONE SUCCESSFULLY!\n")
cat("SPI drought risk outputs saved in:\n")
cat(output_folder, "\n\n")

cat("Main tables:\n")
cat(file.path(output_folder, "Tables", "future_SPI_timeseries.csv"), "\n")
cat(file.path(output_folder, "Tables", "future_SPI_drought_risk_summary.csv"), "\n")
cat(file.path(output_folder, "Tables", "future_SPI_drought_category_summary.csv"), "\n")
cat(file.path(output_folder, "Tables", "future_SPI_drought_events.csv"), "\n\n")

cat("Main figures:\n")
cat(file.path(output_folder, "Figures", "future_SPI_timeseries_SPI1_SPI3_SPI6_SPI12.png"), "\n")
cat(file.path(output_folder, "Figures", "future_SPI_drought_frequency_SPI1_SPI3_SPI6_SPI12.png"), "\n")
cat(file.path(output_folder, "Figures", "future_minimum_SPI_SPI1_SPI3_SPI6_SPI12.png"), "\n")

