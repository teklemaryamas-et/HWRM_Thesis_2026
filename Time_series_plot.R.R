# =====================================================
# MONTHLY ENSEMBLE SPREAD AND VARIANCE FROM MANY GCM TIF FILES
# Historical period: 1994–2014
# Input: many monthly rainfall tif files from different GCMs
# Output: monthly ensemble mean, spread, variance, and range
# =====================================================

# Install and load terra
if (!require(terra)) {
  install.packages("terra")
  library(terra)
}

# =====================================================
# 1. YOUR FOLDERS
# =====================================================

input_folder <- "D:/Desktop/HWRM_Thesis/From Gh/NEX-GDDP-CMIP6"
output_folder <- "D:/Desktop/HWRM_Thesis/From Gh/Ensemble_Results"

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

if (!dir.exists(input_folder)) {
  stop("ERROR: Input folder does not exist. Check your folder path.")
}

# =====================================================
# 2. FIND ALL TIF FILES
# =====================================================

all_tif_files <- list.files(
  path = input_folder,
  pattern = "\\.(tif|tiff)$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

cat("Total tif files found:", length(all_tif_files), "\n")

if (length(all_tif_files) == 0) {
  stop("ERROR: No tif files found.")
}

# =====================================================
# 3. SELECT HISTORICAL PERIOD 1994–2014
# =====================================================

years <- 1994:2014
months <- 1:12

# =====================================================
# 4. LOOP THROUGH EACH YEAR AND MONTH
# =====================================================

for (yr in years) {
  
  for (mn in months) {
    
    month_2digit <- sprintf("%02d", mn)
    
    # -------------------------------------------------
    # This searches filenames containing year and month
    # Examples it can detect:
    # 1994_01, 1994-01, 199401, 1994.01
    # -------------------------------------------------
    
    pattern1 <- paste0(yr, "_", month_2digit)
    pattern2 <- paste0(yr, "-", month_2digit)
    pattern3 <- paste0(yr, month_2digit)
    pattern4 <- paste0(yr, "\\.", month_2digit)
    
    monthly_files <- all_tif_files[
      grepl(pattern1, basename(all_tif_files), ignore.case = TRUE) |
        grepl(pattern2, basename(all_tif_files), ignore.case = TRUE) |
        grepl(pattern3, basename(all_tif_files), ignore.case = TRUE) |
        grepl(pattern4, basename(all_tif_files), ignore.case = TRUE)
    ]
    
    cat("\nProcessing:", yr, "-", month_2digit, "\n")
    cat("Number of GCM files found for this month:", length(monthly_files), "\n")
    
    if (length(monthly_files) < 2) {
      cat("Skipped: fewer than 2 files found for this month.\n")
      next
    }
    
    # Optional warning
    if (length(monthly_files) != 12) {
      warning(
        paste(
          "For", yr, month_2digit,
          "expected 12 GCM files but found",
          length(monthly_files)
        )
      )
    }
    
    # =================================================
    # 5. READ MONTHLY FILES AS STACK
    # =================================================
    
    gcm_stack <- rast(monthly_files)
    
    # =================================================
    # 6. CALCULATE ENSEMBLE STATISTICS
    # =================================================
    
    ensemble_mean <- mean(gcm_stack, na.rm = TRUE)
    
    ensemble_spread <- app(gcm_stack, fun = sd, na.rm = TRUE)
    
    ensemble_variance <- app(gcm_stack, fun = var, na.rm = TRUE)
    
    ensemble_range <- max(gcm_stack, na.rm = TRUE) - min(gcm_stack, na.rm = TRUE)
    
    # =================================================
    # 7. SAVE OUTPUT FILES
    # =================================================
    
    writeRaster(
      ensemble_mean,
      filename = file.path(
        output_folder,
        paste0("ensemble_mean_", yr, "_", month_2digit, ".tif")
      ),
      overwrite = TRUE
    )
    
    writeRaster(
      ensemble_spread,
      filename = file.path(
        output_folder,
        paste0("ensemble_spread_std_", yr, "_", month_2digit, ".tif")
      ),
      overwrite = TRUE
    )
    
    writeRaster(
      ensemble_variance,
      filename = file.path(
        output_folder,
        paste0("ensemble_variance_", yr, "_", month_2digit, ".tif")
      ),
      overwrite = TRUE
    )
    
    writeRaster(
      ensemble_range,
      filename = file.path(
        output_folder,
        paste0("ensemble_range_", yr, "_", month_2digit, ".tif")
      ),
      overwrite = TRUE
    )
    
    cat("Finished:", yr, "-", month_2digit, "\n")
  }
}

cat("\nALL DONE SUCCESSFULLY!\n")
cat("Results saved in:\n")
cat(output_folder, "\n")
# =====================================================
# CHIRPS OBSERVED VS ENSEMBLE MEAN TIME SERIES
# Period: 1994–2014
# Output: CSV table + PNG time series plot
# =====================================================

# Install packages if missing
if (!require(terra)) {
  install.packages("terra")
  library(terra)
}

if (!require(ggplot2)) {
  install.packages("ggplot2")
  library(ggplot2)
}

if (!require(dplyr)) {
  install.packages("dplyr")
  library(dplyr)
}

# =====================================================
# 1. FOLDER PATHS
# =====================================================

# PASTE YOUR CHIRPS OBSERVED TIF FOLDER HERE
chirps_folder <- "D:/Desktop/HWRM_Thesis/From Gh/TM proposal/data/CHIRPS_Monthly_Abay/XGBOOST_corrected/resampled_0_25deg"

# Your ensemble mean folder
ensemble_folder <- "D:/Desktop/HWRM_Thesis/From Gh/Ensemble_Results"

# Output folder for graph and CSV
output_folder <- "D:/Desktop/HWRM_Thesis/From Gh/TimeSeries_Result"

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# =====================================================
# 2. CHECK FOLDERS
# =====================================================

if (!dir.exists(chirps_folder)) {
  stop("ERROR: CHIRPS folder does not exist. Please check chirps_folder path.")
}

if (!dir.exists(ensemble_folder)) {
  stop("ERROR: Ensemble folder does not exist. Please check ensemble_folder path.")
}

# =====================================================
# 3. FIND CHIRPS AND ENSEMBLE MEAN TIF FILES
# =====================================================

chirps_files <- list.files(
  path = chirps_folder,
  pattern = "\\.(tif|tiff)$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

ensemble_files <- list.files(
  path = ensemble_folder,
  pattern = "^ensemble_mean_.*\\.(tif|tiff)$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

cat("CHIRPS files found:", length(chirps_files), "\n")
cat("Ensemble mean files found:", length(ensemble_files), "\n")

if (length(chirps_files) == 0) {
  stop("ERROR: No CHIRPS tif files found.")
}

if (length(ensemble_files) == 0) {
  stop("ERROR: No ensemble_mean tif files found.")
}

# =====================================================
# 4. FUNCTION TO EXTRACT YEAR AND MONTH FROM FILE NAME
# Works with names like:
# 1994_01, 1994-01, 199401, 1994.01
# =====================================================

extract_year_month <- function(file_path) {
  
  fname <- basename(file_path)
  
  # Pattern: 1994_01 or 1994-01 or 1994.01
  match1 <- regmatches(
    fname,
    regexpr("(199[4-9]|200[0-9]|201[0-4])[_\\.-](0[1-9]|1[0-2])", fname)
  )
  
  if (length(match1) > 0 && match1 != "") {
    ym <- gsub("[_\\.-]", "-", match1)
    return(ym)
  }
  
  # Pattern: 199401
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

# =====================================================
# 5. CREATE FILE TABLES
# =====================================================

chirps_table <- data.frame(
  date_id = sapply(chirps_files, extract_year_month),
  chirps_file = chirps_files
)

ensemble_table <- data.frame(
  date_id = sapply(ensemble_files, extract_year_month),
  ensemble_file = ensemble_files
)

# Remove files without year-month
chirps_table <- chirps_table[!is.na(chirps_table$date_id), ]
ensemble_table <- ensemble_table[!is.na(ensemble_table$date_id), ]

# Keep only 1994–2014
chirps_table <- chirps_table[
  substr(chirps_table$date_id, 1, 4) >= "1994" &
    substr(chirps_table$date_id, 1, 4) <= "2014",
]

ensemble_table <- ensemble_table[
  substr(ensemble_table$date_id, 1, 4) >= "1994" &
    substr(ensemble_table$date_id, 1, 4) <= "2014",
]

# Merge CHIRPS and ensemble files by year-month
matched_table <- merge(
  chirps_table,
  ensemble_table,
  by = "date_id"
)

matched_table <- matched_table[order(matched_table$date_id), ]

cat("Matched monthly pairs:", nrow(matched_table), "\n")

if (nrow(matched_table) == 0) {
  stop("ERROR: No matching CHIRPS and ensemble_mean files found by year-month.")
}

# =====================================================
# 6. EXTRACT SPATIAL MEAN RAINFALL FOR EACH MONTH
# =====================================================

results <- data.frame(
  date = as.Date(character()),
  chirps_observed = numeric(),
  ensemble_mean = numeric()
)

for (i in 1:nrow(matched_table)) {
  
  date_id <- matched_table$date_id[i]
  
  cat("Processing:", date_id, "\n")
  
  chirps_r <- rast(matched_table$chirps_file[i])
  ensemble_r <- rast(matched_table$ensemble_file[i])
  
  # If CRS/resolution/extent are different, resample ensemble to CHIRPS grid
  if (!compareGeom(chirps_r, ensemble_r, stopOnError = FALSE)) {
    ensemble_r <- resample(ensemble_r, chirps_r, method = "bilinear")
  }
  
  chirps_mean_value <- global(chirps_r, fun = "mean", na.rm = TRUE)[1, 1]
  ensemble_mean_value <- global(ensemble_r, fun = "mean", na.rm = TRUE)[1, 1]
  
  results <- rbind(
    results,
    data.frame(
      date = as.Date(paste0(date_id, "-01")),
      chirps_observed = chirps_mean_value,
      ensemble_mean = ensemble_mean_value
    )
  )
}

# =====================================================
# 7. SAVE TIME SERIES TABLE
# =====================================================

csv_output <- file.path(
  output_folder,
  "CHIRPS_observed_vs_ensemble_mean_timeseries_1994_2014.csv"
)

write.csv(results, csv_output, row.names = FALSE)

# =====================================================
# 8. PREPARE DATA FOR PLOTTING
# =====================================================

plot_data <- data.frame(
  date = rep(results$date, 2),
  rainfall = c(results$chirps_observed, results$ensemble_mean),
  dataset = c(
    rep("CHIRPS Observed", nrow(results)),
    rep("Ensemble Mean", nrow(results))
  )
)

# =====================================================
# 9. DRAW TIME SERIES PLOT
# =====================================================

p <- ggplot(plot_data, aes(x = date, y = rainfall, color = dataset)) +
  geom_line(linewidth = 0.8) +
  labs(
    title = "CHIRPS Observed vs Multi-Model Ensemble Mean Rainfall",
    subtitle = "Monthly rainfall time series, 1994–2014",
    x = "Year",
    y = "Mean Monthly Rainfall(mm)",
    color = "Dataset"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

print(p)

# =====================================================
# 10. SAVE PLOT
# =====================================================

plot_output <- file.path(
  output_folder,
  "CHIRPS_observed_vs_ensemble_mean_timeseries_1994_2014.png"
)

ggsave(
  filename = plot_output,
  plot = p,
  width = 12,
  height = 6,
  dpi = 300
)

cat("\nDONE SUCCESSFULLY!\n")
cat("CSV saved at:\n", csv_output, "\n")
cat("Plot saved at:\n", plot_output, "\n")
# =====================================================
# CHIRPS OBSERVED VS ENSEMBLE MEAN WITH ENSEMBLE SPREAD
# Graph style similar to your example
# Period: 1994–2014
# =====================================================

library(terra)
library(ggplot2)
library(dplyr)

# =====================================================
# 1. FOLDER PATHS
# =====================================================

chirps_folder <- "D:/Desktop/HWRM_Thesis/From Gh/TM proposal/data/CHIRPS_Monthly_Abay/XGBOOST_corrected/resampled_0_25deg"
ensemble_folder <- "D:/Desktop/HWRM_Thesis/From Gh/Ensemble_Results"
output_folder <- "D:/Desktop/HWRM_Thesis/From Gh/TimeSeries_Result"

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# =====================================================
# 2. FIND FILES
# =====================================================

chirps_files <- list.files(
  chirps_folder,
  pattern = "\\.(tif|tiff)$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

ensemble_mean_files <- list.files(
  ensemble_folder,
  pattern = "^ensemble_mean_.*\\.(tif|tiff)$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

ensemble_spread_files <- list.files(
  ensemble_folder,
  pattern = "^ensemble_spread_std_.*\\.(tif|tiff)$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

cat("CHIRPS files:", length(chirps_files), "\n")
cat("Ensemble mean files:", length(ensemble_mean_files), "\n")
cat("Ensemble spread files:", length(ensemble_spread_files), "\n")

if (length(chirps_files) == 0) stop("No CHIRPS files found.")
if (length(ensemble_mean_files) == 0) stop("No ensemble mean files found.")
if (length(ensemble_spread_files) == 0) stop("No ensemble spread files found.")

# =====================================================
# 3. FUNCTION TO EXTRACT YEAR-MONTH FROM FILE NAME
# =====================================================

extract_year_month <- function(file_path) {
  
  fname <- basename(file_path)
  
  match1 <- regmatches(
    fname,
    regexpr("(199[4-9]|200[0-9]|201[0-4])[_\\.-](0[1-9]|1[0-2])", fname)
  )
  
  if (length(match1) > 0 && match1 != "") {
    ym <- gsub("[_\\.-]", "-", match1)
    return(ym)
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

# =====================================================
# 4. CREATE TABLES
# =====================================================

chirps_table <- data.frame(
  date_id = sapply(chirps_files, extract_year_month),
  chirps_file = chirps_files
)

mean_table <- data.frame(
  date_id = sapply(ensemble_mean_files, extract_year_month),
  ensemble_mean_file = ensemble_mean_files
)

spread_table <- data.frame(
  date_id = sapply(ensemble_spread_files, extract_year_month),
  ensemble_spread_file = ensemble_spread_files
)

chirps_table <- chirps_table[!is.na(chirps_table$date_id), ]
mean_table <- mean_table[!is.na(mean_table$date_id), ]
spread_table <- spread_table[!is.na(spread_table$date_id), ]

# Merge all by year-month
matched_table <- chirps_table %>%
  inner_join(mean_table, by = "date_id") %>%
  inner_join(spread_table, by = "date_id") %>%
  arrange(date_id)

cat("Matched months:", nrow(matched_table), "\n")

if (nrow(matched_table) == 0) {
  stop("No matching CHIRPS, ensemble mean, and ensemble spread files found.")
}

# =====================================================
# 5. EXTRACT AREA-AVERAGE VALUES
# =====================================================

results <- data.frame(
  date = as.Date(character()),
  chirps_observed = numeric(),
  ensemble_mean = numeric(),
  ensemble_spread = numeric()
)

for (i in 1:nrow(matched_table)) {
  
  date_id <- matched_table$date_id[i]
  cat("Processing:", date_id, "\n")
  
  chirps_r <- rast(matched_table$chirps_file[i])
  mean_r <- rast(matched_table$ensemble_mean_file[i])
  spread_r <- rast(matched_table$ensemble_spread_file[i])
  
  # Resample if grids are different
  if (!compareGeom(chirps_r, mean_r, stopOnError = FALSE)) {
    mean_r <- resample(mean_r, chirps_r, method = "bilinear")
  }
  
  if (!compareGeom(chirps_r, spread_r, stopOnError = FALSE)) {
    spread_r <- resample(spread_r, chirps_r, method = "bilinear")
  }
  
  chirps_value <- global(chirps_r, fun = "mean", na.rm = TRUE)[1, 1]
  mean_value <- global(mean_r, fun = "mean", na.rm = TRUE)[1, 1]
  spread_value <- global(spread_r, fun = "mean", na.rm = TRUE)[1, 1]
  
  results <- rbind(
    results,
    data.frame(
      date = as.Date(paste0(date_id, "-01")),
      chirps_observed = chirps_value,
      ensemble_mean = mean_value,
      ensemble_spread = spread_value
    )
  )
}

# =====================================================
# 6. KEEP ONLY 1994–2014
# =====================================================

results <- results[
  results$date >= as.Date("1994-01-01") &
    results$date <= as.Date("2014-12-31"),
]

# Create upper and lower ensemble spread bands
results$ensemble_lower <- results$ensemble_mean - results$ensemble_spread
results$ensemble_upper <- results$ensemble_mean + results$ensemble_spread

# Do not allow negative rainfall
results$ensemble_lower[results$ensemble_lower < 0] <- 0

# Save CSV
write.csv(
  results,
  file.path(output_folder, "CHIRPS_vs_ensemble_mean_with_spread_1994_2014.csv"),
  row.names = FALSE
)

# =====================================================
# 7. PLOT LIKE YOUR EXAMPLE
# =====================================================

p <- ggplot(results, aes(x = date)) +
  
  # Light blue uncertainty/spread area
  geom_ribbon(
    aes(ymin = ensemble_lower, ymax = ensemble_upper),
    fill = "royalblue",
    alpha = 0.25
  ) +
  
  # CHIRPS observed black line
  geom_line(
    aes(y = chirps_observed, color = "CHIRPS Observation"),
    linewidth = 0.7
  ) +
  
  # Ensemble mean blue line
  geom_line(
    aes(y = ensemble_mean, color = "Ensemble Mean"),
    linewidth = 0.8
  ) +
  
  scale_color_manual(
    values = c(
      "CHIRPS Observation" = "red",
      "Ensemble Mean" = "blue"
    )
  ) +
  
  scale_x_date(
    limits = c(as.Date("1994-01-01"), as.Date("2014-12-31")),
    breaks = seq(
      as.Date("1994-01-01"),
      as.Date("2014-01-01"),
      by = "5 years"
    ),
    date_labels = "%Y",
    expand = c(0, 0)
  ) +
  
  labs(
    x = "Date (Year-Month)",
    y = "Perciptation (mm)",
    color = NULL
  ) +
  
  theme_bw() +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_blank(),
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_line(color = "grey92"),
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 11),
    plot.margin = margin(10, 10, 10, 10)
  )

print(p)

# =====================================================
# 8. SAVE PLOT
# =====================================================

plot_output <- file.path(
  output_folder,
  "CHIRPS_observed_vs_ensemble_mean_with_spread_1994_2014.png"
)

ggsave(
  filename = plot_output,
  plot = p,
  width = 12,
  height = 5,
  dpi = 300
)

cat("\nDone successfully!\n")
cat("Plot saved at:\n", plot_output, "\n")

