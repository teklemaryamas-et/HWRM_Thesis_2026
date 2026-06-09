# ============================================================
# FAST BCCA-PT DOWNSCALING USING MONTHLY PT PARAMETERS
# Reference: XGBoost-corrected CHIRPS
# GCM: Resampled MIROC6
# Historical period: 1994–2014
# Future period: 2030–2080
# Scenarios: SSP2-4.5 and SSP5-8.5
# ============================================================

# ============================================================
# 1. LOAD PACKAGES
# ============================================================

packages <- c("terra", "dplyr")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# ============================================================
# 2. INPUT FOLDERS
# ============================================================

chirps_folder <- "D:/Desktop/HWRM_Thesis/From Gh/Gebres/TM proposal/data/CHIRPS_Monthly_Abay/XGBOOST_corrected/resampled_0_25deg"

historical_gcm_folder <- "D:/Desktop/HWRM_Thesis/From Gh/Resampled_MIROC6_to_CHIRPS/Historical_MIROC6_1994_2014"

future_gcm_folder <- "D:/Desktop/HWRM_Thesis/From Gh/Resampled_MIROC6_to_CHIRPS/Future_MIROC6_2030_2080"

basin_folder <- "D:/Desktop/HWRM_Thesis/From Gh/Gebres/TM proposal/data/CHIRPS_Monthly_Abay/XGBOOST_corrected/abay"

output_folder <- "D:/Desktop/HWRM_Thesis/From Gh/BCCA_PT_Downscaled_MIROC6_FAST"

dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_folder, "01_PT_Parameters"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_folder, "02_Historical_Downscaled"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_folder, "03_Future_Downscaled"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_folder, "04_Tables"), recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 3. SETTINGS
# ============================================================

selected_gcm <- "MIROC6"
ssp_scenarios <- c("ssp245", "ssp585")

historical_start <- 1994
historical_end <- 2014

future_start <- 2030
future_end <- 2080

eps <- 0.001

# ============================================================
# 4. LOAD BASIN SHAPEFILE
# ============================================================

basin_shp_files <- list.files(
  path = basin_folder,
  pattern = "\\.shp$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

if (length(basin_shp_files) == 0) {
  warning("No basin shapefile found. Code will continue without masking.")
  basin_vect <- NULL
} else {
  basin_shapefile <- basin_shp_files[1]
  cat("Using basin shapefile:\n", basin_shapefile, "\n")
  basin_vect <- terra::vect(basin_shapefile)
}

# ============================================================
# 5. FUNCTIONS
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
  
  if (length(files) == 0) {
    stop(paste("ERROR: No tif files found in:", folder))
  }
  
  inv <- data.frame(
    dataset = dataset_name,
    date_id = sapply(files, extract_year_month),
    file = files
  )
  
  inv <- inv[!is.na(inv$date_id), ]
  inv$year <- as.numeric(substr(inv$date_id, 1, 4))
  inv$month <- as.numeric(substr(inv$date_id, 6, 7))
  inv <- inv[order(inv$date_id), ]
  
  return(inv)
}

apply_pt <- function(r, a, b, eps = 0.001) {
  out <- a * ((r + eps) ^ b)
  out[out < 0] <- 0
  return(out)
}

# ============================================================
# 6. INVENTORIES
# ============================================================

chirps_inv <- make_inventory(chirps_folder, "CHIRPS_XGBoost")

chirps_inv <- chirps_inv[
  chirps_inv$year >= historical_start &
    chirps_inv$year <= historical_end,
]

hist_gcm_inv <- make_inventory(historical_gcm_folder, "MIROC6_resampled_historical")

hist_gcm_inv <- hist_gcm_inv[
  hist_gcm_inv$year >= historical_start &
    hist_gcm_inv$year <= historical_end &
    grepl(selected_gcm, hist_gcm_inv$file, ignore.case = TRUE),
]

cat("CHIRPS historical files:", nrow(chirps_inv), "\n")
cat("MIROC6 historical resampled files:", nrow(hist_gcm_inv), "\n")

matched_hist <- inner_join(
  chirps_inv[, c("date_id", "year", "month", "file")],
  hist_gcm_inv[, c("date_id", "file")],
  by = "date_id",
  suffix = c("_chirps", "_gcm")
)

matched_hist <- matched_hist[order(matched_hist$date_id), ]

cat("Matched historical months:", nrow(matched_hist), "\n")

if (nrow(matched_hist) < 30) {
  stop("ERROR: Too few matched historical months.")
}

write.csv(
  matched_hist,
  file.path(output_folder, "04_Tables", "matched_CHIRPS_MIROC6.csv"),
  row.names = FALSE
)

# ============================================================
# 7. PREPARE REFERENCE RASTER AND BASIN
# ============================================================

ref_raster <- rast(matched_hist$file_chirps[1])

if (!is.null(basin_vect)) {
  basin_vect <- terra::project(basin_vect, terra::crs(ref_raster))
}

# ============================================================
# 8. CALIBRATE MONTHLY PT PARAMETERS
# ============================================================

pt_params <- data.frame()

for (m in 1:12) {
  
  cat("Calibrating PT for month:", m, "\n")
  
  month_data <- matched_hist[matched_hist$month == m, ]
  
  obs_all <- c()
  raw_all <- c()
  
  for (i in 1:nrow(month_data)) {
    
    obs_r <- rast(month_data$file_chirps[i])
    raw_r <- rast(month_data$file_gcm[i])
    
    if (!compareGeom(obs_r, raw_r, stopOnError = FALSE)) {
      stop("ERROR: CHIRPS and MIROC6 are not aligned. Check resampling.")
    }
    
    if (!is.null(basin_vect)) {
      obs_r <- mask(crop(obs_r, basin_vect), basin_vect)
      raw_r <- mask(crop(raw_r, basin_vect), basin_vect)
    }
    
    obs_v <- values(obs_r, mat = FALSE)
    raw_v <- values(raw_r, mat = FALSE)
    
    good <- is.finite(obs_v) & is.finite(raw_v) & obs_v >= 0 & raw_v >= 0
    
    obs_all <- c(obs_all, obs_v[good])
    raw_all <- c(raw_all, raw_v[good])
  }
  
  if (length(obs_all) < 30) {
    warning(paste("Too few valid values for month", m, "- using a=1 and b=1"))
    a <- 1
    b <- 1
  } else {
    fit <- lm(log(obs_all + eps) ~ log(raw_all + eps))
    b <- as.numeric(coef(fit)[2])
    a <- exp(as.numeric(coef(fit)[1]))
    
    if (!is.finite(a) || !is.finite(b)) {
      a <- 1
      b <- 1
    }
  }
  
  pt_params <- rbind(
    pt_params,
    data.frame(
      month = m,
      a = a,
      b = b
    )
  )
}

write.csv(
  pt_params,
  file.path(output_folder, "01_PT_Parameters", "monthly_PT_parameters.csv"),
  row.names = FALSE
)

cat("\nMonthly PT parameters:\n")
print(pt_params)

# ============================================================
# 9. APPLY PT TO HISTORICAL MIROC6 FOR VALIDATION
# ============================================================

hist_out_folder <- file.path(output_folder, "02_Historical_Downscaled", selected_gcm)
dir.create(hist_out_folder, recursive = TRUE, showWarnings = FALSE)

validation_table <- data.frame()

for (i in 1:nrow(matched_hist)) {
  
  date_id <- matched_hist$date_id[i]
  m <- matched_hist$month[i]
  
  a <- pt_params$a[pt_params$month == m]
  b <- pt_params$b[pt_params$month == m]
  
  obs_r <- rast(matched_hist$file_chirps[i])
  raw_r <- rast(matched_hist$file_gcm[i])
  
  if (!is.null(basin_vect)) {
    obs_r <- mask(crop(obs_r, basin_vect), basin_vect)
    raw_r <- mask(crop(raw_r, basin_vect), basin_vect)
  }
  
  corrected_r <- apply_pt(raw_r, a, b, eps)
  
  out_file <- file.path(
    hist_out_folder,
    paste0(selected_gcm, "_historical_BCCA_PT_", date_id, ".tif")
  )
  
  writeRaster(
    corrected_r,
    out_file,
    overwrite = TRUE,
    wopt = list(datatype = "FLT4S", gdal = "COMPRESS=LZW")
  )
  
  validation_table <- rbind(
    validation_table,
    data.frame(
      date_id = date_id,
      year = matched_hist$year[i],
      month = m,
      CHIRPS_XGBoost = global(obs_r, "mean", na.rm = TRUE)[1, 1],
      MIROC6_resampled_raw = global(raw_r, "mean", na.rm = TRUE)[1, 1],
      MIROC6_BCCA_PT = global(corrected_r, "mean", na.rm = TRUE)[1, 1],
      downscaled_file = out_file
    )
  )
  
  if (i %% 25 == 0) {
    cat("Historical BCCA-PT completed:", i, "of", nrow(matched_hist), "\n")
  }
}

write.csv(
  validation_table,
  file.path(output_folder, "04_Tables", "historical_BCCA_PT_validation_timeseries.csv"),
  row.names = FALSE
)

# ============================================================
# 10. PERFORMANCE METRICS
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
  Method = c("Resampled_Raw_MIROC6", "BCCA_PT_MIROC6"),
  NSE = c(
    NSE(validation_table$CHIRPS_XGBoost, validation_table$MIROC6_resampled_raw),
    NSE(validation_table$CHIRPS_XGBoost, validation_table$MIROC6_BCCA_PT)
  ),
  RMSE = c(
    RMSE(validation_table$CHIRPS_XGBoost, validation_table$MIROC6_resampled_raw),
    RMSE(validation_table$CHIRPS_XGBoost, validation_table$MIROC6_BCCA_PT)
  ),
  MAE = c(
    MAE(validation_table$CHIRPS_XGBoost, validation_table$MIROC6_resampled_raw),
    MAE(validation_table$CHIRPS_XGBoost, validation_table$MIROC6_BCCA_PT)
  ),
  PBIAS = c(
    PBIAS(validation_table$CHIRPS_XGBoost, validation_table$MIROC6_resampled_raw),
    PBIAS(validation_table$CHIRPS_XGBoost, validation_table$MIROC6_BCCA_PT)
  ),
  r = c(
    R_value(validation_table$CHIRPS_XGBoost, validation_table$MIROC6_resampled_raw),
    R_value(validation_table$CHIRPS_XGBoost, validation_table$MIROC6_BCCA_PT)
  )
)

write.csv(
  metrics,
  file.path(output_folder, "04_Tables", "historical_BCCA_PT_performance_metrics.csv"),
  row.names = FALSE
)

cat("\nPerformance metrics:\n")
print(metrics)

# ============================================================
# 11. APPLY BCCA-PT TO FUTURE MIROC6
# ============================================================

future_inv <- make_inventory(future_gcm_folder, "MIROC6_resampled_future")

future_inv <- future_inv[
  future_inv$year >= future_start &
    future_inv$year <= future_end &
    grepl(selected_gcm, future_inv$file, ignore.case = TRUE),
]

cat("\nFuture resampled MIROC6 files:", nrow(future_inv), "\n")

if (nrow(future_inv) == 0) {
  stop("ERROR: No future resampled MIROC6 files found.")
}

future_mean_table <- data.frame()

for (ssp in ssp_scenarios) {
  
  ssp_inv <- future_inv[
    grepl(ssp, future_inv$file, ignore.case = TRUE),
  ]
  
  if (nrow(ssp_inv) == 0) {
    warning(paste("No files found for", ssp))
    next
  }
  
  ssp_inv <- ssp_inv[order(ssp_inv$date_id), ]
  
  ssp_out_folder <- file.path(
    output_folder,
    "03_Future_Downscaled",
    selected_gcm,
    ssp
  )
  
  dir.create(ssp_out_folder, recursive = TRUE, showWarnings = FALSE)
  
  cat("\nApplying BCCA-PT to future scenario:", ssp, "\n")
  
  for (i in 1:nrow(ssp_inv)) {
    
    date_id <- ssp_inv$date_id[i]
    m <- ssp_inv$month[i]
    
    a <- pt_params$a[pt_params$month == m]
    b <- pt_params$b[pt_params$month == m]
    
    raw_r <- rast(ssp_inv$file[i])
    
    if (!is.null(basin_vect)) {
      raw_r <- mask(crop(raw_r, basin_vect), basin_vect)
    }
    
    corrected_r <- apply_pt(raw_r, a, b, eps)
    
    out_file <- file.path(
      ssp_out_folder,
      paste0(selected_gcm, "_", ssp, "_BCCA_PT_", date_id, ".tif")
    )
    
    writeRaster(
      corrected_r,
      out_file,
      overwrite = TRUE,
      wopt = list(datatype = "FLT4S", gdal = "COMPRESS=LZW")
    )
    
    future_mean_table <- rbind(
      future_mean_table,
      data.frame(
        GCM = selected_gcm,
        SSP = ssp,
        date_id = date_id,
        year = ssp_inv$year[i],
        month = m,
        resampled_raw_future_mean = global(raw_r, "mean", na.rm = TRUE)[1, 1],
        BCCA_PT_future_mean = global(corrected_r, "mean", na.rm = TRUE)[1, 1],
        downscaled_file = out_file
      )
    )
    
    if (i %% 25 == 0) {
      cat(ssp, "future BCCA-PT completed:", i, "of", nrow(ssp_inv), "\n")
    }
  }
}

write.csv(
  future_mean_table,
  file.path(output_folder, "04_Tables", "future_BCCA_PT_MIROC6_basin_mean_timeseries.csv"),
  row.names = FALSE
)

# ============================================================
# 12. SAVE METHOD SUMMARY
# ============================================================

summary_text <- paste0(
  "BCCA-PT downscaling was applied after historical and future MIROC6 rainfall were ",
  "resampled to the XGBoost-corrected CHIRPS grid. The PT method was selected from ",
  "the comparison of GEQM, GQM, LS, and PT bias correction methods. Monthly PT ",
  "parameters were calibrated using historical resampled MIROC6 and XGBoost-corrected ",
  "CHIRPS rainfall for 1994-2014. The calibrated monthly PT parameters were then applied ",
  "to future MIROC6 rainfall under SSP2-4.5 and SSP5-8.5 for 2030-2080."
)

writeLines(
  summary_text,
  file.path(output_folder, "04_Tables", "BCCA_PT_downscaling_method_summary.txt")
)

# ============================================================
# 13. FINAL MESSAGE
# ============================================================

cat("\nDONE SUCCESSFULLY!\n")
cat("BCCA-PT downscaled outputs saved in:\n")
cat(output_folder, "\n\n")

cat("Use this folder for rainfall projection and SPI:\n")
cat(file.path(output_folder, "03_Future_Downscaled"), "\n\n")

cat("Validation metrics saved at:\n")
cat(file.path(output_folder, "04_Tables", "historical_BCCA_PT_performance_metrics.csv"), "\n\n")

cat("Future downscaled basin mean time series saved at:\n")
cat(file.path(output_folder, "04_Tables", "future_BCCA_PT_MIROC6_basin_mean_timeseries.csv"), "\n")

