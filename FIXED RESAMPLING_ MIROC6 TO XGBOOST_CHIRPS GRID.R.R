# ============================================================
# FIXED RESAMPLING: MIROC6 TO XGBOOST-CHIRPS GRID
# Uses complete and valid future MIROC6 folder
# Keeps full CHIRPS grid size: 113 x 120
# Does NOT crop output grid
# ============================================================

library(terra)
library(dplyr)

# ============================================================
# 1. INPUT FOLDERS
# ============================================================

# Corrected XGBoost-CHIRPS grid
chirps_folder <- "D:/Desktop/HWRM_Thesis/From Gh/XGBOOST_corrected"

# Raw historical MIROC6 folder
historical_gcm_folder <- "D:/Desktop/HWRM_Thesis/From Gh/NEX-GDDP-CMIP6"

# Complete and valid future MIROC6 folder after interpolation
future_gcm_folder <- "D:/Desktop/HWRM_Thesis/From Gh/NEX_GDDP_CMIP6_Future_FULL_REDONE/MIROC6"

# Basin shapefile folder
basin_folder <- "D:/Desktop/HWRM_Thesis/From Gh/Gebres/TM proposal/data/CHIRPS_Monthly_Abay/XGBOOST_corrected/abay"

# Output folder for fixed resampled data
output_folder <- "D:/Desktop/HWRM_Thesis/From Gh/Resampled_MIROC6_to_CHIRPS_FIXED"

historical_output_folder <- file.path(output_folder, "Historical_MIROC6_1994_2014")
future_output_folder <- file.path(output_folder, "Future_MIROC6_2030_2080")
table_output_folder <- file.path(output_folder, "Tables")

dir.create(historical_output_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(future_output_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(table_output_folder, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 2. SETTINGS
# ============================================================

selected_gcm <- "MIROC6"

ssp_scenarios <- c("ssp245", "ssp585")

historical_start <- 1994
historical_end <- 2014

future_start <- 2030
future_end <- 2080

# ============================================================
# 3. FUNCTIONS
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
  
  # Remove backup/log files so they are not counted as climate data
  files <- files[
    !grepl("ZERO_BACKUP", files, ignore.case = TRUE) &
      !grepl("Copy_Check_Log", files, ignore.case = TRUE) &
      !grepl("Interpolated_9_Files_Log", files, ignore.case = TRUE) &
      !grepl("Check_Result", files, ignore.case = TRUE)
  ]
  
  if (length(files) == 0) {
    stop(paste("ERROR: No valid tif files found in:", folder))
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

# ============================================================
# 4. PREPARE CHIRPS REFERENCE GRID
# ============================================================

chirps_inv <- make_inventory(chirps_folder, "XGBoost_CHIRPS")

chirps_inv <- chirps_inv[
  chirps_inv$year >= historical_start &
    chirps_inv$year <= historical_end,
]

if (nrow(chirps_inv) == 0) {
  stop("ERROR: No CHIRPS files found for 1994–2014.")
}

ref_raster <- rast(chirps_inv$file[1])

cat("CHIRPS reference grid:\n")
print(ref_raster)

# ============================================================
# 5. LOAD BASIN AND CREATE MASK ON CHIRPS GRID
# ============================================================

basin_shp <- list.files(
  basin_folder,
  pattern = "\\.shp$",
  full.names = TRUE,
  recursive = TRUE
)[1]

if (is.na(basin_shp)) {
  stop("ERROR: No basin shapefile found.")
}

basin <- vect(basin_shp)
basin <- project(basin, crs(ref_raster))

# Rasterize basin to exact CHIRPS grid
basin_mask <- rasterize(
  basin,
  ref_raster,
  field = 1,
  touches = TRUE
)

writeRaster(
  basin_mask,
  file.path(output_folder, "basin_mask_on_CHIRPS_grid.tif"),
  overwrite = TRUE
)

cat("\nBasin mask grid:\n")
print(basin_mask)

# ============================================================
# 6. FIXED RESAMPLING FUNCTION
# ============================================================

resample_keep_chirps_grid <- function(input_file, output_file, ref_raster, basin_mask) {
  
  r <- rast(input_file)
  
  # Reproject if CRS is different
  if (!same.crs(r, ref_raster)) {
    r <- project(r, ref_raster, method = "bilinear")
  }
  
  # Resample to exact CHIRPS grid
  r_resampled <- resample(r, ref_raster, method = "bilinear")
  
  # Apply basin mask WITHOUT cropping raster extent
  r_masked <- mask(r_resampled, basin_mask)
  
  # Remove negative rainfall values
  r_masked[r_masked < 0] <- 0
  
  # Confirm geometry remains same as CHIRPS
  if (!compareGeom(r_masked, ref_raster, stopOnError = FALSE)) {
    stop("ERROR: Output raster geometry changed.")
  }
  
  writeRaster(
    r_masked,
    output_file,
    overwrite = TRUE,
    wopt = list(datatype = "FLT4S", gdal = "COMPRESS=LZW")
  )
  
  valid_cells <- global(!is.na(r_masked), "sum", na.rm = TRUE)[1, 1]
  basin_cells <- global(!is.na(basin_mask), "sum", na.rm = TRUE)[1, 1]
  valid_percent <- 100 * valid_cells / basin_cells
  
  return(data.frame(
    output_file = output_file,
    valid_cells = valid_cells,
    basin_cells = basin_cells,
    valid_percent_inside_basin = valid_percent
  ))
}

# ============================================================
# 7. RESAMPLE HISTORICAL MIROC6
# ============================================================

hist_inv_all <- make_inventory(historical_gcm_folder, "Historical_MIROC6")

hist_inv <- hist_inv_all[
  grepl(selected_gcm, hist_inv_all$file, ignore.case = TRUE) &
    hist_inv_all$year >= historical_start &
    hist_inv_all$year <= historical_end,
]

cat("\nHistorical MIROC6 files found:", nrow(hist_inv), "\n")

if (nrow(hist_inv) == 0) {
  stop("ERROR: No historical MIROC6 files found.")
}

hist_table <- data.frame()

for (i in 1:nrow(hist_inv)) {
  
  date_id <- hist_inv$date_id[i]
  
  out_file <- file.path(
    historical_output_folder,
    paste0(selected_gcm, "_historical_", date_id, "_resampled_CHIRPS.tif")
  )
  
  check <- resample_keep_chirps_grid(
    input_file = hist_inv$file[i],
    output_file = out_file,
    ref_raster = ref_raster,
    basin_mask = basin_mask
  )
  
  hist_table <- rbind(
    hist_table,
    data.frame(
      GCM = selected_gcm,
      date_id = date_id,
      year = hist_inv$year[i],
      month = hist_inv$month[i],
      raw_file = hist_inv$file[i],
      resampled_file = out_file,
      valid_cells = check$valid_cells,
      basin_cells = check$basin_cells,
      valid_percent_inside_basin = check$valid_percent_inside_basin
    )
  )
  
  if (i %% 25 == 0) {
    cat("Historical completed:", i, "of", nrow(hist_inv), "\n")
  }
}

write.csv(
  hist_table,
  file.path(table_output_folder, "Historical_MIROC6_resampled_inventory_FIXED.csv"),
  row.names = FALSE
)

# ============================================================
# 8. RESAMPLE FUTURE MIROC6
# ============================================================

future_inv_all <- make_inventory(future_gcm_folder, "Future_MIROC6")

future_inv <- future_inv_all[
  grepl(selected_gcm, future_inv_all$file, ignore.case = TRUE) &
    future_inv_all$year >= future_start &
    future_inv_all$year <= future_end,
]

cat("\nFuture MIROC6 files found:", nrow(future_inv), "\n")

if (nrow(future_inv) == 0) {
  stop("ERROR: No future MIROC6 files found.")
}

# Check expected number
cat("SSP245 files:", sum(grepl("ssp245", future_inv$file, ignore.case = TRUE)), "\n")
cat("SSP585 files:", sum(grepl("ssp585", future_inv$file, ignore.case = TRUE)), "\n")

if (nrow(future_inv) != 1224) {
  warning("Expected 1224 future files. Check future inventory table.")
}

future_table <- data.frame()

for (ssp in ssp_scenarios) {
  
  ssp_inv <- future_inv[
    grepl(ssp, future_inv$file, ignore.case = TRUE),
  ]
  
  if (nrow(ssp_inv) == 0) {
    warning(paste("No future files found for", ssp))
    next
  }
  
  ssp_out_folder <- file.path(
    future_output_folder,
    selected_gcm,
    ssp
  )
  
  dir.create(ssp_out_folder, recursive = TRUE, showWarnings = FALSE)
  
  cat("\nProcessing:", ssp, "\n")
  cat("Files:", nrow(ssp_inv), "\n")
  
  for (i in 1:nrow(ssp_inv)) {
    
    date_id <- ssp_inv$date_id[i]
    
    out_file <- file.path(
      ssp_out_folder,
      paste0(selected_gcm, "_", ssp, "_", date_id, "_resampled_CHIRPS.tif")
    )
    
    check <- resample_keep_chirps_grid(
      input_file = ssp_inv$file[i],
      output_file = out_file,
      ref_raster = ref_raster,
      basin_mask = basin_mask
    )
    
    future_table <- rbind(
      future_table,
      data.frame(
        GCM = selected_gcm,
        SSP = ssp,
        date_id = date_id,
        year = ssp_inv$year[i],
        month = ssp_inv$month[i],
        raw_file = ssp_inv$file[i],
        resampled_file = out_file,
        valid_cells = check$valid_cells,
        basin_cells = check$basin_cells,
        valid_percent_inside_basin = check$valid_percent_inside_basin
      )
    )
    
    if (i %% 25 == 0) {
      cat(ssp, "completed:", i, "of", nrow(ssp_inv), "\n")
    }
  }
}

write.csv(
  future_table,
  file.path(table_output_folder, "Future_MIROC6_resampled_inventory_FIXED.csv"),
  row.names = FALSE
)

# ============================================================
# 9. FINAL CHECK
# ============================================================

test_hist <- rast(hist_table$resampled_file[1])
test_future <- rast(future_table$resampled_file[1])

cat("\n================================================\n")
cat("FINAL CHECK\n")
cat("================================================\n")

cat("CHIRPS grid:\n")
print(ref_raster)

cat("\nHistorical output grid:\n")
print(test_hist)

cat("\nFuture output grid:\n")
print(test_future)

cat("\nHistorical matches CHIRPS:\n")
print(compareGeom(test_hist, ref_raster, stopOnError = FALSE))

cat("\nFuture matches CHIRPS:\n")
print(compareGeom(test_future, ref_raster, stopOnError = FALSE))

cat("\nMinimum valid percent in historical files:\n")
print(min(hist_table$valid_percent_inside_basin, na.rm = TRUE))

cat("\nMinimum valid percent in future files:\n")
print(min(future_table$valid_percent_inside_basin, na.rm = TRUE))

cat("\nFuture file count after resampling:\n")
print(table(future_table$SSP))

cat("\nDONE SUCCESSFULLY!\n")
cat("Fixed resampled files saved in:\n")
cat(output_folder, "\n")