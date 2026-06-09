# ============================================================
# RESAMPLE HISTORICAL AND FUTURE MIROC6 TO CHIRPS GRID
# Study area: Upper Blue Nile Basin
# Reference grid: XGBoost-corrected CHIRPS
# GCM: MIROC6
# Historical: 1994–2014
# Future: 2030–2080
# Scenarios: SSP2-4.5 and SSP5-8.5
# ============================================================

# ============================================================
# 1. INSTALL AND LOAD PACKAGES
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

# XGBoost-corrected CHIRPS reference grid
chirps_folder <- "D:/Desktop/HWRM_Thesis/From Gh/XGBOOST_corrected"

# Raw historical MIROC6 folder
historical_gcm_folder <- "D:/Desktop/HWRM_Thesis/From Gh/NEX-GDDP-CMIP6"

# Raw future MIROC6 folder
future_gcm_folder <- "D:/Desktop/HWRM_Thesis/From Gh/NEX_GDDP_CMIP6_Future/MIROC6"

# Basin boundary folder
basin_folder <- "D:/Desktop/HWRM_Thesis/From Gh/Gebres/TM proposal/data/CHIRPS_Monthly_Abay/XGBOOST_corrected/abay"

# Output folder
output_folder <- "D:/Desktop/HWRM_Thesis/From Gh/Resampled_MIROC6_to_CHIRPS"

historical_output_folder <- file.path(output_folder, "Historical_MIROC6_1994_2014")
future_output_folder <- file.path(output_folder, "Future_MIROC6_2030_2080")
table_output_folder <- file.path(output_folder, "Tables")

dir.create(historical_output_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(future_output_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(table_output_folder, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 3. SETTINGS
# ============================================================

selected_gcm <- "MIROC6"

ssp_scenarios <- c("ssp245", "ssp585")

historical_start <- 1994
historical_end <- 2014

future_start <- 2030
future_end <- 2080

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
  warning("No basin shapefile found. Resampling will continue without basin masking.")
  basin_vect <- NULL
} else {
  basin_shapefile <- basin_shp_files[1]
  cat("Using basin shapefile:\n", basin_shapefile, "\n")
  basin_vect <- terra::vect(basin_shapefile)
}

# ============================================================
# 5. EXTRACT YEAR-MONTH FROM FILE NAME
# Works with: 1994_01, 1994-01, 1994.01, 199401
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
# 6. MAKE INVENTORY FUNCTION
# ============================================================

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

# ============================================================
# 7. PREPARE CHIRPS REFERENCE GRID
# ============================================================

chirps_inv <- make_inventory(chirps_folder, "CHIRPS_XGBoost")

chirps_inv <- chirps_inv[
  chirps_inv$year >= historical_start &
    chirps_inv$year <= historical_end,
]

if (nrow(chirps_inv) == 0) {
  stop("ERROR: No CHIRPS files found for 1994–2014.")
}

cat("CHIRPS files found:", nrow(chirps_inv), "\n")

# Use the first CHIRPS raster as reference grid
ref_raster <- terra::rast(chirps_inv$file[1])

# Project basin to CHIRPS CRS
if (!is.null(basin_vect)) {
  basin_vect <- terra::project(basin_vect, terra::crs(ref_raster))
  ref_raster <- terra::mask(terra::crop(ref_raster, basin_vect), basin_vect)
}

terra::writeRaster(
  ref_raster,
  file.path(output_folder, "CHIRPS_reference_grid.tif"),
  overwrite = TRUE
)

cat("\nReference CHIRPS grid:\n")
print(ref_raster)

# ============================================================
# 8. RESAMPLE FUNCTION
# ============================================================

resample_to_chirps <- function(input_file, output_file, ref_raster, basin_vect = NULL) {
  
  r <- terra::rast(input_file)
  
  # Reproject if CRS is different
  if (!terra::same.crs(r, ref_raster)) {
    r <- terra::project(r, ref_raster, method = "bilinear")
  }
  
  # Resample to CHIRPS grid
  if (!terra::compareGeom(r, ref_raster, stopOnError = FALSE)) {
    r <- terra::resample(r, ref_raster, method = "bilinear")
  }
  
  # Crop and mask to basin if available
  if (!is.null(basin_vect)) {
    r <- terra::mask(terra::crop(r, basin_vect), basin_vect)
  }
  
  # Final alignment check
  if (!terra::compareGeom(r, ref_raster, stopOnError = FALSE)) {
    r <- terra::resample(r, ref_raster, method = "bilinear")
  }
  
  # Remove negative rainfall values
  r[r < 0] <- 0
  
  terra::writeRaster(
    r,
    output_file,
    overwrite = TRUE,
    wopt = list(datatype = "FLT4S", gdal = "COMPRESS=LZW")
  )
  
  return(output_file)
}

# ============================================================
# 9. RESAMPLE HISTORICAL MIROC6
# ============================================================

hist_inv_all <- make_inventory(historical_gcm_folder, "Historical_MIROC6")

hist_inv <- hist_inv_all[
  grepl(selected_gcm, hist_inv_all$file, ignore.case = TRUE) &
    hist_inv_all$year >= historical_start &
    hist_inv_all$year <= historical_end,
]

cat("\nHistorical MIROC6 files found:", nrow(hist_inv), "\n")

if (nrow(hist_inv) == 0) {
  stop("ERROR: No historical MIROC6 files found for 1994–2014.")
}

hist_resampled_table <- data.frame()

for (i in 1:nrow(hist_inv)) {
  
  date_id <- hist_inv$date_id[i]
  
  output_file <- file.path(
    historical_output_folder,
    paste0(selected_gcm, "_historical_", date_id, "_resampled_CHIRPS.tif")
  )
  
  if (!file.exists(output_file)) {
    resample_to_chirps(
      input_file = hist_inv$file[i],
      output_file = output_file,
      ref_raster = ref_raster,
      basin_vect = basin_vect
    )
  }
  
  hist_resampled_table <- rbind(
    hist_resampled_table,
    data.frame(
      GCM = selected_gcm,
      date_id = date_id,
      year = hist_inv$year[i],
      month = hist_inv$month[i],
      raw_file = hist_inv$file[i],
      resampled_file = output_file
    )
  )
  
  if (i %% 25 == 0) {
    cat("Historical resampled:", i, "of", nrow(hist_inv), "\n")
  }
}

write.csv(
  hist_resampled_table,
  file.path(table_output_folder, "Historical_MIROC6_resampled_inventory_1994_2014.csv"),
  row.names = FALSE
)

# ============================================================
# 10. RESAMPLE FUTURE MIROC6
# ============================================================

future_inv_all <- make_inventory(future_gcm_folder, "Future_MIROC6")

future_inv <- future_inv_all[
  grepl(selected_gcm, future_inv_all$file, ignore.case = TRUE) &
    future_inv_all$year >= future_start &
    future_inv_all$year <= future_end,
]

cat("\nFuture MIROC6 files found:", nrow(future_inv), "\n")

if (nrow(future_inv) == 0) {
  stop("ERROR: No future MIROC6 files found for 2030–2080.")
}

future_resampled_table <- data.frame()

for (ssp in ssp_scenarios) {
  
  ssp_inv <- future_inv[
    grepl(ssp, future_inv$file, ignore.case = TRUE),
  ]
  
  if (nrow(ssp_inv) == 0) {
    warning(paste("No future files found for", ssp))
    next
  }
  
  ssp_inv <- ssp_inv[order(ssp_inv$date_id), ]
  
  ssp_out_folder <- file.path(
    future_output_folder,
    selected_gcm,
    ssp
  )
  
  dir.create(ssp_out_folder, recursive = TRUE, showWarnings = FALSE)
  
  cat("\nResampling future scenario:", ssp, "\n")
  cat("Files:", nrow(ssp_inv), "\n")
  
  for (i in 1:nrow(ssp_inv)) {
    
    date_id <- ssp_inv$date_id[i]
    
    output_file <- file.path(
      ssp_out_folder,
      paste0(selected_gcm, "_", ssp, "_", date_id, "_resampled_CHIRPS.tif")
    )
    
    if (!file.exists(output_file)) {
      resample_to_chirps(
        input_file = ssp_inv$file[i],
        output_file = output_file,
        ref_raster = ref_raster,
        basin_vect = basin_vect
      )
    }
    
    future_resampled_table <- rbind(
      future_resampled_table,
      data.frame(
        GCM = selected_gcm,
        SSP = ssp,
        date_id = date_id,
        year = ssp_inv$year[i],
        month = ssp_inv$month[i],
        raw_file = ssp_inv$file[i],
        resampled_file = output_file
      )
    )
    
    if (i %% 25 == 0) {
      cat("Future", ssp, "resampled:", i, "of", nrow(ssp_inv), "\n")
    }
  }
}

write.csv(
  future_resampled_table,
  file.path(table_output_folder, "Future_MIROC6_resampled_inventory_2030_2080.csv"),
  row.names = FALSE
)

# ============================================================
# 11. FINAL GEOMETRY CHECK
# ============================================================

test_hist <- terra::rast(hist_resampled_table$resampled_file[1])
test_future <- terra::rast(future_resampled_table$resampled_file[1])

hist_match <- terra::compareGeom(test_hist, ref_raster, stopOnError = FALSE)
future_match <- terra::compareGeom(test_future, ref_raster, stopOnError = FALSE)

cat("\n================================================\n")
cat("FINAL RESAMPLING SUMMARY\n")
cat("================================================\n")

cat("Historical MIROC6 resampled files:", nrow(hist_resampled_table), "\n")
cat("Future MIROC6 resampled files:", nrow(future_resampled_table), "\n\n")

cat("Historical output matches CHIRPS grid:", hist_match, "\n")
cat("Future output matches CHIRPS grid:", future_match, "\n\n")

cat("Reference CHIRPS grid saved at:\n")
cat(file.path(output_folder, "CHIRPS_reference_grid.tif"), "\n\n")

cat("Historical resampled files saved at:\n")
cat(historical_output_folder, "\n\n")

cat("Future resampled files saved at:\n")
cat(future_output_folder, "\n\n")

cat("Inventory tables saved at:\n")
cat(table_output_folder, "\n")

if (hist_match == TRUE && future_match == TRUE) {
  cat("\nRESULT:\n")
  cat("Historical and future MIROC6 were successfully resampled to the CHIRPS grid.\n")
  cat("Now use these resampled folders for BCCA-PT downscaling.\n")
} else {
  cat("\nWARNING:\n")
  cat("Some outputs still do not match CHIRPS. Check CRS, extent, or basin shapefile alignment.\n")
}

