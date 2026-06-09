# ============================================================
# REPLACE 9 ZERO-VALUE MIROC6 FUTURE FILES BY TEMPORAL INTERPOLATION
# Use after Drive export produced copied files with mean = 0
# ============================================================

library(terra)
library(dplyr)

future_folder <- "D:/Desktop/HWRM_Thesis/From Gh/NEX_GDDP_CMIP6_Future_FULL_REDONE/MIROC6"

log_folder <- file.path(future_folder, "Interpolated_9_Files_Log")
dir.create(log_folder, recursive = TRUE, showWarnings = FALSE)

missing_months <- data.frame(
  scenario = c(
    "ssp245", "ssp245", "ssp245", "ssp245", "ssp245", "ssp245",
    "ssp585", "ssp585", "ssp585"
  ),
  year = c(
    2031, 2039, 2040, 2067, 2069, 2078,
    2032, 2059, 2071
  ),
  month = c(
    1, 12, 2, 12, 12, 2,
    1, 2, 2
  )
)

find_file <- function(folder, scenario, year, month) {
  
  pattern <- paste0(
    "MIROC6_",
    scenario,
    "_",
    year,
    "_",
    sprintf("%02d", month),
    "\\.tif$"
  )
  
  files <- list.files(
    path = folder,
    pattern = pattern,
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  
  if (length(files) == 0) {
    return(NA)
  }
  
  return(files[1])
}

get_prev_next <- function(year, month) {
  
  prev_year <- year
  prev_month <- month - 1
  
  if (prev_month == 0) {
    prev_month <- 12
    prev_year <- year - 1
  }
  
  next_year <- year
  next_month <- month + 1
  
  if (next_month == 13) {
    next_month <- 1
    next_year <- year + 1
  }
  
  return(list(
    prev_year = prev_year,
    prev_month = prev_month,
    next_year = next_year,
    next_month = next_month
  ))
}

interpolation_log <- data.frame()

for (i in 1:nrow(missing_months)) {
  
  scenario <- missing_months$scenario[i]
  year <- missing_months$year[i]
  month <- missing_months$month[i]
  
  target_file <- find_file(future_folder, scenario, year, month)
  
  if (is.na(target_file)) {
    stop(paste("Target file not found:", scenario, year, month))
  }
  
  pn <- get_prev_next(year, month)
  
  prev_file <- find_file(
    future_folder,
    scenario,
    pn$prev_year,
    pn$prev_month
  )
  
  next_file <- find_file(
    future_folder,
    scenario,
    pn$next_year,
    pn$next_month
  )
  
  if (is.na(prev_file) | is.na(next_file)) {
    stop(paste(
      "Neighbor file missing for",
      scenario,
      year,
      sprintf("%02d", month)
    ))
  }
  
  r_prev <- rast(prev_file)
  r_next <- rast(next_file)
  
  prev_mean <- global(r_prev, "mean", na.rm = TRUE)[1, 1]
  next_mean <- global(r_next, "mean", na.rm = TRUE)[1, 1]
  
  if (!is.finite(prev_mean) | !is.finite(next_mean)) {
    stop(paste(
      "Neighbor file contains invalid values for",
      scenario,
      year,
      sprintf("%02d", month)
    ))
  }
  
  # Temporal interpolation
  r_fill <- (r_prev + r_next) / 2
  
  # Backup the zero file before replacing
  backup_file <- paste0(
    tools::file_path_sans_ext(target_file),
    "_ZERO_BACKUP.tif"
  )
  
  file.copy(
    from = target_file,
    to = backup_file,
    overwrite = TRUE
  )
  
  # Replace target file
  writeRaster(
    r_fill,
    target_file,
    overwrite = TRUE,
    wopt = list(datatype = "FLT4S", gdal = "COMPRESS=LZW")
  )
  
  fill_mean <- global(r_fill, "mean", na.rm = TRUE)[1, 1]
  fill_min <- global(r_fill, "min", na.rm = TRUE)[1, 1]
  fill_max <- global(r_fill, "max", na.rm = TRUE)[1, 1]
  fill_valid <- global(!is.na(r_fill), "sum", na.rm = TRUE)[1, 1]
  
  interpolation_log <- rbind(
    interpolation_log,
    data.frame(
      scenario = scenario,
      year = year,
      month = month,
      target_file = target_file,
      backup_file = backup_file,
      previous_file = prev_file,
      next_file = next_file,
      previous_mean = prev_mean,
      next_mean = next_mean,
      interpolated_mean = fill_mean,
      interpolated_min = fill_min,
      interpolated_max = fill_max,
      valid_cells = fill_valid
    )
  )
  
  cat(
    "Interpolated:",
    scenario,
    year,
    sprintf("%02d", month),
    "mean =",
    fill_mean,
    "\n"
  )
}

write.csv(
  interpolation_log,
  file.path(log_folder, "Interpolated_9_MIROC6_Future_Files_Log.csv"),
  row.names = FALSE
)

cat("\nDONE. The 9 zero files were replaced by temporal interpolation.\n")
cat("Log saved at:\n")
cat(file.path(log_folder, "Interpolated_9_MIROC6_Future_Files_Log.csv"), "\n")

