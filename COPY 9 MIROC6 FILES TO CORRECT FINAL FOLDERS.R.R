# ============================================================
# COPY 9 MIROC6 FILES TO CORRECT FINAL FOLDERS
# Source folder:
# D:/Desktop/MIROC6_missing_9_all/MIROC6_missing_9_drive_export
# ============================================================

library(terra)
library(dplyr)

source_folder <- "D:/Desktop/MIROC6_missing_9_all/MIROC6_missing_9_drive_export"

destination_root <- "D:/Desktop/HWRM_Thesis/From Gh/NEX_GDDP_CMIP6_Future_FULL_REDONE/MIROC6"

source_files <- list.files(
  path = source_folder,
  pattern = "\\.(tif|tiff)$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

cat("TIF files found in source folder:", length(source_files), "\n")

if (length(source_files) == 0) {
  stop("ERROR: No tif files found in the source folder. Check the path.")
}

expected_files <- data.frame(
  scenario = c("ssp245","ssp245","ssp245","ssp245","ssp245","ssp245","ssp585","ssp585","ssp585"),
  year = c(2031,2039,2040,2067,2069,2078,2032,2059,2071),
  month = c(1,12,2,12,12,2,1,2,2)
)

expected_files$file_name <- paste0(
  "MIROC6_",
  expected_files$scenario,
  "_",
  expected_files$year,
  "_",
  sprintf("%02d", expected_files$month),
  ".tif"
)

copy_log <- data.frame()

for (i in 1:nrow(expected_files)) {
  
  fname <- expected_files$file_name[i]
  scenario <- expected_files$scenario[i]
  year <- expected_files$year[i]
  
  match_file <- source_files[basename(source_files) == fname]
  
  if (length(match_file) == 0) {
    copy_log <- rbind(
      copy_log,
      data.frame(
        file_name = fname,
        status = "MISSING_IN_SOURCE_FOLDER",
        source_file = NA,
        destination_file = NA,
        valid_count = NA,
        mean_value = NA
      )
    )
    next
  }
  
  match_file <- match_file[1]
  
  if (year >= 2030 && year <= 2050) {
    period_folder <- "Near_future_2030_2050"
  } else {
    period_folder <- "Far_future_2051_2080"
  }
  
  destination_folder <- file.path(
    destination_root,
    scenario,
    period_folder
  )
  
  dir.create(destination_folder, recursive = TRUE, showWarnings = FALSE)
  
  destination_file <- file.path(destination_folder, fname)
  
  file.copy(
    from = match_file,
    to = destination_file,
    overwrite = TRUE
  )
  
  r <- rast(destination_file)
  
  valid_count <- global(!is.na(r), "sum", na.rm = TRUE)[1, 1]
  mean_value <- global(r, "mean", na.rm = TRUE)[1, 1]
  
  status <- ifelse(
    valid_count > 0 && is.finite(mean_value),
    "COPIED_AND_VALID",
    "COPIED_BUT_INVALID"
  )
  
  copy_log <- rbind(
    copy_log,
    data.frame(
      file_name = fname,
      status = status,
      source_file = match_file,
      destination_file = destination_file,
      valid_count = valid_count,
      mean_value = mean_value
    )
  )
}

log_folder <- file.path(destination_root, "Copy_Check_Log")
dir.create(log_folder, recursive = TRUE, showWarnings = FALSE)

write.csv(
  copy_log,
  file.path(log_folder, "MIROC6_missing_9_copy_log.csv"),
  row.names = FALSE
)

cat("\nCopy summary:\n")
print(copy_log[, c("file_name", "status", "valid_count", "mean_value")])

if (all(copy_log$status == "COPIED_AND_VALID")) {
  cat("\nRESULT: All 9 files were copied correctly and are valid.\n")
} else {
  cat("\nWARNING: Some files are missing or invalid. Check the copy log:\n")
  cat(file.path(log_folder, "MIROC6_missing_9_copy_log.csv"), "\n")
}

